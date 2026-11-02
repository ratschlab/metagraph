#include <json/json.h>
#include <tsl/hopscotch_map.h>
#include <server_http.hpp>
#include <csignal>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "query.hpp"
#include "align.hpp"
#include "server_utils.hpp"
#include "server_cache.hpp"
#include "server_request.hpp"
#include "cli/load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;


Json::Value process_search_request(const SearchRequest &request,
                                   const graph::AnnotatedDBG &anno_graph) {
    const auto &config = request.config;

    // Throw client an error if they try to query coordinates/kmer-counts on unsupported indexes
    if ((config.query_mode == COUNTS || config.query_mode == COUNTS_SUM)
            && !dynamic_cast<const annot::matrix::IntMatrix *>(
                            &anno_graph.get_annotator().get_matrix())) {
        throw std::invalid_argument("Annotation does not support k-mer count queries");
    }

    if (config.query_mode == COORDS
            && !dynamic_cast<const annot::matrix::MultiIntMatrix *>(
                            &anno_graph.get_annotator().get_matrix())) {
        throw std::invalid_argument("Annotation does not support k-mer coordinate queries");
    }

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (request.align) {
        aligner_config.reset(new align::DBGAlignerConfig(
            initialize_aligner_config(config, anno_graph.get_graph())
        ));
    }

    // Need mutex while appending to vector
    std::vector<SeqSearchResult> search_results;
    std::mutex result_mutex;

    // writing to temporary file in order to reuse query code. This is not optimal and
    // may turn out to be an issue in production. However, adapting FastaParser to
    // work on strings seems non-trivial. An alternative would be to use
    // read_fasta_from_string.
    utils::TempFile tf(config.tmp_dir);
    tf.ofstream() << request.fasta;
    tf.ofstream().close();

    // Query sequences and callback by appending result to vector with mutex for thread safety
    query_fasta(tf.name(),
        [&](const SeqSearchResult &result) {
            std::lock_guard<std::mutex> lock(result_mutex);
            search_results.emplace_back(std::move(result));
        },
        config, anno_graph, aligner_config.get()
    );

    // Ensure JSON results are sorted by their ID
    std::sort(search_results.begin(), search_results.end(),
              [](const SeqSearchResult &lhs, const SeqSearchResult &rhs) {
                  return lhs.get_sequence().id < rhs.get_sequence().id;
              });

    // Create full JSON object
    Json::Value search_response(Json::arrayValue);
    for (const auto &seq_result : search_results) {
        search_response.append(seq_result.to_json(config.verbose_output, anno_graph.get_graph().get_k()));
    }

    return search_response;
}

// TODO: implement alignment_result.to_json as in process_search_request
Json::Value process_align_request(const std::string &received_message,
                                  const graph::DeBruijnGraph &graph,
                                  const Config &config_orig) {
    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];

    Json::Value root = Json::Value(Json::arrayValue);

    Config config(config_orig);

    config.alignment_num_alternative_paths = json.get(
        "max_alternative_alignments",
        (uint64_t)config.alignment_num_alternative_paths).asInt();

    if (!config.alignment_num_alternative_paths) {
        // TODO: better throw an exception and send an error response to the client
        logger->warn("[Server] Got invalid value of alignment_num_alternative_paths = {}."
                     " The default value of 1 will be used instead...", config.alignment_num_alternative_paths);
        config.alignment_num_alternative_paths = 1;
    }

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::invalid_argument(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    align::DBGAligner aligner(graph, initialize_aligner_config(config, graph));
    const align::DBGAlignerConfig &aligner_config = aligner.get_config();

    // TODO: make parallel?
    seq_io::read_fasta_from_string(fasta.asString(),
                                   [&](seq_io::kseq_t *read_stream) {
        Json::Value align_entry;
        align_entry[SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD] = read_stream->name.s;

        // not supporting reverse complement yet
        Json::Value alignments = Json::Value(Json::arrayValue);

        for (const auto &path : aligner.align(read_stream->seq.s)) {
            Json::Value a;
            a[SeqSearchResult::SCORE_JSON_FIELD] = path.get_score();
            a[SeqSearchResult::MAX_SCORE_JSON_FIELD] = aligner_config.match_score(read_stream->seq.s)
                + aligner_config.left_end_bonus + aligner_config.right_end_bonus;
            aligner.get_config().match_score(read_stream->seq.s);
            a[SeqSearchResult::SEQUENCE_JSON_FIELD] = std::string(path.get_sequence());
            a[SeqSearchResult::CIGAR_JSON_FIELD] = path.get_cigar().to_string();
            a[SeqSearchResult::ORIENTATION_JSON_FIELD] = path.get_orientation();

            alignments.append(a);
        }

        align_entry[SeqSearchResult::ALIGNMENT_JSON_FIELD] = alignments;

        root.append(align_entry);
    });

    return root;
}

std::thread start_server(HttpServer &server_startup, Config &config, size_t num_threads,
                         std::function<void(unsigned short)> on_start) {
    server_startup.config.thread_pool_size = num_threads;

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;
    server_startup.config.timeout_request = 30;    // 30 sec to finish headers
    server_startup.config.timeout_content = 900;   // 15 minutes for body/compute (per request) max

    logger->info("[Server] Will listen on {} port {}",
                 server_startup.config.address, server_startup.config.port);
    logger->info("[Server] Maximum connections: {}", num_threads);
    return std::thread([&server_startup, on_start]() { server_startup.start(on_start); });
}

std::vector<std::string> filter_graphs_from_list(
        const tsl::hopscotch_map<std::string, std::vector<std::pair<std::string, std::string>>> &indexes,
        const Json::Value &content_json,
        size_t request_id,
        size_t max_names_without_filtering = 10) {
    std::vector<std::string> graphs_to_query;

    if (content_json.isMember("graphs") && content_json["graphs"].isArray()) {
        for (const auto &item : content_json["graphs"]) {
            graphs_to_query.push_back(item.asString());
            if (!indexes.count(graphs_to_query.back()))
                throw std::invalid_argument("Request with an uninitialized graph " + graphs_to_query.back());
        }
        // deduplicate
        std::sort(graphs_to_query.begin(), graphs_to_query.end());
        graphs_to_query.erase(std::unique(graphs_to_query.begin(), graphs_to_query.end()), graphs_to_query.end());
    } else {
        if (indexes.size() > max_names_without_filtering) {
            throw std::invalid_argument(
                fmt::format("Bad request: requests without names (no \"graphs\" field) are "
                            "only supported for small indexes (<={} names)",
                            max_names_without_filtering));
        }
        // query all graphs from list `config->fnames`
        for (const auto &[name, _] : indexes) {
            graphs_to_query.push_back(name);
        }
    }
    return graphs_to_query;
}


int run_server(Config *config) {
    assert(config);
    // Register before loading so termination also works when running as PID 1.
    asio::io_service shutdown_service;
    asio::signal_set shutdown_signals(shutdown_service, SIGINT, SIGTERM);
    std::atomic<size_t> num_requests = 0;

    std::shared_future<std::unique_ptr<AnnotatedDBG>> anno_graph;
    // Destroy the loader before the future holding its result.
    ThreadPool graph_loader(1, 1);

    tsl::hopscotch_map<std::string, std::vector<std::pair<std::string, std::string>>> indexes;

    ThreadPool graphs_pool(get_num_threads(), 1000 /* max_num_tasks */);
    size_t num_server_threads = std::max(1u, get_num_threads());
    set_num_threads(std::max(1u, config->parallel_each));
    config->parallel_each = 1;  // query one batch at a time
    logger->info("[Server] Threads per graph: {}", get_num_threads());

    VectorMap<std::pair<std::string, std::string>, std::unique_ptr<AnnotatedDBG>> graphs_cache;
    bool loaded_with_mmap = utils::with_mmap();

    if (config->infbase_annotators.size() == 1) {
        assert(config->fnames.empty());
        anno_graph = graph_loader.enqueue([&]() {
            logger->info("[Server] Loading graph and annotation in parallel...");
            auto anno_graph = initialize_annotated_dbg(*config);
            logger->info("[Server] Annotated graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);
            return anno_graph;
        });
    } else {
        assert(config->fnames.size() == 1);

        std::ifstream file(config->fnames[0]);

        if (!file.is_open()) {
            logger->error("[Server] Could not open file {} for reading", config->fnames[0]);
            std::exit(1);
        }

        size_t num_indexes = 0;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty())
                continue; // skip empty lines

            std::stringstream ss(line);

            std::string name;
            std::string graph_path;
            std::string annotation_path;

            std::getline(ss, name, ',');
            std::getline(ss, graph_path, ',');
            if (ss.eof()) {
                logger->error("[Server] Invalid line in the csv file: `{}`", line);
                std::exit(1);
            }
            std::getline(ss, annotation_path, ',');

            indexes[name].emplace_back(std::move(graph_path), std::move(annotation_path));
            num_indexes++;
        }
        std::vector<std::string> names;
        for (const auto &[name, _] : indexes) {
            names.push_back(name);
        }
        logger->info("[Server] Loaded a list of {} graphs for {} names: {}",
                     num_indexes, indexes.size(), fmt::join(names, ", "));
        if (loaded_with_mmap) {
            logger->info("[Server] Graphs will be loaded with mmap (--mmap set)."
                         " Make sure they're on a fast drive.");
        } else {
            logger->info("[Server] Graphs will be loaded into RAM (--mmap not set)."
                         " Total memory ≈ sum of all graph+annotation file sizes.");
        }

        // Deduplicate graph paths so the same underlying graph file is loaded only
        // once, even if it appears in multiple (graph, annotation) entries.
        std::vector<std::string> unique_graph_paths;
        tsl::hopscotch_map<std::string, size_t> graph_path_to_idx;
        for (const auto &[name, graphs] : indexes) {
            for (const auto &[graph_fname, anno_fname] : graphs) {
                if (graph_path_to_idx.emplace(graph_fname, unique_graph_paths.size()).second)
                    unique_graph_paths.push_back(graph_fname);
                graphs_cache[{ graph_fname, anno_fname }] = nullptr;
            }
        }
        if (graphs_cache.empty()) {
            logger->error("[Server] No graphs to serve. Exiting.");
            exit(1);
        }
        if (unique_graph_paths.size() < graphs_cache.size()) {
            logger->info("[Server] Deduplicated {} graph references down to {} unique graph files.",
                         graphs_cache.size(), unique_graph_paths.size());
        }

        logger->info("[Server] Loading {} unique graph(s)...", unique_graph_paths.size());
        std::vector<std::shared_ptr<DeBruijnGraph>> loaded_graphs(unique_graph_paths.size());
        #pragma omp parallel for num_threads(get_num_threads() * num_server_threads) schedule(dynamic)
        for (size_t i = 0; i < unique_graph_paths.size(); ++i) {
            loaded_graphs[i] = load_critical_dbg(unique_graph_paths[i]);
        }

        logger->info("[Server] Loading {} annotation(s)...", graphs_cache.size());
        #pragma omp parallel for num_threads(get_num_threads() * num_server_threads) schedule(dynamic)
        for (size_t i = 0; i < graphs_cache.size(); ++i) {
            Config config_copy = *config;
            auto it = graphs_cache.nth(i);
            const auto &[graph_fname, anno_fname] = it.key();
            config_copy.infbase = graph_fname;
            config_copy.infbase_annotators = { anno_fname };
            it.value() = initialize_annotated_dbg(loaded_graphs[graph_path_to_idx.at(graph_fname)],
                                                  config_copy);
        }
        logger->info("[Server] All graphs were loaded ({}). Ready to serve queries.",
                     loaded_with_mmap ? "with mmap" : "into RAM");
        // Dynamic per-request loads (in_ram path) should always be in RAM.
        utils::set_mmap(false);
    }

    size_t memory_all = config->memory_available * 1e9;
    size_t memory_left = memory_all;
    std::atomic<size_t> graphs_being_queried = 0;
    std::condition_variable space_cv;
    std::mutex space_mutex;

    size_t cache_size_bytes = static_cast<size_t>(config->server_cache_size * 1e9);
    ServerQueryCache search_cache(cache_size_bytes);
    if (cache_size_bytes) {
        logger->info("[Server] /search result cache enabled, max size: {} GB",
                     config->server_cache_size);
    } else {
        logger->info("[Server] /search result cache disabled (--cache-size 0)");
    }

    HttpServer server;

    // Keep the cache handle alive until the asynchronous response write completes.
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        auto handle_holder = cache_size_bytes
                ? std::make_shared<ServerQueryCache::Handle>() : nullptr;
        process_request(response, request, request_id,
            [&, handle_holder](const std::string &content) -> std::string {
                if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                    throw CurrentlyInitializingError();  // graph still loading

                Json::Value content_json = parse_json_string(content);
                if (logger->should_log(spdlog::level::trace))
                    logger->trace("[Server] Request {}: {}", request_id, content_json.toStyledString());

                const auto search_request = parse_search_request(content_json, *config);
                // Validate fields omitted from the key before returning a cached result.
                std::vector<std::string> graphs_to_query;
                bool load_in_ram = false;
                if (config->fnames.empty()) {
                    if (content_json.isMember("graphs"))
                        throw std::invalid_argument("Bad request: no support for filtering graphs on this server");
                } else {
                    graphs_to_query = filter_graphs_from_list(indexes, content_json, request_id);
                    load_in_ram = content_json.get("in_ram", false).asBool() && loaded_with_mmap;
                }

                if (handle_holder) {
                    std::string graph_identity = config->infbase;
                    if (!config->fnames.empty()) {
                        std::ostringstream oss;
                        for (const auto &name : graphs_to_query) {
                            for (const auto &[g, a] : indexes[name]) {
                                oss << g.size() << ':' << g << a.size() << ':' << a;
                            }
                        }
                        graph_identity = oss.str();
                    }
                    *handle_holder = search_cache.acquire(
                        make_search_cache_key(search_request, graph_identity));
                    if (!handle_holder->is_miss()) {
                        auto cached = handle_holder->get();
                        logger->info("[Server] Request {}: cache HIT  {:.1f} KB  "
                                     "(occupancy {:.1f}/{:.0f} MB, {} entries)",
                                     request_id, cached->size() / 1e3,
                                     search_cache.size_bytes() / 1e6,
                                     cache_size_bytes / 1e6,
                                     search_cache.entry_count());
                        return *cached;
                    }
                }

                // Cache miss: run the search, publish the result. If the
                // compute throws, propagate the exception to any duplicate
                // waiters via set_exception() before re-throwing.
                Json::Value result;
                try {
                    if (!config->fnames.size()) {
                        logger->trace("Request {}: Started querying graph {}, in total graphs being queried at the moment: {}",
                                      request_id, config->infbase, graphs_being_queried.fetch_add(1) + 1);
                        try {
                            result = process_search_request(search_request, *anno_graph.get());
                        } catch (...) {
                            graphs_being_queried--;
                            throw;
                        }
                        graphs_being_queried--;
                    } else {
                        std::mutex mu;
                        struct PendingQueries {
                            std::vector<std::shared_future<std::exception_ptr>> futures;
                            ~PendingQueries() {
                                // Tasks capture request-local state, including mu and result.
                                for (auto &future : futures)
                                    future.wait();
                            }
                        } pending;
                        auto &futures = pending.futures;
                        size_t num_queries = 0;
                        for (const auto &name : graphs_to_query)
                            num_queries += indexes.at(name).size();
                        // Once a task is launched, recording its future must not allocate.
                        futures.reserve(num_queries);
                        for (const auto &name : graphs_to_query) {
                            for (const auto &[graph_fname, anno_fname] : indexes[name]) {
                                futures.push_back(graphs_pool.enqueue([&,config,graph_fname=graph_fname,anno_fname=anno_fname]() {
                                    logger->trace("Request {}: Started querying graph {}. In total graphs being queried at the moment: {}",
                                                  request_id, graph_fname, graphs_being_queried.fetch_add(1) + 1);
                                    size_t index_size_reserved = 0;
                                    auto release_memory = [&]() {
                                        if (!index_size_reserved)
                                            return;
                                        {
                                            std::unique_lock<std::mutex> lock(space_mutex);
                                            memory_left += index_size_reserved;
                                        }
                                        index_size_reserved = 0;
                                        space_cv.notify_all();
                                    };
                                    try {
                                        std::unique_ptr<AnnotatedDBG> index_loaded;
                                        const AnnotatedDBG *index;
                                        bool in_ram = load_in_ram;
                                        size_t index_size = in_ram ? std::filesystem::file_size(graph_fname)
                                                                        + std::filesystem::file_size(anno_fname)
                                                                   : -1;
                                        if (in_ram && index_size > memory_all) {
                                            logger->warn("Request {}: Graph of size {} GB is too large to fit into "
                                                         "RAM (reserved memory: {} GB). It will be queried with mmap",
                                                         request_id, index_size / 1e9, memory_all / 1e9);
                                            in_ram = false;
                                        }
                                        Timer timer;
                                        if (in_ram) {
                                            {
                                                std::unique_lock<std::mutex> lock(space_mutex);
                                                space_cv.wait(lock, [&]() {
                                                    return memory_left >= index_size;
                                                });
                                                memory_left -= index_size;
                                            }
                                            index_size_reserved = index_size;
                                            logger->trace("Request {}: Loading graph {} of size {} GB to RAM...",
                                                          request_id, graph_fname, index_size / 1e9);
                                            timer.reset();
                                            Config config_copy = *config;
                                            config_copy.infbase = graph_fname;
                                            config_copy.infbase_annotators = { anno_fname };
                                            index_loaded = initialize_annotated_dbg(config_copy);
                                            index = index_loaded.get();
                                        } else {
                                            index = graphs_cache.at({ graph_fname, anno_fname }).get();
                                        }

                                        auto json = process_search_request(search_request, *index);
                                        logger->trace("Request {}: {} graph {} {} in {} sec",
                                                      request_id, in_ram ? "Loaded and searched" : "Searched",
                                                      graph_fname, anno_fname, timer.elapsed());

                                        index_loaded.reset();
                                        release_memory();

                                        std::lock_guard<std::mutex> lock(mu);
                                        if (result.empty()) {
                                            result = std::move(json);
                                        } else {
                                            assert(json.size() == result.size());
                                            for (Json::ArrayIndex i = 0; i < result.size(); ++i) {
                                                if (result[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]
                                                        != json[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]) {
                                                    throw std::logic_error("ERROR: Results for different sequences can't be merged");
                                                }
                                                for (auto&& value : json[i]["results"]) {
                                                    result[i]["results"].append(std::move(value));
                                                }
                                            }
                                        }
                                    } catch (...) {
                                        release_memory();
                                        graphs_being_queried--;
                                        return std::current_exception();
                                    }
                                    graphs_being_queried--;
                                    return std::exception_ptr();
                                }));
                            }
                        }
                        for (auto &future : futures) {
                            if (auto ex = future.get())
                                std::rethrow_exception(ex);
                        }
                    }
                } catch (...) {
                    if (handle_holder)
                        handle_holder->set_exception(std::current_exception());
                    throw;
                }

                // Serialize once: the cached form *is* the response
                // body. Cache hits then skip serialization entirely.
                std::string serialized = Json::writeString(Json::StreamWriterBuilder(), result);
                if (handle_holder) {
                    bool retained = handle_holder->set_result(serialized);
                    logger->info("[Server] Request {}: cache {} {:.1f} KB  "
                                 "(occupancy {:.1f}/{:.0f} MB, {} entries)",
                                 request_id, retained ? "STORE" : "BYPASS", serialized.size() / 1e3,
                                 search_cache.size_bytes() / 1e6,
                                 cache_size_bytes / 1e6,
                                 search_cache.entry_count());
                }
                return serialized;
            },
            [handle_holder, request_id](const SimpleWeb::error_code &ec) {
                if (!handle_holder || !*handle_holder)
                    return;  // never acquired (e.g. early throw before acquire)
                if (ec) {
                    logger->info("[Server] Request {}: delivery FAILED ({})", request_id, ec.message());
                    handle_holder->mark_protected();
                } else {
                    handle_holder->mark_delivered();
                }
            });
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            if (config->fnames.size()) {
                throw std::invalid_argument("Bad request: alignment requests are not yet supported for "
                                            "servers with multiple graphs");
            }
            return Json::writeString(Json::StreamWriterBuilder(),
                process_align_request(content, anno_graph.get()->get_graph(), *config));
        });
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            Json::Value root(Json::arrayValue);
            if (!config->fnames.size()) {
                auto labels = anno_graph.get()->get_annotator().get_label_encoder().get_labels();
                for (const std::string &label : labels) {
                    root.append(label);
                }
            } else {
                for (const auto &[name, graphs] : indexes) {
                    for (const auto &[graph_fname, anno_fname] : graphs) {
                        const auto &labels = graphs_cache.at({ graph_fname, anno_fname })->get_annotator().get_label_encoder().get_labels();
                        for (const std::string &label : labels) {
                            root.append(label);
                        }
                    }
                }
            }
            return Json::writeString(Json::StreamWriterBuilder(), root);
        });
    };

    server.resource["^/stats"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                            shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            auto get_num_labels = [](const AnnotatedDBG &anno_dbg) {
                uint64_t num_labels = 0;
                if (const auto *coord_to_header = anno_dbg.get_coord_to_header()) {
                    for (uint64_t col = 0; col < coord_to_header->num_columns(); ++col) {
                        num_labels += coord_to_header->num_sequences(col);
                    }
                } else {
                    num_labels = anno_dbg.get_annotator().num_labels();
                }
                return num_labels;
            };

            Json::Value root;
            if (config->fnames.size()) {
                // for scenarios with multiple graphs
                const auto &graph = graphs_cache.begin()->second->get_graph();
                size_t k = graph.get_k();
                bool is_consistent_k = true;
                bool is_canonical = (graph.get_mode() == graph::DeBruijnGraph::CANONICAL);
                bool is_consistent_canonical = true;
                for (const auto &[graph_anno, anno_dbg] : graphs_cache) {
                    const auto &graph = anno_dbg->get_graph();
                    if (k != graph.get_k())
                        is_consistent_k = false;
                    if (is_canonical != (graph.get_mode() == graph::DeBruijnGraph::CANONICAL))
                        is_consistent_canonical = false;
                }
                if (is_consistent_k)
                    root["graph"]["k"] = static_cast<uint64_t>(k);
                if (is_consistent_canonical)
                    root["graph"]["is_canonical_mode"] = is_canonical;
                uint64_t num_labels = 0;
                for (const auto &[name, graphs] : indexes) {
                    for (const auto &[graph_fname, anno_fname] : graphs) {
                        num_labels += get_num_labels(*graphs_cache.at({ graph_fname, anno_fname }));
                    }
                }
                root["annotation"]["labels"] = num_labels;
            } else {
                root["graph"]["filename"] = std::filesystem::path(config->infbase).filename().string();
                root["graph"]["k"] = static_cast<uint64_t>(anno_graph.get()->get_graph().get_k());
                root["graph"]["nodes"] = anno_graph.get()->get_graph().num_nodes();
                root["graph"]["is_canonical_mode"] = (anno_graph.get()->get_graph().get_mode()
                                                        == graph::DeBruijnGraph::CANONICAL);
                const auto &annotation = anno_graph.get()->get_annotator();
                root["annotation"]["filename"] = std::filesystem::path(config->infbase_annotators.front()).filename().string();
                root["annotation"]["labels"] = get_num_labels(*anno_graph.get());
                root["annotation"]["objects"] = static_cast<uint64_t>(annotation.num_objects());
            }
            return Json::writeString(Json::StreamWriterBuilder(), root);
        });
    };

    server.default_resource["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                        shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->warn("[Server] Not found {} for {} request {} from {}",
                     request->path, request->method, request_id,
                     request->remote_endpoint().address().to_string());
        response->write(SimpleWeb::StatusCode::client_error_not_found,
                        "Could not find path " + request->path);
    };
    server.default_resource["POST"] = server.default_resource["GET"];

    server.on_error = [](shared_ptr<HttpServer::Request> /*request*/,
                         const SimpleWeb::error_code &ec) {
        // Handle errors here, ignoring a few trivial ones.
        if (ec.value() != asio::stream_errc::eof
                && ec.value() != asio::error::operation_aborted) {
            logger->warn("[Server] Got error {} {} {}",
                         ec.message(), ec.category().name(), ec.value());
        }
    };

    std::promise<void> listening;
    std::thread server_thread = start_server(server, *config, num_server_threads,
                                            [&](unsigned short port) {
        logger->info("[Server] Listening on port {}", port);
        listening.set_value();
    });
    listening.get_future().wait();
    shutdown_signals.async_wait([&](const asio::error_code &error, int) {
        if (!error) {
            logger->info("[Server] Shutting down...");
            // Only stop the HTTP workers here. Closing the acceptor and the
            // connections from this thread would race with the accept and read
            // handlers the workers are running (and could crash in asio).
            server.io_service->stop();
        }
    });
    shutdown_service.run();
    server_thread.join();
    // No HTTP worker is left: close the acceptor and the open connections.
    server.stop();

    return 0;
}

} // namespace cli
} // namespace mtg
