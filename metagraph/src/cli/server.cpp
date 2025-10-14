#include <json/json.h>
#include <tsl/hopscotch_map.h>
#include <server_http.hpp>

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
#include "cli/load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;


Json::Value process_search_request(const Json::Value &json,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const Config &config_orig) {
    const auto &fasta = json["FASTA"];
    if (fasta.isNull())
        throw std::domain_error("No input sequences received from client");

    Config config(config_orig);
    // discovery_fraction a proxy of 1 - %similarity
    config.discovery_fraction
            = json.get("discovery_fraction", config.discovery_fraction).asDouble();

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    if (config.discovery_fraction < 0.0 || config.discovery_fraction > 1.0) {
        throw std::domain_error(
                "Discovery fraction should be within [0, 1.0]. Instead got "
                + std::to_string(config.discovery_fraction));
    }

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::domain_error(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.num_top_labels = json.get("top_labels", config.num_top_labels).asInt();

    if (json.get("query_coords", false).asBool()) {
        config.query_mode = COORDS;
    } else if (json.get("query_counts", false).asBool()) {
        config.query_mode = COUNTS;
    } else if (json.get("with_signature", false).asBool()) {
        config.query_mode = SIGNATURE;
    } else if (json.get("abundance_sum", false).asBool()) {
        config.query_mode = COUNTS_SUM;
    } else {
        config.query_mode = MATCHES;
    }

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
    if (json.get("align", false).asBool()) {
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
    tf.ofstream() << fasta.asString();
    tf.ofstream().close();

    // dummy pool doing everything in the caller thread
    ThreadPool dummy_pool(0);
    QueryExecutor engine(config, anno_graph, std::move(aligner_config), dummy_pool);

    // Query sequences and callback by appending result to vector with mutex for thread safety
    engine.query_fasta(tf.name(),
        [&](const SeqSearchResult &result) {
            std::lock_guard<std::mutex> lock(result_mutex);
            search_results.emplace_back(std::move(result));
        }
    );

    // Ensure JSON results are sorted by their ID
    std::sort(search_results.begin(), search_results.end(),
              [](const SeqSearchResult &lhs, const SeqSearchResult &rhs) {
                  return lhs.get_sequence().id < rhs.get_sequence().id;
              });

    // Create full JSON object
    Json::Value search_response(Json::arrayValue);
    for (const auto &seq_result : search_results) {
        search_response.append(seq_result.to_json(config.verbose_output, anno_graph));
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
        throw std::domain_error(
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

std::thread start_server(HttpServer &server_startup, Config &config) {
    server_startup.config.thread_pool_size = std::max(1u, get_num_threads());

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;
    server_startup.config.timeout_request = 30;    // 30 sec to finish headers
    server_startup.config.timeout_content = 1800;  // 30 minutes for body/compute (per request) max

    logger->info("[Server] Will listen on {} port {}",
                 server_startup.config.address, server_startup.config.port);
    logger->info("[Server] Maximum connections: {}", get_num_threads());
    return std::thread([&server_startup]() { server_startup.start(); });
}

template<typename T>
bool check_data_ready(std::shared_future<T> &data, shared_ptr<HttpServer::Response> response) {
    if (data.wait_for(0s) != std::future_status::ready) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        response->write(SimpleWeb::StatusCode::server_error_service_unavailable,
                        "Server is currently initializing, please come back later.");
        return false;
    }

    return true;
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
    std::atomic<size_t> num_requests = 0;

    ThreadPool graph_loader(1, 1);
    std::shared_future<std::unique_ptr<AnnotatedDBG>> anno_graph;

    tsl::hopscotch_map<std::string, std::vector<std::pair<std::string, std::string>>> indexes;

    if (config->infbase_annotators.size() == 1) {
        assert(config->fnames.empty());
        anno_graph = graph_loader.enqueue([&]() {
            logger->info("[Server] Loading graph...");
            auto graph = load_critical_dbg(config->infbase);
            logger->info("[Server] Graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);

            auto anno_graph = initialize_annotated_dbg(graph, *config);
            logger->info("[Server] Annotated graph loaded too. Current mem usage: {} MiB", get_curr_RSS() >> 20);
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
        logger->info("[Server] Loaded paths for {} graphs for {} names: {}",
                     num_indexes, indexes.size(), fmt::join(names, ", "));
        if (!utils::with_mmap()) {
            logger->warn("[Server] --mmap wasn't passed but all indexes will be loaded with mmap."
                         " Make sure they're on a fast disk.");
            utils::with_mmap(true);
        }
    }

    ThreadPool graphs_pool(get_num_threads());

    logger->info("Collecting graph stats...");
    tsl::hopscotch_map<std::string, std::vector<std::string>> name_labels;
    for (const auto &[name, graphs] : indexes) {
        for (const auto &[graph_fname, anno_fname] : graphs) {
            auto &out = name_labels[name];
            const auto &labels = read_labels(anno_fname);
            out.insert(out.end(), labels.begin(), labels.end());
        }
    }

    logger->info("All graphs were loaded and stats collected. Ready to serve queries.");

    // the actual server
    HttpServer server;
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->info("[Server] {} request {} from {}", request->path, request_id,
                     request->remote_endpoint().address().to_string());

        if (!config->fnames.size() && !check_data_ready(anno_graph, response))
            return;  // the index is not loaded yet, so we can't process the request

        process_request(response, request, [&](const std::string &content) {
            Json::Value content_json = parse_json_string(content);
            logger->info("Request {}: {}", request_id, content_json.toStyledString());
            Json::Value result;

            // simple case with a single graph pair
            if (!config->fnames.size()) {
                if (content_json.isMember("graphs"))
                    throw std::invalid_argument("Bad request: no support for filtering graphs on this server");
                result = process_search_request(content_json, *anno_graph.get(), *config);
            } else {
                std::vector<std::string> graphs_to_query
                        = filter_graphs_from_list(indexes, content_json, request_id);
                std::mutex mu;
                std::vector<std::shared_future<void>> futures;
                for (const auto &name : graphs_to_query) {
                    for (const auto &[graph_fname, anno_fname] : indexes[name]) {
                        Config config_copy = *config;
                        config_copy.infbase = graph_fname;
                        config_copy.infbase_annotators = { anno_fname };
                        futures.push_back(graphs_pool.enqueue([config_copy{std::move(config_copy)},&content_json,&result,&mu]() {
                            auto index = initialize_annotated_dbg(config_copy);
                            auto json = process_search_request(content_json, *index, config_copy);
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
                        }));
                    }
                }
                for (auto &future : futures) {
                    future.wait();
                }
            }
            logger->info("Request {} finished", request_id);
            return result;
        });
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->info("[Server] {} request {} from {}", request->path, request_id,
                     request->remote_endpoint().address().to_string());

        if (!config->fnames.size() && !check_data_ready(anno_graph, response))
            return;  // the index is not loaded yet, so we can't process the request

        process_request(response, request, [&](const std::string &content) {
            if (!config->fnames.size())
                return process_align_request(content, anno_graph.get()->get_graph(), *config);

            throw std::invalid_argument("Bad request: alignment requests are not yet supported for "
                                        "servers with multiple graphs");
        });
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->info("[Server] {} request {} from {}", request->path, request_id,
                     request->remote_endpoint().address().to_string());

        if (!config->fnames.size() && !check_data_ready(anno_graph, response))
            return;  // the index is not loaded yet, so we can't process the request

        process_request(response, request, [&](const std::string &) {
            Json::Value root(Json::arrayValue);
            if (!config->fnames.size()) {
                auto labels = anno_graph.get()->get_annotator().get_label_encoder().get_labels();
                for (const std::string &label : labels) {
                    root.append(label);
                }
            } else {
                for (const auto &[name, labels] : name_labels) {
                    for (const std::string &label : labels) {
                        root.append(label);
                    }
                }
            }
            return root;
        });
    };

    server.resource["^/stats"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                            shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->info("[Server] {} request {} from {}", request->path, request_id,
                     request->remote_endpoint().address().to_string());

        if (!config->fnames.size() && !check_data_ready(anno_graph, response))
            return;  // the index is not loaded yet, so we can't process the request

        process_request(response, request, [&](const std::string &) {
            Json::Value root;
            if (config->fnames.size()) {
                // for scenarios with multiple graphs
                uint64_t num_labels = 0;
                for (const auto &[name, labels] : name_labels) {
                    num_labels += labels.size();
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
                root["annotation"]["labels"] = static_cast<uint64_t>(annotation.num_labels());
                root["annotation"]["objects"] = static_cast<uint64_t>(annotation.num_objects());
                root["annotation"]["relations"] = static_cast<uint64_t>(annotation.num_relations());
            }
            return root;
        });
    };

    server.default_resource["GET"] = [](shared_ptr<HttpServer::Response> response,
                                        shared_ptr<HttpServer::Request> request) {
        logger->info("Not found " + request->path);
        response->write(SimpleWeb::StatusCode::client_error_not_found,
                        "Could not find path " + request->path);
    };
    server.default_resource["POST"] = server.default_resource["GET"];

    server.on_error = [](shared_ptr<HttpServer::Request> /*request*/,
                         const SimpleWeb::error_code &ec) {
        // Handle errors here, ignoring a few trivial ones.
        if (ec.value() != asio::stream_errc::eof
                && ec.value() != asio::error::operation_aborted) {
            logger->info("[Server] Got error {} {} {}",
                         ec.message(), ec.category().name(), ec.value());
        }
    };

    std::thread server_thread = start_server(server, *config);
    server_thread.join();

    return 0;
}

} // namespace cli
} // namespace mtg
