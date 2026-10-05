#include <map>
#include <mutex>

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
#include "traverse.hpp"
#include "traverse_attempts.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "server_utils.hpp"
#include "cli/load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

// How long the HTTP server gives a request after its header was read (the body, the handler
// and sending the response) before it shuts the connection: also the cap on the duration
// bound of a /traverse attempt (traverse_attempts.hpp), less a second for the response
constexpr long kContentTimeoutS = 900;


Json::Value process_search_request(const Json::Value &json,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const Config &config_orig) {
    const auto &fasta = json["FASTA"];
    if (fasta.isNull())
        throw std::invalid_argument("No input sequences received from client");

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
        throw std::invalid_argument(
                "Discovery fraction should be within [0, 1.0]. Instead got "
                + std::to_string(config.discovery_fraction));
    }

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::invalid_argument(
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

std::thread start_server(HttpServer &server_startup, Config &config, size_t num_threads) {
    server_startup.config.thread_pool_size = num_threads;

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;
    server_startup.config.timeout_request = 30;    // 30 sec to finish headers
    server_startup.config.timeout_content = kContentTimeoutS;   // 15 minutes for body/compute (per request) max

    logger->info("[Server] Will listen on {} port {}",
                 server_startup.config.address, server_startup.config.port);
    logger->info("[Server] Maximum connections: {}", num_threads);
    return std::thread([&server_startup]() { server_startup.start(); });
}

using GraphPair = std::pair<std::string, std::string>;
using GraphIndexes = tsl::hopscotch_map<std::string, std::vector<GraphPair>>;

/**
 * Address exactly one physical (graph, annotation) pair for a traversal request or the probe.
 * A traversal must stay inside one graph, so a name covering several shards is rejected
 * unless the request disambiguates it with |graph_path|; a graph_path naming one graph with
 * several annotations under the name is rejected too — a traversal reads ONE annotation, and
 * answering from the first of them would leave the others unread without saying so.
 */
const GraphPair& select_traverse_pair(const std::string &name, const std::string *graph_path,
                                      const GraphIndexes &indexes) {
    auto it = indexes.find(name);
    if (it == indexes.end())
        throw InvalidRequest("Bad request: unknown graph '" + name + "'");

    // the same pair listed twice under the name is one index (/search keeps its list as is)
    std::vector<const GraphPair*> pairs;
    for (const auto &pair : it->second) {
        if (std::none_of(pairs.begin(), pairs.end(),
                         [&](const GraphPair *p) { return *p == pair; })) {
            pairs.push_back(&pair);
        }
    }
    if (pairs.size() == 1)
        return *pairs[0];

    // a name whose pairs share one graph has several annotations, which graph_path cannot
    // choose between (review of pass 5: it was told to pass graph_path, then refused for it)
    auto several_annotations = [&](const std::string &graph) {
        std::string annotations;
        for (const GraphPair *pair : pairs) {
            if (pair->first == graph)
                annotations += (annotations.empty() ? "" : ", ") + pair->second;
        }
        return InvalidRequest("Bad request: " + (graph_path
                                  ? "graph_path '" + graph + "' names a graph"
                                  : "index '" + name + "' lists one graph (" + graph + ")")
                              + " with several annotations" + (graph_path ? " under index '"
                                  + name + "'" : std::string()) + " (" + annotations
                              + "): a traversal reads one annotation, which graph_path cannot "
                              "choose; list each annotation under a name of its own");
    };
    if (graph_path) {
        const GraphPair *found = nullptr;
        size_t matching = 0;
        for (const GraphPair *pair : pairs) {
            if (pair->first != *graph_path)
                continue;
            matching++;
            if (!found)
                found = pair;
        }
        if (!found) {
            throw InvalidRequest("Bad request: graph_path '" + *graph_path + "' is not part of "
                                 "index '" + name + "'");
        }
        if (matching > 1)
            throw several_annotations(*graph_path);
        return *found;
    }
    std::vector<std::string> graphs;
    for (const GraphPair *pair : pairs) {
        if (std::find(graphs.begin(), graphs.end(), pair->first) == graphs.end())
            graphs.push_back(pair->first);
    }
    if (graphs.size() == 1)
        throw several_annotations(graphs[0]);
    std::string candidates;
    for (const std::string &graph : graphs) {
        candidates += (candidates.empty() ? "" : ", ") + graph;
    }
    throw InvalidRequest("Bad request: index '" + name + "' spans several graphs (" + candidates
                         + "); pass 'graph_path' to pick one");
}

/**
 * The (graph, annotation) pair of a POST /resolve or /traverse request: the single index, or in
 * multi-graph mode the one its "graph" (and "graph_path") fields select.
 */
const graph::AnnotatedDBG&
resolve_traverse_index(const Json::Value &json, const Config &config,
                       const std::shared_future<std::unique_ptr<graph::AnnotatedDBG>> &anno_graph,
                       const GraphIndexes &indexes,
                       const VectorMap<GraphPair, std::unique_ptr<graph::AnnotatedDBG>> &graphs_cache) {
    if (config.fnames.empty()) {
        if (json.isMember("graph") || json.isMember("graph_path")) {
            throw InvalidRequest("Bad request: this server hosts a single graph; "
                                 "remove the 'graph' / 'graph_path' field");
        }
        return *anno_graph.get();
    }
    if (!json.isMember("graph") || !json["graph"].isString())
        throw InvalidRequest("Bad request: 'graph' (index name) is required in multi-graph mode");
    std::string wanted;
    const bool has_path = json.isMember("graph_path") && json["graph_path"].isString();
    if (has_path)
        wanted = json["graph_path"].asString();
    return *graphs_cache.at(select_traverse_pair(json["graph"].asString(),
                                                 has_path ? &wanted : nullptr, indexes));
}

std::vector<std::string> filter_graphs_from_list(
        const GraphIndexes &indexes,
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


// The reverse index of the sequence headers (CoordToHeader::find_header), built while the
// index loads rather than by the first request naming a header: on refseq33m (33M headers)
// that request spent seconds in it, inside its time budget and outside every read timer (R10:
// 9.3 s of a 1 s budget)
static void build_header_index(const graph::AnnotatedDBG &index) {
    const auto *coord_to_header = index.get_coord_to_header();
    if (!coord_to_header)
        return;
    Timer timer;
    const size_t headers = coord_to_header->build_header_index();
    logger->info("[Server] Sequence header index built: {} headers in {:.1f} s", headers,
                 timer.elapsed());
}

int run_server(Config *config) {
    assert(config);
    std::atomic<size_t> num_requests = 0;

    ThreadPool graph_loader(1, 1);
    std::shared_future<std::unique_ptr<AnnotatedDBG>> anno_graph;

    GraphIndexes indexes;

    ThreadPool graphs_pool(get_num_threads(), 1000 /* max_num_tasks */);
    size_t num_server_threads = std::max(1u, get_num_threads());
    set_num_threads(std::max(1u, config->parallel_each));
    config->parallel_each = 1;  // query one batch at a time
    logger->info("[Server] Threads per graph: {}", get_num_threads());

    VectorMap<std::pair<std::string, std::string>, std::unique_ptr<AnnotatedDBG>> graphs_cache;
    bool loaded_with_mmap = utils::with_mmap();

    // The index identity stated by /traverse and /resolve (DESIGN-traverse-graphlet.md
    // §3.1). Single-index mode: --index-name and --index-manifest, checked against the
    // loaded files before the index is served; the meta fingerprint once the index is
    // loaded. Multi-index mode: per (graph, annotation) pair, the name and the manifest of the
    // graph list's optional columns, checked before loading, and the meta fingerprint computed
    // on first use (a pair without them states null, as before).
    IndexIdentity single_identity;
    std::mutex identities_mutex;
    std::map<const AnnotatedDBG*, IndexIdentity> identities;
    auto identity_of = [&](const AnnotatedDBG &index) -> IndexIdentity {
        if (config->fnames.empty())
            return single_identity;
        std::lock_guard<std::mutex> lock(identities_mutex);
        IndexIdentity &id = identities[&index];
        if (id.meta_fp.empty())
            id.meta_fp = index_meta_fingerprint(graph::traversal::LabelOracle(index));
        return id;
    };

    if (config->infbase_annotators.size() == 1) {
        assert(config->fnames.empty());
        // a manifest of another bundle is refused before hours of loading, not after
        try {
            single_identity.name = config->index_name;
            if (!config->index_manifest.empty()) {
                single_identity.fp = index_manifest_fingerprint(
                        config->index_manifest,
                        index_bundle_files(config->infbase, config->infbase_annotators[0],
                                           !config->no_coord_mapping),
                        nullptr,
                        index_unloaded_optional_files(config->infbase,
                                                      config->infbase_annotators[0],
                                                      !config->no_coord_mapping));
                logger->info("[Server] Index manifest {}: index_fp {}", config->index_manifest,
                             single_identity.fp);
            }
        } catch (const std::exception &e) {
            logger->error("[Server] {}", e.what());
            std::exit(1);
        }
        anno_graph = graph_loader.enqueue([&]() {
            logger->info("[Server] Loading graph and annotation in parallel...");
            auto anno_graph = initialize_annotated_dbg(*config);
            logger->info("[Server] Annotated graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);
            // set before the future is ready, so every request sees it complete
            single_identity.meta_fp = index_meta_fingerprint(graph::traversal::LabelOracle(*anno_graph));
            build_header_index(*anno_graph);
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
        std::vector<GraphListEntry> entries;
        for (size_t line_no = 1; std::getline(file, line); ++line_no) {
            if (line.empty())
                continue; // skip empty lines
            try {
                entries.push_back(parse_graph_list_line(line, line_no));
            } catch (const std::exception &e) {
                logger->error("[Server] Invalid line in the csv file: {}", e.what());
                std::exit(1);
            }
            const GraphListEntry &e = entries.back();
            indexes[e.name].emplace_back(e.graph_path, e.annotation_path);
            num_indexes++;
        }
        // The per-graph identity (DESIGN-traverse-graphlet.md §21): every listed manifest is
        // checked against the files its pair loads — sizes and the digest of its list, no
        // re-hashing, as --index-manifest — before hours of loading, and a mismatch, or two
        // lines stating different identities for one pair, refuses to start
        std::map<GraphPair, std::pair<std::string, std::string>> pair_identity;
        try {
            // what each pair loads (the loader dependency inventory of its listed spellings)
            auto inventory = [&](const GraphListEntry &e) {
                return index_bundle_files(e.graph_path, e.annotation_path,
                                          !config->no_coord_mapping);
            };
            pair_identity = graph_list_identities(entries, [&](const GraphListEntry &e) {
                std::string stated;
                std::string fp;
                try {
                    fp = index_manifest_fingerprint(
                            e.manifest_path, inventory(e), &stated,
                            index_unloaded_optional_files(e.graph_path, e.annotation_path,
                                                          !config->no_coord_mapping));
                } catch (const std::exception &ex) {
                    throw std::invalid_argument("line " + std::to_string(e.line) + " of the "
                                                "graph list: " + ex.what());
                }
                if (!stated.empty() && stated != e.index_ns) {
                    logger->warn("[Server] Index manifest {} names the index '{}'; line {} of "
                                 "the graph list names it '{}', which its responses state",
                                 e.manifest_path, stated, e.line, e.index_ns);
                }
                return fp;
            }, inventory);
        } catch (const std::exception &e) {
            logger->error("[Server] {}", e.what());
            std::exit(1);
        }
        for (const auto &[pair, id] : pair_identity) {
            if (!id.first.empty() || !id.second.empty()) {
                logger->info("[Server] Index ({}, {}): index_ns {}, index_fp {}", pair.first,
                             pair.second, id.first.empty() ? "null" : id.first,
                             id.second.empty() ? "null (no manifest)" : id.second);
            }
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
        for (const auto &[pair, index] : graphs_cache) {
            build_header_index(*index);
        }
        // every response computed on a pair states its identity
        for (const auto &[pair, index] : graphs_cache) {
            const auto &[name, fp] = pair_identity.at(pair);
            IndexIdentity &id = identities[index.get()];
            id.name = name;
            id.fp = fp;
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

    // the actual server
    HttpServer server;
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string& content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();  // the index is not loaded yet, so we can't process the request

            Json::Value content_json = parse_json_string(content);
            logger->info("[Server] Request {}: {}", request_id, content_json.toStyledString());
            Json::Value result;

            // simple case with a single graph pair
            if (!config->fnames.size()) {
                if (content_json.isMember("graphs"))
                    throw std::invalid_argument("Bad request: no support for filtering graphs on this server");
                logger->trace("Request {}: Started querying graph {}, in total graphs being queried at the moment: {}",
                              request_id, config->infbase, graphs_being_queried.fetch_add(1) + 1);
                try {
                    result = process_search_request(content_json, *anno_graph.get(), *config);
                } catch (...) {
                    graphs_being_queried--;
                    throw;
                }
                graphs_being_queried--;
            } else {
                std::vector<std::string> graphs_to_query
                        = filter_graphs_from_list(indexes, content_json, request_id);
                std::mutex mu;
                std::vector<std::shared_future<std::exception_ptr>> futures;
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
                                bool in_ram = content_json.isMember("in_ram") && content_json["in_ram"].asBool();
                                if (in_ram && !loaded_with_mmap)
                                    in_ram = false;  // already in RAM, no need to re-load
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

                                auto json = process_search_request(content_json, *index, *config);
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
            return result;
        });
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            if (!config->fnames.size())
                return process_align_request(content, anno_graph.get()->get_graph(), *config);

            throw std::invalid_argument("Bad request: alignment requests are not yet supported for "
                                        "servers with multiple graphs");
        });
    };

    // The ledger-managed /traverse attempts of this process (requests with attempt_id): the
    // backend half of stage 4 of DESIGN-traverse-graphlet.md §14 (traverse_attempts.hpp)
    AttemptSettings attempt_settings;
    attempt_settings.allowance_ms = config->traverse_attempt_allowance_ms;
    attempt_settings.hard_cap_ms = kContentTimeoutS * 1000.0 - 1000;
    attempt_settings.retention_s = config->traverse_attempt_retention_s;
    attempt_settings.retention_count = config->traverse_attempt_retention;
    attempt_settings.tombstone_max_s = config->traverse_attempt_tombstone_max_s;
    attempt_settings.clock_skew_ms = config->traverse_clock_skew_ms;
    attempt_settings.content_timeout_s = kContentTimeoutS;
    attempt_settings.delivery_compress_mbps = config->traverse_delivery_compress_mbps;
    attempt_settings.delivery_build_mbps = config->traverse_delivery_build_mbps;
    // until measured longer, the walk is taken to end at most a chunk of a read and 950 ms (the
    // heads between two readings of the clock, the stopped seed's finalisation) after its
    // walk-until (calibrated in the efficiency pass: 352-1,001 ms measured on SRA, 1,699 once
    // under load; 200 assumed before)
    attempt_settings.delivery_stop_ms = static_cast<double>(config->traverse_chunk_target_ms) + 950;
    AttemptRegistry attempts(attempt_settings);
    logger->info("[Server] Traverse attempts: server_instance {}, allowance {} ms, {}",
                 attempts.server_instance(), attempt_settings.allowance_ms,
                 attempts.retention_text());

    // The traversal routes' transport: compact JSON, compressed (when the client accepts it)
    // at a faster zlib level than the other routes' 9 — the time to build a response is what
    // an attempt's delivery window bounds, and level 1 writes 3-4 times faster than 9 at
    // about 1.8 times the bytes (measured on real responses); the decompressed bytes are the
    // same at every level
    ResponseControl traversal_io;
    traversal_io.compression_level = config->traverse_compression_level;

    // Report where a query is supported and which labels carry which blocks, and
    // optionally freeze seeds for /traverse. No graph traversal.
    server.resource["^/resolve$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();

            Json::Value json = parse_json_string(content);
            const auto &index = resolve_traverse_index(json, *config, anno_graph,
                                                       indexes, graphs_cache);
            const IndexIdentity identity = identity_of(index);
            // a client that is gone is not answered: abandoned between the request's phases
            try {
                return process_resolve_request(json, index, config->index_release,
                                               config->resolve_max_query_bp, &identity,
                                               [&request]() { return client_gone(*request); });
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            }
        }, /* compact */ true, &traversal_io);
    };

    // Extend frozen seeds along consistent annotation labels.
    server.resource["^/traverse$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                 shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        // Every traversal is stopped when its client is gone: its walk polls the client's
        // connection, and nothing is written. A request with attempt_id is also an attempt of
        // the service's ledger: registered, cancellable by id, bounded in duration by the
        // server itself, and every response to it states its usage (traverse_attempts.hpp)
        // (the client check holds the request weakly: a finished attempt is retained for its
        // state, not with the request's body)
        auto attempt = std::make_shared<Attempt>(
                request_id, request->header_read_time, attempts.settings(),
                attempts.server_instance(), /* enforced */ true,
                [weak = std::weak_ptr<HttpServer::Request>(request)]() {
                    auto r = weak.lock();
                    return !r || client_gone(*r);
                });
        // what this server measured of its deliveries replaces the configured rates and ratios
        // in the attempt's reserve (the slowest, the smallest, of its recent responses)
        attempt->set_measured(attempts.measured());
        bool registered = false;
        // what an error states besides its message once the attempt is registered: its usage
        auto with_usage = [&](int status, const std::string &what, const std::string &reason) {
            Json::Value body;
            body["error"] = what;
            body["usage"] = attempt->usage_json(reason);
            return HttpError(status, std::move(body));
        };
        ResponseControl control;
        control.compression_level = config->traverse_compression_level;
        control.on_compressed = [&](size_t text_bytes, double seconds) {
            if (text_bytes >= kMeasuredTextBytes && seconds > 0)
                attempts.note_compress_rate(static_cast<double>(text_bytes) / seconds / 1e6);
        };
        // each seed's result is written as text once built (its tree freed at once, its bytes
        // known to the attempt's delivery reserve); the response is assembled from them, byte
        // for byte the text of the whole tree
        ResultTexts texts;
        control.write = [&texts, attempt](const Json::Value &envelope,
                                          const std::function<void()> &check) {
            // the longest stretch between two checks (deadline_check, finding 6)
            double gap = 0;
            std::string text = texts.active
                ? assemble_traverse_response(envelope, texts.texts, check, &gap)
                : json_text(envelope, true, check, &gap);
            attempt->note_delivery_gap_ms(gap);
            return text;
        };
        // the attempt's client and bound while the response is written and compressed
        control.check = [&]() {
            try {
                attempt->check_delivery();
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            }
        };
        control.on_written = [&](int status, std::optional<size_t> bytes) {
            // the longest single annotation read of any /traverse (deadline_check), and the
            // slowest build rate measured on its large seeds (the delivery reserve)
            attempts.note_uninterruptible(attempt->max_read_ms());
            attempts.note_uninterruptible(attempt->max_delivery_gap_ms());
            attempts.note_build_rate(attempt->own_build_mbps());
            attempts.note_account_per_text_byte(attempt->delivery_detail(),
                                                attempt->own_account_per_text_byte());
            attempts.note_stop_latency(attempt->own_stop_ms());
            if (!registered) {
                if (!status) {
                    logger->info("[Server] Request {}: client gone, {}; no response written",
                                 request_id, attempt->stop_summary());
                }
                return;
            }
            const std::string reason = !status ? "client_gone"
                                     : status == 503 ? "deadline"
                                     : status >= 400 ? "error"
                                     : attempt->walk_reason();
            attempts.finish(attempt, reason, status, bytes);
            logger->info("[Server] Attempt {} (request {}) finished ({}): {}; {}",
                         attempt->ids().attempt_id, request_id, reason, attempt->stop_summary(),
                         status ? fmt::format("response {}, {} bytes, {:.0f} ms", status,
                                              bytes.value_or(0), attempt->elapsed_ms())
                                : std::string("no response written (the client is gone)"));
        };
        process_request(response, request, request_id, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();

            Json::Value json = parse_json_string(content);
            // a malformed id or not_after_ms is refused before anything is registered (400, no
            // usage)
            attempt->set_ids(attempt_ids(json));
            if (attempt->managed()) {
                // an attempt runs once per server process: a second request with a running or
                // retained id is refused, without usage (it would be reconciled against the
                // other attempt); and one whose not_after_ms has passed is not started at all
                // (its ledger may already have released it), also without usage
                if (auto refused = attempts.start(attempt)) {
                    if (refused->instance_mismatch) {
                        logger->info("[Server] Attempt {} (request {}): not started, {}",
                                     attempt->ids().attempt_id, request_id,
                                     refused->body["error"].asString());
                        throw HttpError(409, std::move(refused->body));
                    }
                    if (refused->expired) {
                        logger->info("[Server] Attempt {} (request {}): not started, {}",
                                     attempt->ids().attempt_id, request_id,
                                     refused->body["error"].asString());
                        throw HttpError(409, std::move(refused->body));
                    }
                    Json::Value body;
                    if (refused->body["tombstone"].asBool()) {
                        // its attempt object states the suppression, judged against this
                        // request's not_after_ms (the tombstone now covers it when it can)
                        body["error"] = "attempt_id '" + attempt->ids().attempt_id + "' was "
                                        "cancelled before this request arrived (tombstoned: "
                                      + attempts.retention_text() + "): it is not run";
                    } else {
                        body["error"] = "attempt_id '" + attempt->ids().attempt_id + "' is "
                                        "running or was used on this server within the "
                                        "retention period (" + attempts.retention_text()
                                      + "): an attempt runs once";
                    }
                    body["attempt"] = std::move(refused->body);
                    throw HttpError(409, std::move(body));
                }
                registered = true;
                logger->info("[Server] Attempt {} (request {}): registered",
                             attempt->ids().attempt_id, request_id);
            } else if (const uint64_t now = attempts.now_ms();
                       not_after_passed(attempt->ids(), now)) {
                // not_after_ms without attempt_id is honoured the same way: refused at
                // handler start once passed, nothing run
                Json::Value body = expired_json(attempt->ids(), now, attempts.server_instance());
                logger->info("[Server] Request {}: not started, {}", request_id,
                             body["error"].asString());
                throw HttpError(409, std::move(body));
            }
            try {
                const auto &index = resolve_traverse_index(json, *config, anno_graph,
                                                           indexes, graphs_cache);
                TraverseLimits limits;
                limits.max_time_ms = config->traverse_max_time_ms;
                limits.max_seeds = config->traverse_max_seeds;
                limits.max_seed_bp = config->traverse_max_seed_bp;
                limits.max_seed_labels = config->traverse_max_seed_labels;
                limits.max_memory_mb = config->traverse_max_memory_mb;
                limits.max_work_units = config->traverse_max_work_units;
                limits.chunk_target_ms = static_cast<double>(config->traverse_chunk_target_ms);
                limits.path_cache_bytes = config->traverse_path_cache_mb << 20;
                const IndexIdentity identity = identity_of(index);
                return process_traverse_request(json, index, config->index_release, limits,
                                                &identity, attempt.get(), &texts);
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            } catch (const std::exception &e) {
                if (!attempt->managed())
                    throw;
                throw with_usage(400, e.what(), "error");
            } catch (...) {
                if (!attempt->managed())
                    throw;
                throw with_usage(500, "Internal server error", "error");
            }
        }, /* compact */ true, &control);
    };

    // Cancel a running attempt by its id (POST {"attempt_id", "wait_ms"?, "not_after_ms"?}):
    // 200 when it was asked to stop (state stopping, or finished within wait_ms), 404 when it
    // has finished or is unknown (an unknown id tombstoned, held to the not_after_ms given +
    // the clock skew allowance), 429 when an unknown id was not tombstoned. Not refused while
    // the index loads: an attempt may be cancelled at any time
    server.resource["^/traverse/cancel$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                        shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string &content) {
            Json::Value json = parse_json_string(content);
            if (!json.isObject())
                throw InvalidRequest("request: expected an object");
            for (const std::string &name : json.getMemberNames()) {
                if (name != "attempt_id" && name != "wait_ms" && name != "not_after_ms")
                    throw InvalidRequest("request: unknown field '" + name + "'");
            }
            if (!json["attempt_id"].isString() || !valid_attempt_id(json["attempt_id"].asString())) {
                throw InvalidRequest("request.attempt_id: expected a string of 1 to 128 "
                                     "characters from [A-Za-z0-9._:-]");
            }
            uint64_t wait_ms = 0;
            if (json.isMember("wait_ms")) {
                const Json::Value &w = json["wait_ms"];
                if (!w.isIntegral() || (w.isInt64() && w.asInt64() < 0) || w.asUInt64() > 10'000)
                    throw InvalidRequest("request.wait_ms: expected an integer in [0, 10000]");
                wait_ms = w.asUInt64();
            }
            // the not_after_ms of the request being cancelled (the same rule as /traverse's): a
            // tombstone of an unknown id is held until it has passed + the skew allowance
            const std::optional<uint64_t> not_after_ms = read_not_after_ms(json, "request");
            const std::string id = json["attempt_id"].asString();
            auto [status, body] = attempts.cancel(id, wait_ms, not_after_ms);
            logger->info("[Server] Attempt {}: cancel requested (request {}): {} {}", id,
                         request_id, status, body.get("state", "").asString());
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, /* compact */ true, &traversal_io);
    };

    // The state of an attempt (running | stopping | finished, with the reason and when it
    // stopped), kept for the retention period after it finished; 404 after that
    server.resource["^/traverse/attempt/([^/]+)$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                               shared_ptr<HttpServer::Request> request) {
        const std::string id = request->path_match[1].str();
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!valid_attempt_id(id)) {
                throw InvalidRequest("attempt_id: expected 1 to 128 characters from "
                                     "[A-Za-z0-9._:-]");
            }
            auto [status, body] = attempts.state(id);
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, /* compact */ true, &traversal_io);
    };

    // The content encodings of the traversal routes (compact JSON, Accept-Encoding honoured:
    // gzip preferred, deflate accepted)
    auto encodings_json = []() {
        Json::Value encodings(Json::arrayValue);
        encodings.append("gzip");
        encodings.append("deflate");
        return encodings;
    };

    // How a deadline reaches the walk (both capabilities routes): the time-sized chunks of
    // the annotation reads it may fall into, what stays uninterruptible, and the longest single
    // piece seen
    auto deadline_check_json = [&]() {
        Json::Value d;
        d["chunk_target_ms"] = static_cast<Json::UInt64>(config->traverse_chunk_target_ms);
        // no bound on one row's decode exists before stage 3c (selected-label decoding): a
        // row-diff row reads its whole dependency path, however wide
        d["max_uninterruptible_ms"] = Json::Value();
        d["observed_max_uninterruptible_ms"]
            = static_cast<Json::UInt64>(attempts.observed_max_uninterruptible_ms());
        const graph::traversal::DecodePacer pacer;
        d["rule"] = fmt::format(
            "under a deadline — the seed's bounds.time_budget_ms (at depth > 0 and in a "
            "derivation) and, with attempt_id, the attempt's walk-until — every /traverse "
            "annotation read (a level's fetch, the lookahead, a seed's validation, a "
            "derivation's window) that the deadline may fall into is decoded in chunks, the "
            "deadline checked before each: a read is one piece when, at the slowest per-row "
            "time the request has seen, it would take less than 1/{} of the time left (rows "
            "more than {} times slower than any seen before can make such a read overrun); "
            "else its first chunk is at most {} rows, each next one at most 4 times the "
            "previous and sized at the rate the previous measured to take chunk_target_ms (or "
            "the time left), taken in the walk's order (the rows of one path share their "
            "row-diff decoding), and the rest is one piece once predicted at that rate to take "
            "less than 1/{} of the time left. A read far from its deadline is thus one piece, "
            "as before, and a cancel or a gone client is seen after it. A seed's validation is "
            "stopped by the attempt only, never by its own time budget. One chunk, at least "
            "one row, is uninterruptible, and no bound on one row exists before stage 3c "
            "(max_uninterruptible_ms: null); observed_max_uninterruptible_ms is the longest "
            "single piece of this process, a whole read far from its deadline included, an "
            "observation, not a bound. A chunked read returns exactly what one read would (the "
            "same rows, caches and counters) unless the deadline stops it, and a stopped read "
            "censors the walk at the read (the unchunked walk ran it to its end, past the "
            "deadline, before its next check). The text of each seed's result and of the "
            "response is written under the attempt's delivery check every 64 KiB, a larger "
            "piece copied in pieces up to the next check; the preparation of one token (one "
            "JSON value, e.g. a graphlet string of many MB, is escaped whole before it is "
            "copied) is not interrupted, and the longest time between two such checks is part "
            "of observed_max_uninterruptible_ms. Not chunked: /resolve, the mapping of a "
            "seed's k-mers, a head's processing (checked every work_check_interval units), a "
            "seed's finalisation and summary, the building of a seed's JSON tree between the "
            "attempt's delivery checks (every 4096 objects), and the transport; "
            "chunk_target_ms 0: one piece per read",
            pacer.far_factor, pacer.far_factor, pacer.first_rows, pacer.rest_factor);
        return d;
    };

    // What one index of this deployment supports, so a client can pick a strategy before
    // asking: GET /traverse/capabilities (per graph in multi-graph mode)
    auto probe_json = [&](const AnnotatedDBG &index, const IndexIdentity &identity) {
        graph::traversal::LabelOracle oracle(index);
        Json::Value caps = capabilities_to_json(oracle, config->index_release, &identity);
        caps["max_time_ms"] = config->traverse_max_time_ms;
        caps["max_seeds"] = static_cast<Json::UInt64>(config->traverse_max_seeds);
        caps["max_seed_bp"] = static_cast<Json::UInt64>(config->traverse_max_seed_bp);
        caps["max_seed_labels"] = static_cast<Json::UInt64>(config->traverse_max_seed_labels);
        caps["max_query_bp"] = static_cast<Json::UInt64>(config->resolve_max_query_bp);
        // the request budgets of DESIGN-traverse-graphlet.md §14 (bounds.max_memory_mb,
        // bounds.max_work_units) and W, the interval in charged work units at which the
        // walker reads the clock at the latest; stated here, not in every response, where
        // the per-request capabilities stay as they were. How far a work stop can exceed
        // its budget is stated with what bounds it, not as a fixed maximum: a fetch
        // call's rows are decoded and charged whole (GPT review of stage 2, finding 2),
        // and each stop states the most its seed charged between two comparisons (the
        // review of the stage-2 fixes, F7: no fixed kind of charge bounds them all)
        Json::Value budgets(Json::arrayValue);
        budgets.append("max_memory_mb");
        budgets.append("max_work_units");
        caps["budgets"] = budgets;
        // the server's maxima of those budgets (feature level 4, R16; 0: off): a larger
        // budget is lowered to it and an omitted one set to it, echoed in strategy.clamped
        caps["max_memory_mb"] = static_cast<Json::UInt64>(config->traverse_max_memory_mb);
        caps["max_work_units"] = static_cast<Json::UInt64>(config->traverse_max_work_units);
        caps["work_check_interval"]
            = static_cast<Json::UInt64>(graph::traversal::kWorkCheckInterval);
        // Work is deterministic LOGICAL work, not measured decode effort (review of stage 3,
        // answer 1): the physical decode counters are in each response's timing
        caps["work_bound"] = "bounds.max_work_units counts deterministic logical work, not "
            "measured decode effort: 4 per successor enumeration; per annotation row a fetch "
            "returns 8 per key and 1 per entry and coordinate, and on a budget-aware "
            "(row-diff) annotation 8 per row-diff dependency row and 1 per entry it stores, "
            "whatever the decode shared or cached; 1 per pair evaluation, refusal-scan "
            "entry, edge-reuse probe and step. The walk compares the budget after every "
            "charge, so a stop exceeds it by at most what was charged since the previous "
            "comparison, one indivisible charge (a fetch call's rows, sized from the budget "
            "left down to one key; a label-state scan; or the roots' rows with the end of the "
            "seed phase), and the stop states the most its seed charged between two "
            "comparisons; the seed phase is compared every work_check_interval units and "
            "fails at a comparison finding it at least that much over budget; the deadline "
            "is read before every head and at least every work_check_interval units; the "
            "physical decode counters are in timing";
        caps["memory_bound"] = "soft";
        // the row-diff path cache of the reads (feature level 4, the efficiency pass)
        Json::Value decode_cache;
        decode_cache["path_cache_mb"] = static_cast<Json::UInt64>(config->traverse_path_cache_mb);
        decode_cache["rule"] = "on a row-diff annotation rows a /traverse request's reads "
            "reconstruct are kept in a cache of at most path_cache_mb MiB, so that a later read's "
            "row-diff path stops at a cached row instead of decoding to its anchor again (in two "
            "generations: the older is dropped when the current one fills half the bound): the "
            "rows asked for, the 8 rows after each on its path, every row whose distance to its "
            "anchor is a multiple of 16 (anchors included) and every row whose copy holds less "
            "than 4096 bytes — not every row of every path, which copied each wide row of a long "
            "path and made first reads slower than without the cache (tuple rows are kept flat: "
            "columns, ends, coordinates). What a "
            "read returns, the work units charged (each row with its whole row-diff path) and "
            "the memory admissions (each row by the demand of its whole path) do not depend on "
            "it; the decode time and the physical counters in timing do. Without a memory budget "
            "the cache is the request's, kept from seed to seed; under bounds.max_memory_mb it is "
            "each seed's, off during the seed phase and an annotate root's read, then within what "
            "the label cache leaves of its allotment (min(budget / 4, 64 MiB), held by the "
            "account), and the lookahead's reads are admitted as without it (each also charges "
            "what the cache spared it). 0: off";
        caps["decode_cache"] = std::move(decode_cache);
        // which walk a response gives (every /traverse response names it too)
        caps["algorithm_version"] = kTraverseAlgorithmVersion;
        // the ledger-managed attempts (requests with attempt_id): how they are named,
        // cancelled, queried, kept and bounded (traverse_attempts.hpp)
        caps["attempts"] = attempts.capabilities_json();
        // transport: the traversal routes write compact JSON and honour Accept-Encoding
        // (gzip preferred, deflate accepted), at this zlib level
        caps["content_encodings"] = encodings_json();
        caps["compression_level"] = config->traverse_compression_level;
        caps["deadline_check"] = deadline_check_json();
        return caps;
    };

    // What this deployment supports, so a client can pick a strategy before asking.
    server.resource["^/traverse/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                            shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();
            const auto query = request->parse_query_string();
            if (config->fnames.empty()) {
                // a single index: the parameters that would select one are refused (they name
                // what this server does not have), any other is ignored, as it always was
                for (const auto &[key, value] : query) {
                    if (key == "graph" || key == "graph_path") {
                        throw InvalidRequest("Bad request: this server hosts a single graph; "
                                             "remove the 'graph' / 'graph_path' parameter");
                    }
                }
                return probe_json(*anno_graph.get(), single_identity);
            }
            // Multi-graph mode: the probe of one (graph, annotation) pair, selected by the
            // request fields' rules (?graph=<name>, and graph_path=<path> when the name spans
            // several graphs); it states which pair it describes
            std::optional<std::string> name, graph_path;
            for (const auto &[key, value] : query) {
                std::optional<std::string> *slot = key == "graph" ? &name
                                                 : key == "graph_path" ? &graph_path : nullptr;
                if (!slot) {
                    throw InvalidRequest("Bad request: unknown parameter '" + key + "' of GET "
                                         "/traverse/capabilities (graph, graph_path)");
                }
                if (*slot)
                    throw InvalidRequest("Bad request: the parameter '" + key + "' is given twice");
                *slot = value;
            }
            if (!name) {
                throw InvalidRequest("Bad request: in multi-graph mode GET /traverse/capabilities "
                                     "needs ?graph=<name> (and graph_path=<path> when the name "
                                     "spans several graphs); GET /capabilities lists the graphs");
            }
            const GraphPair &pair = select_traverse_pair(*name, graph_path ? &*graph_path : nullptr,
                                                         indexes);
            const AnnotatedDBG &index = *graphs_cache.at(pair);
            Json::Value caps = probe_json(index, identity_of(index));
            caps["graph"] = *name;
            caps["graph_path"] = pair.first;
            return caps;
        }, /* compact */ true, &traversal_io);
    };

    // The server-wide capabilities (DESIGN-traverse-graphlet.md §21, the owner's note): the
    // routes and features this server offers, its mode and graphs, the attempts and how
    // deadlines are checked — answered while the single index loads (ready: false), so that a
    // service can learn the server's instance and contract before it routes anything to it
    server.resource["^/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                   shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            const bool multi = !config->fnames.empty();
            Json::Value c;
            c["algorithm_version"] = kTraverseAlgorithmVersion;
            c["attempts"] = attempts.capabilities_json();
            c["compression_level"] = config->traverse_compression_level;
            c["content_encodings"] = encodings_json();
            c["deadline_check"] = deadline_check_json();
            c["feature_level"] = kTraverseFeatureLevel;
            Json::Value features(Json::arrayValue);
            Json::Value routes;
            routes["capabilities"] = "GET /capabilities";
            routes["search"] = "POST /search";
            features.append("search");
            if (!multi) {
                // a multi-graph server answers /align with 400
                routes["align"] = "POST /align";
                features.append("align");
            }
            routes["resolve"] = "POST /resolve";
            routes["traverse"] = "POST /traverse";
            routes["traverse_capabilities"] = multi
                ? "GET /traverse/capabilities?graph={name}[&graph_path={path}]"
                : "GET /traverse/capabilities";
            routes["cancel"] = "POST /traverse/cancel";
            routes["attempt"] = "GET /traverse/attempt/{attempt_id}";
            routes["column_labels"] = "GET /column_labels";
            routes["stats"] = "GET /stats";
            for (const char *f : { "resolve", "traverse", "attempts" }) {
                features.append(f);
            }
            c["features"] = std::move(features);
            if (multi) {
                std::vector<std::string> names;
                for (const auto &[name, _] : indexes) {
                    names.push_back(name);
                }
                std::sort(names.begin(), names.end());
                Json::Value graphs(Json::arrayValue);
                for (const std::string &name : names) {
                    graphs.append(name);
                }
                c["graphs"] = std::move(graphs);
            } else {
                c["graphs"] = Json::Value();
            }
            c["mode"] = multi ? "multi" : "single";
            // a multi-graph server loads every index before it listens
            c["ready"] = multi || anno_graph.wait_for(0s) == std::future_status::ready;
            c["release"] = config->index_release;
            c["routes"] = std::move(routes);
            c["schema_version"] = 1;
            c["server_instance"] = attempts.server_instance();
            return c;
        }, /* compact */ true, &traversal_io);
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
            return root;
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
            return root;
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

    std::thread server_thread = start_server(server, *config, num_server_threads);
    server_thread.join();

    return 0;
}

} // namespace cli
} // namespace mtg
