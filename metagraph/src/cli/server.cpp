#include <json/json.h>
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


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;


std::string process_search_request(const std::string &received_message,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const Config &config_orig) {
    Json::Value json = parse_json_string(received_message);

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

    config.count_labels = true;
    config.num_top_labels = json.get("num_labels", config.num_top_labels).asInt();
    config.fast = json.get("fast", config.fast).asBool();
    config.print_signature = json.get("with_signature", config.print_signature).asBool();
    config.query_coords = json.get("query_coords", config.query_coords).asBool();
    config.count_kmers = json.get("abundance_sum", config.count_kmers).asBool();
    config.query_counts = json.get("query_counts", config.query_counts).asBool();

    // Throw client an error if they try to query coordinates/kmer-counts on unsupported indexes
    if ((config.count_kmers || config.query_counts)
            && !(dynamic_cast<const annot::matrix::IntMatrix *>(
                            &anno_graph.get_annotator().get_matrix()))) {
        throw std::invalid_argument("Annotation does not support k-mer count queries");
    }

    if (config.query_coords) {
        if (!dynamic_cast<const annot::matrix::MultiIntMatrix *>(
                        &anno_graph.get_annotator().get_matrix())) {
            throw std::invalid_argument("Annotation does not support k-mer coordinate queries");
        }

        if (config.fast) {
            config.fast = false;
            logger->warn("Attempted to query k-mer coordinates in batch mode. "
                         "Defaulted to basic query mode.");
        }
    }

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (json.get("align", false).asBool())
        aligner_config.reset(new align::DBGAlignerConfig(initialize_aligner_config(config)));

    // Need mutex while appending to vector
    std::vector<SeqSearchResult> search_results;
    std::mutex result_mutex;

    // writing to temporary file in order to reuse query code. This is not optimal and
    // may turn out to be an issue in production. However, adapting FastaParser to
    // work on strings seems non-trivial. An alternative would be to use
    // read_fasta_from_string for non fast queries.
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

    // Return JSON string
    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, search_response);
}

// TODO: implement alignment_result.to_json as in process_search_request
std::string process_align_request(const std::string &received_message,
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

    align::DBGAligner aligner(graph, initialize_aligner_config(config));

    // TODO: make parallel?
    seq_io::read_fasta_from_string(fasta.asString(),
                                   [&](seq_io::kseq_t *read_stream) {
        Json::Value align_entry;
        align_entry[SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD] = read_stream->name.s;

        // not supporting reverse complement yet
        Json::Value alignments = Json::Value(Json::arrayValue);

        // Only align if we have a non-empty sequence
        // TODO: Investigate why calling aligner->align on empty sequence fails
        std::string_view seq = read_stream->seq.s;
        if (!seq.empty()) {
            const auto paths = aligner.align(seq);

            for (const auto &path : paths.data()) {
                Json::Value a;
                a[SeqSearchResult::SCORE_JSON_FIELD] = path.get_score();
                a[SeqSearchResult::SEQUENCE_JSON_FIELD] = path.get_sequence();
                a[SeqSearchResult::CIGAR_JSON_FIELD] = path.get_cigar().to_string();

                alignments.append(a);
            };
        }

        align_entry[SeqSearchResult::ALIGNMENT_JSON_FIELD] = alignments;

        root.append(align_entry);
    });

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::string process_column_label_request(const graph::AnnotatedDBG &anno_graph) {
    auto labels = anno_graph.get_annotator().get_all_labels();

    Json::Value root = Json::Value(Json::arrayValue);

    for (const std::string &label : labels) {
        Json::Value entry = label;
        root.append(entry);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::string process_stats_request(const graph::AnnotatedDBG &anno_graph,
                                  const std::string &graph_filename,
                                  const std::string &annotation_filename) {
    Json::Value root;

    Json::Value graph_stats;
    graph_stats["filename"] = std::filesystem::path(graph_filename).filename().string();
    graph_stats["k"] = static_cast<uint64_t>(anno_graph.get_graph().get_k());
    graph_stats["nodes"] = anno_graph.get_graph().num_nodes();
    graph_stats["is_canonical_mode"] = (anno_graph.get_graph().get_mode()
                                            == graph::DeBruijnGraph::CANONICAL);
    root["graph"] = graph_stats;

    Json::Value annotation_stats;
    const auto &annotation = anno_graph.get_annotator();
    annotation_stats["filename"] = std::filesystem::path(annotation_filename).filename().string();
    annotation_stats["labels"] = static_cast<uint64_t>(annotation.num_labels());
    annotation_stats["objects"] = static_cast<uint64_t>(annotation.num_objects());
    annotation_stats["relations"] = static_cast<uint64_t>(annotation.num_relations());

    root["annotation"] = annotation_stats;

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::thread start_server(HttpServer &server_startup, Config &config) {
    server_startup.config.thread_pool_size = std::max(1u, get_num_threads());

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;

    logger->info("[Server] Will listen on {} port {}", config.host_address,
                 server_startup.config.port);
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

int run_server(Config *config) {
    assert(config);

    assert(config->infbase_annotators.size() == 1);

    ThreadPool graph_loader(1, 1);

    logger->info("[Server] Loading graph...");

    auto anno_graph = graph_loader.enqueue([&]() {
        auto graph = load_critical_dbg(config->infbase);
        logger->info("[Server] Graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);

        auto anno_graph = initialize_annotated_dbg(graph, *config);
        logger->info("[Server] Annotated graph loaded too. Current mem usage: {} MiB", get_curr_RSS() >> 20);
        return anno_graph;
    });

    // defaults for the server
    config->num_top_labels = 10000;
    config->fast = true;

    // the actual server
    HttpServer server;
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &content) {
                return process_search_request(content, *anno_graph.get(), *config);
            });
        }
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &content) {
                return process_align_request(content, anno_graph.get()->get_graph(), *config);
            });
        }
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &) {
                return process_column_label_request(*anno_graph.get());
            });
        }
    };

    server.resource["^/stats"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                            shared_ptr<HttpServer::Request> request) {
        if (check_data_ready(anno_graph, response)) {
            process_request(response, request, [&](const std::string &) {
                return process_stats_request(*anno_graph.get(), config->infbase,
                                             config->infbase_annotators.front());
            });
        }
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
