#include <json/json.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "query.hpp"
#include "align.hpp"

#include "asio.hpp"

#include "server_http.hpp"

using mg::common::logger;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

Json::Value adjust_for_types(const string &v) {
    if(v == "nan") {
        return Json::nullValue;
    }

    try {
        return Json::Value(std::stof(v));
    } catch(...) {};

    return Json::Value(v);
}

string convert_to_json(const string &ret_str) {
    string no_endl = utils::remove_suffix(ret_str, "\n");
    vector<string> parts = utils::split_string(no_endl, "\t");

    Json::Value root;
    // only supporting single query at the moment
    root["query"] =  parts[0];

    root["results"] = Json::Value(Json::arrayValue);

    for (auto it = ++parts.begin(); it != parts.end(); ++it) {
        Json::Value sampleEntry;

        vector<string> entries = utils::split_string(*it, ":");

        vector<string> labels = utils::split_string(entries[0].substr(1, entries[0].size()-2), ";");

        sampleEntry["sampleName"] = labels[0];

        Json::Value properties = Json::objectValue;

        for (auto lit = ++labels.begin(); lit != labels.end(); ++lit) {
            vector<string> keyValue = utils::split_string(*lit, "=");
            properties[keyValue[0]] = adjust_for_types(keyValue[1]);
        }

        sampleEntry["properties"] = properties;
        sampleEntry["sampleCount"] = (int) atoi(entries[1].c_str());

        root["results"].append(sampleEntry);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

std::string form_client_reply(const std::string &received_message,
                              const AnnotatedDBG &anno_graph,
                              const Config &config,
                              IDBGAligner *aligner = nullptr) {
    // TODO: query with alignment to graph
    // TODO: fast query
    std::ignore = aligner;

    Json::Value json;

    {
        Json::CharReaderBuilder rbuilder;
        std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
        std::string errors;

        if (!reader->parse(received_message.data(),
                           received_message.data() + received_message.size(),
                           &json,
                           &errors)) {
            logger->error("Bad json file:\n{}", errors);
            throw std::domain_error("Bad json received: " + errors);
        }
    }

    const auto &fasta = json["FASTA"];

    // discovery_fraction a proxy of 1 - %similarity
    auto discovery_fraction = json.get("discovery_fraction",
                                       config.discovery_fraction).asDouble();
    auto count_labels = json.get("count_labels", config.count_labels).asBool();
    auto print_signature = json.get("print_signature", config.print_signature).asBool();
    auto num_top_labels = json.get("num_labels", config.num_top_labels).asInt();

    std::ostringstream oss;

    // query callback shared by FASTA and sequence modes
    auto execute_server_query = [&](const std::string &name, const std::string &sequence) {
        oss << QueryExecutor::execute_query(name, sequence, count_labels, print_signature,
                                            config.suppress_unlabeled, num_top_labels,
                                            discovery_fraction,
                                            config.anno_labels_delimiter, anno_graph);
    };

    if (!fasta.isNull()) {
        // input is a FASTA sequence
        read_fasta_from_string(
            fasta.asString(),
            [&](kseq_t *read_stream) {
                execute_query(read_stream->name.s,
                              read_stream->seq.s,
                              count_labels,
                              print_signature,
                              config.suppress_unlabeled,
                              num_top_labels,
                              discovery_fraction,
                              config.anno_labels_delimiter,
                              anno_graph,
                              oss);
            }
        );
    } else {
        logger->error("No input sequences received from client");
        throw std::domain_error("No input sequences received from client");
    }

    return convert_to_json(oss.str());
}

std::string form_column_labels_reply(const AnnotatedDBG &anno_graph) {
    auto labels = anno_graph.get_annotation().get_all_labels();

    Json::Value root = Json::Value(Json::arrayValue);

    for(const auto &label : labels) {
        Json::Value entry = label;
        root.append(entry);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

string json_str_with_error_msg(const string &msg) {
    Json::Value root;
    root["error"] = msg;
    return Json::writeString(Json::StreamWriterBuilder(), root);
}

int run_server(Config *config) {
    assert(config);

    assert(config->infbase_annotators.size() == 1);

    logger->info("Loading graph...");

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    logger->info("Graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);

    std::unique_ptr<IDBGAligner> aligner;
    if (config->align_sequences && !config->fast)
        aligner.reset(build_aligner(*graph, *config).release());

    const size_t num_threads = std::max(1u, get_num_threads());

    HttpServer server;
    server.config.thread_pool_size = num_threads;
    if(config->host_address != "") {
        logger->info("Will listen on interface {}", config->host_address);
        server.config.address = config->host_address;
    }

    server.config.port = config->port;

    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
        // Retrieve string:
        auto content = request->content.string();

        logger->info("Got search request from " + request->remote_endpoint().address().to_string());

        try {
            auto ret = form_client_reply(content, *anno_graph, *config, aligner.get());
            response->write(SimpleWeb::StatusCode::success_ok, ret);
        }  catch(const std::exception &e) {
            logger->info("Error on request " + string(e.what()));
            response->write(SimpleWeb::StatusCode::client_error_bad_request, json_str_with_error_msg(e.what()));
        }
        catch(...) {
            logger->info("Error on request ");
            response->write(SimpleWeb::StatusCode::server_error_internal_server_error, json_str_with_error_msg("Internal server error"));
        }
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
        logger->info("Got column_labels request from " + request->remote_endpoint().address().to_string());

        try {
            auto ret = form_column_labels_reply(*anno_graph);
            response->write(SimpleWeb::StatusCode::success_ok, ret);
        }  catch(const std::exception &e) {
            logger->info("Error on request " + string(e.what()));
            response->write(SimpleWeb::StatusCode::client_error_bad_request, json_str_with_error_msg(e.what()));
        }
        catch(...) {
            logger->info("Error on request ");
            response->write(SimpleWeb::StatusCode::server_error_internal_server_error, json_str_with_error_msg("Internal server error"));
        }
    };

    server.default_resource["GET"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
        logger->info("Not found " + request->path);
        response->write(SimpleWeb::StatusCode::client_error_not_found, "Could not find path " + request->path);
    };
    server.default_resource["POST"] = server.default_resource["GET"];

    server.on_error = [](shared_ptr<HttpServer::Request> /*request*/, const SimpleWeb::error_code &ec) {
        // Handle errors here
        // Note that connection timeouts will also call this handle with ec set to SimpleWeb::errc::operation_canceled
        if(ec.value() != asio::stream_errc::eof) {
            logger->info("Got error " + ec.message() + " " + ec.category().name() + " " + std::to_string(ec.value()));
        }
    };

    promise<unsigned short> server_port;
    thread server_thread([&server, &server_port]() {
        // Start server
        server.start([&server_port](unsigned short port) {
            server_port.set_value(port);
        });
    });

    logger->info("Initializing a HTTP server with {} threads"
                 ", listening on port {}", num_threads, server_port.get_future().get());

    server_thread.join();

    return 0;
}
