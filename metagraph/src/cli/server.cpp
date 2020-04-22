#include <json/json.h>
#include <zlib.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
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

const std::string SEQ_DESCRIPTION_JSON_FIELD = "seq_description";
const std::string SCORE_JSON_FIELD = "score";
const std::string SEQUENCE_JSON_FIELD = "sequence";

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
    vector<string> queries = utils::split_string(ret_str, "\n");

    std::map<int, Json::Value> query_results;

    for (auto qit = queries.begin(); qit != queries.end(); ++qit) {
        vector<string> parts = utils::split_string(*qit, "\t");

        // no sequences found
        if(parts.size() <= 2) {
            continue;
        }

        Json::Value res_obj;
        vector<string> query_desc_parts = utils::split_string(parts[1], ":");
        res_obj[SEQ_DESCRIPTION_JSON_FIELD] = query_desc_parts[0].substr(1);

        if(query_desc_parts.size() > 1) {
            // we aligned first, so extracting aligned sequence and score:

            res_obj[SEQUENCE_JSON_FIELD] = query_desc_parts[1];
            res_obj[SCORE_JSON_FIELD] = (int) atoi(query_desc_parts[2].c_str());
        }

        res_obj["results"] = Json::Value(Json::arrayValue);

        for (auto it = ++(++parts.begin()); it != parts.end(); ++it) {
            Json::Value sampleEntry;

            vector<string> entries = utils::split_string(*it, ":");

            vector<string> labels = utils::split_string(entries[0].substr(1, entries[0].size()-2), ";");

            sampleEntry["sampleName"] = labels[0];

            Json::Value properties = Json::objectValue;

            for (auto lit = ++labels.begin(); lit != labels.end(); ++lit) {
                vector<string> keyValue = utils::split_string(*lit, "=");
                properties[keyValue[0]] = adjust_for_types(keyValue[1]);
            }

            if(properties.size() > 0) {
                sampleEntry["properties"] = properties;
            }
            sampleEntry["sampleCount"] = (int) atoi(entries[1].c_str());

            res_obj["results"].append(sampleEntry);
        }

        query_results[(int) atoi(parts[0].c_str())] = res_obj;
    }

    Json::Value root = Json::Value(Json::arrayValue);
    // output by query id
    for(const auto it : query_results) {
        root.append(it.second);
    }

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
}

Json::Value parse_json_string(const std::string &msg) {
    Json::Value json;

    {
        Json::CharReaderBuilder rbuilder;
        std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
        std::string errors;

        if (!reader->parse(msg.data(),
                           msg.data() + msg.size(),
                           &json,
                           &errors)) {
            logger->error("Bad json file:\n{}", errors);
            throw std::domain_error("Bad json received: " + errors);
        }
    }

    return json;
}


std::string form_client_reply(const std::string &received_message,
                              const AnnotatedDBG &anno_graph,
                              const Config &config_orig,
                              IDBGAligner *aligner,
                              ThreadPool &thread_pool) {
    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];

    Config config(config_orig);

    // discovery_fraction a proxy of 1 - %similarity
    config.discovery_fraction = json.get("discovery_fraction", config.discovery_fraction).asDouble();

    config.count_labels = true;
    config.num_top_labels = json.get("num_labels", config.num_top_labels).asInt();
    config.fast = json.get("fast", config.fast).asBool();

    bool do_alignment = json.get("align", false).asBool();

    std::ostringstream oss;

    if (!fasta.isNull()) {
        // writing to temporary file in order to reuse query code. This is not optimal and
        // may turn out to be an issue in production. However, adapting FastaParser to work
        // on strings seems non-trivial. An alternative would be to use read_fasta_from_string
        // for non fast queries.
        utils::TempFile tf(config.tmp_dir);
        tf.ofstream() << fasta.asString();
        tf.ofstream().close();

        std::mutex oss_mutex;

        QueryExecutor qe(config, anno_graph,
                do_alignment ? aligner : nullptr,
                thread_pool);
        qe.query_fasta(tf.name(), [&](const std::string res){
          std::lock_guard<std::mutex> lock(oss_mutex);
          oss << res;
        });
    }
    else {
        logger->error("No input sequences received from client");
        throw std::domain_error("No input sequences received from client");
    }

    return convert_to_json(oss.str());
}

std::string form_align_reply(const std::string &received_message,
                              const Config &config,
                              IDBGAligner *aligner = nullptr) {
    std::ignore = config;

    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];

    Json::Value root = Json::Value(Json::arrayValue);

    read_fasta_from_string(fasta.asString(), [&](kseq_t *read_stream) {
        auto name = read_stream->name.s;
        auto query = read_stream->seq.s;

        const auto paths = aligner->align(query);

        Json::Value align_entry;
        align_entry[SEQ_DESCRIPTION_JSON_FIELD] = name;

        // not supporting reverse complement yet
        if(!paths.empty()) {
            auto path = paths.front();

            align_entry[SCORE_JSON_FIELD] = path.get_score();
            align_entry[SEQUENCE_JSON_FIELD] = path.get_sequence();
        } else {
            align_entry[SEQUENCE_JSON_FIELD] = "";
        }

        root.append(align_entry);
    });

    Json::StreamWriterBuilder builder;
    return Json::writeString(builder, root);
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

//https://panthema.net/2007/0328-ZLibString.html
/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
std::string compress_string(const std::string& str,
                            int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}

void write_compressed_if_possible(SimpleWeb::StatusCode status, const string& msg,
        shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
    auto encoding_header = request->header.find("Accept-Encoding");

    if (encoding_header != request->header.end() && encoding_header->second.find("deflate") != std::string::npos) {
        auto compressed = compress_string(msg);
        auto header = SimpleWeb::CaseInsensitiveMultimap({{"Content-Encoding", "deflate"},
                                                     {"Content-Length",   std::to_string(compressed.size())}});
        response->write(status, compressed, header);
    } else {
        response->write(status, msg);
    }
}


int run_server(Config *config) {
    assert(config);

    assert(config->infbase_annotators.size() == 1);

    logger->info("Loading graph...");

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    logger->info("Graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);

    std::unique_ptr<IDBGAligner> aligner;
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
            ThreadPool query_thread_pool(num_threads); // TODO: smarter way to determine how many threads to use?
            auto ret = form_client_reply(content, *anno_graph, *config, aligner.get(), query_thread_pool);
            write_compressed_if_possible(SimpleWeb::StatusCode::success_ok, ret, response, request);
        }  catch(const std::exception &e) {
            logger->info("Error on request " + string(e.what()));
            response->write(SimpleWeb::StatusCode::client_error_bad_request, json_str_with_error_msg(e.what()));
        }
        catch(...) {
            logger->info("Error on request ");
            response->write(SimpleWeb::StatusCode::server_error_internal_server_error, json_str_with_error_msg("Internal server error"));
        }
    };

    // TODO: refactor request business in lambda, use switch?
    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
        // Retrieve string:
        auto content = request->content.string();

        logger->info("Got align request from " + request->remote_endpoint().address().to_string());

        try {
            auto ret = form_align_reply(content, *config, aligner.get());
            write_compressed_if_possible(SimpleWeb::StatusCode::success_ok, ret, response, request);
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
