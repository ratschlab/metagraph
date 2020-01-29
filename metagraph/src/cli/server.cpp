#include "server.hpp"

#include <json/json.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "common/network/server.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "query.hpp"
#include "align.hpp"

using mg::common::logger;


std::string form_client_reply(const std::string &received_message,
                              const AnnotatedDBG &anno_graph,
                              const Config &config,
                              IDBGAligner *aligner = nullptr) {
    // TODO: incorporate aligner
    // TODO: fast query
    std::ignore = aligner;

    try {
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
                //TODO: send error message back in a json file
                throw std::domain_error("bad json received");
            }
        }

        const auto &fasta = json["FASTA"];
        const auto &seq = json["SEQ"];

        // discovery_fraction a proxy of 1 - %similarity
        auto discovery_fraction = json.get("discovery_fraction",
                                           config.discovery_fraction).asDouble();
        auto count_labels = json.get("count_labels", config.count_labels).asBool();
        auto print_signature = json.get("print_signature", config.print_signature).asBool();
        auto num_top_labels = json.get("num_labels", config.num_top_labels).asInt();

        std::ostringstream oss;

        // query callback shared by FASTA and sequence modes
        auto execute_server_query = [&](const std::string &name,
                                        const std::string &sequence) {
            execute_query(name,
                          sequence,
                          count_labels,
                          print_signature,
                          config.suppress_unlabeled,
                          num_top_labels,
                          discovery_fraction,
                          config.anno_labels_delimiter,
                          anno_graph,
                          oss);
        };

        if (!seq.isNull()) {
            // input is plain sequence
            execute_server_query(seq.asString(), seq.asString());
        } else if (!fasta.isNull()) {
            // input is a FASTA sequence
            read_fasta_from_string(
                fasta.asString(),
                [&](kseq_t *read_stream) {
                    execute_server_query(read_stream->name.s,
                                         read_stream->seq.s);
                }
            );
        } else {
            logger->error("No input sequences received from client");
            // TODO: no input sequences -> form an error message for the client
            throw std::domain_error("No input sequences");
        }

        return oss.str();

    } catch (const Json::LogicError &e) {
        logger->error("Bad json file: {}", e.what());
        //TODO: send errors in a json file
        throw;
    } catch (const std::exception &e) {
        logger->error("Processing request error: {}", e.what());
        //TODO: send errors in a json file
        throw;
    } catch (...) {
        logger->error("Processing request error");
        //TODO: send errors in a json file
        throw;
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
    if (config->align_sequences && !config->fast)
        aligner.reset(build_aligner(*graph, *config).release());

    const size_t num_threads = std::max(1u, get_num_threads());

    logger->info("Initializing a TCP service with {} threads"
                 ", listening on port {}", num_threads, config->port);

    try {
        asio::io_context io_context;

        asio::signal_set signals(io_context, SIGINT, SIGTERM);
        signals.async_wait([&](auto, auto) { io_context.stop(); });

        Server server(io_context, config->port,
            [&](const std::string &received_message) {
                return form_client_reply(
                    received_message,
                    *anno_graph,
                    *config,
                    aligner.get()
                );
            }
        );

        std::vector<std::thread> workers;
        for (size_t i = 0; i < std::max(1u, get_num_threads()); ++i) {
            workers.emplace_back([&io_context]() { io_context.run(); });
        }
        for (auto &thread : workers) {
            thread.join();
        }
    } catch (const std::exception &e) {
        logger->error("Exception: {}", e.what());
    } catch (...) {
        logger->error("Unknown exception");
    }

    return 0;
}
