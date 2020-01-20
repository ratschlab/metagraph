#include <fmt/format.h>
#include <json/json.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/algorithms.hpp"
#include "common/threads/threading.hpp"
#include "graph/alignment/aligner_methods.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "common/taxid_mapper.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_graph.hpp"
#include "cli/load/load_annotation.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/build.hpp"
#include "cli/annotate.hpp"
#include "cli/stats.hpp"
#include "cli/augment.hpp"
#include "cli/clean.hpp"
#include "cli/merge.hpp"
#include "cli/align.hpp"
#include "cli/query.hpp"
#include "cli/assemble.hpp"
#include "cli/server.hpp"
#include "cli/transform_graph.hpp"
#include "cli/transform_annotation.hpp"

using mg::common::logger;


int main(int argc, char *argv[]) {
    auto config = std::make_unique<Config>(argc, argv);

    logger->set_level(utils::get_verbose() ? spdlog::level::trace : spdlog::level::info);
    //logger->set_pattern("%^date %x....%$  %v");
    //spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
    //console_sink->set_color(spdlog::level::trace, "\033[37m");
    spdlog::flush_every(std::chrono::seconds(1));

    logger->trace("Metagraph started");

    switch (config->identity) {
        case Config::EXPERIMENT: {
            break;
        }
        case Config::BUILD: {
            return build_graph(config.get());
        }
        case Config::EXTEND: {
            return augment_graph(config.get());
        }
        case Config::ANNOTATE: {
            return annotate_graph(config.get());
        }
        case Config::ANNOTATE_COORDINATES: {
            return annotate_graph_with_genome_coordinates(config.get());
        }
        case Config::MERGE_ANNOTATIONS: {
            return merge_annotation(config.get());
        }
        case Config::QUERY: {
            return query_graph(config.get());
        }
        case Config::SERVER_QUERY: {
            return run_server(config.get());
        }
        case Config::COMPARE: {
            return compare(config.get());
        }
        case Config::CONCATENATE: {
            return concatenate_graph_chunks(config.get());
        }
        case Config::MERGE: {
            return merge_graph(config.get());
        }
        case Config::CLEAN: {
            return clean_graph(config.get());
        }
        case Config::STATS: {
            return print_stats(config.get());
        }
        case Config::TRANSFORM_ANNOTATION: {
            return transform_annotation(config.get());
        }
        case Config::TRANSFORM: {
            return transform_graph(config.get());
        }
        case Config::ASSEMBLE: {
            return assemble(config.get());
        }
        case Config::RELAX_BRWT: {
            return relax_multi_brwt(config.get());
        }
        case Config::ALIGN: {
            return align_to_graph(config.get());
        }
        case Config::CALL_VARIANTS: {
            assert(config->infbase_annotators.size() == 1);

            std::unique_ptr<TaxIDMapper> taxid_mapper;
            if (config->taxonomy_map.length()) {
                taxid_mapper.reset(new TaxIDMapper());
                std::ifstream taxid_mapper_in(config->taxonomy_map, std::ios::binary);
                if (!taxid_mapper->load(taxid_mapper_in)) {
                    logger->error("Failed to read accession->taxid map");
                    exit(1);
                }
            }

            auto anno_graph = initialize_annotated_dbg(*config);
            auto masked_graph = mask_graph(*anno_graph, config.get());

            logger->trace("Filter out: {}", fmt::join(config->label_filter, " "));

            std::ostream *outstream = config->outfbase.size()
                ? new std::ofstream(config->outfbase)
                : &std::cout;

            std::unique_ptr<Json::StreamWriter> json_writer;
            if (config->output_json) {
                Json::StreamWriterBuilder builder;
                builder["indentation"] = "";
                json_writer.reset(builder.newStreamWriter());
            } else {
                *outstream << "Index"
                           << "\t" << "Ref"
                           << "\t" << "Var"
                           << "\t" << "Label";

                if (taxid_mapper.get())
                    *outstream << "\t" << "TaxID";

                *outstream << std::endl;
            }

            std::sort(config->label_filter.begin(), config->label_filter.end());

            ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1);
            std::mutex print_label_mutex;
            std::atomic_uint64_t num_calls = 0;

            auto mask_in_labels = utils::join_strings(config->label_mask_in, ",");

            auto print_variant =
                [&](auto&& alignment, const std::string &query, auto&& vlabels) {
                    // filter out labels
                    std::sort(vlabels.begin(), vlabels.end());

                    auto it = config->label_filter.begin();
                    for (const auto &label : vlabels) {
                        while (it != config->label_filter.end() && *it < label)
                            ++it;

                        if (it == config->label_filter.end())
                            break;

                        if (*it == label)
                            return;
                    }

                    num_calls++;

                    auto label = utils::join_strings(vlabels, ",");

                    // map labels to Taxonomy IDs
                    if (taxid_mapper.get()) {
                        label += "\t" + std::accumulate(
                            vlabels.begin(),
                            vlabels.end(),
                            std::string(),
                            [&](std::string &taxids, const std::string &label) {
                                return std::move(taxids) + ","
                                    + std::to_string(taxid_mapper->gb_to_taxid(label));
                            }
                        );
                    }

                    // print labels
                    std::lock_guard<std::mutex> lock(print_label_mutex);
                    if (config->output_json) {
                        json_writer->write(
                            alignment.to_json(
                                query,
                                masked_graph->get_graph(),
                                false,
                                mask_in_labels + ":" + std::to_string(num_calls),
                                label
                            ),
                            outstream
                        );

                        *outstream << std::endl;
                    } else {
                        logger->info("{}\t{}\t{}\t{}", alignment.front(), query,
                                     alignment.get_sequence(), label);
                    }
                };

            if (config->call_bubbles) {
                annotated_graph_algorithm::call_bubbles(
                    *masked_graph,
                    *anno_graph,
                    print_variant,
                    &thread_pool
                );
            } else if (config->call_breakpoints) {
                annotated_graph_algorithm::call_breakpoints(
                    *masked_graph,
                    *anno_graph,
                    print_variant,
                    &thread_pool
                );
            } else {
                logger->error("No variant calling mode selected. Exiting");
                exit(1);
            }

            thread_pool.join();

            logger->trace("# nodes checked: {}", masked_graph->num_nodes());
            logger->trace("# called: {}", num_calls);

            return 0;
        }
        case Config::PARSE_TAXONOMY: {
            TaxIDMapper taxid_mapper;
            if (config->accession2taxid.length()
                && !taxid_mapper.parse_accession2taxid(config->accession2taxid)) {
                logger->error("Failed to read accession->taxid file");
                exit(1);
            }

            if (config->taxonomy_nodes.length()
                && !taxid_mapper.parse_nodes(config->taxonomy_nodes)) {
                logger->error("Failed to read nodes.dmp file");
                exit(1);
            }

            std::ofstream out(config->outfbase + ".taxonomy.map", std::ios::binary);
            taxid_mapper.serialize(out);
            return 0;
        }
        case Config::NO_IDENTITY: {
            assert(false);
            break;
        }
    }

    return 0;
}
