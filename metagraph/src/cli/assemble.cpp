#include "assemble.hpp"

#include <fmt/format.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

int assemble(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    if (config->infbase_annotators.size()) {
        assert(config->label_mask_file.size());
        auto anno_graph = initialize_annotated_dbg(graph, *config);

        logger->trace("Generating masked graphs...");

        std::filesystem::remove(
            utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz"
        );

        std::mutex file_open_mutex;
        std::mutex write_mutex;

        size_t num_threads_per_traversal = get_num_threads() / config->parallel_assemblies;

        call_masked_graphs(*anno_graph, config,
            [&](const graph::MaskedDeBruijnGraph &graph, const std::string &header) {
                std::lock_guard<std::mutex> file_lock(file_open_mutex);
                seq_io::FastaWriter writer(config->outfbase, header,
                                           config->enumerate_out_sequences,
                                           get_num_threads() > 1, /* async write */
                                           true /* append */);

                if (config->unitigs || config->min_tip_size > 1) {
                    graph.call_unitigs([&](const auto &unitig, auto&&) {
                                           std::lock_guard<std::mutex> lock(write_mutex);
                                           writer.write(unitig);
                                       },
                                       num_threads_per_traversal,
                                       config->min_tip_size,
                                       config->kmers_in_single_form);
                } else {
                    graph.call_sequences([&](const auto &seq, auto&&) {
                                             std::lock_guard<std::mutex> lock(write_mutex);
                                             writer.write(seq);
                                         },
                                         num_threads_per_traversal,
                                         config->kmers_in_single_form);
                }
            },
            config->parallel_assemblies,
            num_threads_per_traversal
        );

        return 0;
    }

    logger->trace("Extracting sequences from graph...");

    timer.reset();

    if (config->to_gfa) {
        if (!config->unitigs) {
            logger->error("Flag '--unitigs' must be set for GFA output");
            exit(1);
        }

        logger->trace("Writing graph to GFA...");

        std::ofstream gfa_file(utils::remove_suffix(config->outfbase, ".gfa") + ".gfa");
        std::mutex str_mutex;

        gfa_file << "H\tVN:Z:1.0" << std::endl;
        size_t k = graph->get_k();
        size_t overlap = k - 1;
        graph->call_unitigs(
            [&](const auto &unitig, const auto &path) {
                std::ostringstream ostr;
                if (config->output_compacted) {
                    ostr << fmt::format("S\t{}\t{}\n", path.back(), unitig);
                    graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                        ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                            node, path.back(), overlap);
                    });
                } else {
                    for (size_t i = 0; i < path.size(); ++i) {
                        ostr << fmt::format("S\t{}\t{}\n", path[i],
                                            std::string_view(unitig.data() + i, k));
                        if (i)
                            ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                                path[i - 1], path[i], overlap);
                    }
                    graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                        ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                            node, path.front(), overlap);
                    });
                }
                std::lock_guard<std::mutex> lock(str_mutex);
                gfa_file << ostr.str();
            },
            get_num_threads(),
            config->min_tip_size
        );

        return 0;
    }

    seq_io::FastaWriter writer(config->outfbase, config->header,
                               config->enumerate_out_sequences,
                               get_num_threads() > 1);
    std::mutex write_mutex;

    if (config->unitigs || config->min_tip_size > 1) {
        graph->call_unitigs([&](const auto &unitig, auto&&) {
                                std::lock_guard<std::mutex> lock(write_mutex);
                                writer.write(unitig);
                            },
                            get_num_threads(),
                            config->min_tip_size,
                            config->kmers_in_single_form);
    } else {
        graph->call_sequences([&](const auto &contig, auto&&) {
                                  std::lock_guard<std::mutex> lock(write_mutex);
                                  writer.write(contig);
                              },
                              get_num_threads(),
                              config->kmers_in_single_form);
    }

    logger->trace("Sequences extracted in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
