#include "assemble.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"

using mg::common::logger;


int assemble(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    std::unique_ptr<AnnotatedDBG> anno_graph;
    if (config->infbase_annotators.size()) {
        anno_graph = initialize_annotated_dbg(graph, *config);

        logger->trace("Masking graph...");

        graph = mask_graph(*anno_graph, config);

        logger->trace("Masked in {} sec", timer.elapsed());
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

        gfa_file << "H\tVN:Z:1.0" << std::endl;
        graph->call_unitigs(
            [&](const auto &unitig, const auto &path) {
                gfa_file << "S\t" << path.back() << "\t" << unitig << std::endl;
                graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                    gfa_file << "L\t" << node << "\t+\t" << path.back() << "\t+\t0M" << std::endl;
                });
            },
            config->min_tip_size
        );
    }

    FastaWriter writer(config->outfbase, config->header,
                       config->enumerate_out_sequences,
                       get_num_threads() > 1);

    if (config->unitigs || config->min_tip_size > 1) {
        graph->call_unitigs([&](const auto &unitig, auto&&) { writer.write(unitig); },
                            config->min_tip_size,
                            config->kmers_in_single_form);
    } else {
        graph->call_sequences([&](const auto &contig, auto&&) { writer.write(contig); },
                              config->kmers_in_single_form);
    }

    logger->trace("Sequences extracted in {} sec", timer.elapsed());

    return 0;
}
