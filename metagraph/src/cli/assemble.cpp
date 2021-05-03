#include "assemble.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/masked_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::graph::DeBruijnGraph;
using mtg::graph::MaskedDeBruijnGraph;
using mtg::graph::AnnotatedDBG;


std::unique_ptr<MaskedDeBruijnGraph>
mask_graph(const AnnotatedDBG &anno_graph, Config *config) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    logger->trace("Masked in: {}", fmt::join(config->label_mask_in, " "));
    logger->trace("Masked out: {}", fmt::join(config->label_mask_out, " "));

    if (!config->filter_by_kmer) {
        return std::make_unique<MaskedDeBruijnGraph>(
            graph,
            mask_nodes_by_unitig_labels(
                anno_graph,
                config->label_mask_in,
                config->label_mask_out,
                std::max(1u, get_num_threads()),
                config->label_mask_in_fraction,
                config->label_mask_out_fraction,
                config->label_other_fraction
            )
        );
    }

    return std::make_unique<MaskedDeBruijnGraph>(
        graph,
        mask_nodes_by_node_label(
            anno_graph,
            config->label_mask_in,
            config->label_mask_out,
            [config,&anno_graph](auto index,
                                 auto get_num_in_labels,
                                 auto get_num_out_labels) {
                assert(index != DeBruijnGraph::npos);

                size_t num_in_labels = get_num_in_labels();

                if (num_in_labels < config->label_mask_in_fraction
                                        * config->label_mask_in.size())
                    return false;

                size_t num_out_labels = get_num_out_labels();

                if (num_out_labels < config->label_mask_out_fraction
                                        * config->label_mask_out.size())
                    return false;

                size_t num_total_labels = anno_graph.get_labels(index).size();

                return num_total_labels - num_in_labels - num_out_labels
                            <= config->label_other_fraction * num_total_labels;
            },
            std::max(1u, get_num_threads())
        )
    );
}


int assemble(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    std::unique_ptr<graph::AnnotatedDBG> anno_graph;
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
