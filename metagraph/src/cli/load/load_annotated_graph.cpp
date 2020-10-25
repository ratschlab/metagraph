#include "load_annotated_graph.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/logger.hpp"
#include "load_graph.hpp"
#include "load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;

using mtg::common::logger;


std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config) {
    // TODO: check and wrap into canonical only if the graph is primary
    if (config.canonical && !graph->is_canonical_mode())
        graph = std::make_shared<CanonicalDBG>(graph, true);

    uint64_t max_index = graph->max_index();
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(graph.get()))
        max_index = canonical->get_graph().max_index();

    auto annotation_temp = config.infbase_annotators.size()
            ? initialize_annotation(config.infbase_annotators.at(0), config, 0)
            : initialize_annotation(config.anno_type, config, max_index);

    if (config.infbase_annotators.size()) {
        if (!annotation_temp->load(config.infbase_annotators.at(0))) {
            logger->error("Cannot load annotations for graph {}, file corrupted",
                          config.infbase);
            exit(1);
        }
        const Config::AnnotationType input_anno_type
                = parse_annotation_type(config.infbase_annotators.at(0));
        // row_diff annotation is special, as it must know the graph structure
        if (input_anno_type == Config::AnnotationType::RowDiff
            || input_anno_type == Config::AnnotationType::RowDiffBRWT
            || input_anno_type == Config::AnnotationType::RowDiffRowSparse) {
            auto dbg_graph = dynamic_cast<const DBGSuccinct *>(graph.get());

            if (!dbg_graph) {
                if (auto *canonical = dynamic_cast<CanonicalDBG *>(graph.get()))
                    dbg_graph = dynamic_cast<const DBGSuccinct *>(&canonical->get_graph());
            }

            if (!dbg_graph) {
                logger->error(
                        "Only succinct de Bruijn graph representations are "
                        "supported for row-diff annotations");
                std::exit(1);
            }
            // this is really ugly, but the alternative is to add a set_graph method to
            // all annotations
            using namespace annot::binmat;
            BinaryMatrix &matrix = const_cast<BinaryMatrix &>(annotation_temp->get_matrix());
            if (input_anno_type == Config::AnnotationType::RowDiff) {
                dynamic_cast<RowDiff<ColumnMajor> &>(matrix).set_graph(dbg_graph);
                std::string anchor_fname = config.infbase + kRowDiffAnchorExt;
                dynamic_cast<RowDiff<ColumnMajor> &>(matrix).load_anchor(anchor_fname);
            } else if (input_anno_type == Config::AnnotationType::RowDiffRowSparse) {
                dynamic_cast<RowDiff<RowSparse> &>(matrix).set_graph(dbg_graph);
            } else {
                dynamic_cast<RowDiff<BRWT> &>(matrix).set_graph(dbg_graph);
            }
        }
    }

    // load graph
    auto anno_graph
            = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation_temp));

    if (!anno_graph->check_compatibility()) {
        logger->error("Graph and annotation are not compatible");
        exit(1);
    }

    return anno_graph;
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    return initialize_annotated_dbg(load_critical_dbg(config.infbase), config);
}

std::unique_ptr<MaskedDeBruijnGraph>
mask_graph(const AnnotatedDBG &anno_graph, Config *config) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    // Remove non-present labels
    config->label_mask_in.erase(
        std::remove_if(config->label_mask_in.begin(),
                       config->label_mask_in.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists)
                               logger->trace("Removing mask-in label {}", label);

                           return !exists;
                       }),
        config->label_mask_in.end()
    );

    config->label_mask_out.erase(
        std::remove_if(config->label_mask_out.begin(),
                       config->label_mask_out.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists)
                               logger->trace("Removing mask-out label {}", label);

                           return !exists;
                       }),
        config->label_mask_out.end()
    );

    logger->trace("Masked in: {}", fmt::join(config->label_mask_in, " "));
    logger->trace("Masked out: {}", fmt::join(config->label_mask_out, " "));

    if (!config->filter_by_kmer) {
        return std::make_unique<MaskedDeBruijnGraph>(
            graph,
            mask_nodes_by_unitig_labels(
                anno_graph,
                config->label_mask_in,
                config->label_mask_out,
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
            }
        )
    );
}

} // namespace cli
} // namespace mtg
