#include "load_annotated_graph.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/logger.hpp"
#include "cli/config/config.hpp"
#include "load_graph.hpp"
#include "load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using mtg::common::logger;


std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config,
                                                       size_t max_chunks_open) {
    uint64_t max_index = graph->max_index();
    const auto *dbg_graph = dynamic_cast<const DBGSuccinct*>(graph.get());

    if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
        graph = std::make_shared<CanonicalDBG>(graph);
        logger->trace("Primary graph wrapped into canonical");
    }

    auto annotation_temp = config.infbase_annotators.size()
            ? initialize_annotation(config.infbase_annotators.at(0), config, 0, max_chunks_open)
            : initialize_annotation(config.anno_type, config, max_index, max_chunks_open);

    if (config.infbase_annotators.size()) {
        bool loaded = false;
        if (auto *cc = dynamic_cast<annot::ColumnCompressed<>*>(annotation_temp.get())) {
            loaded = cc->merge_load(config.infbase_annotators);
        } else {
            if (config.infbase_annotators.size() > 1) {
                logger->warn("Cannot merge annotations of this type. Only the first"
                             " file {} will be loaded.", config.infbase_annotators.at(0));
            }
            loaded = annotation_temp->load(config.infbase_annotators.at(0));
        }
        if (!loaded) {
            logger->error("Cannot load annotations for graph {}, file corrupted",
                          config.infbase);
            exit(1);
        }

        // row_diff annotation is special, as it must know the graph structure
        using namespace annot::matrix;
        BinaryMatrix &matrix = const_cast<BinaryMatrix &>(annotation_temp->get_matrix());
        if (IRowDiff *row_diff = dynamic_cast<IRowDiff*>(&matrix)) {
            if (!dbg_graph) {
                logger->error("Only succinct de Bruijn graph representations"
                              " are supported for row-diff annotations");
                std::exit(1);
            }

            row_diff->set_graph(dbg_graph);

            if (auto *row_diff_column = dynamic_cast<RowDiff<ColumnMajor> *>(&matrix)) {
                row_diff_column->load_anchor(config.infbase + kRowDiffAnchorExt);
                row_diff_column->load_fork_succ(config.infbase + kRowDiffForkSuccExt);
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


} // namespace cli
} // namespace mtg
