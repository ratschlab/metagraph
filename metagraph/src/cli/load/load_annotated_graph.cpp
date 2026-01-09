#include "load_annotated_graph.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/coord_to_header.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/logger.hpp"
#include "common/utils/string_utils.hpp"
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

    auto base_graph = graph;
    if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
        graph = std::make_shared<CanonicalDBG>(graph);
        logger->trace("Primary graph wrapped into canonical");
    }

    auto annotation = config.infbase_annotators.size()
            ? initialize_annotation(config.infbase_annotators.at(0), config, 0, max_chunks_open)
            : initialize_annotation(config.anno_type, config, max_index, max_chunks_open);

    std::unique_ptr<annot::CoordToHeader> coord_to_header;

    if (config.infbase_annotators.size()) {
        bool loaded = false;
        if (auto *cc = dynamic_cast<annot::ColumnCompressed<>*>(annotation.get())) {
            loaded = cc->merge_load(config.infbase_annotators);
        } else {
            if (config.infbase_annotators.size() > 1) {
                logger->warn("Cannot merge annotations of this type. Only the first"
                             " file {} will be loaded.", config.infbase_annotators.at(0));
            }
            loaded = annotation->load(config.infbase_annotators.at(0));
        }
        if (!loaded) {
            logger->error("Cannot load annotations for graph {}, file corrupted",
                          config.infbase);
            exit(1);
        }

        // row_diff annotation is special, as it must know the graph structure
        using namespace annot::matrix;
        BinaryMatrix &matrix = const_cast<BinaryMatrix &>(annotation->get_matrix());
        if (IRowDiff *row_diff = dynamic_cast<IRowDiff*>(&matrix)) {
            row_diff->set_graph(base_graph.get());

            if (auto *row_diff_column = dynamic_cast<RowDiff<ColumnMajor> *>(&matrix)) {
                row_diff_column->load_anchor(config.infbase + kRowDiffAnchorExt);
                row_diff_column->load_fork_succ(config.infbase + kRowDiffForkSuccExt);
            }
        }

        if (dynamic_cast<const MultiIntMatrix *>(&annotation->get_matrix())
                && !config.no_coord_mapping && config.identity != Config::ANNOTATE) {
            // Load CoordToHeader mapping if exists
            auto cth_fname = utils::remove_suffix(config.infbase_annotators.at(0),
                                                  annotation->file_extension())
                                            + annot::CoordToHeader::kExtension;
            if (std::filesystem::exists(cth_fname)) {
                coord_to_header = std::make_unique<annot::CoordToHeader>();
                if (!coord_to_header->load(cth_fname)) {
                    logger->error("Failed loading the CoordToHeader mapping from {}", cth_fname);
                    exit(1);
                }
                logger->trace("CoordToHeader mapping loaded successfully from {}. "
                              "All queries will be performed against individual sequences.",
                              cth_fname);
            } else {
                logger->warn("CoordToHeader mapping file not found at {}. Querying will be done "
                             "against original labels indexed in annotation {}. Results will show "
                             "file-based coordinates (e.g., '<file_37.fa>:0-1086-1090') instead of "
                             "sequence-based coordinates (e.g., '<seq_9>:0-1-5'). For sequence-"
                             "based queries, create a mapping with '--index-header-coords' during "
                             "annotation, or pass '--no-coord-mapping' to suppress this warning.",
                             cth_fname, config.infbase_annotators.at(0));
            }
        }
    }

    // load graph
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation),
                                                     false, std::move(coord_to_header));

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
