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
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "cli/config/config.hpp"
#include "load_graph.hpp"
#include "load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using mtg::common::logger;

namespace {

// Kick off graph loading on a worker thread.
std::shared_future<std::shared_ptr<DeBruijnGraph>>
async_load_critical_dbg(const Config &config) {
    return std::async(std::launch::async, [&config]() -> std::shared_ptr<DeBruijnGraph> {
        return load_critical_dbg(config.infbase);
    }).share();
}

// Load (or initialize, if no annotators are configured) the annotation matrix.
// When annotators are configured, metadata init runs while the graph is still
// loading on another thread; the graph is only awaited at the point where it
// is needed (to set the base graph on a row_diff annotator).
std::unique_ptr<annot::MultiLabelAnnotation<std::string>>
load_annotation_async(std::shared_future<std::shared_ptr<DeBruijnGraph>> graph_future,
                      const Config &config,
                      size_t max_chunks_open) {
    if (!config.infbase_annotators.size()) {
        // No annotators configured: build an empty annotation sized to the graph.
        auto graph = graph_future.get();
        return initialize_annotation(config.anno_type, config,
                                     graph->max_index(), max_chunks_open);
    }

    // Annotation metadata init runs in parallel with graph load.
    auto annotation = initialize_annotation(config.infbase_annotators.at(0),
                                            config, 0, max_chunks_open);
    auto graph = graph_future.get();

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
    if (!loaded)
        exit(1);

    // row_diff annotation needs the underlying base graph (unwrapped).
    // If a CanonicalDBG was passed in, walk one level down.
    using namespace annot::matrix;
    BinaryMatrix &matrix = const_cast<BinaryMatrix &>(annotation->get_matrix());
    if (IRowDiff *row_diff = dynamic_cast<IRowDiff*>(&matrix)) {
        const DeBruijnGraph *base_graph = graph.get();
        if (auto *canonical = dynamic_cast<const CanonicalDBG *>(graph.get()))
            base_graph = &canonical->get_graph();
        row_diff->set_graph(base_graph);

        if (auto *row_diff_column = dynamic_cast<RowDiff<ColumnMajor> *>(&matrix)) {
            row_diff_column->load_anchor(config.infbase + kRowDiffAnchorExt);
            row_diff_column->load_fork_succ(config.infbase + kRowDiffForkSuccExt);
        }
    }

    return annotation;
}

std::unique_ptr<annot::CoordToHeader>
load_coord_to_header(const annot::MultiLabelAnnotation<std::string> &annotation,
                     const Config &config) {
    std::unique_ptr<annot::CoordToHeader> coord_to_header;

    if (dynamic_cast<const annot::matrix::MultiIntMatrix *>(&annotation.get_matrix())
            && !config.no_coord_mapping && config.identity != Config::ANNOTATE) {
        // Load CoordToHeader mapping if exists
        auto cth_fname = utils::remove_suffix(config.infbase_annotators.at(0),
                                              annotation.file_extension())
                                        + annot::CoordToHeader::kExtension;
        if (std::filesystem::exists(cth_fname)) {
            coord_to_header = std::make_unique<annot::CoordToHeader>();
            if (!coord_to_header->load(cth_fname)) {
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

    return coord_to_header;
}

// Build the AnnotatedDBG from a graph future: wait for the annotation and
// graph, wrap PRIMARY into CanonicalDBG, attach the CoordToHeader index,
// and verify compatibility.
std::unique_ptr<AnnotatedDBG>
build_annotated_dbg(std::shared_future<std::shared_ptr<DeBruijnGraph>> graph_future,
                    const Config &config,
                    size_t max_chunks_open) {
    auto annotation = load_annotation_async(graph_future, config, max_chunks_open);
    auto graph = graph_future.get();
    if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
        graph = std::make_shared<CanonicalDBG>(graph);
        logger->trace("Primary graph wrapped into canonical");
    }
    auto coord_to_header = load_coord_to_header(*annotation, config);
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation),
                                                     false, std::move(coord_to_header));
    if (!anno_graph->check_compatibility()) {
        logger->error("Graph and annotation are not compatible");
        exit(1);
    }
    return anno_graph;
}

} // namespace

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config,
                                                       size_t max_chunks_open) {
    std::promise<std::shared_ptr<DeBruijnGraph>> ready;
    ready.set_value(std::move(graph));
    return build_annotated_dbg(ready.get_future().share(), config, max_chunks_open);
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    return build_annotated_dbg(async_load_critical_dbg(config), config, kDefaultMaxChunksOpen);
}

std::pair<std::shared_ptr<DeBruijnGraph>,
          std::future<std::unique_ptr<AnnotatedDBG>>>
load_graph_with_async_annotation(const Config &config) {
    auto graph_future = async_load_critical_dbg(config);
    std::future<std::unique_ptr<AnnotatedDBG>> anno_dbg_future;
    if (config.infbase_annotators.size()) {
        anno_dbg_future = std::async(std::launch::async, [graph_future, &config] {
            return build_annotated_dbg(graph_future, config, kDefaultMaxChunksOpen);
        });
    }
    return std::make_pair(graph_future.get(), std::move(anno_dbg_future));
}


} // namespace cli
} // namespace mtg
