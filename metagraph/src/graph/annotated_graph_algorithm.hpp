#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "annotated_dbg.hpp"
#include "bitmap.hpp"
#include "dbg_aligner.hpp"
#include "threading.hpp"
#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "masked_graph.hpp"


typedef std::function<uint64_t()> UInt64Callback;

namespace annotated_graph_algorithm {

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    std::function<bool(uint64_t, uint64_t)> keep_node,
                    double lazy_evaluation_density_cutoff = 0.05);

// Allows for lazy evaluation of the counts passed to keep_node
std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    std::function<bool(UInt64Callback, UInt64Callback)> keep_node,
                    double lazy_evaluation_density_cutoff = 0.05);

template <class Index, typename... Args>
using VariantCallback = std::function<void(Alignment<Index>&&,
                                           const std::string&, // query sequence
                                           Args&&...)>;

typedef Alignment<DeBruijnGraph::node_index> DBGAlignment;
typedef VariantCallback<DeBruijnGraph::node_index,
                        std::vector<AnnotatedDBG::Annotator::Label>&&> VariantLabelCallback;


// These functions will callback nothing if graph is equal to the graph stored
// in anno_graph

void call_bubbles(const DeBruijnGraph &graph,
                  const AnnotatedDBG &anno_graph,
                  const VariantLabelCallback &callback,
                  ThreadPool *thread_pool = nullptr,
                  const std::function<bool()> &terminate = []() { return false; });

void call_breakpoints(const DeBruijnGraph &graph,
                      const AnnotatedDBG &anno_graph,
                      const VariantLabelCallback &callback,
                      ThreadPool *thread_pool = nullptr,
                      const std::function<bool()> &terminate = []() { return false; });

// Given a seed, construct a masked graph that only includes the labels of the seed node,
// then extend on that.
Extender<DeBruijnGraph::node_index>
build_masked_graph_extender(const AnnotatedDBG &anno_graph,
                            double seed_label_discovery_fraction = 1.0,
                            Extender<DeBruijnGraph::node_index>&& extender
                                = default_extender<DeBruijnGraph::node_index>);

// Seed until the first breakpoint, relative to the given foreground graph
MapExtendSeederBuilder<DeBruijnGraph::node_index>
build_breakpoint_seeder_builder(const DeBruijnGraph &foreground);

// Given a seed, only extend on the background graph, relative to the
// given foreground graph
Extender<DeBruijnGraph::node_index>
build_background_graph_extender(const DeBruijnGraph &foreground,
                                Extender<DeBruijnGraph::node_index>&& extender
                                    = default_extender<DeBruijnGraph::node_index>);

template <class AlignmentCompare = std::less<Alignment<DeBruijnGraph::node_index>>,
          class ColumnCompare = std::less<typename Alignment<DeBruijnGraph::node_index>::Column>>
void call_variants(const DeBruijnGraph &foreground,
                   const AnnotatedDBG &anno_graph,
                   const VariantLabelCallback &callback,
                   const DBGAlignerConfig &config,
                   const Extender<DeBruijnGraph::node_index> &extender
                       = default_extender<DeBruijnGraph::node_index>,
                   ThreadPool *thread_pool = nullptr,
                   const std::function<bool()> &terminate = []() { return false; }) {
    assert(dynamic_cast<const DeBruijnGraph*>(anno_graph.get_graph_ptr().get()));

    const auto &background = dynamic_cast<const DeBruijnGraph&>(anno_graph.get_graph());

    if (&foreground == &background)
        return;

    auto process_path = [&](std::string sequence) {
        if (terminate())
            return;

        DBGAligner<AlignmentCompare, ColumnCompare> variant_aligner(
            foreground,
            config,
            suffix_seeder<DeBruijnGraph::node_index>, // extend_mapping ignores this anyway
            build_background_graph_extender(background,
                                            Extender<DeBruijnGraph::node_index>(extender))
        );

        auto paths = variant_aligner.extend_mapping_forward_and_reverse_complement(
            sequence,
            0,
            build_breakpoint_seeder_builder(background)
        );

        assert(paths.get_query() == sequence);

        for (auto&& path : paths) {
            if (path.is_exact_match())
                continue;

            auto labels = anno_graph.get_labels(path.get_sequence(), 1.0);

            if (labels.empty())
                continue;

            callback(std::move(path),
                     path.get_orientation()
                         ? paths.get_query_reverse_complement()
                         : paths.get_query(),
                     std::move(labels));
        }
    };

    if (thread_pool) {
        foreground.call_unitigs(
            [&](const auto &sequence) {
                // don't add to the thread pool if this'll terminatate immediately
                if (terminate())
                    return;

                thread_pool->enqueue(process_path, sequence);
            }
        );

        thread_pool->join();
    } else {
        foreground.call_unitigs(process_path);
    }
}


} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
