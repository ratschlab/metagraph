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


typedef std::function<size_t()> LabelCountCallback;

namespace annotated_graph_algorithm {


// Given an AnnotatedDBG and vectors of foreground (labels_in) and
// background (labels_out) labels, construct a bitmap of length
// anno_graph.num_nodes() + 1. An index i is set to 1 if
//
// is_node_in_mask(get # of labels of node i in labels_in,
//                 get # of labels of node i in labels_out)
//
// where the getters for the numbers of overlapping labels are callbacks of
// type LabelCountCallback.
//
// For labels which appear in more than
// num_nodes() * lazy_evaluation_label_frequency_cutoff
// nodes, precompute the numbers of overlapping labels.
std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &labels_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &labels_out,
                    std::function<bool(LabelCountCallback /* get_num_labels_in */,
                                       LabelCountCallback /* get_num_labels_out */)> is_node_in_mask,
                    double lazy_evaluation_label_frequency_cutoff = 0.05);

std::unique_ptr<bitmap>
mask_nodes_by_unitig_label(const AnnotatedDBG &anno_graph,
                           const std::vector<AnnotatedDBG::Annotator::Label> &labels_in,
                           const std::vector<AnnotatedDBG::Annotator::Label> &labels_out,
                           double unitig_labels_in_admixture = 1.0,
                           double unitig_labels_out_admixture = 0.0,
                           double lazy_evaluation_label_frequency_cutoff = 0.05);

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
                            Extender<DeBruijnGraph::node_index>&& extender = default_extender<>);

// Seed until the first breakpoint, relative to the given foreground graph
template <typename NodeType = typename DeBruijnGraph::node_index>
class BreakpointSeeder : public MEMSeeder<NodeType> {
  public:
    BreakpointSeeder(const DeBruijnGraph &foreground, const DBGAlignerConfig &config)
          : MEMSeeder<NodeType>(
                foreground,
                config,
                std::make_unique<bitmap_lazy>(
                    [&, background = dynamic_cast<const MaskedDeBruijnGraph*>(&foreground)](auto i) {
                        assert(background);
                        assert(i != DeBruijnGraph::npos);
                        return background->get_graph().outdegree(i)
                            > foreground.outdegree(i);
                    },
                    foreground.num_nodes() + 1
                )
            ) {}
};

// Given a seed, only extend on the background graph, relative to the
// given foreground graph
Extender<DeBruijnGraph::node_index>
build_background_graph_extender(const DeBruijnGraph &foreground,
                                Extender<DeBruijnGraph::node_index>&& extender = default_extender<>);

template <class AlignmentCompare = std::less<Alignment<DeBruijnGraph::node_index>>,
          class ColumnCompare = std::less<typename Alignment<DeBruijnGraph::node_index>::Column>>
void call_variants(const DeBruijnGraph &foreground,
                   const AnnotatedDBG &anno_graph,
                   const VariantLabelCallback &callback,
                   const DBGAlignerConfig &config,
                   const Extender<DeBruijnGraph::node_index> &extender = default_extender<>,
                   ThreadPool *thread_pool = nullptr,
                   const std::function<bool()> &terminate = []() { return false; }) {
    assert(dynamic_cast<const DeBruijnGraph*>(anno_graph.get_graph_ptr().get()));
    assert(dynamic_cast<const MaskedDeBruijnGraph*>(&foreground));

    const auto &background = dynamic_cast<const DeBruijnGraph&>(anno_graph.get_graph());

    auto process_path = [&](std::string sequence) {
        if (terminate())
            return;

        DBGAligner<BreakpointSeeder<>, AlignmentCompare, ColumnCompare> variant_aligner(
            foreground,
            config,
            build_background_graph_extender(background,
                                            Extender<DeBruijnGraph::node_index>(extender))
        );

        auto paths = variant_aligner.align(sequence, 0);

        for (auto&& path : paths) {
            assert(path.is_valid(background));

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
