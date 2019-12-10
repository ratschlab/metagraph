#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "common/threading.hpp"
#include "graph/annotated_dbg.hpp"
#include "utils/bit_vectors/bitmap.hpp"

typedef std::function<size_t()> LabelCountCallback;

template <typename NodeType>
class Alignment;

namespace annotated_graph_algorithm {

typedef std::function<bool(const std::string&,
                           const std::vector<DeBruijnGraph::node_index>&)> KeepUnitigPath;

// Given a DeBruijnGraph and a bool-returning string callback, return a bitmap of
// length graph.max_index() + 1. An index is set to 1 if it is contained in a
// unitig satisfying keep_unitig(unitig).
std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig(const DeBruijnGraph &graph,
                     const KeepUnitigPath &keep_unitig);

// Given an AnnotatedDBG and sets of foreground (in) and background (out) labels,
// return a bitmap of length anno_graph.get_graph().max_index() + 1. An index i
// is set to 1 if there is a unitig containing i such that
// at least (label_mask_in_fraction * 100)% of the total possible number of in labels is present,
// at most (label_mask_out_fraction * 100)% of the total possible number of out labels is present,
// and (label_other_fraction * 100)% of the labels are neither in or out masked.
std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig_labels(const AnnotatedDBG &anno_graph,
                            const std::vector<AnnotatedDBG::Annotator::Label> &labels_in,
                            const std::vector<AnnotatedDBG::Annotator::Label> &labels_out,
                            double label_mask_in_fraction = 1.0,
                            double label_mask_out_fraction = 0.0,
                            double label_other_fraction = 1.0);

// Given an AnnotatedDBG and vectors of foreground (labels_in) and
// background (labels_out) labels, construct a bitmap of length
// anno_graph.get_graph().max_index() + 1. An index i is set to 1 if
//
// is_node_in_mask(i,
//                 get # of labels of node i in labels_in,
//                 get # of labels of node i in labels_out)
//
// where the getters for the numbers of overlapping labels are callbacks of
// type LabelCountCallback.
//
// For labels which appear in more than
// num_nodes() * lazy_evaluation_label_frequency_cutoff
// nodes, precompute the numbers of overlapping labels.
std::unique_ptr<bitmap>
mask_nodes_by_node_label(const AnnotatedDBG &anno_graph,
                         const std::vector<AnnotatedDBG::Annotator::Label> &labels_in,
                         const std::vector<AnnotatedDBG::Annotator::Label> &labels_out,
                         std::function<bool(DeBruijnGraph::node_index,
                                            LabelCountCallback, /* get_num_labels_in */
                                            LabelCountCallback /* get_num_labels_out */)> is_node_in_mask,
                         double min_frequency_for_frequent_label = 0.05);

template <class Index, typename... Args>
using VariantCallback = std::function<void(Alignment<Index>&&,
                                           const std::string&, // query sequence
                                           Args&&...)>;

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

} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
