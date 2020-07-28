#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "common/threads/threading.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/masked_graph.hpp"


namespace mtg {
namespace graph {

typedef std::function<size_t()> LabelCountCallback;

template <typename NodeType>
class Alignment;

typedef std::function<bool(const std::string&,
                           const std::vector<DeBruijnGraph::node_index>&)> KeepUnitigPath;

typedef AnnotatedDBG::Annotator::Label Label;

// Given a DeBruijnGraph and a bool-returning string callback, return a bitmap of
// length graph.max_index() + 1. An index is set to 1 if it is contained in a
// unitig satisfying keep_unitig(unitig).
void mask_nodes_by_unitig(const DeBruijnGraph &graph,
                          const KeepUnitigPath &keep_unitig,
                          bitmap_vector *mask_to_update,
                          bool unset,
                          size_t num_threads = 1);


// Given an AnnotatedDBG and sets of foreground (in) and background (out) labels,
// return a bitmap of length anno_graph.get_graph().max_index() + 1. An index i
// is set to 1 if there is a unitig containing i such that
// at least (label_mask_in_fraction * 100)% of the total possible number of in labels is present,
// at most (label_mask_out_fraction * 100)% of the total possible number of out labels is present,
// and (label_other_fraction * 100)% of the labels are neither in or out masked.
MaskedDeBruijnGraph
make_masked_graph_by_unitig_labels(const AnnotatedDBG &anno_graph,
                                   const std::vector<Label> &labels_in,
                                   const std::vector<Label> &labels_out,
                                   size_t num_threads = 1,
                                   double label_mask_in_fraction = 1.0,
                                   double label_mask_out_fraction = 0.0,
                                   double label_other_fraction = 1.0,
                                   bool add_complement = false);

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
                         const std::vector<Label> &labels_in,
                         const std::vector<Label> &labels_out,
                         const std::function<bool(DeBruijnGraph::node_index,
                                                  const LabelCountCallback & /* get_num_labels_in */,
                                                  const LabelCountCallback & /* get_num_labels_out */)> &is_node_in_mask,
                         size_t num_threads = 1,
                         double min_frequency_for_frequent_label = 0.05);

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
