#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "graph/annotated_dbg.hpp"
#include "common/vectors/bitmap.hpp"


namespace mtg {
namespace graph {


class MaskedDeBruijnGraph;

/**
 * A container for differential assembly parameters. These are:
 * label_mask_in_kmer_fraction: minimum fraction of in-labels which should be
 *                              present for a node to be kept
 * label_mask_in_unitig_fraction: minimum fraction of nodes within a unitig which
 *                                must satisfy the in-label criteria
 * label_mask_out_kmer_fraction: maximum fraction of out-labels which can be
 *                               present in a kept node
 * label_mask_out_unitig_fraction: maximum fraction of nodes within a unitig which
 *                                 can satisfy the out-label criteria
 * label_mask_other_unitig_fraction: maximum fraction of nodes within a unitig which
 *                                   can contain labels not in the in- and out-groups
 * add_complement: also consider the reverse complements of nodes
 */
struct DifferentialAssemblyConfig {
    double label_mask_in_unitig_fraction = 1.0;
    double label_mask_in_kmer_fraction = 1.0;
    double label_mask_out_unitig_fraction = 0.0;
    double label_mask_out_kmer_fraction = 0.0;
    double label_mask_other_unitig_fraction = 1.0;
    bool add_complement = false;
};

/**
 * Return an int_vector<>, bitmap pair, each of length anno_graph.get_graph().max_index().
 * For an index i, the int_vector will contain a packed integer representing the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with. The least significant half of each integer represents the count
 * from labels_in, while the most significant half represents the count from
 * labels_out.
 * The returned bitmap is a binarization of the int_vector
 */
std::pair<sdsl::int_vector<>, std::unique_ptr<bitmap>>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<AnnotatedDBG::Annotator::Label> &labels_in,
                  const std::vector<AnnotatedDBG::Annotator::Label> &labels_out,
                  size_t num_threads);


/**
 * Given an AnnotatedDBG and sets of foreground (in) and background (out) labels,
 * return a MaskedDeBruijnGraph with the nodes of anno_graph masked according to
 * the parameters specified by config.
 */
MaskedDeBruijnGraph
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<typename AnnotatedDBG::Annotator::Label> &labels_in,
                    const std::vector<typename AnnotatedDBG::Annotator::Label> &labels_out,
                    const std::vector<typename AnnotatedDBG::Annotator::Label> &labels_in_post,
                    const std::vector<typename AnnotatedDBG::Annotator::Label> &labels_out_post,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads = 1);

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
