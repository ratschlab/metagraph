#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>

#include "graph/annotated_dbg.hpp"


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
