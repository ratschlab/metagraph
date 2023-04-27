#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <tsl/hopscotch_set.h>

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
    bool add_complement = true;
    bool count_kmers = true;
    double pvalue = 0.05; // Differential assembly statistical test cutoff
    double family_wise_error_rate = 0; // The family wise rate for bonferonni multiple testing correction. If 0, no multiple testing is used.
    bool test_by_unitig = false;
    bool evaluate_assembly = false; // Myrthe temporary : evaluate an assembly
};

/**
 * Given an AnnotatedDBG and sets of foreground (in) and background (out) labels,
 * return a MaskedDeBruijnGraph with the nodes of anno_graph masked according to
 * the parameters specified by config. In round 1, a masked graph is constructed
 * containing the nodes corresponding to labels_in and labels_out, with a vector
 * counting how many of labels_in and labels_out correspond to each node,
 * respectively.
 * In round 2, if any of the sequences of the initial masked graph have the
 * labels in labels_in_round2 or labels_out_round2, or if they have other labels,
 * the count
 * vector is updated accordingly.
 * In round 3, the config rules are used to discard (i.e., mask out) nodes from
 * the masked graph. The resulting graph is returned.
 */
std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_in,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_out,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_in_round2,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_out_round2,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads = 1,
                    size_t num_parallel_files = std::numeric_limits<size_t>::max());

std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    const std::vector<std::string> &files,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_in,
                    const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_out,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads = 1,
                    size_t num_parallel_files = std::numeric_limits<size_t>::max());

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
