#ifndef __TAXONOMIC_CLASSIFIER_HPP__
#define __TAXONOMIC_CLASSIFIER_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "annotation/representation/base/annotation.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

/**
 * TaxoClassifier imports the data from TaxonomicDB and does the taxonomic
 * classification for a given sequence.
 */
class TaxoClassifier {
  public:
    typedef std::uint64_t TaxId;
    typedef annot::MultiLabelEncoded<std::string> Annotator;
    using KmerId = Annotator::Index;
    using DeBruijnGraph = mtg::graph::DeBruijnGraph;

    /**
     * Construct a TaxoClassifier
     *
     * @param [input] filepath to the file exported by TaxonomyDB.
     */
    TaxoClassifier(const std::string &filepath);
    TaxoClassifier(){};

    /**
     * Assign a LCA taxid to a given sequence.
     * Consider matches[node] = number of kmers in 'sequence' for which the taxonomic_map (LCA) points to 'node'.
     *          weight[node] = matches[node] / #(kmers in sequence). (To obtain values in [0, 1])
     *          score[node] = sum(weight[node*]) where node* is iterating over node's subtree + node's ancestors (Obtain values in [0, 1])
     * The assigned taxid is the farthest node to the root with score[node] >= lca_coverage_threshold (unique).
      */
    TaxId assign_class(const mtg::graph::DeBruijnGraph &graph,
                       const std::string &sequence,
                       const double &lca_coverage_threshold);

  private:
    TaxId root_node;

    /**
     * node_parent[node] returns the taxid of node's parent.
    */
    tsl::hopscotch_map<TaxId, TaxId> node_parent;

    /**
     * taxonomic_map returns the taxid LCA for a given kmer.
     */
    sdsl::int_vector<> taxonomic_map;

    /**
     * Import 'this->taxonomic_map' and the taxonomic tree (as parent list)
     * from the given filepath (created by TaxonomyDB in annotation cli).
     */
    void import_taxonomy(const std::string &filepath);

    /**
     * Update current node scores and best LCA taxid found.
     *
     * @param [input] 'start_node' starting node to update ancestors and descendants in the taxonomic tree
     * @param [input] 'num_kmers_per_node[taxid]' the number of kmers that point to 'taxid' according to the 'taxonomic_map'.
     * @param [input] 'desired_number_kmers' represents the threshold score that a node have to exceed in order to be considered the solution.
     * @param [modified] 'node_scores' the current scores for each node in the tree.
     * @param [modified] 'nodes_already_propagated' list of nodes that were already considered as 'start_node'.
     * @param [modified] 'best_lca' the node furthest to the root that exceeds the `desired_number_kmers` threshold.
     * @param [modified] 'best_lca_dist_to_root' represents the distance to the root for the current best lca taxid.
     */
     void update_scores_and_lca(const TaxId start_node,
                                const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                                const uint64_t &desired_number_kmers,
                                tsl::hopscotch_map<TaxId, uint64_t> &node_scores,
                                tsl::hopscotch_set<TaxId> &nodes_already_propagated,
                                TaxId &best_lca,
                                uint64_t &best_lca_dist_to_root);
};

} // namespace annot
} // namespace mtg
#endif // __TAXONOMIC_CLASSIFIER_HPP__
