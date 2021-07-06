#ifndef __TAXONOMIC_CLASSIFIER_HPP__
#define __TAXONOMIC_CLASSIFIER_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "annotation/representation/base/annotation.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include <sdsl/int_vector.hpp>
#include <sdsl/dac_vector.hpp>

namespace mtg {
namespace annot {

/**
 * TaxClassifier imports data from TaxonomicDB and implements a method for taxonomic
 * classification a given DNA sequencing read.
 */
class TaxClassifier {
  public:
    using TaxId = std::uint64_t;
    using Annotator = annot::MultiLabelEncoded<std::string>;
    using KmerId = Annotator::Index;
    using DeBruijnGraph = mtg::graph::DeBruijnGraph;

    TaxClassifier(const std::string &taxdb_filepath,
                  const double lca_coverage_rate,
                  const double kmers_discovery_rate);
    TaxClassifier(){};

    /**
     * Assign a LCA taxid to a given DNA sequencing read.
     * Consider matches[node] = number of kmers (in the given read) that are pointing to 'node'.
     *          weight[node] = matches[node] / (the total number of kmers). (To obtain values in [0, 1])
     *          score[node] = sum(weight[node*]) where 'node*' is any of node's descendants or node's ancestors (Obtain values in [0, 1])
     * The assigned taxid is the farthest node to the root with score[node] >= lca_coverage_threshold (the node with highest score in case of equal distances).
      */
    TaxId assign_class(const mtg::graph::DeBruijnGraph &graph,
                       const std::string &sequence) const;

    static TaxId find_lca(const TaxId a,
                          const TaxId b,
                          const TaxId root_node,
                          const tsl::hopscotch_map<TaxId, TaxId> &node_parent);

    /**
     * Update the current node_scores and best_lca by taking into account the kmers in start_node and all it's ancestors.
     *
     * @param [input] 'start_node' starting node to update ancestors and descendants in the taxonomic tree
     * @param [input] 'num_kmers_per_node[taxid]' the number of kmers that point to 'taxid' according to the 'taxonomic_map'.
     * @param [input] 'desired_number_kmers' is the threshold score that a node has to exceed in order to be considered as a valid solution.
     * @param [input] 'root_node' = a copy of "this->root_node" for a static function
     * @param [input] 'node_parent' = a copy of "this->node_parent" for a static function
     * @param [modified] 'node_scores' the current scores for each node in the tree.
     * @param [modified] 'nodes_already_propagated' list of nodes that were already processed.
     * @param [modified] 'best_lca' the node furthest to the root that exceeds the `desired_number_kmers` threshold (in case of equal distances - the one with the highest score).
     * @param [modified] 'best_lca_dist_to_root' represents the distance to the root for the current best lca taxid.
     */
     static void update_scores_and_lca(const TaxId start_node,
                                       const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                                       const uint64_t desired_number_kmers,
                                       const TaxId root_node,
                                       const tsl::hopscotch_map<TaxId, TaxId> &node_parent,
                                       tsl::hopscotch_map<TaxId, uint64_t> *node_scores,
                                       tsl::hopscotch_set<TaxId> *nodes_already_propagated,
                                       TaxId *best_lca,
                                       uint64_t *best_lca_dist_to_root);
private:
    /**
     * Import 'this->node_parents', 'this->code_to_taxid', 'this->code_to_taxid'
     * from the given taxDB filepath (created by './metagraph transform_anno_tax').
     *
     * -> "this->node_parents" represents the taxonomic tree given as list of parents
     * -> The taxonomic_map[k] (kmer to lca taxid) can be composed as code_to_taxid[code[k]]
     */
    void import_taxonomy(const std::string &taxdb_filepath);

    TaxId root_node;

    /**
     * node_parent[node] returns the taxid of node's parent.
    */
    tsl::hopscotch_map<TaxId, TaxId> node_parent;

    /**
     * The code_to_taxid[code[i]] (taxonomic_map[i]) returns the taxid LCA for a given kmer 'i'.
     */
    sdsl::dac_vector_dp<sdsl::rrr_vector<>> code;
    sdsl::int_vector<> code_to_taxid;

    double lca_coverage_rate;
    double kmers_discovery_rate;
};

} // namespace annot
} // namespace mtg

#endif // __TAXONOMIC_CLASSIFIER_HPP__
