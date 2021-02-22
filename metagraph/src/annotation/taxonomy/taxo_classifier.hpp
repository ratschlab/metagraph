#ifndef __TAXONOMIC_CLASSIFIER_HPP__
#define __TAXONOMIC_CLASSIFIER_HPP__

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

  public:
    /**
     * Construct a TaxoClassifier
     *
     * @param [input] filepath to the file exported by TaxonomyDB.
     */
    TaxoClassifier(const std::string &filepath);

    /**
     * Assign a LCA taxid to a given sequence.
     * Consider matches[node] = number of kmers in 'sequence' for which the taxonomic_map points to 'node'.
     *          weight[node] = matches[node] / #(kmers in sequence). (Values in [0, 1])
     *          score[node] = sum(weight[node*]) where node* is iterating over node's subtree + node's ancestors (Values in [0, 1])
     * The assigned taxid is the farthest node to the root with score[node] >= lca_coverage_threshold (unique).
      */
    TaxId assign_class(const mtg::graph::DeBruijnGraph &graph,
                       const std::string &sequence,
                       const double &lca_coverage_threshold);
};

}
}
#endif // __TAXONOMIC_CLASSIFIER_HPP__
