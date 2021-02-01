#ifndef __TAXONOMIC_CLASSIFIER_HPP__
#define __TAXONOMIC_CLASSIFIER_HPP__

#include "annotation/representation/base/annotation.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

const double DEFAULT_LCA_COVERAGE_THRESHOLD = 0.66;
/**
 * TaxoClassifier imports the data from TaxonomicDB and does the taxonomic
 * classification for a given sequence.
 */
class TaxoClassifier {
  public:
    typedef std::uint64_t NormalizedTaxId;
    typedef annot::MultiLabelEncoded<std::string> Annotator;
    using KmerId = Annotator::Index;
    using DeBruijnGraph = mtg::graph::DeBruijnGraph;


  private:
    NormalizedTaxId root_node;

    /**
     * node_parent[node] returns the taxid of the corresponding parent.
    */
    std::vector<NormalizedTaxId> node_parent;

    /**
     * node_to_acc_version maps normalized taxid to accession version.
     */
    std::vector<std::string> node_to_acc_version;

    /**
     * taxonomic_map returns the taxid LCA for a given kmer.
     */
    tsl::hopscotch_map<KmerId, NormalizedTaxId> taxonomic_map;

    /**
     * Import 'this->taxonomic_map', 'this->node_to_acc_version' and the taxonomic tree (as parent list)
     * from the given filepath (created by TaxonomyDB in annotation cli).
     */
    void import_taxonomy(const std::string &filepath);

  public:
    /**
     * Construct a TaxoClassifier
     *
     * @param [input] filepath path to file exported by TaxonomyDB.
     */
    TaxoClassifier(const std::string &filepath);

    /**
     * Assign a LCA accession version to a given sequence.
     * Consider matches[node] = number of kmers in 'sequence' for which taxonomic_map points to 'node'.
     *          weight[node] = matches[node] / #(kmers in sequence). (Values in [0, 1])
     *          score[node] = sum(weight[node*]) where node* is in node's subtree or on root->node path. (Values in [0, 1])
     * The assigned LCA accession version is the farthest node to the root with score[node] >= lca_coverage_threshold (unique).
      */
    std::string assign_class(const mtg::graph::DeBruijnGraph &graph,
                             const std::string &sequence,
                             const double &lca_coverage_threshold = DEFAULT_LCA_COVERAGE_THRESHOLD);
};

}
}
#endif // __TAXONOMIC_CLASSIFIER_HPP__