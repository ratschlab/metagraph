#ifndef __TAX_CLASSIFIER_HPP__
#define __TAX_CLASSIFIER_HPP__

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace annot {

using TaxId = std::uint32_t;
using ChildrenList = tsl::hopscotch_map<TaxId, std::vector<TaxId>>;

class TaxonomyBase {
  public:
    using KmerId = annot::MultiLabelEncoded<std::string>::Index;
    using node_index = graph::SequenceGraph::node_index;

    enum LabelType {
        GEN_BANK, // e.g. ">gi|1070643132|ref|NC_031224.1| Arthrobacter phage Mudcat, complete genome"
        TAXID, // e.g. ">kraken:taxid|2016032|NC_047834.1 Alteromonas virus vB_AspP-H4/4, complete genome"
    };

    TaxonomyBase() {}
    TaxonomyBase(double lca_coverage_rate, double kmers_discovery_rate)
        : lca_coverage_rate_(lca_coverage_rate),
          kmers_discovery_rate_(kmers_discovery_rate) {}
    virtual ~TaxonomyBase() {}

  protected:
    std::string get_accession_version_from_label(const std::string &label) const;

    /** Reads the accession version to taxid lookup table.
     * If 'anno_matrix' is not NULL, only the labels that exist in the given annotation matrix will be stored.
     * If 'anno_matrix' is NULL, the entire content of 'filepath' will be read and stored.
     *
     * @param [input] filepath -> a ".accession2taxid" file.
     * @param [optional input] anno_matrix -> pointer to the annotation matrix.
     */
    void read_accversion_to_taxid_map(const std::string &filepath,
                                      const graph::AnnotatedDBG *anno_matrix = NULL);

    LabelType label_type_;

    /**
     * node_depth_ returns the depth for each node in the taxonomic tree.
     * The root is the unique node with maximal depth and all the leaves have depth equal to 1.
     */
    tsl::hopscotch_map<TaxId, uint32_t> node_depth_;

    TaxId root_node_;

    /**
     *  node_parent_ stores a taxonomic tree representation as a taxid to taxid parent list.
     */
    tsl::hopscotch_map<TaxId, TaxId> node_parent_;

    tsl::hopscotch_map<std::string, TaxId> accversion_to_taxid_map_;

    double lca_coverage_rate_;
    double kmers_discovery_rate_;
};

class TaxonomyClsAnno : public TaxonomyBase {
  public:
    /**
     * TaxonomyCls constructor
     *
     * @param [input] anno -> the annotation matrix
     * @param [input] lca_coverage_rate -> threshold used for taxonomic classification.
     * @param [input] kmers_discovery_rate -> threshold used for taxonomic classification.
     * @param [input] tax_tree_filepath ->  path to a taxonomic tree ("nodes.dmp" file).
     * @param [optional input] label_taxid_map_filepath ->  path to the acc_version to
     * taxid lookup table (".accession2taxid"). Mandatory if the taxid is not mentioned in
     * the label string.
     */
    TaxonomyClsAnno(const graph::AnnotatedDBG &anno,
                    const std::string &tax_tree_filepath,
                    double lca_coverage_rate = 0,
                    double kmers_discovery_rate = 0,
                    const std::string &label_taxid_map_filepath = "");
    TaxonomyClsAnno() {}

  private:
    /**
     * Reads and returns the taxonomic tree as a list of children.
     *
     * @param [input] tax_tree_filepath -> path to a "nodes.dmp" file.
     * @param [output] tree -> tree stored as a list of children.
     */
    void read_tree(const std::string &tax_tree_filepath,
                   ChildrenList *tree);

    /**
     * rmq_preprocessing computes 'this->rmq_data' field.
     *
     * @param [input] tree_linearization -> the linearization of the taxonomic tree.
     */
    void rmq_preprocessing(const std::vector<TaxId> &tree_linearization);

    /**
     * dfs_statistics method calculates the following fields:
     *      + tree_linearization;
     *      + this->node_depth;
     *      + this->node_to_linearization_idx.
     *
     * @param [input] node -> the node that is currently processed.
     * @param [input] tree -> the taxonomic tree stored as a list of children.
     * @param [output] tree_linearization -> the linearization of the received tree.
     */
    void dfs_statistics(TaxId node,
                        const ChildrenList &tree,
                        std::vector<TaxId> *tree_linearization);

    /**
     * rmq_data_[0] contains the taxonomic tree linearization
     *          (e.g. for root 1 and edges={1-2; 1-3}, the linearization is "1 2 1 3 1").
     * rmq_data_[l][x] returns the node with the maximal depth among positions [x, x+2^l-1] in the linearization
     *          (e.g. rmq_data_[3][6] return the node with max depth in [6, 13]).
     */
    std::vector<std::vector<TaxId>> rmq_data_;

    /**
     * node_to_linearization_idx_[node] returns the index of the first occurrence of node
     * in the tree linearization order. This array will be further used inside a RMQ query.
     */
    tsl::hopscotch_map<TaxId, uint32_t> node_to_linearization_idx_;

    const graph::AnnotatedDBG *anno_matrix_ = NULL;
};

} // namespace annot
} // namespace mtg

#endif // __TAX_CLASSIFIER_HPP__
