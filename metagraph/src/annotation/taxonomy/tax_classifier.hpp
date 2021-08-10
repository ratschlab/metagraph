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
        UNASSIGNED,
        GEN_BANK, // e.g. ">gi|1070643132|ref|NC_031224.1| Arthrobacter phage Mudcat, complete genome"
        TAXID, // e.g. ">kraken:taxid|2016032|NC_047834.1 Alteromonas virus vB_AspP-H4/4, complete genome"
    };

    TaxonomyBase() {};
    TaxonomyBase(const double lca_coverage_rate, const double kmers_discovery_rate)
        : lca_coverage_rate_(lca_coverage_rate),
          kmers_discovery_rate_(kmers_discovery_rate) {};

    TaxId assign_class(const std::string &sequence) const;

  protected:
    virtual TaxId find_lca(const std::vector<TaxId> &taxids) const = 0;

    std::string get_accession_version_from_label(const std::string &label) const;

    bool get_taxid_from_label(const std::string &label, TaxId *taxid) const;

    /** Reads the accession version to taxid lookup table.
     * If 'anno_matrix' is mentioned, only the labels that exist in the given annotation matrix will be stored.
     * If 'anno_matrix' is not mentioned, the entire content of 'filepath' will be read and stored.
     *
     * @param [input] filepath -> a ".accession2taxid" file.
     * @param [optional input] anno_matrix -> pointer to the annotation matrix.
     */
    void read_accversion_to_taxid_map(const std::string &filepath, const graph::AnnotatedDBG *anno_matrix);

    /**
     * Update the current node_scores and best_lca by taking into account the weight of the start_node and
     * all its ancestors.
     *
     * @param [input] 'start_node' -> the starting node to update 'node_scores'.
     * @param [input] 'num_kmers_per_node[taxid]' -> the number of kmers 'k' with taxonomic_map[k]=taxid.
     * @param [input] 'desired_number_kmers' -> the threshold score that a node has to exceed in order to be
     *                                          considered as a valid solution.
     * @param [modified] 'node_scores' -> the current score for each node in the tree.
     * @param [modified] 'nodes_already_propagated' -> the set of nodes that were previously processed.
     * @param [modified] 'best_lca' -> the current classification prediction (node that exceeds the
     *                                 `desired_number_kmers` threshold and is placed as close as possible to the leaves).
     * @param [modified] 'best_lca_dist_to_root' -> the distance to the root for the current classification prediction.
     */
    void update_scores_and_lca(const TaxId start_node,
                               const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                               const uint64_t desired_number_kmers,
                               tsl::hopscotch_map<TaxId, uint64_t> *node_scores,
                               tsl::hopscotch_set<TaxId> *nodes_already_propagated,
                               TaxId *best_lca,
                               uint32_t *best_lca_dist_to_root) const;

    /**
     * Get the list of LCA taxid for each kmer in a given sequences.
     * The sequence can be given in forward or in reversed orientation.
     */
    virtual std::vector<TaxId> get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const = 0;

    LabelType label_type_ = LabelType::UNASSIGNED;

    /**
     * node_depth returns the depth for each node in the taxonomic tree.
     * The root is the unique node with maximal depth and all the leaves have depth equal to 1.
     */
    tsl::hopscotch_map<TaxId, uint32_t> node_depth_;

    TaxId root_node_;

    /**
     *  node_parent stores a taxonomic tree representation as a taxid to taxid parent list.
     */
    tsl::hopscotch_map<TaxId, TaxId> node_parent_;

    tsl::hopscotch_map<std::string, TaxId> accversion_to_taxid_map_;

    double lca_coverage_rate_;
    double kmers_discovery_rate_;
};

class TaxonomyClsImportDB : public TaxonomyBase {
  public:
    // todo implement
    TaxonomyClsImportDB(const std::string &taxdb_filepath,
                        const double lca_coverage_rate,
                        const double kmers_discovery_rate);

  private:
    std::vector<TaxId> get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const;
    TaxId find_lca(const std::vector<TaxId> &taxids) const;
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
     * @param [optional input] label_taxid_map_filepath ->  path to the acc_version to taxid lookup table (".accession2taxid").
     *                                                      Mandatory if the taxid is not mentioned in the label string.
     */
    TaxonomyClsAnno(const graph::AnnotatedDBG &anno,
                    const double lca_coverage_rate,
                    const double kmers_discovery_rate,
                    const std::string &tax_tree_filepath,
                    const std::string &label_taxid_map_filepath = "");
    TaxonomyClsAnno() {};
    virtual ~TaxonomyClsAnno() {};

    // todo implement
    void export_taxdb(const std::string &filepath) const;

    /** assign_class_toplabels is a faster and less precise taxonomic classification approach.
     * It takes into account only the labels with the highest frequency among the kmers in one read.
     *
     * Given a set of kmers taken from one read, a query on the annotation matrix will return all the labels that
     * correspond to at least $label_fraction$ percent of the kmers.
     * Then, returns the LCA of the taxid nodes where these frequent labels are mapping to.
     *
     * @param [input] sequence -> the given sequence to classify.
     * @param [input] label_fraction -> threshold used for taxonomic classification.
     * @return -> the classification result: a taxid node
     */
    TaxId assign_class_toplabels(const std::string &sequence, const double label_fraction) const;

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
    void dfs_statistics(const TaxId node,
                        const ChildrenList &tree,
                        std::vector<TaxId> *tree_linearization);

    TaxId find_lca(const std::vector<TaxId> &taxids) const;

    std::vector<TaxId> get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const;

    /**
     * rmq_data[0] contains the taxonomic tree linearization
     *          (e.g. for root 1 and edges={1-2; 1-3}, the linearization is "1 2 1 3 1").
     * rmq_data[l][x] returns the node with the maximal depth among positions [x, x+2^l-1] in the linearization
     *          (e.g. rmq_data[3][6] return the node with max depth in [6, 13]).
     */
    std::vector<std::vector<TaxId>> rmq_data_;

    /**
     * node_to_linearization_idx[node] returns the index of the first occurrence of node
     * in the tree linearization order. This array will be further used inside a RMQ query.
     */
    tsl::hopscotch_map<TaxId, uint32_t> node_to_linearization_idx_;

    const graph::AnnotatedDBG *anno_matrix_ = NULL;
};

} // namespace annot
} // namespace mtg

#endif // __TAX_CLASSIFIER_HPP__
