#ifndef __TAXONOMIC_DB_HPP__
#define __TAXONOMIC_DB_HPP__

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/base/annotation.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace annot {

/**
 * TaxonomyDB constructs a taxonomic map (kmer to taxid) in a similar way
 * to how Kraken2 works. The file exported by this object will be further
 * used by taxonomic DNA sequencing read classifier in './metagraph tax_class'.
 */
class TaxonomyDB {
  public:
    using AccessionVersion = std::string;
    using TaxId = std::uint64_t;
    using NormalizedTaxId = std::uint64_t;
    using ChildrenList = std::vector<std::vector<NormalizedTaxId>>;
    using Annotator = annot::MultiLabelEncoded<AccessionVersion>;
    using KmerId = Annotator::Index;
    using node_index = graph::SequenceGraph::node_index;

    /**
     * TaxonomyDB constructor
     *
     * @param [input] tax_tree_filepath path to a "nodes.dmp" file.
     * @param [input] label_taxid_map_filepath path to a ".accession2taxid" file.
     * @param [input] input_accessions contains all the accession version in the annotation matrix input.
     */
    TaxonomyDB(const std::string &tax_tree_filepath,
               const std::string &label_taxid_map_filepath,
               const tsl::hopscotch_set<AccessionVersion> &input_accessions);

    /**
     * Iterate the received annotation matrix (kmer-OX labels-OY) for updating the
     * LCA taxid per kmer. The updated data is stored in "this->taxonomic_map".
     *
     * @param [input] annot - the annotation matrix object.
     */
    void kmer_to_taxid_map_update(const annot::MultiLabelEncoded<std::string> &annot);

    /**
     * Exports 'taxonomic_map' and the taxonomic tree (as parent list)
     * to the given filepath.
     */
    void export_to_file(const std::string &filepath);

    /**
     * Find LCA for a set of nodes in the tree.
     */
    NormalizedTaxId find_normalized_lca(const std::vector<NormalizedTaxId> &taxids) const;

    bool get_normalized_taxid(const std::string accession_version, NormalizedTaxId *taxid) const;
    static std::string get_accession_version_from_label(const std::string &label);
    bool get_normalized_taxid_from_label(const std::string &label, NormalizedTaxId *taxid) const;

    TaxId assign_class_getrows(const graph::AnnotatedDBG &anno,
                               const std::string &sequence,
                               const double lca_coverage_rate,
                               const double kmers_discovery_rate) const;

    TaxId assign_class_toplabels(const graph::AnnotatedDBG &anno,
                                 const std::string &sequence,
                                 const double label_fraction) const;

  private:
    /**
     * Reads and returns the taxonomic tree as a children list.
     *
     * @param [input] tax_tree_filepath path to a "nodes.dmp" file.
     * @param [output] tree -> tree stored as a children list.
     */
    void read_tree(const std::string &tax_tree_filepath,
                   ChildrenList *tree);

    /**
     * Reads and saves to 'this->label_taxid_map[]' the accession version to taxid map corresponding to the received set of used accession versions.
     *
     * @param [input] label_taxid_map_filepath path to a ".accession2taxid" file.
     * @param [input] input_accessions contains all the accession version in the input.
     */
    void read_label_taxid_map(const std::string &label_taxid_map_filepath,
                              const tsl::hopscotch_set<AccessionVersion> &input_accessions);

    /**
     * Computes the rmq data in this->rmq_data. Beside this, calculates
     *      the fast tables: this->precalc_log and this->precalc_pow2.
     *
     * @param [input] tree_linearization -> the linearization of the received tree.
     */
    void rmq_preprocessing(const std::vector<NormalizedTaxId> &tree_linearization);

    /**
     * dfs_statistics calculates a tree_linearization, this->node_depth and this->node_to_linearization_idx.
     *
     * @param [input] node -> the node that is currently processed.
     * @param [input] tree -> tree stored as list of children.
     * @param [output] tree_linearization -> the linearization of the received tree.
     */
    void dfs_statistics(const NormalizedTaxId node,
                        const ChildrenList &tree,
                        std::vector<NormalizedTaxId> *tree_linearization);

    /**
     * node_depth returns the depth for each node in the taxonomic tree.
     * The root is the unique node with maximal depth and all the leaves have depth 1.
     */
    std::vector<uint64_t> node_depth;
    TaxId root_node;
    tsl::hopscotch_map<TaxId, TaxId> node_parent;

    /**
     * rmq_data[0] contains the taxonomic tree linearization
     *          (e.g. for root 1 and edges={1-2; 1-3}, the linearization is "1 2 1 3 1")
     * rmq_data[l][x] returns the node with the maximal depth among positions [x, x+2^l-1]
     * in the linearization.
     *          (e.g. rmq_data[3][6] return the node with max depth in [6, 13])
     */
    std::vector<std::vector<NormalizedTaxId>> rmq_data;

    /**
     * node_to_linearization_idx[node] returns the index of the first occurrence of node
     * in the linearization.
     */
    std::vector<uint64_t> node_to_linearization_idx;

    /**
     * label_taxid_map maps accession version to taxid.
     */
    tsl::hopscotch_map<AccessionVersion, TaxId> label_taxid_map;

    /**
     * fast_log2 is a table for a fast compute of log2(x).
     */
    std::vector<uint64_t> fast_log2;

    /**
     * fast_pow2 is a table for a fast compute of pow2(x).
     */
    std::vector<uint64_t> fast_pow2;

    /**
     * Maps taxid to its normalized index. Used for optimizing the runtime performance.
     */
    tsl::hopscotch_map<TaxId, NormalizedTaxId> normalized_taxid;

    /**
     * Maps normalized taxid to its denormalized counterpart.
     */
    std::vector<TaxId> denormalized_taxid;

    /**
    * taxonomic_map[kmer] returns the taxid LCA for the given kmer.
    */
    sdsl::int_vector<> taxonomic_map;

    // num_external_get_taxid_calls and num_external_get_taxid_calls_failed used only for logging purposes.
    static uint64_t num_get_taxid_calls;
    static uint64_t num_get_taxid_calls_failed;
};

} // namespace annot
} // namespace mtg

#endif // __TAXONOMIC_DB_HPP__
