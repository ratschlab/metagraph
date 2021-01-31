#ifndef __TAXONOMIC_DB_HPP__
#define __TAXONOMIC_DB_HPP__

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace annot {

/**
 * TaxonomyDB constructs a taxonomic map (kmer to taxid) in a similar way
 * to how Kraken2 works. The file exported by this object will be further
 * used for taxonomic sequence classifier by query metagraph cli cmd.
 */
class TaxonomyDB {
  public:
    typedef std::string AccessionVersion;
    typedef std::uint64_t TaxId;
    typedef std::uint64_t NormalizedTaxId;
    typedef std::vector<std::vector<NormalizedTaxId>> ChildrenList;
    typedef annot::MultiLabelEncoded<AccessionVersion> Annotator;
    using KmerId = Annotator::Index;

  private:
    /**
     * node_depth returns the depth for each node.
     * The root is the unique node with maximal depth and all the leaves have depth 0.
     */
    std::vector<uint64_t> node_depth; // The root has the maximal depth;
    /**
     * rmq_data[0] contains the taxo tree linearization
     *         (e.g. for root 1 and 1-2 and 1-3 edges, the linearization is 1 2 1 3 1)
     * rmq_data[l][x] returns the node with maximal depth among the interval [x, x+2^l - 1]
     * in the linearization. (e.g. rmq_data[3][6] return the node with max depth in [6, 13])
     */
    std::vector<std::vector<NormalizedTaxId>> rmq_data;
    /**
     * node_to_linearization_idx[node] returns the index of the first occurence of node
     * in the linearization.
     */
    std::vector<uint64_t> node_to_linearization_idx;
    /**
     * lookup_table maps accession version to normalized taxids.
     */
    tsl::hopscotch_map<AccessionVersion, NormalizedTaxId> lookup_table;
    /**
     * node_to_acc_version maps normalized taxids to accession version.
     */
    std::vector<AccessionVersion> node_to_acc_version;
    /**
     * precalc_log2 is a table for faster compute of log2.
     */
    std::vector<uint64_t> precalc_log2;
    /**
     * precalc_pow2 is a table for faster compute of pow2.
     */
    std::vector<uint64_t> precalc_pow2;
    /**
     * taxonomic_map returns the taxid LCA for a given kmer.
     */
    tsl::hopscotch_map<KmerId, NormalizedTaxId> taxonomic_map;

    /**
     * Reads the taxonomic tree
     * @param [input] taxo_tree_filepath path to a "nodes.dmp" file.
     * @param [input] reversed_lookup_table maps taxid (before normalizing) to accession version.
     * @param [output] tree -> tree stored as list of children.
     * @param [output] root_node -> normalized id of the root of the tree.
     */
    void read_tree(const std::string &taxo_tree_filepath,
                   const tsl::hopscotch_map<TaxId, AccessionVersion> &reversed_lookup_table,
                   ChildrenList &tree,
                   NormalizedTaxId &root_node);
    /**
     * Reads the lookup table
     * @param [input] lookup_table_filepath path to a ".accession2taxid" file.
     * @param [input] input_accessions set containing all the accession version in the fasta input.
     * @param [output] reversed_lookup_table maps taxid (before normalizing) to accession version.
     */
    void read_lookup_table(const std::string &lookup_table_filepath,
                           const tsl::hopscotch_set<AccessionVersion> &input_accessions,
                           tsl::hopscotch_map<TaxId, AccessionVersion> &reversed_lookup_table);
    /**
     * Reads the fasta headers
     * @param [input] fasta_headers_filepath path to a ".fasta.fai" file.
     * @param [output] input_accessions set containing all the accession version in the fasta input.
     */
    void get_input_accessions(const std::string &fasta_headers_filepath,
                              tsl::hopscotch_set<AccessionVersion> &input_accessions);
    /**
     * Stores the node depths in "this->node_depth" by doing a dfs over the tree.
     * @param [input] node - the node that is currently processed.
     * @param [input] tree -> tree stored as list of children.
     */
    void dfs_node_depth(const NormalizedTaxId &node, const ChildrenList &tree);
    /**
     * Stores the rmq data in "this->rmq_data" by doing a linearization of the tree.
     * @param [input] tree -> tree stored as list of children.
     * @param [input] root_node -> normalized id of the root of the tree.
     */
    void rmq_preprocessing(const ChildrenList &tree, const NormalizedTaxId &root_node);
    /**
     * Stores the rmq data in "this->rmq_data" by doing a linearization of the tree.
     * @param [input] node - the node that is currently processed.
     * @param [input] tree -> tree stored as list of children.
     * @param [output] tree_linearization -> the linearization of the received tree.
     */
    void dfs_linearization(const NormalizedTaxId &node,
                           const ChildrenList &tree,
                           std::vector<NormalizedTaxId> &tree_linearization);
    NormalizedTaxId find_lca(const NormalizedTaxId &taxid1, const NormalizedTaxId &taxid2);
    NormalizedTaxId find_lca(const std::vector<NormalizedTaxId> &taxids);

  public:
    /**
     * Constructs a TaxonomyDB
     * @param [input] taxo_tree_filepath path to a "nodes.dmp" file.
     * @param [input] lookup_table_filepath path to a ".accession2taxid" file.
     * @param [input] fasta_headers_filepath path to a ".fasta.fai" file.
     */
    TaxonomyDB(const std::string &taxo_tree_filepath,
               const std::string &lookup_table_filepath,
               const std::string &fasta_headers_filepath);
    /**
     * Update the corresponding "this->taxonomic_map" for each received kmer to consider
     * the new data: taxonomic_map[kmer] = find_lca(taxonomic_map[kmer], lca);
     * @param [input] kmers - list of kmers found in an accession version.
     * @param [input] lca - taxid corresponding to an accession version that contains the kmers.
     */
    void update_taxonomic_map(const std::vector<KmerId> &kmers,
                              const NormalizedTaxId &lca);
    /**
     * Exports 'this->taxonomic_map', 'this->index_to_label' and the parent list for the
     * stored taxonomic tree to the given filepath.
     */
    void export_to_file(const std::string &filepath);
    bool find_lca(const std::vector<AccessionVersion> &accessions,
                  NormalizedTaxId &lca);

  private:
    // num_external_lca_calls counts the number of find_lca() calls only for logging purposes.
    uint64_t num_external_lca_calls;
    uint64_t num_external_lca_calls_failed;
};

}
}

#endif // __TAXONOMIC_HPP__