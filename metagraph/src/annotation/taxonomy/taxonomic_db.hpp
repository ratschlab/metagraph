#ifndef __TAXONOMIC_DB_HPP__
#define __TAXONOMIC_DB_HPP__

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/base/annotation.hpp"
#include "graph/annotated_dbg.hpp"
#include "cli/config/config.hpp"


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
     * The root is the unique node with maximal depth and all the leaves have depth 1.
     */
    std::vector<uint64_t> node_depth;

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
     * lookup_table maps accession version to taxid.
     */
    tsl::hopscotch_map<AccessionVersion, TaxId> lookup_table;

    /**
     * precalc_log2 is a table for faster compute of log2.
     */
    std::vector<uint64_t> precalc_log2;

    /**
     * precalc_pow2 is a table for faster compute of pow2.
     */
    std::vector<uint64_t> precalc_pow2;

    /**
     * Maps taxid to its normalized index. Used for optimizing the runtime performance.
     */
    tsl::hopscotch_map<TaxId, NormalizedTaxId> normalized_taxid;
    std::vector<TaxId> denormalized_taxid;

    /**
    * taxonomic_map returns the taxid LCA for a given kmer.
    */
    sdsl::int_vector<> taxonomic_map;

    /*  *
     * Reads and returns the taxonomic tree
     *
     * @param [input] taxo_tree_filepath path to a "nodes.dmp" file.
     * @param [output] tree -> tree stored as list of children.
     * @param [output] root_node -> normalized id of the root of the tree.
     */
    void read_tree(const std::string &taxo_tree_filepath,
                   ChildrenList &tree, NormalizedTaxId &root_node);

    /**
     * Reads and returns the lookup table corresponding to the received accession versions.
     *
     * @param [input] lookup_table_filepath path to a ".accession2taxid" file.
     * @param [input] input_accessions contains all the accession version in the input.
     */
    void read_lookup_table(const std::string &lookup_table_filepath,
                           const tsl::hopscotch_set<AccessionVersion> &input_accessions);

    /**
     * Computes the rmq data in "this->rmq_data". Beside this, calculates
     *      the fast tables: 'this->precalc_log' and 'this->precalc_pow2'.
     *
     * @param [input] tree_linearization -> the linearization of the received tree.
     */
    void rmq_preprocessing(const std::vector<NormalizedTaxId> &tree_linearization);

    /**
     * dfs_statistics calculates a tree_linearization, "this->node_depth" and
     *      'this->node_to_linearization_idx'.
     *
     * @param [input] node - the node that is currently processed.
     * @param [input] tree -> tree stored as list of children.
     * @param [output] tree_linearization -> the linearization of the received tree.
     */
    void dfs_statistics(const NormalizedTaxId &node, const ChildrenList &tree,
                        std::vector<NormalizedTaxId> &tree_linearization);

  public:
    /**
     * Constructs a TaxonomyDB
     *
     * @param [input] taxo_tree_filepath path to a "nodes.dmp" file.
     * @param [input] lookup_table_filepath path to a ".accession2taxid" file.
     * @param [input] input_accessions contains all the accession version in the annotation matrix input.
     */
    TaxonomyDB(const std::string &taxo_tree_filepath,
               const std::string &lookup_table_filepath,
               const tsl::hopscotch_set<AccessionVersion> &input_accessions);

    /**
     * Iterate the annotation matrix (kmer-OX labels-OY) from the received set of files
     * for updating the LCA taxid per kmer. The new data is stored in "this->taxonomic_map".
     *
     * @param [input] anno_graph - the annotation matrix object.
     */
    void kmer_to_taxid_map_update(const std::vector<std::string> &files, cli::Config *config);

    /**
     * Exports 'taxonomic_map' and the taxonomic tree (as parent list)
     * to the given filepath.
     */
    void export_to_file(const std::string &filepath);
    NormalizedTaxId find_lca(const std::vector<NormalizedTaxId> &taxids) const;

    bool get_normalized_taxid(const std::string accession_version, NormalizedTaxId &taxid) const;
    static std::string get_accession_version_from_label(const std::string &label);

  private:
    // num_external_get_taxid_calls and num_external_get_taxid_calls_failed used only for logging purposes.
    static uint64_t num_external_get_taxid_calls;
    static uint64_t num_external_get_taxid_calls_failed;
};

}
}

#endif // __TAXONOMIC_DB_HPP__
