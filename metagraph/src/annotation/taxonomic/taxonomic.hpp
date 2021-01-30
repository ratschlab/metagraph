#ifndef __TAXONOMIC_HPP__
#define __TAXONOMIC_HPP__

#include "annotation/representation/base/annotation.hpp"
#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>


namespace mtg {
namespace annot {


// TODO add some description.
class Taxonomy {
  public:
    typedef std::string AccessionVersion; // Maybe better to point this to the real structure.
    typedef std::uint64_t TaxId;
    typedef std::uint64_t TaxNormalizedId;
    typedef std::vector<std::vector<TaxNormalizedId>> ChildrenList;
    typedef annot::MultiLabelEncoded<AccessionVersion> Annotator;
    using KmerId = Annotator::Index;

  private:
    TaxNormalizedId root_node;
    uint64_t num_external_lca_calls;
    uint64_t num_external_lca_calls_failed;

//   TODO what is this rmq. https://stackoverflow.com/questions/15488470/syntax-for-dynamically-allocating-a-2d-array-of-smart-pointers
    std::vector<std::vector<TaxNormalizedId>> rmq_data;
    std::vector<uint64_t> node_depth; // The root has the maximal depth;
    std::vector<uint64_t> linearization_idx;
    tsl::hopscotch_map<TaxId, AccessionVersion> reversed_lookup_table;
    tsl::hopscotch_map<TaxId, TaxNormalizedId> normalized_tax_it;
    tsl::hopscotch_map<AccessionVersion, TaxNormalizedId> lookup_table;
    std::vector<AccessionVersion> index_to_label;
    std::vector<uint64_t> precalc_log;
    std::vector<uint64_t> precalc_pow2;
    tsl::hopscotch_map<KmerId, TaxNormalizedId> taxonomic_map;

    void read_tree(const std::string &tree_filepath,
                   ChildrenList &tree);
    void parse_lookup_table(const std::string &lookup_table_filepath,
                            const tsl::hopscotch_set<AccessionVersion> &input_accessions);
    void get_input_accessions(const std::string &fasta_fai,
                              tsl::hopscotch_set<AccessionVersion> &input_accessions);
    void calculate_node_depth(const TaxNormalizedId &node, const ChildrenList &tree);
    void rmq_preprocessing(const ChildrenList &tree);
    void dfs_linearization(const TaxNormalizedId &node,
                           const ChildrenList &tree,
                           std::vector<TaxNormalizedId> &tree_linearization);
    TaxNormalizedId find_lca(const TaxNormalizedId &label1, const TaxNormalizedId &label2);
    TaxNormalizedId find_lca(const std::vector<TaxNormalizedId> &labels);
  public:
    Taxonomy(const std::string &tree_filepath,
             const std::string &lookup_table_filepath,
             const std::string &fasta_fai);
    void update_kmers_lca(const std::vector<KmerId> &indices,
                          const TaxNormalizedId &lca);
    void export_to_file(const std::string &filepath);
    bool find_lca(const std::vector<AccessionVersion> &accessions,
                             TaxNormalizedId &lca);
};

}
}

#endif // __TAXONOMIC_HPP__
