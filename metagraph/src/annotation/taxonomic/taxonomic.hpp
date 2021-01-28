#ifndef __TAXONOMIC_HPP__
#define __TAXONOMIC_HPP__

#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace annot {


// TODO add some description.
class Taxonomy {
  public:
    typedef std::string Label; // Maybe better to point this to the real structure.
    typedef std::uint64_t TaxoLabel;
    typedef std::vector<std::vector<TaxoLabel>> ChildrenList;
    typedef annot::MultiLabelEncoded<Label> Annotator;
    using row_index = Annotator::Index;

  private:
    TaxoLabel root_node;

//   TODO what is this rmq. https://stackoverflow.com/questions/15488470/syntax-for-dynamically-allocating-a-2d-array-of-smart-pointers
    std::vector<std::vector<TaxoLabel>> rmq_data;
    std::vector<uint64_t> node_depth; // The root has the maximal depth;
    std::vector<uint64_t> linearization_idx;
    tsl::hopscotch_map<Label, TaxoLabel> label_to_index;
    std::vector<Label> index_to_label;
    std::vector<uint64_t> precalc_log;
    std::vector<uint64_t> precalc_pow2;
    tsl::hopscotch_map<row_index, TaxoLabel> taxonomic_map;

    void read_tree(const std::string &tree_filepath, ChildrenList &tree);
    void calculate_node_depth(const TaxoLabel &node, const ChildrenList &tree);
    void rmq_preprocessing(const ChildrenList &tree);
    void dfs_linearization(const TaxoLabel &node,
                           const ChildrenList &tree,
                           std::vector<TaxoLabel> &tree_linearization);
    TaxoLabel find_lca(const TaxoLabel &label1, const TaxoLabel &label2);
    TaxoLabel find_lca(const std::vector<Label> &labels);
    TaxoLabel find_lca(const std::vector<TaxoLabel> &labels);
  public:
    Taxonomy(const std::string &tree_filepath);
    void update_row_indices(const std::vector<row_index> &indices,
                            const std::vector<Label> &labels);
    void export_to_file(const std::string &filepath);
};

}
}

#endif // __TAXONOMIC_HPP__
