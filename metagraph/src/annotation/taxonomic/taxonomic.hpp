#ifndef __TAXONOMIC_HPP__
#define __TAXONOMIC_HPP__

#include <cassert>
#include <tsl/hopscotch_map.h>


namespace mtg {
namespace annot {


// TODO add some description.
class Taxonomy {
  public:
    typedef std::string Label;
    typedef std::uint64_t TaxoLabel;
    typedef std::vector<std::vector<TaxoLabel>> ChildrenList;

  private:
    TaxoLabel root_node;

//   TODO what is this rmq. https://stackoverflow.com/questions/15488470/syntax-for-dynamically-allocating-a-2d-array-of-smart-pointers
    std::unique_ptr<TaxoLabel[]> *rmq_data;
    std::vector<uint64_t> node_depth; // The root has the maximal depth;
    std::vector<uint64_t> linerization_idx;
    tsl::hopscotch_map<Label, TaxoLabel> label_to_index;

    void read_tree(const std::string &tree_filepath, ChildrenList &tree);
    void calculate_node_depth(const TaxoLabel &node, const ChildrenList &tree);
    void rmq_preprocessing(const ChildrenList &tree);
    void dfs_linearization(
        const TaxoLabel &node,
        const ChildrenList &tree,
        std::vector<TaxoLabel> &tree_linearization
    );
  public:
    Taxonomy(const std::string &tree_filepath);
};

}
}

#endif // __TAXONOMIC_HPP__
