#ifndef __TAXONOMIC_CLASSIFIER_HPP__
#define __TAXONOMIC_CLASSIFIER_HPP__

#include "annotation/representation/base/annotation.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {


// TODO add some description.
class Classifier {
  public:
//    typedef std::string Label; // Maybe better to point this to the real structure.
    typedef std::uint64_t TaxoLabel;
    typedef annot::MultiLabelEncoded<std::string> Annotator;
    using row_index = Annotator::Index;
    using DeBruijnGraph = mtg::graph::DeBruijnGraph;


  private:
    TaxoLabel root_node;
    uint64_t num_nodes;
    std::vector<uint64_t> node_depth; // The root has the maximal depth;
    std::vector<uint64_t> node_parent;
    std::vector<std::string> index_to_label;
    tsl::hopscotch_map<row_index, TaxoLabel> taxonomic_map;

    void import_taxonomy(const std::string &filepath, std::vector<TaxoLabel> &linearization);

  public:
    Classifier(const std::string &filepath); // receive the file exported by Taxonomy obj.
    std::string assign_class(const mtg::graph::DeBruijnGraph &graph,
                       const std::string &sequence,
                       const double &lca_coverage_threshold);
};

}
}
#endif // __TAXONOMIC_CLASSIFIER_HPP__
