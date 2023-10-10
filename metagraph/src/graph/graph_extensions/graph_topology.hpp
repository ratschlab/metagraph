#ifndef __GRAPH_TOPOLOGY_HPP__
#define __GRAPH_TOPOLOGY_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace graph {

class GraphTopology : public SequenceGraph::GraphExtension {
  public:
    using Annotator = AnnotatedDBG::Annotator;
    using Column = size_t;
    using Row = Annotator::Index;
    using Coord = std::pair<Column, uint64_t>;
    using TopologyIndex = std::pair<size_t, size_t>;
    using Coords = Vector<std::pair<Coord, TopologyIndex>>;

    GraphTopology(const graph::DeBruijnGraph &graph,
                  std::unique_ptr<Annotator>&& annotator,
                  std::unique_ptr<Annotator>&& unitigs,
                  std::unique_ptr<Annotator>&& clusters);

    std::vector<Coords> get_coords(const std::vector<graph::DeBruijnGraph::node_index> &nodes) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const graph::SequenceGraph &, bool = true) const { return true; }

    const Annotator& get_coord_annotator() const {
        assert(annotator_);
        return *annotator_;
    }

  private:
    const DeBruijnGraph &graph_;
    std::unique_ptr<Annotator> annotator_;
    std::unique_ptr<Annotator> unitig_annotator_;
    std::unique_ptr<Annotator> cluster_annotator_;
};

} // namespace graph
} // namespace mtg

#endif // __GRAPH_TOPOLOGY_HPP__
