#ifndef __GRAPH_TOPOLOGY_HPP__
#define __GRAPH_TOPOLOGY_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/seq_indexed/seq_indexed.hpp"

namespace mtg::graph {

namespace align {
    class AnnotationBuffer;
}

class GraphTopology : public SequenceGraph::GraphExtension {
  public:
    using Annotator = AnnotatedDBG::Annotator;
    using Label = Annotator::Label;
    using Column = size_t;
    using Row = Annotator::Index;
    using Coord = std::pair<Column, uint64_t>;
    using TopologyIndex = std::pair<size_t, size_t>;
    using Coords = Vector<std::pair<Coord, TopologyIndex>>;

    GraphTopology(const graph::DeBruijnGraph &graph,
                  std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator = {});

    GraphTopology(const graph::DeBruijnGraph &graph,
                  std::shared_ptr<align::AnnotationBuffer> buffer);

    std::vector<Coords> get_coords(const std::vector<graph::DeBruijnGraph::node_index> &nodes) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const graph::SequenceGraph &, bool = true) const { return true; }

    const Annotator& get_annotator() const {
        assert(annotator_);
        return *annotator_;
    }

    static std::string cluster_extension() { return kClusterExtension; }

  private:
    const DeBruijnGraph &graph_;
    std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator_;
    std::shared_ptr<align::AnnotationBuffer> buffer_;

    static constexpr auto kClusterExtension = ".clusters";
};

} // namespace mtg::graph

#endif // __GRAPH_TOPOLOGY_HPP__
