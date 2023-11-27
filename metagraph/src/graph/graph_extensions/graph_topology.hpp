#ifndef __GRAPH_TOPOLOGY_HPP__
#define __GRAPH_TOPOLOGY_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/seq_indexed/seq_indexed.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

namespace mtg::graph {

namespace align {
    class AnnotationBuffer;
}

class GraphSeqAnnotator : public SequenceGraph::GraphExtension {
  public:
    using Annotator = AnnotatedDBG::Annotator;
    using Label = Annotator::Label;
    using Column = annot::matrix::BinaryMatrix::Column;
    using RowTuples = annot::matrix::MultiIntMatrix::RowTuples;
    using SeqIds = typename annot::SeqIndexedAnnotator<Label>::SeqIds;

    GraphSeqAnnotator(const graph::DeBruijnGraph &graph,
                      std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator = {});

    GraphSeqAnnotator(const graph::DeBruijnGraph &graph,
                      std::shared_ptr<align::AnnotationBuffer> buffer);

    bool load(const std::string &) {
        throw std::runtime_error("Load not implemented");
    }

    void serialize(const std::string &) const {
        throw std::runtime_error("Serialize not implemented");
    }

    bool is_compatible(const graph::SequenceGraph &, bool = true) const { return true; }

    std::pair<std::vector<SeqIds>, std::vector<RowTuples>>
    get_seq_ids(const std::vector<graph::DeBruijnGraph::node_index> &nodes) const;

    const Annotator& get_annotator() const {
        assert(annotator_);
        return *annotator_;
    }

  protected:
    const DeBruijnGraph &graph_;
    std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator_;
    std::shared_ptr<align::AnnotationBuffer> buffer_;
};

class GraphTopology : public GraphSeqAnnotator {
  public:
    using Row = Annotator::Index;
    using Coord = std::pair<Column, uint64_t>;
    using TopologyIndex = std::pair<size_t, size_t>;
    using Coords = Vector<std::pair<Coord, TopologyIndex>>;

    GraphTopology(const graph::DeBruijnGraph &graph,
                  std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator = {});

    GraphTopology(const graph::DeBruijnGraph &graph,
                  std::shared_ptr<align::AnnotationBuffer> buffer);

    std::vector<Coords> get_coords(const std::vector<graph::DeBruijnGraph::node_index> &nodes) const;

    static std::string seq_extension() { return kSeqExtension; }
    static std::string cluster_extension() { return kClusterExtension; }

  private:
    static constexpr auto kSeqExtension = ".seq";
    static constexpr auto kClusterExtension = ".clusters";
};

} // namespace mtg::graph

#endif // __GRAPH_TOPOLOGY_HPP__
