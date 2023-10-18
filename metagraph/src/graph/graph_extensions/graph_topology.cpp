#include "graph_topology.hpp"

#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
#include "common/logger.hpp"
#include "graph/alignment/annotation_buffer.hpp"
#include "graph/representation/succinct/boss.hpp"

namespace mtg::graph {

using namespace mtg::annot::matrix;
using node_index = DeBruijnGraph::node_index;
using Tuple = MultiIntMatrix::Tuple;
using RowTuples = MultiIntMatrix::RowTuples;

using RowColumnFlatMap = std::vector<std::pair<GraphTopology::Row,
                                               Vector<GraphTopology::Column>>>;

static_assert(std::is_same_v<GraphTopology::Column, BinaryMatrix::Column>);
static_assert(std::is_same_v<GraphTopology::Row, BinaryMatrix::Row>);

GraphTopology::GraphTopology(const graph::DeBruijnGraph &graph,
                             std::shared_ptr<const annot::SeqIndexedAnnotator<Label>> annotator)
      : graph_(graph), annotator_(annotator),
        buffer_(std::make_shared<align::AnnotationBuffer>(graph_, *annotator_, false)) {}

GraphTopology::GraphTopology(const graph::DeBruijnGraph &graph,
                             std::shared_ptr<align::AnnotationBuffer> buffer)
      : graph_(graph),
        annotator_(std::shared_ptr<const annot::SeqIndexedAnnotator<Label>>{},
                   dynamic_cast<const annot::SeqIndexedAnnotator<Label>*>(&buffer_->get_annotator())),
        buffer_(buffer) {
    if (!annotator_) {
        common::logger->error("Annotator does not support graph topology.");
        exit(1);
    }
}

auto GraphTopology::get_coords(const std::vector<node_index> &nodes) const
        -> std::vector<Coords> {
    assert(annotator_);
    buffer_->queue_path(std::vector<node_index>(nodes));
    buffer_->fetch_queued_annotations();

    std::vector<RowTuples> tuples;
    tuples.reserve(nodes.size());
    for (node_index node : nodes) {
        auto &row_tuple = tuples.emplace_back();
        const auto [cols, coords] = buffer_->get_labels_and_coords(node);
        assert(cols);
        assert(coords);
        assert(cols->size());
        assert(cols->size() == coords->size());

        for (size_t i = 0; i < cols->size(); ++i) {
            if ((*coords)[i].size())
                row_tuple.emplace_back((*cols)[i], Tuple((*coords)[i].begin(), (*coords)[i].end()));
        }
    }
    assert(tuples.size() == nodes.size());

    std::vector<Coords> result(tuples.size());
    auto row_seq_ids = annotator_->get_seq_ids(tuples);
    assert(row_seq_ids.size() == nodes.size());

    for (size_t i = 0; i < row_seq_ids.size(); ++i) {
        assert(row_seq_ids[i].size() == tuples[i].size());
        for (size_t j = 0; j < row_seq_ids[i].size(); ++j) {
            auto &[c, seq_ids] = row_seq_ids[i][j];

            assert(c == tuples[i][j].first);
            auto &coords = tuples[i][j].second;

            assert(seq_ids.size() == 2);
            assert(seq_ids[0].size() == seq_ids[1].size());
            assert(coords.size() == seq_ids[0].size());

            for (size_t k = 0; k < seq_ids[0].size(); ++k) {
                result[i].emplace_back(Coord(c, coords[k]),
                                       TopologyIndex(seq_ids[0][k], seq_ids[1][k]));
            }
        }
    }

    return result;
}

bool GraphTopology::load(const std::string &filename_base) {
    std::ignore = filename_base;
    // {
    //     std::string fname = filename_base + kUnitigExtension + unitig_annotator_->file_extension();
    //     if (!unitig_annotator_->load(fname)) {
    //         common::logger->error("Failed to load unitig indicator from {}", fname);
    //         return false;
    //     }
    // }

    // {
    //     std::string fname = filename_base + kClusterExtension + cluster_annotator_->file_extension();
    //     if (!cluster_annotator_->load(fname)) {
    //         common::logger->error("Failed to load cluster indicator from {}", fname);
    //         return false;
    //     }
    // }

    // assert(annotator_->num_labels() == unitig_annotator_->num_labels());
    // assert(annotator_->num_labels() == cluster_annotator_->num_labels());

    return true;
}

void GraphTopology::serialize(const std::string &filename_base) const {
    std::ignore = filename_base;
    // unitig_annotator_->serialize(filename_base + kUnitigExtension + unitig_annotator_->file_extension());
    // cluster_annotator_->serialize(filename_base + kClusterExtension + cluster_annotator_->file_extension());
}

} // namespace mtg::graph