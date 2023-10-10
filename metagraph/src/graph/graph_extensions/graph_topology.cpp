#include "graph_topology.hpp"

#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
#include "common/logger.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "graph/alignment/annotation_buffer.hpp"

namespace mtg::graph {

using namespace mtg::annot::matrix;
using node_index = DeBruijnGraph::node_index;
using Tuple = MultiIntMatrix::Tuple;
using RowTuples = Vector<std::pair<GraphTopology::Column,
                                   align::Alignment::Tuple>>;

using RowColumnFlatMap = std::vector<std::pair<GraphTopology::Row,
                                               Vector<GraphTopology::Column>>>;

static_assert(std::is_same_v<GraphTopology::Column, BinaryMatrix::Column>);
static_assert(std::is_same_v<GraphTopology::Row, BinaryMatrix::Row>);

GraphTopology::GraphTopology(const graph::DeBruijnGraph &graph,
                             std::shared_ptr<const Annotator> annotator,
                             std::unique_ptr<Annotator>&& unitigs,
                             std::unique_ptr<Annotator>&& clusters)
      : graph_(graph), annotator_(annotator),
        buffer_(std::make_shared<align::AnnotationBuffer>(graph_, *annotator_)),
        unitig_annotator_(std::move(unitigs)), cluster_annotator_(std::move(clusters)) {
    assert(annotator_);
    assert(unitig_annotator_);
    assert(cluster_annotator_);
    assert(annotator_->num_labels() == unitig_annotator_->num_labels());
    assert(annotator_->num_labels() == cluster_annotator_->num_labels());
}

GraphTopology::GraphTopology(const graph::DeBruijnGraph &graph,
                             std::shared_ptr<align::AnnotationBuffer> buffer,
                             std::unique_ptr<Annotator>&& unitigs,
                             std::unique_ptr<Annotator>&& clusters)
      : graph_(graph),
        annotator_(std::shared_ptr<const Annotator>{}, &buffer_->get_annotator()),
        buffer_(buffer), unitig_annotator_(std::move(unitigs)),
        cluster_annotator_(std::move(clusters)) {
    assert(annotator_);
    assert(unitig_annotator_);
    assert(cluster_annotator_);
    assert(annotator_->num_labels() == unitig_annotator_->num_labels());
    assert(annotator_->num_labels() == cluster_annotator_->num_labels());
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
            row_tuple.emplace_back(std::move((*cols)[i]), std::move((*coords)[i]));
        }
    }

    assert(unitig_annotator_);
    assert(cluster_annotator_);

    std::vector<Coords> result(tuples.size());

    std::vector<RowColumnFlatMap> coords;
    coords.reserve(tuples.size());
    for (size_t i = 0; i < tuples.size(); ++i) {
        const auto &row_tuples = tuples[i];
        VectorMap<Row, Vector<Column>> row_coords;

        for (const auto &[c, tuple] : row_tuples) {
            assert(tuple.size() <= 2);
            for (auto coord : tuple) {
                row_coords[coord].emplace_back(c);
                result[i].emplace_back(Coord(c, coord), TopologyIndex());
            }
        }

        assert(row_coords.size());
        coords.emplace_back(const_cast<RowColumnFlatMap&&>(row_coords.values_container()));
    }
    assert(coords.size() == tuples.size());

    {
        auto unitigs = unitig_annotator_->get_matrix().get_ranks(coords);
        assert(unitigs.size() == nodes.size());
        for (size_t i = 0; i < unitigs.size(); ++i) {
            const auto &row_unitigs = unitigs[i];
            assert(row_unitigs.size());
            assert(row_unitigs.size() == result[i].size());
            assert(result[i].size() == row_unitigs.size());

            for (size_t j = 0; j < row_unitigs.size(); ++j) {
                const auto &[c, unitig_id] = row_unitigs[j];
                auto &[coord, topo_index] = result[i][j];
                assert(coord.first == c);
                topo_index.first = unitig_id;
            }
        }
    }

    {
        auto clusters = cluster_annotator_->get_matrix().get_ranks(coords);
        assert(clusters.size() == nodes.size());
        for (size_t i = 0; i < clusters.size(); ++i) {
            const auto &row_clusters = clusters[i];
            assert(row_clusters.size());
            assert(row_clusters.size() == result[i].size());
            assert(result[i].size() == row_clusters.size());

            for (size_t j = 0; j < row_clusters.size(); ++j) {
                const auto &[c, cluster_id] = row_clusters[j];
                auto &[coord, topo_index] = result[i][j];
                assert(coord.first == c);
                topo_index.second = cluster_id;
            }
        }
    }

    assert(result.size() == nodes.size());
    return result;
}

bool GraphTopology::load(const std::string &filename_base) {
    {
        std::string fname = filename_base + unitig_annotator_->file_extension();
        if (!unitig_annotator_->load(fname)) {
            common::logger->error("Failed to load unitig indicator from {}", fname);
            return false;
        }
    }

    {
        std::string fname = filename_base + cluster_annotator_->file_extension();
        if (!cluster_annotator_->load(fname)) {
            common::logger->error("Failed to load cluster indicator from {}", fname);
            return false;
        }
    }

    assert(annotator_->num_labels() == unitig_annotator_->num_labels());
    assert(annotator_->num_labels() == cluster_annotator_->num_labels());

    return true;
}

void GraphTopology::serialize(const std::string &filename_base) const {
    unitig_annotator_->serialize(filename_base + unitig_annotator_->file_extension());
    cluster_annotator_->serialize(filename_base + cluster_annotator_->file_extension());
}

} // namespace mtg::graph