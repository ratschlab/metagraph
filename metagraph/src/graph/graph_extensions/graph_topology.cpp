#include "graph_topology.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"

namespace mtg::graph {

auto GraphTopology::get_coords(const std::vector<annot::matrix::BinaryMatrix::Row> &rows)
        -> std::vector<Coords> {
    std::vector<Coords> result(rows.size());

    assert(coords_);
    assert(unitig_indicator_);
    assert(cluster_indicator_);

    {
        auto coords = coords_->get_row_tuples(rows);
        assert(coords.size() == rows.size());
        for (size_t i = 0; i < coords.size(); ++i) {
            const auto &row_tuples = coords[i];
            assert(result[i].empty() || result[i].size() == row_tuples.size());

            if (result[i].empty())
                result[i].resize(row_tuples.size());

            for (size_t j = 0; j < row_tuples.size(); ++j) {
                const auto &[c, tuple] = row_tuples[j];
                assert(tuple.size() == 1);
                auto &[rcol, ru, rsb, rcoord] = result[i][j];
                rcol = c;
                rcoord = tuple[0];
            }
        }
    }

    {
        auto unitigs = unitig_indicator_->get_column_ranks(rows);
        assert(unitigs.size() == rows.size());
        for (size_t i = 0; i < unitigs.size(); ++i) {
            const auto &row_unitigs = unitigs[i];
            assert(result[i].size() == row_unitigs.size());

            for (size_t j = 0; j < row_unitigs.size(); ++j) {
                const auto &[c, unitig_id] = row_unitigs[j];
                auto &[rcol, ru, rsb, rcoord] = result[i][j];
                assert(rcol == c);
                ru = unitig_id;
            }
        }
    }

    {
        auto clusters = cluster_indicator_->get_column_ranks(rows);
        assert(clusters.size() == rows.size());
        for (size_t i = 0; i < clusters.size(); ++i) {
            const auto &row_clusters = clusters[i];
            assert(result[i].size() == row_clusters.size());

            for (size_t j = 0; j < row_clusters.size(); ++j) {
                const auto &[c, cluster_id] = row_clusters[j];
                auto &[rcol, ru, rsb, rcoord] = result[i][j];
                assert(rcol == c);
                rsb = cluster_id;
            }
        }
    }

    return result;
}

} // namespace mtg::graph