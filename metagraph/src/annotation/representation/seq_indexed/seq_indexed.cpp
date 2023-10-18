#include "seq_indexed.hpp"

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg::annot {

using namespace mtg::annot::matrix;

template <typename Label>
auto SeqIndexedAnnotator<Label>
::get_seq_ids(const std::vector<RowTuples> &row_tuples) const -> std::vector<SeqIds> {
    assert(std::all_of(row_tuples.begin(), row_tuples.end(), [](const auto &row_tuple) {
        for (const auto &[c, tuple] : row_tuple) {
            if (tuple.empty())
                return false;
        }

        return true;
    }));

    using Row = matrix::BinaryMatrix::Row;
    using Column = matrix::BinaryMatrix::Column;

    std::vector<std::vector<std::pair<Row, Vector<Column>>>> row_coords;
    row_coords.reserve(row_tuples.size());
    for (const auto &row_tuple : row_tuples) {
        using RowColumnFlatMap = std::vector<std::pair<Row, Vector<Column>>>;
        VectorMap<Row, Vector<Column>> cur_row_coords;

        for (const auto &[c, tuple] : row_tuple) {
            for (auto coord : tuple) {
                cur_row_coords[coord].emplace_back(c);
            }
        }

        row_coords.emplace_back(const_cast<RowColumnFlatMap&&>(
            cur_row_coords.values_container()
        ));
    }

    std::vector<std::vector<Vector<std::pair<Column, uint64_t>>>> matrix_ranks;
    matrix_ranks.reserve(seq_indexes_.size());
    for (const auto &seq_index : seq_indexes_) {
        matrix_ranks.emplace_back(seq_index->get_matrix().get_ranks(row_coords));
        assert(matrix_ranks.back().size() == row_tuples.size());
    }

    std::vector<SeqIds> results;
    results.reserve(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        tsl::hopscotch_map<Column, Ids> rearrange;
        for (size_t j = 0; j < matrix_ranks.size(); ++j) {
            const auto &ranks = matrix_ranks[j][i];
            for (const auto &[c, r] : ranks) {
                auto &bucket = rearrange[c];
                bucket.resize(matrix_ranks.size());
                bucket[j].emplace_back(r);
            }
        }

        assert(rearrange.size() == row_tuples[i].size());

        auto &result = results.emplace_back();
        result.reserve(rearrange.size());
        for (const auto &[c, seq_ids] : rearrange) {
            result.emplace_back(c, seq_ids);
        }
    }

    return results;
}

template class SeqIndexedAnnotator<>;

} // namespace mtg::annot