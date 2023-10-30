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

        return std::is_sorted(row_tuple.begin(), row_tuple.end());
    }));

    using Row = matrix::BinaryMatrix::Row;
    using Column = matrix::BinaryMatrix::Column;

    std::vector<std::vector<std::pair<Row, Vector<Column>>>> row_coords;
    row_coords.reserve(row_tuples.size());
    for (const auto &row_tuple : row_tuples) {
        VectorMap<Row, Vector<Column>> cur_row_coords;

        for (const auto &[c, tuple] : row_tuple) {
            static_assert(std::is_same_v<matrix::MultiIntMatrix::Tuple::value_type, uint64_t>);
            for (int64_t coord : tuple) {
                if (coord > 0)
                    cur_row_coords[coord - 1].emplace_back(c);
            }
        }

        row_coords.emplace_back(const_cast<std::vector<std::pair<Row, Vector<Column>>>&&>(
            cur_row_coords.values_container()
        ));
    }

    std::vector<std::vector<Vector<std::pair<Column, uint64_t>>>> matrix_ranks;
    matrix_ranks.reserve(seq_indexes_.size());
    for (const auto &seq_index : seq_indexes_) {
        auto &merged_cur_matrix_ranks = matrix_ranks.emplace_back();
        merged_cur_matrix_ranks.reserve(row_tuples.size());

        auto cur_matrix_ranks = seq_index->get_matrix().get_ranks(row_coords);
        assert(cur_matrix_ranks.size() == row_tuples.size());
        for (size_t i = 0; i < row_tuples.size(); ++i) {
            auto &cur_merged_row = merged_cur_matrix_ranks.emplace_back();
            auto it = cur_matrix_ranks[i].begin();
            assert(std::is_sorted(cur_matrix_ranks[i].begin(),
                                  cur_matrix_ranks[i].end()));
            for (const auto &[c, tuple] : row_tuples[i]) {
                static_assert(std::is_same_v<matrix::MultiIntMatrix::Tuple::value_type, uint64_t>);
                for (int64_t coord : tuple) {
                    if (coord > 0) {
                        assert(it != cur_matrix_ranks[i].end());
                        assert(c == it->first);
                        cur_merged_row.emplace_back(c, it->second);
                        ++it;
                    } else {
                        cur_merged_row.emplace_back(c, 0);
                    }
                }
            }
            assert(it == cur_matrix_ranks[i].end());
        }
    }

    std::vector<SeqIds> results;
    results.reserve(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        VectorMap<Column, Ids> rearrange;
        for (size_t j = 0; j < matrix_ranks.size(); ++j) {
            const auto &ranks = matrix_ranks[j][i];
            for (const auto &[c, r] : ranks) {
                auto &bucket = rearrange[c];
                bucket.resize(matrix_ranks.size());
                bucket[j].emplace_back(r);
            }
        }

        assert(rearrange.size() == row_tuples[i].size());
        assert(std::equal(rearrange.begin(), rearrange.end(),
                          row_tuples[i].begin(), row_tuples[i].end(),
                          [&](const auto &a, const auto &b) {
            return a.first == b.first;
        }));

        auto &result = results.emplace_back();
        result.reserve(rearrange.size());
        auto it = row_tuples[i].begin();
        for (const auto &[c, seq_ids] : rearrange) {
            assert(it->first == c);
            assert(std::all_of(seq_ids.begin(), seq_ids.end(), [&](const auto &s) {
                return s.size() == it->second.size();
            }));
            result.emplace_back(c, seq_ids);
            ++it;
        }
        assert(it == row_tuples[i].end());
    }

    return results;
}

template class SeqIndexedAnnotator<>;

} // namespace mtg::annot