#include "seq_indexed.hpp"

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg::annot {

using namespace mtg::annot::matrix;

template <class BinaryMatrix, typename Label>
auto SeqIndexedAnnotator<BinaryMatrix, Label>
::get_seq_ids(const std::vector<RowTuples> &row_tuples) const -> std::vector<SeqIds> {
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
        matrix_ranks.emplace_back(seq_index->get_ranks(row_coords));
    }

    std::vector<SeqIds> results;
    results.reserve(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        const auto &tuple = row_tuples[i];
        auto &result = results.emplace_back();
        result.reserve(tuple.size());
        for (size_t j = 0; j < tuple.size(); ++i) {
            const auto &[c, coords] = tuple[j];
            auto &seq_ids = result.emplace_back(c, Ids{}).second;
            seq_ids.reserve(matrix_ranks.size());
            for (const auto &ranks : matrix_ranks) {
                assert(ranks.size() == row_tuples.size());
                assert(ranks[i].size() == tuple.size());
                assert(ranks[i][j].first == c);
                seq_ids.emplace_back(ranks[i][j].second);
            }
        }
    }

    return results;
}

template class SeqIndexedAnnotator<TupleRowDiff<CoordRowDisk>, std::string>;

template class SeqIndexedAnnotator<TupleCSCMatrix<ColumnMajor>, std::string>;
template class SeqIndexedAnnotator<TupleCSCMatrix<BRWT>, std::string>;

template class SeqIndexedAnnotator<TupleRowDiff<TupleCSCMatrix<ColumnMajor>>, std::string>;
template class SeqIndexedAnnotator<TupleRowDiff<TupleCSCMatrix<BRWT>>, std::string>;

} // namespace mtg::annot