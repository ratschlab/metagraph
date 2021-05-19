#include "test_matrix_helpers.hpp"

#include <filesystem>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "common/vectors/bitmap_mergers.hpp"


namespace mtg {

using namespace mtg::annot::binmat;

namespace annot {
using CallColumn = std::function<void(std::unique_ptr<bit_vector>&&)>;
std::unique_ptr<Rainbow<BRWT>>
convert_to_RainbowBRWT(const std::function<void(const CallColumn &)> &call_columns,
                       size_t max_brwt_arity = 1);
} // namespace annot


namespace test {

typedef std::function<void(const BinaryMatrix::SetBitPositions &)> RowCallback;


template <typename BinMat>
BinMat
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t,
                       uint64_t) {
    return BinMat(generate_rows, num_columns);
}

template <>
RowConcatenated<>
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return RowConcatenated<>(generate_rows, num_columns, num_rows, num_relations);
}

template <>
RowSparse
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return RowSparse(generate_rows, num_columns, num_rows, num_relations);
}

template <>
BinRelWT
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t,
                       uint64_t num_relations) {
    return BinRelWT(generate_rows, num_relations, num_columns);
}

template <>
BinRelWT_sdsl
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t,
                       uint64_t num_relations) {
    return BinRelWT_sdsl(generate_rows, num_relations, num_columns);
}


template <typename BinMat>
BinMat build_matrix_from_columns(const BitVectorPtrArray &columns, uint64_t num_rows) {
    auto num_columns = columns.size();
    auto num_set_bits = 0;
    for (const auto &column : columns) {
        num_set_bits += column->num_set_bits();
    }

    return build_matrix_from_rows<BinMat>(
        [&](auto row_callback) {
            utils::RowsFromColumnsTransformer(columns).call_rows(row_callback);
        },
        num_columns, num_rows, num_set_bits
    );
}
template BinRelWT build_matrix_from_columns<BinRelWT>(const BitVectorPtrArray&, uint64_t);
template BinRelWT_sdsl build_matrix_from_columns<BinRelWT_sdsl>(const BitVectorPtrArray&, uint64_t);
template RowConcatenated<> build_matrix_from_columns<RowConcatenated<>>(const BitVectorPtrArray&, uint64_t);
template RowSparse build_matrix_from_columns<RowSparse>(const BitVectorPtrArray&, uint64_t);
template UniqueRowBinmat build_matrix_from_columns<UniqueRowBinmat>(const BitVectorPtrArray&, uint64_t);
template Rainbowfish build_matrix_from_columns<Rainbowfish>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<1> build_matrix_from_columns<RainbowfishBuffer<1>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<2> build_matrix_from_columns<RainbowfishBuffer<2>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<3> build_matrix_from_columns<RainbowfishBuffer<3>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<4> build_matrix_from_columns<RainbowfishBuffer<4>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<5> build_matrix_from_columns<RainbowfishBuffer<5>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<6> build_matrix_from_columns<RainbowfishBuffer<6>>(const BitVectorPtrArray&, uint64_t);
template <>
ColumnMajor build_matrix_from_columns<ColumnMajor>(const BitVectorPtrArray &columns, uint64_t) {
    BitVectorPtrArray columns_copy;
    for (const auto &col_ptr : columns) {
        columns_copy.emplace_back(new bit_vector_sd(col_ptr->copy_to<bit_vector_sd>()));
    }
    return ColumnMajor(std::move(columns_copy));
}
template <>
BRWT build_matrix_from_columns<BRWT>(const BitVectorPtrArray &columns, uint64_t) {
    BitVectorPtrArray columns_copy;
    for (const auto &col_ptr : columns) {
        columns_copy.push_back(col_ptr->copy());
    }
    BRWT matrix(BRWTBottomUpBuilder::build(std::move(columns_copy)));
    EXPECT_TRUE(matrix.avg_arity() <= 2) << matrix.avg_arity();
    return matrix;
}
template <>
BRWTOptimized build_matrix_from_columns<BRWTOptimized>(const BitVectorPtrArray &columns, uint64_t) {
    BitVectorPtrArray columns_copy;
    for (const auto &col_ptr : columns) {
        columns_copy.push_back(col_ptr->copy());
    }
    return BRWTOptimized(BRWTBottomUpBuilder::build(std::move(columns_copy)));
}
template <>
Rainbow<BRWT> build_matrix_from_columns<Rainbow<BRWT>>(const BitVectorPtrArray &columns, uint64_t) {
    return std::move(*annot::convert_to_RainbowBRWT([&](const auto &callback) {
        for (size_t j = 0; j < columns.size(); ++j) {
            callback(columns[j]->copy());
        }
    }));
}

template <typename TypeParam>
void test_serialization(const TypeParam &matrix) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_bad);

    {
        std::ofstream out(test_dump_basename_vec_good, std::ios::binary);
        matrix.serialize(out);
        out.close();
    }

    {
        TypeParam loaded;
        std::ifstream in(test_dump_basename_vec_bad, std::ios::binary);
        ASSERT_FALSE(loaded.load(in));
    }

    TypeParam loaded;
    {
        std::ifstream in(test_dump_basename_vec_good, std::ios::binary);
        ASSERT_TRUE(loaded.load(in));
        ASSERT_TRUE(in.good());
    }

    ASSERT_EQ(matrix.num_columns(), loaded.num_columns());
    ASSERT_EQ(matrix.num_rows(), loaded.num_rows());
    ASSERT_EQ(matrix.num_relations(), loaded.num_relations());
    for (size_t j = 0; j < loaded.num_columns(); ++j) {
        EXPECT_EQ(matrix.get_column(j), loaded.get_column(j));
    }

    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_bad);
}

template <typename TypeParam>
void test_matrix(const TypeParam &matrix, const BitVectorPtrArray &columns) {
    // check if the number of columns is the same
    ASSERT_EQ(columns.size(), matrix.num_columns());
    // if no columns, check if the number of rows is zero as well
    if (!columns.size()) {
        ASSERT_EQ(0u, matrix.num_rows());
        return;
    }

    // check if the number of rows is the same
    ASSERT_EQ(columns.at(0)->size(), matrix.num_rows());

    // check get_column
    for (size_t j = 0; j < matrix.num_columns(); ++j) {
        const auto &col = *columns[j];

        assert(col.size() == matrix.num_rows());

        auto col_set_bits = matrix.get_column(j);

        // make sure all returned indexes are unique
        ASSERT_EQ(col_set_bits.size(), convert_to_set(col_set_bits).size());

        ASSERT_EQ(col.num_set_bits(), col_set_bits.size())
            << "Column: " << j;

        for (auto i : col_set_bits) {
            ASSERT_TRUE(i < matrix.num_rows());
            ASSERT_TRUE(col[i]);
        }
    }

    // check call_columns
    for (size_t m : { size_t(0),
                      size_t(matrix.num_columns() / 2),
                      size_t(matrix.num_columns()) }) {
        std::vector<uint64_t> indices(m);
        std::iota(indices.begin(), indices.end(), 0);

        std::vector<std::vector<BinaryMatrix::Row>> column_map(m);
        matrix.call_columns(indices, [&](auto j, const bitmap &rows) {
            rows.call_ones([&](auto i) { column_map[j].push_back(i); });
        });
        std::vector<uint64_t> slice;
        for (size_t j = 0; j < column_map.size(); ++j) {
            slice.insert(slice.end(), column_map[j].begin(), column_map[j].end());
            slice.push_back(std::numeric_limits<BinaryMatrix::Row>::max());
        }

        ASSERT_GE(slice.size(), indices.size());

        auto column_begin = slice.begin();

        for (size_t j = 0; j < indices.size(); ++j) {
            // every row in `slice` ends with `-1`
            auto column_end = std::find(column_begin, slice.end(),
                                        std::numeric_limits<BinaryMatrix::Row>::max());
            std::vector<BinaryMatrix::Row> column_set_bits(column_begin, column_end);
            column_begin = column_end + 1;

            // make sure all returned indexes are unique
            ASSERT_EQ(column_set_bits.size(), convert_to_set(column_set_bits).size());

            for (auto i : column_set_bits) {
                ASSERT_TRUE(i < matrix.num_rows());
                EXPECT_TRUE((*columns[j])[i]);
            }

            auto set_bits = convert_to_set(column_set_bits);
            for (size_t i = 0; i < matrix.num_rows(); ++i) {
                EXPECT_EQ((*columns[j])[i], set_bits.count(i));
            }
        }
    }

    // check get_row
    for (size_t i = 0, n_rows = matrix.num_rows(); i < n_rows; ++i) {
        auto row_set_bits = matrix.get_row(i);

        // make sure all returned indexes are unique
        ASSERT_EQ(row_set_bits.size(), convert_to_set(row_set_bits).size());

        for (auto j : row_set_bits) {
            ASSERT_TRUE(j < matrix.num_columns());
            EXPECT_TRUE((*columns[j])[i]);
        }

        auto set_bits = convert_to_set(row_set_bits);
        for (size_t j = 0; j < columns.size(); ++j) {
            EXPECT_EQ((*columns[j])[i], set_bits.count(j));
        }
    }

    // check get_rows, query first |n| rows
    for (size_t n : { size_t(0),
                      size_t(matrix.num_rows() / 2),
                      size_t(matrix.num_rows()) }) {
        std::vector<uint64_t> indices(n);
        std::iota(indices.begin(), indices.end(), 0);

        auto rows = matrix.get_rows(indices);

        ASSERT_EQ(n, rows.size());

        for (size_t i = 0; i < rows.size(); ++i) {
            const auto &row_set_bits = rows[i];

            // make sure all returned indexes are unique
            ASSERT_EQ(row_set_bits.size(), convert_to_set(row_set_bits).size());

            for (auto j : row_set_bits) {
                ASSERT_TRUE(j < matrix.num_columns());
                EXPECT_TRUE((*columns[j])[i]);
            }

            auto set_bits = convert_to_set(row_set_bits);
            for (size_t j = 0; j < columns.size(); ++j) {
                EXPECT_EQ((*columns[j])[i], set_bits.count(j));
            }
        }
    }

    // check slice_rows, query first |n| rows
    for (size_t n : { size_t(0),
                      size_t(matrix.num_rows() / 2),
                      size_t(matrix.num_rows()) }) {
        std::vector<uint64_t> indices(n);
        std::iota(indices.begin(), indices.end(), 0);

        auto slice = matrix.slice_rows(indices);

        ASSERT_TRUE(slice.size() >= indices.size());

        auto row_begin = slice.begin();

        for (size_t i = 0; i < indices.size(); ++i) {
            // every row in `slice` ends with `-1`
            auto row_end = std::find(row_begin, slice.end(),
                                     std::numeric_limits<BinaryMatrix::Column>::max());
            std::vector<BinaryMatrix::Column> row_set_bits(row_begin, row_end);
            row_begin = row_end + 1;

            // make sure all returned indexes are unique
            ASSERT_EQ(row_set_bits.size(), convert_to_set(row_set_bits).size());

            for (auto j : row_set_bits) {
                ASSERT_TRUE(j < matrix.num_columns());
                EXPECT_TRUE((*columns[j])[i]);
            }

            auto set_bits = convert_to_set(row_set_bits);
            for (size_t j = 0; j < columns.size(); ++j) {
                EXPECT_EQ((*columns[j])[i], set_bits.count(j));
            }
        }
    }

    // check get
    for (size_t i = 0, n_rows = matrix.num_rows(); i < n_rows; ++i) {
        for (size_t j = 0; j < matrix.num_columns(); ++j) {
            EXPECT_EQ(columns[j]->operator[](i), matrix.get(i, j))
                << i << " " << j;
        }
    }

    // check serialization
    test_serialization(matrix);
}

template void test_matrix<BRWT>(const BRWT&, const BitVectorPtrArray &);
template void test_matrix<BRWTOptimized>(const BRWTOptimized&, const BitVectorPtrArray &);
template void test_matrix<ColumnMajor>(const ColumnMajor&, const BitVectorPtrArray &);
template void test_matrix<BinRelWT>(const BinRelWT&, const BitVectorPtrArray &);
template void test_matrix<BinRelWT_sdsl>(const BinRelWT_sdsl&, const BitVectorPtrArray &);
template void test_matrix<RowConcatenated<>>(const RowConcatenated<>&, const BitVectorPtrArray &);
template void test_matrix<RowSparse>(const RowSparse&, const BitVectorPtrArray &);
template void test_matrix<UniqueRowBinmat>(const UniqueRowBinmat&, const BitVectorPtrArray &);
template void test_matrix<Rainbow<BRWT>>(const Rainbow<BRWT>&, const BitVectorPtrArray &);
template void test_matrix<Rainbowfish>(const Rainbowfish&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<1>>(const RainbowfishBuffer<1>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<2>>(const RainbowfishBuffer<2>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<3>>(const RainbowfishBuffer<3>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<4>>(const RainbowfishBuffer<4>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<5>>(const RainbowfishBuffer<5>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<6>>(const RainbowfishBuffer<6>&, const BitVectorPtrArray &);

} // namespace test
} // namespace mtg
