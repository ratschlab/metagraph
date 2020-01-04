#include "test_matrix_helpers.hpp"

#include <filesystem>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/binary_matrix/multi_brwt/BRWT.hpp"
#include "annotation/binary_matrix/multi_brwt/BRWT_builders.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "common/vectors/bitmap_mergers.hpp"

typedef BinaryMatrix::SetBitPositions RowSetBits;
typedef std::function<void(const RowSetBits &)> RowCallback;


// Used to generate a set of columns from a row generator
template <typename BinMat>
BinMat
build_matrix_from_columns(const std::function<void(const RowCallback &)> &generate_rows,
                          uint64_t num_columns,
                          uint64_t num_rows,
                          uint64_t) {
    BitVectorPtrArray columns;
    while (num_columns--) {
        columns.emplace_back(new bit_vector_stat(num_rows, false));
    }

    uint64_t cur_row = 0;
    generate_rows([&](auto row) {
        for (const auto &column : row) {
            columns.at(column)->set(cur_row, true);
        }
        cur_row++;
    });

    return build_matrix_from_columns<BinMat>(std::move(columns), num_rows);
}

template <typename BinMat>
BinMat build_matrix_from_columns(BitVectorPtrArray&& columns, uint64_t) {
    return BinMat(std::move(columns));
}

template BinRelWT build_matrix_from_columns<BinRelWT>(BitVectorPtrArray&&, uint64_t);

template <>
Rainbowfish build_matrix_from_columns(BitVectorPtrArray&& columns, uint64_t num_rows) {
    return build_matrix_from_rows<Rainbowfish>(std::move(columns), num_rows);
}

#define RBFBufferCol(n) \
template <> \
RainbowfishBuffer<n> build_matrix_from_columns(BitVectorPtrArray&& columns, uint64_t num_rows) { \
    return build_matrix_from_rows<RainbowfishBuffer<n>>(std::move(columns), num_rows); \
}
RBFBufferCol(1)
RBFBufferCol(2)
RBFBufferCol(3)
RBFBufferCol(4)
RBFBufferCol(5)
RBFBufferCol(6)

template <>
RowConcatenated<> build_matrix_from_columns(BitVectorPtrArray&& columns, uint64_t num_rows) {
    return build_matrix_from_rows<RowConcatenated<>>(std::move(columns), num_rows);
}

template <>
BinRelWT_sdsl build_matrix_from_columns(BitVectorPtrArray&& columns, uint64_t) {
    return build_matrix_from_rows<BinRelWT_sdsl>(std::move(columns));
}

template <>
BRWT build_matrix_from_columns<BRWT>(BitVectorPtrArray&& columns, uint64_t) {
    BRWT matrix(BRWTBottomUpBuilder::build(std::move(columns)));
    EXPECT_TRUE(matrix.avg_arity() <= 2) << matrix.avg_arity();
    return matrix;
}

template <>
BRWTOptimized build_matrix_from_columns<BRWTOptimized>(BitVectorPtrArray&& columns, uint64_t) {
    return BRWTOptimized(BRWTBottomUpBuilder::build(std::move(columns)));
}

template <>
ColumnMajor build_matrix_from_columns<ColumnMajor>(BitVectorPtrArray&& columns, uint64_t) {
    std::vector<std::unique_ptr<bit_vector>> columns_sd;
    columns_sd.reserve(columns.size());
    for (auto&& column : columns) {
        columns_sd.emplace_back(new bit_vector_sd(column->convert_to<bit_vector_sd>()));
    }
    return ColumnMajor(std::move(columns_sd));
}


template <typename BinMat>
BinMat build_matrix_from_columns(const BitVectorPtrArray &columns, uint64_t) {
    return BinMat(columns);
}

#define RBFBufferColConst(n) \
template <> \
RainbowfishBuffer<n> build_matrix_from_columns(const BitVectorPtrArray &columns, uint64_t num_rows) { \
    return build_matrix_from_rows<RainbowfishBuffer<n>>(columns, num_rows); \
}
RBFBufferColConst(1)
RBFBufferColConst(2)
RBFBufferColConst(3)
RBFBufferColConst(4)
RBFBufferColConst(5)
RBFBufferColConst(6)



template <typename BinMat>
BinMat
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return BinMat(generate_rows, num_columns, num_rows, num_relations);
}

template RowConcatenated<> build_matrix_from_rows<RowConcatenated<>>(const std::function<void(const RowCallback &)> &, uint64_t, uint64_t, uint64_t);

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

template <>
Rainbowfish
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t,
                       uint64_t) {
    return Rainbowfish(generate_rows, num_columns);
}

#define RBFBufferRow(n) \
template <> \
RainbowfishBuffer<n> \
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows, \
                       uint64_t num_columns, \
                       uint64_t, \
                       uint64_t) { \
    return RainbowfishBuffer<n>(generate_rows, num_columns); \
}
RBFBufferRow(1)
RBFBufferRow(2)
RBFBufferRow(3)
RBFBufferRow(4)
RBFBufferRow(5)
RBFBufferRow(6)

template <>
BRWT
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return build_matrix_from_columns<BRWT>(generate_rows,
                                           num_columns,
                                           num_rows,
                                           num_relations);
}

template <>
BRWTOptimized
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return build_matrix_from_columns<BRWTOptimized>(generate_rows,
                                                    num_columns,
                                                    num_rows,
                                                    num_relations);
}

template <>
ColumnMajor
build_matrix_from_rows(const std::function<void(const RowCallback &)> &generate_rows,
                       uint64_t num_columns,
                       uint64_t num_rows,
                       uint64_t num_relations) {
    return build_matrix_from_columns<ColumnMajor>(generate_rows,
                                                  num_columns,
                                                  num_rows,
                                                  num_relations);
}

template <typename BinMat>
BinMat build_matrix_from_rows(BitVectorPtrArray&& columns, uint64_t num_rows) {
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
template BRWT build_matrix_from_rows<BRWT>(BitVectorPtrArray&&, uint64_t);
template BRWTOptimized build_matrix_from_rows<BRWTOptimized>(BitVectorPtrArray&&, uint64_t);
template ColumnMajor build_matrix_from_rows<ColumnMajor>(BitVectorPtrArray&&, uint64_t);
template BinRelWT build_matrix_from_rows<BinRelWT>(BitVectorPtrArray&&, uint64_t);
template BinRelWT_sdsl build_matrix_from_rows<BinRelWT_sdsl>(BitVectorPtrArray&&, uint64_t);
template RowConcatenated<> build_matrix_from_rows<RowConcatenated<>>(BitVectorPtrArray&&, uint64_t);
template Rainbowfish build_matrix_from_rows<Rainbowfish>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<1> build_matrix_from_rows<RainbowfishBuffer<1>>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<2> build_matrix_from_rows<RainbowfishBuffer<2>>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<3> build_matrix_from_rows<RainbowfishBuffer<3>>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<4> build_matrix_from_rows<RainbowfishBuffer<4>>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<5> build_matrix_from_rows<RainbowfishBuffer<5>>(BitVectorPtrArray&&, uint64_t);
template RainbowfishBuffer<6> build_matrix_from_rows<RainbowfishBuffer<6>>(BitVectorPtrArray&&, uint64_t);

template <typename BinMat>
BinMat build_matrix_from_rows(const BitVectorPtrArray &columns, uint64_t num_rows) {
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
template Rainbowfish build_matrix_from_rows<Rainbowfish>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<1> build_matrix_from_rows<RainbowfishBuffer<1>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<2> build_matrix_from_rows<RainbowfishBuffer<2>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<3> build_matrix_from_rows<RainbowfishBuffer<3>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<4> build_matrix_from_rows<RainbowfishBuffer<4>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<5> build_matrix_from_rows<RainbowfishBuffer<5>>(const BitVectorPtrArray&, uint64_t);
template RainbowfishBuffer<6> build_matrix_from_rows<RainbowfishBuffer<6>>(const BitVectorPtrArray&, uint64_t);


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
    for (size_t j = 0; j < loaded.num_columns(); ++j) {
        EXPECT_EQ(matrix.get_column(j), loaded.get_column(j));
    }
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

    // check get_row
    for (size_t i = 0; i < matrix.num_rows(); ++i) {
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

    // check get
    for (size_t i = 0; i < matrix.num_rows(); ++i) {
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
template void test_matrix<Rainbowfish>(const Rainbowfish&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<1>>(const RainbowfishBuffer<1>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<2>>(const RainbowfishBuffer<2>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<3>>(const RainbowfishBuffer<3>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<4>>(const RainbowfishBuffer<4>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<5>>(const RainbowfishBuffer<5>&, const BitVectorPtrArray &);
template void test_matrix<RainbowfishBuffer<6>>(const RainbowfishBuffer<6>&, const BitVectorPtrArray &);
