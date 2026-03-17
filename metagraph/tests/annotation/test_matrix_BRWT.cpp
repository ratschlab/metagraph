#include <random>

#include "gtest/gtest.h"

#include "test_matrix_helpers.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "common/vectors/bit_vector_sd.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;
using namespace mtg::annot::matrix;

template <typename BinMat>
class BinaryMatrixBRWTTest : public ::testing::Test { };
typedef ::testing::Types<BRWT, BRWTOptimized> BinMatColumnTypes;
TYPED_TEST_SUITE(BinaryMatrixBRWTTest, BinMatColumnTypes);

TYPED_TEST(BinaryMatrixBRWTTest, ArityEmpty) {
    TypeParam matrix;
    // empty root
    EXPECT_EQ(1u, matrix.num_nodes());
    ASSERT_EQ(0u, matrix.avg_arity());
}

TYPED_TEST(BinaryMatrixBRWTTest, ArityOneCol) {
    {
        BitVectorPtrArray columns;
        columns.emplace_back(new bit_vector_stat(10, true));

        auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
    {
        BitVectorPtrArray columns;
        columns.emplace_back(new bit_vector_stat(10, false));

        auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
}

TYPED_TEST(BinaryMatrixBRWTTest, ArityTwoCol) {
    {
        BitVectorPtrArray columns;
        columns.emplace_back(new bit_vector_stat(10, true));
        columns.emplace_back(new bit_vector_stat(10, false));

        auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
    {
        BitVectorPtrArray columns;
        columns.emplace_back(new bit_vector_stat(10, false));
        columns.emplace_back(new bit_vector_stat(10, true));

        auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
}

TYPED_TEST(BinaryMatrixBRWTTest, GetRowsManyThreadsPreservesMembershipWithDuplicates) {
    BitVectorPtrArray columns;
    columns.emplace_back(new bit_vector_sd({ 1, 0, 0, 0, 1, 0, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 1, 0, 0, 1, 0, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 1, 0, 0, 1, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 1, 0, 1, 0, 0, 0, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 0, 0, 1, 0, 1, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 0, 0, 0, 0, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 1, 1, 1, 0, 0, 0, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 0, 0, 0, 1, 1, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 0, 1, 0, 0, 0, 1 }));
    columns.emplace_back(new bit_vector_sd({ 0, 0, 1, 0, 0, 0, 1, 0 }));
    columns.emplace_back(new bit_vector_sd({ 1, 0, 0, 0, 0, 1, 0, 0 }));
    columns.emplace_back(new bit_vector_sd({ 0, 1, 0, 0, 0, 0, 1, 0 }));

    auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));
    std::vector<BinaryMatrix::Row> row_ids = { 0, 2, 2, 3, 5, 7, 0 };

    auto single_thread_rows = matrix.get_rows(row_ids);
    auto many_thread_rows = matrix.get_rows(row_ids, 4);

    ASSERT_EQ(single_thread_rows.size(), many_thread_rows.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        std::sort(single_thread_rows[i].begin(), single_thread_rows[i].end());
        std::sort(many_thread_rows[i].begin(), many_thread_rows[i].end());
        EXPECT_EQ(single_thread_rows[i], many_thread_rows[i]) << "row index " << row_ids[i];
    }
}

} // namespace
