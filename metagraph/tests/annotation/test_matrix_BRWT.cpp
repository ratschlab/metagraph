#include <random>

#include "gtest/gtest.h"

#include "test_matrix_helpers.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;
using namespace mtg::annot::binmat;

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

} // namespace
