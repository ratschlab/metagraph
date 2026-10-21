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

TYPED_TEST(BinaryMatrixBRWTTest, RowViewsAcrossChunkBoundaries) {
    constexpr size_t num_rows = 131;
    constexpr size_t num_columns = 24;
    std::vector<BinaryMatrix::SetBitPositions> expected_rows(num_rows);
    std::vector<Vector<std::pair<BinaryMatrix::Column, uint64_t>>> expected_ranks(num_rows);
    BitVectorPtrArray columns;
    for (size_t col = 0; col < num_columns; ++col) {
        sdsl::bit_vector bits(num_rows, false);
        uint64_t rank = 0;
        for (size_t row = 0; row < num_rows; ++row) {
            if (row % 13 && (row + col * 7) % 11 < 5) {
                bits[row] = true;
                expected_rows[row].push_back(col);
                expected_ranks[row].emplace_back(col, ++rank);
            }
        }
        columns.emplace_back(new bit_vector_stat(std::move(bits)));
    }
    auto matrix = build_matrix_from_columns<TypeParam>(std::move(columns));

    // Non-owning slices start inside the buffer; IDs are unsorted and repeated.
    // Include empty rows and duplicate IDs on either side of chunk boundaries.
    std::vector<BinaryMatrix::Row> storage(20'005);
    for (size_t i = 0; i < storage.size(); ++i)
        storage[i] = i * 37 % num_rows;
    storage[10'000] = storage[10'001] = 64;
    storage[20'000] = storage[20'001] = 13;

    struct RestoreMode {
        bool previous = get_one_pass_brwt();
        ~RestoreMode() { set_one_pass_brwt(previous); }
    } restore_mode;

    for (bool one_pass : { false, true }) {
        set_one_pass_brwt(one_pass);
        for (size_t threads : { 1, 4 }) {
            for (size_t count : { 0, 9'999, 10'000, 10'001, 20'003 }) {
                SCOPED_TRACE(::testing::Message() << "one-pass=" << one_pass
                             << ", threads=" << threads << ", count=" << count);
                auto row_ids = std::span<const BinaryMatrix::Row>(storage).subspan(1, count);
                auto rows = matrix.get_rows(row_ids, threads);
                auto ranks = matrix.get_column_ranks(row_ids, threads);
                ASSERT_EQ(rows.size(), count);
                ASSERT_EQ(ranks.size(), count);
                for (size_t i = 0; i < count; ++i) {
                    std::sort(rows[i].begin(), rows[i].end());
                    std::sort(ranks[i].begin(), ranks[i].end());
                    ASSERT_EQ(rows[i], expected_rows[row_ids[i]]) << "query position " << i;
                    ASSERT_EQ(ranks[i], expected_ranks[row_ids[i]]) << "query position " << i;
                }
            }
        }
    }
}

} // namespace
