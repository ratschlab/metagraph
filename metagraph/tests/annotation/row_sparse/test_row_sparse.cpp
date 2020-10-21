#include <cstdint>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <graph/representation/succinct/dbg_succinct.hpp>

#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "common/utils/file_utils.hpp"

namespace {
using namespace mtg;
using namespace testing;
using ::testing::_;

TEST(RowSparse, Default) {
    annot::binmat::RowSparse rowdiff;
    EXPECT_EQ(0, rowdiff.num_relations());
    EXPECT_EQ(0, rowdiff.num_columns());
    EXPECT_EQ(0, rowdiff.num_rows());
}

TEST(RowSparse, InitEmpty) {
    Vector<Vector<uint64_t>> set_bits = { {}, {}, {}, {} };
    using RowCallback = annot::binmat::BinaryMatrix::RowCallback;
    annot::binmat::RowSparse under_test([&](const RowCallback& callback) {
        for (auto& row : set_bits) {
            callback(row);
        }
    }, 4, 4, 0);
    EXPECT_EQ(0, under_test.num_relations());
    EXPECT_EQ(4, under_test.num_columns());
    EXPECT_EQ(4, under_test.num_rows());
}

TEST(RowSparse, InitRows) {
    Vector<Vector<uint64_t>> set_bits = { {0}, {1}, {2}, {3} };
    using RowCallback = annot::binmat::BinaryMatrix::RowCallback;
    annot::binmat::RowSparse under_test([&](const RowCallback& callback) {
      for (auto& row : set_bits) {
          callback(row);
      }
    }, 4, 4, 4);
    EXPECT_EQ(4, under_test.num_relations());
    EXPECT_EQ(4, under_test.num_columns());
    EXPECT_EQ(4, under_test.num_rows());
    for (uint32_t i = 0; i < 4; ++i) {
        for(uint32_t j = 0; j < 4; ++j) {
            ASSERT_EQ(under_test.get(i,j), i == j);
        }
    }
}

TEST(RowSparse, Move) {
    Vector<Vector<uint64_t>> set_bits = { {0}, {1}, {2}, {3} };
    using RowCallback = annot::binmat::BinaryMatrix::RowCallback;
    annot::binmat::RowSparse under_test([&](const RowCallback& callback) {
      for (auto& row : set_bits) {
          callback(row);
      }
    }, 4, 4, 4);
    annot::binmat::RowSparse other = std::move(under_test);
    EXPECT_EQ(4, under_test.num_relations());
    EXPECT_EQ(4, under_test.num_columns());
    EXPECT_EQ(4, under_test.num_rows());
    for (uint32_t i = 0; i < 4; ++i) {
        for(uint32_t j = 0; j < 4; ++j) {
            ASSERT_EQ(under_test.get(i,j), i == j);
        }
    }
}

TEST(RowSparse, Serialize) {
    Vector<Vector<uint64_t>> set_bits = { {0}, {1}, {2}, {3} };
    using RowCallback = annot::binmat::BinaryMatrix::RowCallback;
    annot::binmat::RowSparse under_test([&](const RowCallback& callback) {
      for (auto& row : set_bits) {
          callback(row);
      }
    }, 4, 4, 4);

    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    under_test.serialize(out);

    std::ifstream &in = tempfile.ifstream();
    annot::binmat::RowSparse loaded;
    loaded.load(in);

    EXPECT_EQ(4, loaded.num_relations());
    EXPECT_EQ(4, loaded.num_columns());
    EXPECT_EQ(4, loaded.num_rows());
    for (uint32_t i = 0; i < 4; ++i) {
        for(uint32_t j = 0; j < 4; ++j) {
            ASSERT_EQ(loaded.get(i,j), i == j);
        }
    }
}

} // namespace
