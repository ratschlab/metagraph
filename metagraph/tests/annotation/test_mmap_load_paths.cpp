// Loads through a mapped stream (utils::open_ifstream in mmap mode) for types
// whose mapped load had no other test: view-backed members must answer like
// the owned load, and owned members must stay writable.
#include "gtest/gtest.h"
#include "../test_helpers.hpp"
#include "test_matrix_helpers.hpp"

#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/wavelet_tree.hpp"
#include "graph/graph_extensions/node_weights.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

const std::string test_dump_basename = test_dump_dir() + "/mmap_load_paths";

// Write `prefix` padding bytes, then what `write` writes.
template <class Write>
void write_file(const std::string &filename, size_t prefix, Write write) {
    std::ofstream out(filename, std::ios::binary);
    out << std::string(prefix, 'x');
    write(out);
}

// Sets the global mmap mode and restores it on scope exit, also after a failed assertion.
class MmapMode {
  public:
    explicit MmapMode(bool mapped) : old_(utils::with_mmap()) { utils::set_mmap(mapped); }
    ~MmapMode() { utils::set_mmap(old_); }

  private:
    bool old_;
};

// Load `T` from `filename` at byte `prefix`, mapped or not.
template <class T, class Load>
void load_file(const std::string &filename, size_t prefix, bool mapped, Load load) {
    auto in = utils::open_ifstream(filename, mapped);
    in->seekg(prefix);
    load(*in);
}


template <typename WT>
class MmapLoadWaveletTree : public ::testing::Test {};
typedef ::testing::Types<wavelet_tree_stat, wavelet_tree_fast, wavelet_tree_small> MmapWaveletTreeTypes;
TYPED_TEST_SUITE(MmapLoadWaveletTree, MmapWaveletTreeTypes);

TYPED_TEST(MmapLoadWaveletTree, MappedAndOwnedAgree) {
    const auto filename = test_dump_basename + "_wt";
    std::mt19937 rng(1);
    for (uint64_t size : {0, 1, 100, 5000}) {
        for (uint8_t logsigma : {1, 3, 5}) {
            sdsl::int_vector<> symbols(size, 0, logsigma);
            for (auto &&c : symbols)
                c = rng() % (1u << logsigma);
            TypeParam built(logsigma, symbols);
            for (size_t prefix : {0, 3}) {
                write_file(filename, prefix, [&](auto &out) { built.serialize(out); });
                for (bool mapped : {false, true}) {
                    TypeParam loaded(logsigma);
                    load_file<TypeParam>(filename, prefix, mapped, [&](auto &in) {
                        ASSERT_TRUE(loaded.load(in));
                    });
                    ASSERT_EQ(size, loaded.size());
                    for (uint64_t i = 0; i < size; ++i) {
                        ASSERT_EQ(symbols[i], loaded[i]);
                        ASSERT_EQ(built.rank(symbols[i], i), loaded.rank(symbols[i], i));
                    }
                    EXPECT_EQ(symbols, loaded.to_vector());
                    EXPECT_TRUE(built == loaded);
                    // Conversion takes or copies the symbols; the copy is writable.
                    auto converted = loaded.template convert_to<wavelet_tree_dyn>();
                    EXPECT_EQ(symbols, converted.to_vector());
                    if (size) {
                        converted.set(0, (symbols[0] + 1) % (1u << logsigma));
                        EXPECT_NE(symbols[0], converted[0]);
                    }
                }
            }
        }
    }
}


TEST(MmapLoadCounts, CSCMatrixMappedAndOwnedAgree) {
    using CSC = annot::matrix::CSCMatrix<annot::matrix::BRWT, annot::CountsVector>;
    const auto filename = test_dump_basename + "_csc";
    std::mt19937_64 rng(2);
    const uint64_t num_rows = 3000;
    BitVectorPtrArray columns;
    std::vector<std::vector<uint64_t>> counts;
    std::vector<annot::CountsVector> values;
    for (size_t j = 0; j < 7; ++j) {
        sdsl::bit_vector bits(num_rows, 0);
        std::vector<uint64_t> column_counts;
        for (uint64_t i = 0; i < num_rows; ++i) {
            if (rng() % (j + 2) == 0) {
                bits[i] = 1;
                column_counts.push_back(1 + (rng() >> (rng() % 64)));
            }
        }
        columns.emplace_back(new bit_vector_stat(std::move(bits)));
        values.emplace_back(column_counts);
        counts.push_back(std::move(column_counts));
    }
    CSC built(build_matrix_from_columns<annot::matrix::BRWT>(columns, num_rows),
              std::move(values));
    std::vector<uint64_t> rows(num_rows);
    std::iota(rows.begin(), rows.end(), 0);
    const auto expected = built.get_row_values(rows);
    for (size_t prefix : {0, 5}) {
        write_file(filename, prefix, [&](auto &out) { built.serialize(out); });
        for (bool mapped : {false, true}) {
            CSC loaded;
            load_file<CSC>(filename, prefix, mapped, [&](auto &in) { ASSERT_TRUE(loaded.load(in)); });
            ASSERT_EQ(num_rows, loaded.num_rows());
            EXPECT_EQ(expected, loaded.get_row_values(rows));
        }
    }
}


TEST(MmapLoadNodeWeights, MappedLoadIsWritable) {
    const auto filename_base = test_dump_basename + "_weights";
    sdsl::int_vector<> weights(1000, 0, 7);
    for (size_t i = 0; i < weights.size(); ++i)
        weights[i] = i % 100;
    graph::NodeWeights(sdsl::int_vector<>(weights)).serialize(filename_base);
    for (bool mapped : {false, true}) {
        graph::NodeWeights loaded;
        {
            MmapMode mode(mapped);
            ASSERT_TRUE(loaded.load(filename_base));
        }
        for (size_t i = 0; i < weights.size(); ++i)
            ASSERT_EQ(weights[i], loaded[i]);
        loaded.get_data()[5] = 99;
        EXPECT_EQ(99u, loaded[5]);
    }
    graph::NodeWeights check;
    ASSERT_TRUE(check.load(filename_base));
    EXPECT_EQ(5u, check[5]);
}


TEST(MmapLoadRowSparse, MappedAndOwnedAgree) {
    const auto filename = test_dump_basename + "_row_sparse";
    std::mt19937_64 rng(4);
    const uint64_t num_rows = 2000;
    BitVectorPtrArray columns;
    for (size_t j = 0; j < 9; ++j) {
        sdsl::bit_vector bits(num_rows, 0);
        for (uint64_t i = 0; i < num_rows; ++i)
            bits[i] = rng() % (j + 2) == 0;
        columns.emplace_back(new bit_vector_stat(std::move(bits)));
    }
    auto built = build_matrix_from_columns<annot::matrix::RowSparse>(columns, num_rows);
    std::vector<uint64_t> rows(num_rows);
    std::iota(rows.begin(), rows.end(), 0);
    const auto expected = built.get_rows(rows);
    for (size_t prefix : {0, 1, 5}) {
        write_file(filename, prefix, [&](auto &out) { built.serialize(out); });
        for (bool mapped : {false, true}) {
            annot::matrix::RowSparse loaded;
            load_file<annot::matrix::RowSparse>(filename, prefix, mapped, [&](auto &in) {
                ASSERT_TRUE(loaded.load(in));
            });
            ASSERT_EQ(num_rows, loaded.num_rows());
            EXPECT_EQ(expected, loaded.get_rows(rows));
        }
    }
}


template <typename BV>
class MmapLoadBitVector : public ::testing::Test {};
typedef ::testing::Types<bit_vector_stat, bit_vector_il<>, bit_vector_il<4096>, bit_vector_hyb<>,
                         bit_vector_small, bit_vector_smart> MmapBitVectorTypes;
TYPED_TEST_SUITE(MmapLoadBitVector, MmapBitVectorTypes);

TYPED_TEST(MmapLoadBitVector, MappedAndOwnedAgree) {
    const auto filename = test_dump_basename + "_bv";
    std::mt19937_64 rng(3);
    for (uint64_t size : {0, 1, 64, 5000, 300000}) {
        for (double density : {0.01, 0.5, 0.99}) {
            sdsl::bit_vector bits(size, 0);
            std::bernoulli_distribution bit(density);
            for (uint64_t i = 0; i < size; ++i)
                bits[i] = bit(rng);
            TypeParam built(bits);
            for (size_t prefix : {0, 1, 7}) {
                write_file(filename, prefix, [&](auto &out) { built.serialize(out); });
                for (bool mapped : {false, true}) {
                    TypeParam loaded;
                    load_file<TypeParam>(filename, prefix, mapped, [&](auto &in) {
                        ASSERT_TRUE(loaded.load(in));
                    });
                    ASSERT_EQ(size, loaded.size());
                    ASSERT_EQ(built.num_set_bits(), loaded.num_set_bits());
                    for (uint64_t i = 0; i < size; i += 1 + i / 64) {
                        ASSERT_EQ(bits[i], loaded[i]);
                        ASSERT_EQ(built.rank1(i), loaded.rank1(i));
                    }
                    if constexpr (!std::is_same_v<TypeParam, bit_vector_hyb<>>) {
                        for (uint64_t r = 1; r <= loaded.num_set_bits(); r += 1 + r / 32)
                            ASSERT_EQ(built.select1(r), loaded.select1(r));
                    }
                    EXPECT_EQ(bits, loaded.to_vector());
                }
            }
        }
    }
}

} // namespace
