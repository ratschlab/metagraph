#include "gtest/gtest.h"

#include "test_matrix_helpers.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

template <int BufferSize>
void test_rainbowfish_buffer(const uint64_t num_rows) {
    const std::vector<bit_vector_stat> vectors {
        { 0, 1, 1, 0, 0, 1, 0 },
        { 0, 1, 1, 1, 1, 1, 0 }
    };

    BitVectorPtrArray columns, copy;
    for (size_t j = 0; j < vectors.size(); ++j) {
        sdsl::bit_vector bv(num_rows, false);
        ASSERT_GE(vectors.at(j).size(), num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            if (vectors.at(j)[i])
                bv[i] = true;
        }
        columns.emplace_back(new bit_vector_stat(std::move(bv)));
        copy.push_back(columns.back()->copy());
    }

    test_matrix(
        build_matrix_from_columns<RainbowfishBuffer<BufferSize>>(copy, num_rows),
        columns
    );
}

TEST(BinaryMatrixRainbowfishTest, BufferAllSizes) {
    test_rainbowfish_buffer<1>(2);
    test_rainbowfish_buffer<1>(3);
    test_rainbowfish_buffer<1>(4);
    test_rainbowfish_buffer<1>(5);
    test_rainbowfish_buffer<1>(6);
    test_rainbowfish_buffer<1>(7);

    test_rainbowfish_buffer<2>(3);
    test_rainbowfish_buffer<2>(4);
    test_rainbowfish_buffer<2>(5);
    test_rainbowfish_buffer<2>(6);
    test_rainbowfish_buffer<2>(7);

    test_rainbowfish_buffer<3>(4);
    test_rainbowfish_buffer<3>(5);
    test_rainbowfish_buffer<3>(6);
    test_rainbowfish_buffer<3>(7);

    test_rainbowfish_buffer<4>(5);
    test_rainbowfish_buffer<4>(6);
    test_rainbowfish_buffer<4>(7);

    test_rainbowfish_buffer<5>(6);
    test_rainbowfish_buffer<5>(7);

    test_rainbowfish_buffer<6>(7);
}

} // namespace
