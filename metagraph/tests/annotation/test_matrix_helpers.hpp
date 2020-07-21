#ifndef __TEST_MATRIX_HELPERS_HPP__
#define __TEST_MATRIX_HELPERS_HPP__

#include <vector>
#include <string>
#include <memory>

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "common/vectors/bit_vector.hpp"


namespace mtg {

namespace annot {
namespace binmat {

class BRWTOptimized : public annot::binmat::BRWT {
  public:
    template <typename... Args>
    BRWTOptimized(Args&&... args) : BRWT(std::forward<Args>(args)...) {
        annot::binmat::BRWTOptimizer::relax(this);
    }
};

template <int BufferSize>
class RainbowfishBuffer : public annot::binmat::Rainbowfish {
  public:
    template <typename... Args>
    RainbowfishBuffer(Args&... args)
          : Rainbowfish(args..., BufferSize) {}
    RainbowfishBuffer()
          : Rainbowfish() {}
};

} // namespace binmat
} // namespace annot

namespace test {

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_matrix";

typedef std::vector<std::unique_ptr<bit_vector>> BitVectorPtrArray;

template <typename BinMat>
BinMat build_matrix_from_columns(const BitVectorPtrArray &columns = {},
                                 uint64_t num_rows = 0);

template <typename TypeParam>
void test_matrix(const TypeParam &matrix, const BitVectorPtrArray &columns);

} // namespace test
} // namespace mtg

#endif // __TEST_MATRIX_HELPERS_HPP__
