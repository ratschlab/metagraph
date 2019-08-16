#ifndef __TEST_MATRIX_HELPERS_HPP__
#define __TEST_MATRIX_HELPERS_HPP__

#include <vector>
#include <string>
#include <memory>

#include "binary_matrix.hpp"
#include "BRWT.hpp"
#include "BRWT_builders.hpp"
#include "rainbowfish.hpp"
#include "bit_vector.hpp"


class BRWTOptimized : public BRWT {
  public:
    template <typename... Args>
    BRWTOptimized(Args&&... args)
        : BRWT(std::forward<Args>(args)...) { BRWTOptimizer::relax(this); }
};

template <int BufferSize>
class RainbowfishBuffer : public Rainbowfish {
  public:
    template <typename... Args>
    RainbowfishBuffer(Args&... args)
          : Rainbowfish(args..., BufferSize) {}
    RainbowfishBuffer()
          : Rainbowfish() {}
};

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_matrix";

typedef std::vector<std::unique_ptr<bit_vector>> BitVectorPtrArray;

template <typename BinMat>
BinMat build_matrix_from_columns(BitVectorPtrArray&& columns = {}, uint64_t num_rows = 0);

template <typename BinMat>
BinMat build_matrix_from_rows(BitVectorPtrArray&& columns = {}, uint64_t num_rows = 0);


template <typename TypeParam>
void test_matrix(const TypeParam &matrix, const BitVectorPtrArray &columns);


#endif // __TEST_MATRIX_HELPERS_HPP__
