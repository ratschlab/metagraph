#ifndef __TEST_DBG_HELPERS_HPP__
#define __TEST_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_fast.hpp"
#include "dbg_bitmap.hpp"


template <size_t numerator, size_t denominator>
class DBGSuccinctBloomFPR : public DBGSuccinct {
  public:
    template <typename... Args>
    DBGSuccinctBloomFPR(Args&&... args)
          : DBGSuccinct(std::forward<Args>(args)...) {}
};

template <size_t filter_size, size_t num_hash_functions>
class DBGSuccinctBloom : public DBGSuccinct {
  public:
    template <typename... Args>
    DBGSuccinctBloom(Args&&... args)
          : DBGSuccinct(std::forward<Args>(args)...) {}
};

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences = {},
            bool canonical = false);

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences = {},
                  bool canonical = false);

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_iterative(uint64_t k,
                      std::function<void(std::function<void(const std::string&)>)> generate,
                      bool canonical = false);

template <class Graph>
bool check_graph(const std::string &alphabet, bool canonical, bool check_sequence = false);


template <typename Graph>
class DeBruijnGraphTest : public ::testing::Test { };
typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct,
                         DBGSuccinctBloomFPR<1, 1>,
                         DBGSuccinctBloomFPR<1, 10>,
                         DBGSuccinctBloom<4, 1>,
                         DBGSuccinctBloom<4, 50>> GraphTypes;

// in stable graphs the order of input sequences
// does not change the order of k-mers and their indexes
template <typename Graph>
class StableDeBruijnGraphTest : public ::testing::Test { };
typedef ::testing::Types<DBGBitmap,
                         DBGSuccinct,
                         DBGSuccinctBloomFPR<1, 1>,
                         DBGSuccinctBloomFPR<1, 10>,
                         DBGSuccinctBloom<4, 1>,
                         DBGSuccinctBloom<4, 50>> StableGraphTypes;


#endif // __TEST_DBG_HELPERS_HPP__
