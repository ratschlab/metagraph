#ifndef __TEST_DBG_HELPERS_HPP__
#define __TEST_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph;

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

template <size_t indexed_suffix_length>
class DBGSuccinctIndexed : public DBGSuccinct {
  public:
    template <typename... Args>
    DBGSuccinctIndexed(Args&&... args)
          : DBGSuccinct(std::forward<Args>(args)...) {}
};

enum DBGMode {
    NORMAL, CANONICAL, CANONICAL_WRAPPER, PRIMARY
};

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences = {},
            DBGMode mode = NORMAL);

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences = {},
                  DBGMode mode = NORMAL);

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_iterative(uint64_t k,
                      std::function<void(std::function<void(const std::string&)>)> generate,
                      DBGMode mode = NORMAL) {
    std::vector<std::string> sequences;
    generate([&](const auto &sequence) { sequences.push_back(sequence); });
    return build_graph_batch<Graph>(k, sequences, mode);
}

template <class Graph>
bool check_graph(const std::string &alphabet, DBGMode mode, bool check_sequence);

template <class Graph>
bool check_graph_nodes(const Graph &graph) {
    size_t num_nodes = 0;
    graph.call_nodes([&](auto) { num_nodes++; });
    return num_nodes == graph.num_nodes();
}


template <typename Graph>
class DeBruijnGraphTest : public ::testing::Test { };
typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct,
                         DBGSuccinctIndexed<1>,
                         DBGSuccinctIndexed<2>,
                         DBGSuccinctIndexed<10>,
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
                         DBGSuccinctIndexed<1>,
                         DBGSuccinctIndexed<2>,
                         DBGSuccinctIndexed<10>,
                         DBGSuccinctBloomFPR<1, 1>,
                         DBGSuccinctBloomFPR<1, 10>,
                         DBGSuccinctBloom<4, 1>,
                         DBGSuccinctBloom<4, 50>> StableGraphTypes;

typedef ::testing::Types<DBGHashFast,
                         DBGSuccinct,
                         DBGSuccinctIndexed<10>,
                         DBGSuccinctBloomFPR<1, 10>> FewGraphTypes;

} // namespace test
} // namespace mtg

#endif // __TEST_DBG_HELPERS_HPP__
