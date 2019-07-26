#ifndef __TEST_DBG_HELPERS_HPP__
#define __TEST_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <memory>

#include "sequence_graph.hpp"
#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences = {},
            bool canonical = false,
            bool count_kmers = false);

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
                         DBGSuccinct> GraphTypes;

#endif // __TEST_DBG_HELPERS_HPP__
