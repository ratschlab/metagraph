#ifndef __TEST_DBG_HELPERS_HPP__
#define __TEST_DBG_HELPERS_HPP__

#include <string>
#include <vector>
#include <memory>

#include "sequence_graph.hpp"

template <class Graph>
std::unique_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences = {},
            bool canonical = false);

template <class Graph>
std::unique_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences = {},
                  bool canonical = false);

template <class Graph>
bool check_graph(const std::string &alphabet, bool canonical);

#endif // __TEST_DBG_HELPERS_HPP__
