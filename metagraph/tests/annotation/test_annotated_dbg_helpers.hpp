#ifndef __TEST_ANNOTATED_DBG_HELPERS_HPP__
#define __TEST_ANNOTATED_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <vector>
#include <string>

#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace test {

template <class Graph, class Annotation>
std::unique_ptr<graph::AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences = {},
                 const std::vector<std::string> &labels = {},
                 graph::DeBruijnGraph::Mode mode = graph::DeBruijnGraph::BASIC,
                 bool coordinates = false,
                 bool mark_seq_ends = false,
                 bool mask_dummy_kmers = true);

template <class Annotation>
std::unique_ptr<graph::AnnotatedDBG>
build_anno_graph(std::shared_ptr<graph::DeBruijnGraph> graph,
                 const std::vector<std::string> &sequences = {},
                 const std::vector<std::string> &labels = {},
                 bool coordinates = false,
                 bool mark_seq_ends = false);

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
