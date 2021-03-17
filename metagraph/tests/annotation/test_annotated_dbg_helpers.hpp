#ifndef __TEST_ANNOTATED_DBG_HELPERS_HPP__
#define __TEST_ANNOTATED_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <vector>
#include <string>

#define protected public
#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {

namespace graph {
    class AnnotatedDBG;
}

namespace test {

template <class Graph, class Annotation = annot::ColumnCompressed<>>
std::unique_ptr<graph::AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences,
                 const std::vector<std::string> &labels,
                 graph::DeBruijnGraph::Mode mode = graph::DeBruijnGraph::BASIC);

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
