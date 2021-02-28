#ifndef __TEST_ANNOTATED_DBG_HELPERS_HPP__
#define __TEST_ANNOTATED_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <memory>
#include <vector>
#include <string>

#define protected public
#define private public
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/masked_graph.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph;

template <class Graph, class Annotation = annot::ColumnCompressed<>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels,
                                               DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC);

template <class Annotation = annot::ColumnCompressed<>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(std::shared_ptr<DeBruijnGraph> graph,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels);

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
