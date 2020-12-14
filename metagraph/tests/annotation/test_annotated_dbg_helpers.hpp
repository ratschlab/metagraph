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
                                               const std::vector<std::string> &labels);

typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct> MaskedGraphTypes;

typedef ::testing::Types<DBGBitmap,
                         DBGSuccinct> MaskedStableGraphTypes;

typedef ::testing::Types<std::pair<DBGBitmap, annot::ColumnCompressed<>>,
                         std::pair<DBGHashString, annot::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annot::RowFlatAnnotator>,
                         std::pair<DBGHashString, annot::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annot::RowFlatAnnotator>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>
                        > GraphAnnotationPairTypes;

#if ! _PROTEIN_GRAPH

typedef ::testing::Types<std::pair<DBGBitmap, annot::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annot::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annot::RowFlatAnnotator>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>
                        > GraphNoNAnnotationPairTypes;

#endif

typedef ::testing::Types<std::pair<DBGHashString, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashString, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>
                        > GraphWithNAnnotationPairTypes;

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
