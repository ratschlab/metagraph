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

using namespace mtg::bitmap_graph;

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels);

MaskedDeBruijnGraph build_masked_graph(const AnnotatedDBG &anno_graph,
                                       const std::vector<std::string> &ingroup,
                                       const std::vector<std::string> &outgroup,
                                       double mask_in_label_fraction = 1.0,
                                       double mask_out_label_fraction = 0.0,
                                       double other_label_fraction = 1.0,
                                       double lazy_evaluation_density_cutoff = 0.05);

typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct> MaskedGraphTypes;

typedef ::testing::Types<DBGBitmap,
                         DBGSuccinct> MaskedStableGraphTypes;

typedef ::testing::Types<std::pair<DBGBitmap, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashString, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annotate::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annotate::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashString, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashFast, annotate::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annotate::RowFlatAnnotator>
                        > GraphAnnotationPairTypes;

typedef ::testing::Types<std::pair<DBGBitmap, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annotate::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashFast, annotate::RowFlatAnnotator>
                        > GraphNoNAnnotationPairTypes;

typedef ::testing::Types<std::pair<DBGHashString, annotate::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashString, annotate::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annotate::RowFlatAnnotator>
                        > GraphWithNAnnotationPairTypes;

} // namespace test
} // namespace mtg

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
