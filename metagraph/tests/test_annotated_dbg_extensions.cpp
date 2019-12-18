#include "gtest/gtest.h"

#include "column_analysis.hpp"
#include "dbg_succinct.hpp"
#include "annotated_dbg.hpp"
#include "annotate.hpp"
#include "graph/test_dbg_helpers.hpp"
#include "annotation/test_annotated_dbg_helpers.hpp"


const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column_compressed";


template <typename Graph>
class ColumnAnalysisTest : public ::testing::Test {};

TYPED_TEST_CASE(ColumnAnalysisTest, MaskedGraphTypes);

TYPED_TEST(ColumnAnalysisTest, Density) {
    for (size_t k = 2; k <= 20; ++k) {

        std::string seq_first = std::string(k, 'A') + 'C';
        std::string seq_second = std::string(k * 2, 'G') + 'T';

        auto anno_graph = build_anno_graph<TypeParam, annotate::ColumnCompressed<>>(
            k + 1,
            {seq_first, seq_second},
            {"First", "Second"}
        );

        auto annotation = &anno_graph->get_annotation();
        annotation->serialize(test_dump_basename_vec_good);

        anno_graph->add_extension(std::make_shared<ColumnAnalysis<>>());
        auto column_analysis_ext = anno_graph->template get_extension<ColumnAnalysis<>>();
        column_analysis_ext->load({test_dump_basename_vec_good});

        auto densities = column_analysis_ext->densities();

        EXPECT_EQ(densities["First"], 1/double(3));
        EXPECT_EQ(densities["Second"], 2/double(3));
    }
}

TYPED_TEST(ColumnAnalysisTest, Summarize) {
    for (size_t k = 2; k <= 20; ++k) {

        std::string seq_first = std::string(k, 'A') + 'C';
        std::string seq_first2 = std::string(k + 1, 'A');
        std::string seq_second = std::string(k * 2, 'G') + 'T';
        std::string seq_second2 = std::string(k + 1, 'G');

        auto anno_graph = build_anno_graph<TypeParam, annotate::ColumnCompressed<>>(
            k + 1,
            {seq_first, seq_second, seq_first2, seq_second2},
            {"First", "Second", "Second", "Second"}
        );

        auto annotation = &anno_graph->get_annotation();
        annotation->serialize(test_dump_basename_vec_good);

        anno_graph->add_extension(std::make_shared<ColumnAnalysis<>>());
        auto column_analysis_ext = anno_graph->template get_extension<ColumnAnalysis<>>();
        column_analysis_ext->load({test_dump_basename_vec_good});

        auto col1 = column_analysis_ext->summarize_to_new_column(.5f, {std::make_pair("First", 1)});
        auto col2 = column_analysis_ext->summarize_to_new_column(.5f, {std::make_pair("First", -1), std::make_pair("Second", 1)});
        EXPECT_EQ(col1->num_set_bits(), 1u);
        EXPECT_EQ(col2->num_set_bits(), 3u);
    }
}
