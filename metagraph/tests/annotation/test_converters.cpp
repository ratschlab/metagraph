#include <filesystem>
#include <random>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/annotation_converters.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace {

using namespace mtg;
using namespace mtg::annot;
using namespace ::testing;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_row_compressed_merge = test_dump_basename + "_row_compressed_merge";
const std::string test_dump_basename_rowflat_merge = test_dump_basename + "_rowflat_merge";
const std::string test_dump_basename_row_compressed_to_rowflat = test_dump_basename + "_row_compressed_to_rowflat";


class ConvertFromRowCompressed : public ::testing::Test {
  protected:
    static const uint64_t num_rows = 5;
    std::unique_ptr<RowCompressed<>> initial_annotation;
    std::unique_ptr<MultiLabelEncoded<std::string>> annotation;
    std::unique_ptr<graph::DBGSuccinct> graph;

    virtual void SetUp() {
        initial_annotation = std::make_unique<RowCompressed<>>(num_rows);
        initial_annotation->add_labels({ 0 }, { "Label0", "Label2", "Label8" });
        initial_annotation->add_labels({ 2 }, { "Label1", "Label2" });
        initial_annotation->add_labels({ 3 }, { "Label1", "Label2", "Label8" });
        initial_annotation->add_labels({ 4 }, { "Label2" });

        graph.reset(new graph::DBGSuccinct(3));
        graph->add_sequence("ACGTCAC");
        graph->mask_dummy_kmers(1, false);
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_THAT(annotation->get(0), UnorderedElementsAre("Label0", "Label2", "Label8"));
        EXPECT_THAT(annotation->get(1), UnorderedElementsAre());
        EXPECT_THAT(annotation->get(2), UnorderedElementsAre("Label1", "Label2"));
        EXPECT_THAT(annotation->get(3), UnorderedElementsAre("Label1", "Label2", "Label8"));
        EXPECT_THAT(annotation->get(4), UnorderedElementsAre("Label2"));
    }
};

const uint64_t ConvertFromRowCompressed::num_rows;

class MergeAnnotators : public ::testing::Test {
  protected:
    static const uint64_t num_rows = 5;
    static RowCompressed<> *input_annotation_1;
    static RowCompressed<> *input_annotation_2;
    static RowCompressed<> *merged_annotation_expected;
    static MultiLabelEncoded<std::string> *merged_annotation;

    virtual void SetUp() {

        input_annotation_1 = new RowCompressed<>(num_rows);
        input_annotation_1->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        input_annotation_1->add_labels({ 2 }, {"Label1", "Label2"});
        input_annotation_1->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        input_annotation_1->add_labels({ 4 }, {"Label2"});

        input_annotation_2 = new RowCompressed<>(num_rows);
        input_annotation_2->add_labels({ 1 }, {"Label0", "Label3"});
        input_annotation_2->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        input_annotation_2->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});

        merged_annotation_expected = new RowCompressed<>(num_rows);
        merged_annotation_expected->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        merged_annotation_expected->add_labels({ 1 }, {"Label0", "Label3"});
        merged_annotation_expected->add_labels({ 2 }, {"Label0", "Label2", "Label1", "Label9", "Label7"});
        merged_annotation_expected->add_labels({ 3 }, {"Label2", "Label8", "Label1"});
        merged_annotation_expected->add_labels({ 4 }, {"Label2", "Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(merged_annotation);
        for(uint64_t i = 0; i < num_rows; ++i) {
            auto row_expected = merged_annotation_expected->get(i);
            auto row = merged_annotation->get(i);
            EXPECT_EQ(row_expected, row);
        }

        delete input_annotation_1;
        delete input_annotation_2;
        delete merged_annotation_expected;
        delete merged_annotation;
    }
};

const uint64_t MergeAnnotators::num_rows;
RowCompressed<> *MergeAnnotators::input_annotation_1 = nullptr;
RowCompressed<> *MergeAnnotators::input_annotation_2 = nullptr;
RowCompressed<> *MergeAnnotators::merged_annotation_expected = nullptr;
MultiLabelEncoded<std::string> *MergeAnnotators::merged_annotation = nullptr;


class ConvertFromColumnCompressed : public ::testing::Test {
  protected:
    std::unique_ptr<ColumnCompressed<>> initial_annotation;
    std::unique_ptr<MultiLabelEncoded<std::string>> annotation;
    std::unique_ptr<graph::DBGSuccinct> graph;


    virtual void SetUp() {
        initial_annotation = std::make_unique<ColumnCompressed<>>(5);
        initial_annotation->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        initial_annotation->add_labels({ 2 }, {"Label1", "Label2"});
        initial_annotation->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        initial_annotation->add_labels({ 4 }, {"Label2"});

        graph.reset(new graph::DBGSuccinct(3));
        graph->add_sequence("ACGTCAC");
        graph->mask_dummy_kmers(1, false);
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_THAT(annotation->get(0), UnorderedElementsAre("Label0", "Label2", "Label8"));
        EXPECT_THAT(annotation->get(1), UnorderedElementsAre());
        EXPECT_THAT(annotation->get(2), UnorderedElementsAre("Label1", "Label2"));
        EXPECT_THAT(annotation->get(3), UnorderedElementsAre("Label1", "Label2", "Label8"));
        EXPECT_THAT(annotation->get(4), UnorderedElementsAre("Label2"));
    }

};

std::unique_ptr<graph::DBGSuccinct> create_graph(uint32_t k, std::vector<string> sequences) {
    auto graph = std::make_unique<graph::DBGSuccinct>(k);
    for(const auto& s: sequences) {
        graph->add_sequence(s);
    }
    graph->mask_dummy_kmers(1, false);
    return graph;
}

void clean_column_diff_files(const std::string &graph_name) {
    std::string outfbase = test_dump_basename + graph_name;
    // TODO - update extension if needed
    for (const auto &suf :
         { ".succ", ".pred", ".pred_boundary", ".terminal", "row_diff.column.annodbg" }) {
        std::filesystem::remove(outfbase + suf);
    }
}

TEST(ColumnDiff, succ) {
    std::vector<std::unique_ptr<ColumnCompressed<>>> empty_source(1);
    empty_source[0] = std::make_unique<ColumnCompressed<>>(5);

    std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAG" });

    const std::string outfbase = test_dump_basename + "column.diff.succ";
    const std::string succ_file = outfbase + ".succ";
    const std::string pred_file = outfbase + ".pred";
    const std::string pred_boundary_file = outfbase + ".pred_boundary";
    const std::vector<std::vector<uint64_t>> expected_succ
            = { { 0, 0, 0, 0, 0 }, { 0, 4, 1, 5, 0 }, { 0, 4, 1, 5, 3 } };
    const std::vector<std::vector<uint64_t>> expected_pred
            = { {}, { 3, 2, 4 }, { 3, 5, 2, 4 } };
    const std::vector<std::vector<bool>> expected_boundary
            = { { 1, 1, 1, 1, 1 }, { 0, 1, 1, 1, 0, 1, 0, 1 }, { 0, 1, 1, 0, 1, 0, 1, 0, 1 } };

    for (uint32_t max_depth : { 1, 3, 5 }) {
        clean_column_diff_files("column.diff.succ");

        convert_to_column_diff(*graph, empty_source, outfbase, max_depth);

        ASSERT_TRUE(std::filesystem::exists(succ_file));
        ASSERT_TRUE(std::filesystem::exists(pred_file));
        ASSERT_TRUE(std::filesystem::exists(pred_boundary_file));

        sdsl::int_vector_buffer succ(succ_file, std::ios::in);

        uint32_t idx = max_depth / 2;
        ASSERT_EQ(5, succ.size());
        for (uint32_t i = 0; i < succ.size(); ++i) {
            EXPECT_EQ(expected_succ[idx][i], succ[i]);
        }

        sdsl::int_vector_buffer pred(pred_file, std::ios::in);

        EXPECT_EQ(expected_pred[idx].size(), pred.size());
        for (uint32_t i = 0; i < pred.size(); ++i) {
            EXPECT_EQ(expected_pred[idx][i], pred[i]) << max_depth << " " << i;
        }

        sdsl::rrr_vector boundary;
        std::ifstream fpred_boundary(pred_boundary_file, std::ios::binary);
        boundary.load(fpred_boundary);

        ASSERT_EQ(expected_boundary[idx].size(), boundary.size());
        for(uint32_t i = 0; i < expected_boundary[idx].size(); ++i) {
            EXPECT_EQ(expected_boundary[idx][i], boundary[i]);
        }
    }

    clean_column_diff_files("column.diff.succ");
}

TEST(ColumnDiff, ConvertFromColumnCompressedEmpty) {
    std::vector<std::unique_ptr<ColumnCompressed<>>> empty_source(1);
    empty_source[0] = std::make_unique<ColumnCompressed<>>(5);

    std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAG" });


    clean_column_diff_files("column.diff.empty");
    std::string outfbase = test_dump_basename + "column.diff.empty";
    std::vector<std::unique_ptr<ColumnDiffAnnotator>> result
            = convert_to_column_diff<std::string>(*graph, empty_source, outfbase, 1);

    ASSERT_EQ(0u, result.size());

    clean_column_diff_files("column.diff.empty");
}

TEST(CoumnDiff, ConvertFromColumnCompressedSameLabels) {
    const std::string outfbase = test_dump_basename + "column.diff.convert";
    const std::string succ_file = outfbase + ".succ";
    clean_column_diff_files("column.diff.convert");

    std::vector<std::vector<std::string>> label_groups
            = { { "Label0" }, { "Label1", "Label2" }, { "Label0", "Label2", "Label1" } };
    for (const auto &labels : label_groups) {
        auto initial_annotation = std::make_unique<ColumnCompressed<>>(5);
        initial_annotation->add_labels({ 0, 1, 2, 3, 4 }, labels);

        std::vector<std::unique_ptr<ColumnCompressed<>>> source(1);
        source[0] = std::move(initial_annotation);

        auto graph = std::make_unique<graph::DBGSuccinct>(3);
        graph->add_sequence("ACGTCAC");
        graph->mask_dummy_kmers(1, false);

        const uint32_t expected_relations[] = { 5, 2, 1, 1, 1 };

        for (const uint32_t max_depth : { 1, 2, 3, 4, 5 }) {
            std::vector<std::unique_ptr<ColumnDiffAnnotator>> result
                    = convert_to_column_diff(*graph, source, outfbase, max_depth);
            clean_column_diff_files("column.diff.convert");

            ASSERT_EQ(1U, result.size());

            ASSERT_EQ(labels.size(), result[0]->num_labels());
            ASSERT_EQ(5u, result[0]->num_objects());
            EXPECT_EQ(labels.size() * expected_relations[max_depth - 1],
                      result[0]->num_relations());

            for (uint32 i = 0; i < result[0]->num_objects(); ++i) {
                ASSERT_THAT(result[0]->get(i), ContainerEq(labels));
            }
        }
    }
}

TEST(CoumnDiff, ConvertFromColumnCompressedSameLabelsMultipleColumns) {
    const std::string outfbase = test_dump_basename + "column.diff.convert";
    const std::string succ_file = outfbase + ".succ";
    clean_column_diff_files("column.diff.convert");

    std::vector<std::vector<std::string>> label_groups
            = { { "Label0" }, { "Label1", "Label2" }, { "Label0", "Label2", "Label1" } };
    for (const auto &labels : label_groups) {
        std::vector<std::unique_ptr<ColumnCompressed<>>> sources;
        for (const std::string &label : labels) {
            auto initial_annotation = std::make_unique<ColumnCompressed<>>(5);
            initial_annotation->add_labels({ 0, 1, 2, 3, 4 }, { label });
            sources.push_back(std::move(initial_annotation));
        }

        auto graph = std::make_unique<graph::DBGSuccinct>(3);
        graph->add_sequence("ACGTCAC");
        graph->mask_dummy_kmers(1, false);

        const uint32_t expected_relations[] = { 5, 2, 1, 1, 1 };

        for (const uint32_t max_depth : { 1, 2, 3, 4, 5 }) {
            std::vector<std::unique_ptr<ColumnDiffAnnotator>> result
                    = convert_to_column_diff(*graph, sources, outfbase, max_depth);
            clean_column_diff_files("column.diff.convert");

            ASSERT_EQ(labels.size(), result.size());
            for (uint32_t i = 0; i < result.size(); ++i) {
                ASSERT_EQ(1, result[i]->num_labels());
                ASSERT_EQ(5u, result[i]->num_objects());
                EXPECT_EQ(expected_relations[max_depth - 1],
                          result[i]->num_relations());

                for (uint32 idx = 0; idx < result[i]->num_objects(); ++idx) {
                    ASSERT_THAT(result[i]->get(idx), ElementsAre(labels[i]));
                }
            }
        }
    }
}

void test_column_diff(uint32_t k,
                      uint32_t max_depth,
                      const std::vector<std::string> &sequences,
                      const std::vector<std::vector<std::string>> &annotations,
                      const std::string &prefix) {
    const std::string outfbase = test_dump_basename + prefix;
    clean_column_diff_files(prefix);

    auto graph = std::make_unique<graph::DBGSuccinct>(k);
    for (const auto& seq : sequences) {
        graph->add_sequence(seq);
    }
    graph->mask_dummy_kmers(1, false);

    auto initial_annotation = std::make_unique<ColumnCompressed<>>(graph->num_nodes());
    std::vector<uint32_t> added_idx(graph->num_nodes());
    std::unordered_set<std::string> all_labels;
    for (uint32_t anno_idx = 0; anno_idx < graph->num_nodes(); ++anno_idx) {
        const std::vector<std::string> &labels = annotations[anno_idx];
        initial_annotation->add_labels({anno_idx}, labels);
        std::for_each(labels.begin(), labels.end(), [&](auto l) { all_labels.insert(l); });
    }

    std::vector<std::unique_ptr<ColumnCompressed<>>> source(1);
    source[0] = std::move(initial_annotation);

    std::vector<std::unique_ptr<ColumnDiffAnnotator>> result
            = convert_to_column_diff(*graph, source, outfbase, max_depth);
    clean_column_diff_files(prefix);

    ASSERT_EQ(1U, result.size());
    ASSERT_EQ(all_labels.size(), result[0]->num_labels());
    ASSERT_EQ(graph->num_nodes(), result[0]->num_objects());

    for (uint32_t anno_idx = 0; anno_idx < graph->num_nodes(); ++anno_idx) {
        ASSERT_THAT(result[0]->get(anno_idx),
                    UnorderedElementsAreArray(annotations[anno_idx]));
    }
}

TEST(CoumnDiff, ConvertFromColumnCompressed) {
    test_column_diff(3, 5, { "ACGTCAC" },
                     { { "Label0", "Label2", "Label8" },
                       {},
                       { "Label1", "Label2" },
                       { "Label1", "Label2", "Label8" },
                       { "Label2" } },
                     "column.diff.convert");
}

TEST(ColumnDiff, ConvertFromColumnCompressed4Loops) {
    std::vector<std::string> sequences = { std::string(100, 'A'), std::string(100, 'G'),
                                           std::string(100, 'C'), std::string(100, 'T') };

    test_column_diff(3, 5, sequences,
                     { { "Lb0", "Lb2", "Lb8" }, {}, { "Lb9" }, { "Lb0", "Lb8", "Lb2" } },
                     "column.diff.4loops");
}

TEST(ColumnDiff, ConvertFromColumnCompressed4PathsRandomLabels) {
    std::mt19937 gen(12345);
    std::uniform_int_distribution<> distrib(0, 2);


    std::vector<std::vector<std::string>> options = {{"L1"}, {"L1", "L2"}, {"L3"}};
    std::vector<std::vector<std::string>> annotations;
    for (uint32_t anno_idx = 0; anno_idx < 24; ++anno_idx) {
        annotations.push_back(options[distrib(gen)]);
    }

    test_column_diff(3, 5, { "ATCGGAAGA", "TTTAAACCCGGG", "ATACAGCTCGCT", "AAAAAA" },
                     annotations, "column.diff.4paths");
}

TEST(ColumnDiff, ConvertFromColumnCompressed2BigLoops) {
    std::mt19937 gen(12345);
    std::uniform_int_distribution<> distrib(0, 2);

    std::vector<std::string> sequences
            = { "ATCGGAAGAGCACACGTCTG"
                "AACTCCAGACA"
                "CTAAGGCATCTCGTATGCATCGGAAGAGC",
                "GTGAGGCGTCATGCATGCAT"
                "TGTCTGGAGTT"
                "TCGTAGCGGCGGCTAGTGCGCGTAGTGAGGCGTCA" };


    std::vector<std::vector<std::string>> options = {{"L1"}, {"L1", "L2"}, {"L2"}};
    std::vector<std::vector<std::string>> annotations;
    for (uint32_t anno_idx = 0; anno_idx < 104; ++anno_idx) {
        annotations.push_back(options[distrib(gen)]);
    }

    test_column_diff(10, 3, sequences, annotations, "column.diff.2bigloops");
}


TEST(CoumnDiff, ConvertFromColumnCompressedWithMergesAndBifurcations) {
    const std::string outfbase = test_dump_basename + "column.diff.convert";
    const std::string succ_file = outfbase + ".succ";
    clean_column_diff_files("column.diff.convert");

    auto initial_annotation = std::make_unique<ColumnCompressed<>>(5);
    initial_annotation->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
    initial_annotation->add_labels({ 2 }, {"Label1", "Label2"});
    initial_annotation->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
    initial_annotation->add_labels({ 4 }, {"Label2"});

    std::vector<std::unique_ptr<ColumnCompressed<>>> source(1);
    source[0] = std::move(initial_annotation);

    auto graph = std::make_unique<graph::DBGSuccinct>(3);
    graph->add_sequence("ACGTCAC");
    graph->mask_dummy_kmers(1, false);

    constexpr uint32_t max_depth = 5;
    std::vector<std::unique_ptr<ColumnDiffAnnotator>> result
            = convert_to_column_diff(*graph, source, outfbase, max_depth);
    clean_column_diff_files("column.diff.convert");

    ASSERT_EQ(1U, result.size());

    ASSERT_EQ(4u, result[0]->num_labels());
    ASSERT_EQ(5u, result[0]->num_objects());
    // we added 9 relations, and after sparsification we have 12 (because in our
    // test graph there is no correlation between annotation and vicinity)
    ASSERT_EQ(12u, result[0]->num_relations());

    ASSERT_THAT(result[0]->get(0), UnorderedElementsAre("Label0", "Label2", "Label8"));
    ASSERT_THAT(result[0]->get(1), UnorderedElementsAre());
    ASSERT_THAT(result[0]->get(2), UnorderedElementsAre("Label1", "Label2"));
    ASSERT_THAT(result[0]->get(3), UnorderedElementsAre("Label1", "Label2", "Label8"));
    ASSERT_THAT(result[0]->get(4), UnorderedElementsAre("Label2"));
}

// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT) {
//     ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<BinRelWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT) {
    annotation = convert<BinRelWTAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT_sdsl) {
//     ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<BinRelWT_sdslAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT_sdsl) {
    annotation = convert<BinRelWT_sdslAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromColumnCompressedEmpty, to_RowFlat) {
//     ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RowFlat) {
    annotation = convert<RowFlatAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromColumnCompressedEmpty, to_Rainbowfish) {
//     ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RainbowfishAnnotator) {
    annotation = convert<RainbowfishAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromColumnCompressedEmpty, to_GreedyBRWT) {
//     ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert_to_greedy_BRWT<MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_GreedyBRWT) {
    annotation = convert_to_greedy_BRWT<MultiBRWTAnnotator>(std::move(*initial_annotation));
}


TEST(ConvertFromRowCompressedEmpty, to_BinRelWT) {
    RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = convert<BinRelWTAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT) {
    annotation = convert<BinRelWTAnnotator>(std::move(*initial_annotation));
}

TEST(ConvertFromRowCompressedEmpty, to_BinRelWT_sdsl) {
    RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = convert<BinRelWT_sdslAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT_sdsl) {
    annotation = convert<BinRelWT_sdslAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromRowCompressedEmpty, to_RowFlat) {
//     RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RowFlat) {
    annotation = convert<RowFlatAnnotator>(std::move(*initial_annotation));
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat) {
    initial_annotation->serialize(test_dump_basename_row_compressed_to_rowflat);

    annotation = convert<RowFlatAnnotator>(test_dump_basename_row_compressed_to_rowflat);
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat2) {
    const auto rowflat_filename = test_dump_basename_row_compressed_to_rowflat
                                                + RowCompressed<>::kExtension;

    initial_annotation->serialize(rowflat_filename);

    annotation = convert<RowFlatAnnotator>(rowflat_filename);
}

// TEST(ConvertFromRowCompressedEmpty, to_Rainbowfish) {
//     RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert<RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RainbowfishAnnotator) {
    annotation = convert<RainbowfishAnnotator>(std::move(*initial_annotation));
}

// TEST(ConvertFromRowCompressedEmpty, to_GreedyBRWT) {
//     RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = convert_to_greedy_BRWT<MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

// TEST_F(ConvertFromRowCompressed, to_GreedyBRWT) {
//     annotation = convert_to_greedy_BRWT<MultiBRWTAnnotator>(
//         std::move(*initial_annotation)
//     ).release();
// }

TEST_F(MergeAnnotators, RowCompressed) {
    std::vector<std::string> filenames;
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_1";
        input_annotation_1->serialize(filename);
        filenames.push_back(filename + RowCompressed<>::kExtension);
    }
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_2";
        input_annotation_2->serialize(filename);
        filenames.push_back(filename + RowCompressed<>::kExtension);
    }

    merge<RowCompressed<>, std::string>(
        {}, filenames, test_dump_basename_row_compressed_merge + "_merged"
    );

    merged_annotation = new RowCompressed<>(num_rows);
    merged_annotation->merge_load({ test_dump_basename_row_compressed_merge + "_merged" });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowCompressed) {
    std::vector<std::unique_ptr<MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(convert<RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(convert<RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowcompressed";
    merge<RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new RowCompressed<>(num_rows);
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowFlat) {
    std::vector<std::unique_ptr<MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(convert<RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(convert<RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowflat";
    merge<RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new RowFlatAnnotator();
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, Mixed_to_RowFlat) {
    std::vector<std::unique_ptr<MultiLabelEncoded<std::string>>> annotators;
    std::vector<std::string> filenames;
    {
        auto annotator = convert<RowFlatAnnotator>(
            std::move(*input_annotation_1)
        );
        annotators.push_back(std::move(annotator));
    }
    {
        auto annotator = std::make_unique<ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    //TODO
    {
        auto annotator = std::make_unique<ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    {
        //TODO: move into fixture as input_annotation_3 and make non-overlapping
        const std::string filename = test_dump_basename_row_compressed_merge + "_mixed_2";
        auto annotation = std::make_unique<RowCompressed<> >(num_rows);
        annotation->add_labels({ 0 }, {"Label0"});
        annotation->add_labels({ 2 }, {"Label1"});
        annotation->add_labels({ 3 }, {"Label1"});
        annotation->add_labels({ 4 }, {"Label2"});
        annotation->serialize(filename);
        filenames.push_back(filename + RowCompressed<>::kExtension);
    }

    const auto outfile = test_dump_basename_rowflat_merge + "_mixed_to_rowflat";
    merge<RowFlatAnnotator, std::string>(std::move(annotators), filenames, outfile);

    merged_annotation = new RowFlatAnnotator();
    merged_annotation->merge_load({ outfile });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST(ColumnCompressed, ToRowAnnotator) {
    {
        ColumnCompressed<> annotation(0);

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);
    }
    {
        ColumnCompressed<> annotation(1);

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorParallel) {
    {
        ColumnCompressed<> annotation(0);

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);
    }
    {
        ColumnCompressed<> annotation(1);

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreaming) {
    {
        ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreamingParallel) {
    {
        ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename, 10);

        RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

} // namespace
