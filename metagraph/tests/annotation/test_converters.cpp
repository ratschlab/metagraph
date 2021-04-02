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

TEST(RowDiff, succ) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/"row_diff_succ";
    const std::string graph_fname = dst_dir/(std::string("ACGTCAG") + graph::DBGSuccinct::kExtension);
    const std::string succ_file = graph_fname + ".succ";
    const std::string succ_boundary_file = graph_fname + ".succ_boundary";
    const std::string pred_file = graph_fname + ".pred";
    const std::string pred_boundary_file = graph_fname + ".pred_boundary";

    /**
     * Graph:
     *
     * 1: CAG
     * 2: ACG
     * 3: TCA
     * 4: CGT
     * 5: GTC
     *
     * 2 -> 4 -> 5 -> 3 -> 1
     */

    const std::vector<uint64_t> expected_succ = { 3, 0, 4, 2 };
    const std::vector<bool> expected_succ_boundary = { 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    const std::vector<uint64_t> expected_pred = { 2, 4, 1, 3 };
    const std::vector<bool> expected_pred_boundary = { 0, 1, 1, 0, 1, 0, 1, 0, 1 };

    for (uint32_t max_depth : { 1, 3, 5 }) {
        std::filesystem::remove_all(dst_dir);
        std::filesystem::create_directories(dst_dir);

        std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAG" });
        graph->serialize(graph_fname);

        convert_to_row_diff({}, graph_fname, 1e9, max_depth, dst_dir, dst_dir,
                            RowDiffStage::COMPUTE_REDUCTION);

        ASSERT_TRUE(std::filesystem::exists(succ_file));
        ASSERT_TRUE(std::filesystem::exists(succ_boundary_file));
        ASSERT_TRUE(std::filesystem::exists(pred_file));
        ASSERT_TRUE(std::filesystem::exists(pred_boundary_file));

        sdsl::int_vector_buffer succ(succ_file, std::ios::in);
        ASSERT_EQ(expected_succ.size(), succ.size());
        for (uint32_t i = 0; i < succ.size(); ++i) {
            EXPECT_EQ(expected_succ[i], succ[i]) << max_depth << " " << i;
        }

        sdsl::int_vector_buffer<1> succ_boundary(succ_boundary_file, std::ios::in);
        ASSERT_EQ(expected_succ_boundary.size(), succ_boundary.size());
        for(uint32_t i = 0; i < expected_succ_boundary.size(); ++i) {
            EXPECT_EQ(expected_succ_boundary[i], succ_boundary[i]) << max_depth << " " << i;
        }

        sdsl::int_vector_buffer pred(pred_file, std::ios::in);
        EXPECT_EQ(expected_pred.size(), pred.size());
        for (uint32_t i = 0; i < pred.size(); ++i) {
            EXPECT_EQ(expected_pred[i], pred[i]) << max_depth << " " << i;
        }

        sdsl::int_vector_buffer<1> pred_boundary(pred_boundary_file, std::ios::in);
        ASSERT_EQ(expected_pred_boundary.size(), pred_boundary.size());
        for(uint32_t i = 0; i < expected_pred_boundary.size(); ++i) {
            EXPECT_EQ(expected_pred_boundary[i], pred_boundary[i]) << max_depth << " " << i;
        }
    }

    std::filesystem::remove_all(dst_dir);
}

TEST(RowDiff, ConvertFromColumnCompressedEmpty) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/"row_diff_empty";
    std::filesystem::remove_all(dst_dir);
    std::filesystem::create_directories(dst_dir);

    std::string annot_fname
            = dst_dir/(std::string("ACGTCAG") + ColumnCompressed<>::kExtension);
    ColumnCompressed<>(5).serialize(annot_fname);

    std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAG" });
    const std::string graph_fname = dst_dir/(std::string("ACGTCAG") + graph::DBGSuccinct::kExtension);
    graph->serialize(graph_fname);

    convert_to_row_diff({ annot_fname }, graph_fname, 1e9, 1, dst_dir, dst_dir, RowDiffStage::COMPUTE_REDUCTION);
    convert_to_row_diff({ annot_fname }, graph_fname, 1e9, 1, dst_dir, dst_dir, RowDiffStage::CONVERT);

    const std::string dest_fname = dst_dir/(std::string("ACGTCAG") + RowDiffColumnAnnotator::kExtension);
    ASSERT_TRUE(!std::filesystem::exists(dest_fname));
    std::filesystem::remove_all(dst_dir);
}

TEST(RowDiff, ConvertFromColumnCompressedSameLabels) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/"row_diff_same_labels";
    const std::string graph_fname
            = dst_dir/(std::string("ACGTCAG") + graph::DBGSuccinct::kExtension);
    const std::string annot_fname
            = dst_dir/(std::string("ACGTCAG") + ColumnCompressed<>::kExtension);
    const std::string dest_fname
            = dst_dir/(std::string("ACGTCAG") + RowDiffColumnAnnotator::kExtension);

    std::vector<std::vector<std::string>> label_groups = {
            { "Label0" }, { "Label1", "Label2" }, { "Label0", "Label2", "Label1" }
    };

    const uint32_t expected_relations[] = { 5, 2, 1, 1, 1 };
    for (const auto &labels : label_groups) {
        for (const uint32_t max_depth : { 1, 2, 3, 4, 5 }) {
            std::filesystem::remove_all(dst_dir);
            std::filesystem::create_directories(dst_dir);

            std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAC" });
            graph->serialize(graph_fname);

            ColumnCompressed source_annot(5);
            source_annot.add_labels({ 0, 1, 2, 3, 4 }, labels);
            source_annot.serialize(annot_fname);

            convert_to_row_diff({ annot_fname }, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::COMPUTE_REDUCTION);
            convert_to_row_diff({ annot_fname }, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::CONVERT);

            ASSERT_TRUE(std::filesystem::exists(dest_fname));
            RowDiffColumnAnnotator annotator;
            annotator.load(dest_fname);
            const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                    .set_graph(graph.get());
            const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                    .load_anchor(graph_fname + binmat::kRowDiffAnchorExt);

            ASSERT_EQ(labels.size(), annotator.num_labels());
            ASSERT_EQ(5u, annotator.num_objects());
            EXPECT_EQ(labels.size() * expected_relations[max_depth - 1],
                      annotator.num_relations());

            for (uint32 i = 0; i < annotator.num_objects(); ++i) {
                ASSERT_THAT(annotator.get(i), ContainerEq(labels));
            }
        }
    }
    std::filesystem::remove_all(dst_dir);
}

TEST(RowDiff, ConvertFromColumnCompressedSameLabelsMultipleColumns) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/"row_diff_multi_col";
    const std::string graph_fname
            = dst_dir/(std::string("graph") + graph::DBGSuccinct::kExtension);

    std::vector<std::vector<std::string>> label_groups
            = { { "Label0" }, { "Label1", "Label2" }, { "Label0", "Label2", "Label1" } };

    for (const auto &labels : label_groups) {
        for (const uint32_t max_depth : { 1, 2, 3, 4, 5 }) {
            std::vector<std::string> annot_fnames;
            std::filesystem::remove_all(dst_dir);
            std::filesystem::create_directories(dst_dir);

            std::unique_ptr<graph::DBGSuccinct> graph = create_graph(3, { "ACGTCAC" });
            graph->serialize(graph_fname);

            const uint32_t expected_relations[] = { 5, 2, 1, 1, 1 };

            std::vector<std::string> sources;
            for (const std::string &label : labels) {
                ColumnCompressed source_annot(5);
                source_annot.add_labels({ 0, 1, 2, 3, 4 }, { label });
                const std::string annot_fname
                        = dst_dir/(label + ColumnCompressed<>::kExtension);
                source_annot.serialize(annot_fname);
                annot_fnames.push_back(annot_fname);
            }

            convert_to_row_diff(annot_fnames, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::COMPUTE_REDUCTION);
            convert_to_row_diff(annot_fnames, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::CONVERT);

            for (uint32_t i = 0; i < labels.size(); ++i) {
                std::string rd_anno = dst_dir/(labels[i] + RowDiffColumnAnnotator::kExtension);
                ASSERT_TRUE(std::filesystem::exists(rd_anno));
                RowDiffColumnAnnotator annotator;
                annotator.load(rd_anno);
                const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                        .set_graph(graph.get());
                const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                        .load_anchor(graph_fname + binmat::kRowDiffAnchorExt);

                ASSERT_EQ(1, annotator.num_labels());
                ASSERT_EQ(5u, annotator.num_objects());
                EXPECT_EQ(expected_relations[max_depth - 1], annotator.num_relations());

                for (uint32 idx = 0; idx < annotator.num_objects(); ++idx) {
                    ASSERT_THAT(annotator.get(idx), ElementsAre(labels[i]));
                }
            }
        }
    }
    std::filesystem::remove_all(dst_dir);
}

void test_row_diff(uint32_t k,
                   uint32_t max_depth,
                   const std::vector<std::string> &sequences,
                   const std::vector<std::vector<std::string>> &annotations,
                   const std::string &prefix) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/prefix;
    const std::string graph_fname
            = dst_dir/(std::string("graph") + graph::DBGSuccinct::kExtension);
    const std::string annot_fname
            = dst_dir/(std::string("anno") + ColumnCompressed<>::kExtension);
    const std::string dest_fname
            = dst_dir/(std::string("anno") + RowDiffColumnAnnotator::kExtension);


    std::filesystem::remove_all(dst_dir);
    std::filesystem::create_directories(dst_dir);

    auto graph = std::make_unique<graph::DBGSuccinct>(k);
    for (const auto& seq : sequences) {
        graph->add_sequence(seq);
    }
    graph->mask_dummy_kmers(1, false);
    graph->serialize(graph_fname);

    ColumnCompressed initial_annotation(graph->num_nodes());
    std::unordered_set<std::string> all_labels;
    for (uint32_t anno_idx = 0; anno_idx < graph->num_nodes(); ++anno_idx) {
        const std::vector<std::string> &labels = annotations[anno_idx];
        initial_annotation.add_labels({anno_idx}, labels);
        std::for_each(labels.begin(), labels.end(), [&](auto l) { all_labels.insert(l); });
    }

    initial_annotation.serialize(annot_fname);

    convert_to_row_diff({ annot_fname }, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::COMPUTE_REDUCTION);
    convert_to_row_diff({ annot_fname }, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::CONVERT);

    ASSERT_TRUE(std::filesystem::exists(dest_fname));
    RowDiffColumnAnnotator annotator;
    annotator.load(dest_fname);
    const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
            .set_graph(graph.get());
    const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
            .load_anchor(graph_fname + binmat::kRowDiffAnchorExt);

    ASSERT_EQ(all_labels.size(), annotator.num_labels());
    ASSERT_EQ(graph->num_nodes(), annotator.num_objects());

    for (uint32_t anno_idx = 0; anno_idx < graph->num_nodes(); ++anno_idx) {
        ASSERT_THAT(annotator.get(anno_idx),
                    UnorderedElementsAreArray(annotations[anno_idx]));
    }

    std::filesystem::remove_all(dst_dir);
}

void test_row_diff_separate_columns(uint32_t k,
                                    uint32_t max_depth,
                                    const std::vector<std::string> &sequences,
                                    const std::vector<std::vector<std::string>> &annotations,
                                    const std::string &prefix) {
    const auto dst_dir = std::filesystem::path(test_dump_basename)/prefix;
    const std::string graph_fname
            = dst_dir/(std::string("graph") + graph::DBGSuccinct::kExtension);
    std::vector<std::string> annot_fnames;

    std::filesystem::remove_all(dst_dir);
    std::filesystem::create_directories(dst_dir);

    auto graph = std::make_unique<graph::DBGSuccinct>(k);
    for (const auto& seq : sequences) {
        graph->add_sequence(seq);
    }
    graph->mask_dummy_kmers(1, false);
    graph->serialize(graph_fname);

    std::map<std::string, std::vector<uint64_t>> col_annotations;
    for (uint32_t anno_idx = 0; anno_idx < graph->num_nodes(); ++anno_idx) {
        for(const auto& label : annotations[anno_idx]) {
            col_annotations[label].push_back(anno_idx);
        }
    }

    for (const auto& [label, indices] : col_annotations) {
        ColumnCompressed initial_annotation(graph->num_nodes());
        initial_annotation.add_labels(indices, {label});
        std::string annot_fname
                = dst_dir/("anno_" + label + ColumnCompressed<>::kExtension);
        annot_fnames.push_back(annot_fname);
        initial_annotation.serialize(annot_fname);
    }

    convert_to_row_diff(annot_fnames, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::COMPUTE_REDUCTION);
    convert_to_row_diff(annot_fnames, graph_fname, 1e9, max_depth, dst_dir, dst_dir, RowDiffStage::CONVERT);

    for (const auto& [label, indices] : col_annotations) {
        const std::string dest_fname
                = dst_dir/("anno_" + label + RowDiffColumnAnnotator::kExtension);
        ASSERT_TRUE(std::filesystem::exists(dest_fname));
        RowDiffColumnAnnotator annotator;
        annotator.load(dest_fname);
        const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                .set_graph(graph.get());
        const_cast<binmat::RowDiff<binmat::ColumnMajor> &>(annotator.get_matrix())
                .load_anchor(graph_fname + binmat::kRowDiffAnchorExt);

        ASSERT_EQ(graph->num_nodes(), annotator.num_objects());

        std::vector<uint64_t> actual_indices;
        annotator.call_objects(label,
                               [&](uint64_t index) { actual_indices.push_back(index); });
        ASSERT_THAT(indices, UnorderedElementsAreArray(indices));
    }

    std::filesystem::remove_all(dst_dir);
}

TEST(RowDiff, ConvertFromColumnCompressedOneToFiveLabels) {
    const std::vector<std::vector<std::string>> annotations
            = { { "Label0" },
                { "Label1", "Label2" },
                { "Label1", "Label2", "Label3" },
                { "Label1", "Label2", "Label3", "Label4" },
                { "Label1", "Label2", "Label3", "Label4", "Label5" } };
    test_row_diff(3, 5, { "ACGTCAC" }, annotations, "column.diff.convert");
    test_row_diff_separate_columns(3, 5, { "ACGTCAC" }, annotations,
                                   "column.diff.convert");
}

TEST(RowDiff, ConvertFromColumnCompressed) {
    const std::vector<std::vector<std::string>> annotations
            = { { "Label0", "Label2", "Label8" },
                {},
                { "Label1", "Label2" },
                { "Label1", "Label2", "Label8" },
                { "Label2" } };
    test_row_diff(3, 5, { "ACGTCAC" }, annotations, "column.diff.convert");
    test_row_diff_separate_columns(3, 5, { "ACGTCAC" }, annotations,
                                   "column.diff.convert");
}

TEST(RowDiff, ConvertFromColumnCompressed4Loops) {
    std::vector<std::string> sequences = { std::string(100, 'A'), std::string(100, 'G'),
                                           std::string(100, 'C'), std::string(100, 'T') };
    const std::vector<std::vector<std::string>> annotations
            = { { "Lb0", "Lb2", "Lb8" }, {}, { "Lb9" }, { "Lb0", "Lb8", "Lb2" } };

    test_row_diff(3, 5, sequences, annotations, "column.diff.4loops");
    test_row_diff_separate_columns(3, 5, sequences, annotations, "column.diff.4loops");
}

TEST(RowDiff, ConvertFromColumnCompressed4PathsRandomLabels) {
    std::mt19937 gen(12345);
    std::uniform_int_distribution<> distrib(0, 2);


    std::vector<std::vector<std::string>> options = { { "L1" }, { "L1", "L2" }, { "L3" } };
    std::vector<std::vector<std::string>> annotations;
    for (uint32_t anno_idx = 0; anno_idx < 24; ++anno_idx) {
        annotations.push_back(options[distrib(gen)]);
    }
    const std::vector<std::string> sequences = { "ATCGGAAGA", "TTTAAACCCGGG", "ATACAGCTCGCT", "AAAAAA" };
    test_row_diff(3, 5, sequences, annotations, "column.diff.4paths");
    test_row_diff_separate_columns(3, 5, sequences, annotations, "column.diff.4paths");
}

TEST(RowDiff, ConvertFromColumnCompressed2BigLoops) {
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

    test_row_diff(10, 3, sequences, annotations, "column.diff.2bigloops");
    test_row_diff_separate_columns(10, 3, sequences, annotations, "column.diff.2bigloops");
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
    annotation = convert_to_greedy_BRWT(std::move(*initial_annotation));
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
    merged_annotation->load(test_dump_basename_row_compressed_merge + "_merged");
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
    merged_annotation->load(filename);
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
    merged_annotation->load(filename);
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
    merged_annotation->load(outfile);
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
