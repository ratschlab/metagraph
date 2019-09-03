#include <random>

#include <filesystem>

#include "gtest/gtest.h"

#include "test_annotation.hpp"
#include "../test_helpers.hpp"

#define private public
#include "test_matrix_helpers.hpp"
#include "unix_tools.hpp"


TYPED_TEST(AnnotatorTest, EmptyConstructor) {
    EXPECT_EQ(0u, this->annotation->num_labels());
    EXPECT_EQ(0u, this->annotation->num_objects());
    EXPECT_EQ(0u, this->annotation->num_relations());
}

TYPED_TEST(AnnotatorDynamicTest, EmptyRows) {
    this->annotation.reset(new TypeParam(5));
    EXPECT_EQ(0u, this->annotation->num_labels());
    EXPECT_EQ(5u, this->annotation->num_objects());
    EXPECT_EQ(0u, this->annotation->num_relations());

    EXPECT_EQ(0u, this->annotation->get(0).size());
    EXPECT_EQ(0u, this->annotation->get(1).size());
    EXPECT_EQ(0u, this->annotation->get(2).size());
    EXPECT_EQ(0u, this->annotation->get(3).size());
    EXPECT_EQ(0u, this->annotation->get(4).size());
}

TYPED_TEST(AnnotatorPresetTest, GetLabels) {
    EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(std::vector<std::string>({}),
              this->annotation->get_labels(1));
    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorTest, GetLabelsOneRow) {
    annotate::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels(0, {"Label 0", "Label 1", "Label 2", "Label 3", "Label 4"});
    this->set(std::move(column_annotator));
    EXPECT_EQ(convert_to_set({"Label 0", "Label 1", "Label 2", "Label 3", "Label 4"}),
              convert_to_set(this->annotation->get_labels(0)));
}

TYPED_TEST(AnnotatorPresetTest, NumObjects) {
    EXPECT_EQ(5u, this->annotation->num_objects());
}

TYPED_TEST(AnnotatorStaticTest, NumLabelsOneEmptyRow) {
    annotate::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels(0, {});
    this->set(std::move(column_annotator));
    EXPECT_EQ(0u, this->annotation->num_labels());
    EXPECT_EQ(0u, this->annotation->num_objects());
    EXPECT_EQ(0u, this->annotation->num_relations());
}

TYPED_TEST(AnnotatorTest, NumLabelsOneRow) {
    annotate::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels(0, {"Label 0", "Label 1", "Label 2",
                                    "Label 3", "Label 4"});
    this->set(std::move(column_annotator));
    EXPECT_EQ(5u, this->annotation->num_labels());
    EXPECT_EQ(5u, this->annotation->num_objects());
    EXPECT_EQ(5u, this->annotation->num_relations());
}

TYPED_TEST(AnnotatorPresetTest, NumLabels) {
    EXPECT_EQ(4u, this->annotation->num_labels());
}

TYPED_TEST(AnnotatorPresetTest, HasLabel) {
    EXPECT_TRUE(this->annotation->has_label(0, "Label0"));
    EXPECT_TRUE(this->annotation->has_label(0, "Label2"));
    EXPECT_TRUE(this->annotation->has_label(0, "Label8"));
    EXPECT_FALSE(this->annotation->has_label(0, "Label1"));

    EXPECT_FALSE(this->annotation->has_label(1, "Label0"));

    EXPECT_TRUE(this->annotation->has_label(2, "Label1"));
    EXPECT_TRUE(this->annotation->has_label(2, "Label2"));
    EXPECT_FALSE(this->annotation->has_label(2, "Label0"));

    EXPECT_TRUE(this->annotation->has_label(3, "Label1"));
    EXPECT_TRUE(this->annotation->has_label(3, "Label2"));
    EXPECT_TRUE(this->annotation->has_label(3, "Label8"));
    EXPECT_FALSE(this->annotation->has_label(3, "Label0"));

    EXPECT_TRUE(this->annotation->has_label(4, "Label2"))
        << this->annotation->get_labels(4)[0];
    EXPECT_FALSE(this->annotation->has_label(4, "Label0"));
    EXPECT_FALSE(this->annotation->has_label(4, "Label1"))
        << this->annotation->get_labels(4)[0];
    EXPECT_FALSE(this->annotation->has_label(4, "Label8"));
}

TYPED_TEST(AnnotatorPresetTest, HasLabels) {
    EXPECT_TRUE(this->annotation->has_labels(0, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(this->annotation->has_labels(0, { "Label0", "Label2", "Label8", "Label1" }));

    EXPECT_FALSE(this->annotation->has_labels(1, { "Label0" }));

    EXPECT_TRUE(this->annotation->has_labels(2, { "Label1", "Label2" }));
    EXPECT_FALSE(this->annotation->has_labels(2, { "Label1", "Label2", "Label0" }));

    EXPECT_TRUE(this->annotation->has_labels(3, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(this->annotation->has_labels(3, { "Label1", "Label2", "Label8", "Label0" }));

    EXPECT_TRUE(this->annotation->has_labels(4, { "Label2" }))
        << this->annotation->get_labels(4)[0];
    EXPECT_FALSE(this->annotation->has_labels(4, { "Label2", "Label0" }));
    EXPECT_FALSE(this->annotation->has_labels(4, { "Label2", "Label0", "Label1" }));
    EXPECT_FALSE(this->annotation->has_labels(4, { "Label2", "Label0", "Label1", "Label8" }));
}


TYPED_TEST(AnnotatorPreset2Test, has_labels2) {
    EXPECT_FALSE(this->annotation->has_labels(0, { "Label0", "Label1",
                                                   "Label2", "Label4",
                                                   "Label5", "Label8" }));

    EXPECT_FALSE(this->annotation->has_labels(0, { "Label0",
                                                   "Label2", "Label4",
                                                   "Label8" }));

    EXPECT_TRUE(this->annotation->has_labels(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(this->annotation->has_labels(0, { "Label0", "Label8" }));
    EXPECT_TRUE(this->annotation->has_labels(0, { "Label2" }));
    EXPECT_TRUE(this->annotation->has_labels(0, {}));

    EXPECT_FALSE(this->annotation->has_labels(1, { "Label0", "Label1",
                                                   "Label2", "Label4",
                                                   "Label5", "Label8" }));

    EXPECT_FALSE(this->annotation->has_labels(1, { "Label0",
                                                   "Label2", "Label4",
                                                   "Label8" }));

    EXPECT_FALSE(this->annotation->has_labels(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(this->annotation->has_labels(1, { "Label0", "Label8" }));
    EXPECT_FALSE(this->annotation->has_labels(1, { "Label2" }));
    EXPECT_TRUE(this->annotation->has_labels(1, {}));

    EXPECT_FALSE(this->annotation->has_labels(2, { "Label0", "Label1",
                                                   "Label2", "Label4",
                                                   "Label5", "Label8" }));

    EXPECT_FALSE(this->annotation->has_labels(2, { "Label0",
                                                   "Label2", "Label4",
                                                   "Label8" }));

    EXPECT_FALSE(this->annotation->has_labels(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(this->annotation->has_labels(2, { "Label1", "Label8" }));
    EXPECT_TRUE(this->annotation->has_labels(2, { "Label1", "Label2" }));
    EXPECT_TRUE(this->annotation->has_labels(2, { "Label2" }));
    EXPECT_TRUE(this->annotation->has_labels(2, {}));
}

TYPED_TEST(AnnotatorTest, LabelExists) {
    EXPECT_FALSE(this->annotation->label_exists("Label0"));
    EXPECT_FALSE(this->annotation->label_exists("Label1"));
    EXPECT_FALSE(this->annotation->label_exists("Label2"));
    EXPECT_FALSE(this->annotation->label_exists("Label3"));
    EXPECT_FALSE(this->annotation->label_exists("Label8"));
}

TYPED_TEST(AnnotatorPresetTest, LabelExists) {
    EXPECT_TRUE(this->annotation->label_exists("Label0"));
    EXPECT_TRUE(this->annotation->label_exists("Label1"));
    EXPECT_TRUE(this->annotation->label_exists("Label2"));
    EXPECT_FALSE(this->annotation->label_exists("Label3"));
    EXPECT_TRUE(this->annotation->label_exists("Label8"));
}

TYPED_TEST(AnnotatorPreset2Test, LabelExists) {
    EXPECT_TRUE(this->annotation->label_exists("Label0"));
    EXPECT_TRUE(this->annotation->label_exists("Label1"));
    EXPECT_TRUE(this->annotation->label_exists("Label2"));
    EXPECT_FALSE(this->annotation->label_exists("Label3"));
    EXPECT_TRUE(this->annotation->label_exists("Label8"));
}

TYPED_TEST(AnnotatorPreset3Test, LabelExists) {
    EXPECT_TRUE(this->annotation->label_exists("Label0"));
    EXPECT_TRUE(this->annotation->label_exists("Label1"));
    EXPECT_TRUE(this->annotation->label_exists("Label2"));
    EXPECT_FALSE(this->annotation->label_exists("Label3"));
    EXPECT_TRUE(this->annotation->label_exists("Label8"));
}

TYPED_TEST(AnnotatorPresetTest, Sparsity) {
    EXPECT_EQ(1 - static_cast<double>(this->annotation->num_relations())
                                          / this->annotation->num_objects()
                                          / this->annotation->num_labels(), 0.55);
}

TYPED_TEST(AnnotatorPresetTest, CallIndices) {
    std::vector<typename TypeParam::Index> indices;

    this->annotation->call_objects("Label0",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(1u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label1",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(2u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({2u, 3u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label2",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(4u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 2u, 3u, 4u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label3",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(0u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label8",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(2u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 3u}),
              convert_to_set(indices));
    indices.clear();
}

TYPED_TEST(AnnotatorPreset2Test, CallIndices) {
    std::vector<typename TypeParam::Index> indices;

    this->annotation->call_objects("Label0",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(1u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label1",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(1u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({2u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label2",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(2u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 2u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label3",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(0u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label8",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(2u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 4u}),
              convert_to_set(indices));
    indices.clear();
}

TYPED_TEST(AnnotatorPreset3Test, CallIndices) {
    std::vector<typename TypeParam::Index> indices;

    this->annotation->call_objects("Label0",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(1u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label1",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(2u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({2u, 3u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label2",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(4u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 2u, 3u, 4u}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label3",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(0u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({}),
              convert_to_set(indices));
    indices.clear();

    this->annotation->call_objects("Label8",
                                   [&](const auto &i) { indices.push_back(i); });

    EXPECT_EQ(3u, indices.size());
    EXPECT_EQ(convert_to_set<typename TypeParam::Index>({0u, 3u, 4u}),
              convert_to_set(indices));
    indices.clear();
}

TYPED_TEST(AnnotatorDynamicTest, add_label) {
    TypeParam annotation(5);
    annotation.add_label(0, "0");
    annotation.add_label(1, "0");
    annotation.add_label(2, "1");
    annotation.add_label(1, "2");

    EXPECT_EQ(convert_to_set({"0"}), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({"0", "2"}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({"1"}), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(4)));
}

TYPED_TEST(AnnotatorDynamicTest, add_label_long_sparse) {
    TypeParam annotation(200'000'000);
    std::vector<uint64_t> indices(5'000);
    std::iota(indices.begin(), indices.end(), 0);
    annotation.add_labels(indices, { "0" });

    std::iota(indices.begin(), indices.end(), 4'999);
    annotation.add_labels(indices, { "1" });

    std::iota(indices.begin(), indices.end(), 1'000'000);
    annotation.add_labels(indices, { "2" });

    EXPECT_EQ(convert_to_set({"0"}), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({"0", "1"}), convert_to_set(annotation.get(4'999)));
    EXPECT_EQ(convert_to_set({"2"}), convert_to_set(annotation.get(1'000'001)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1'500'000)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1'900'000)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_first_column_to_empty_annotation) {
    this->annotation.reset(new TypeParam(0));

    ASSERT_EQ(0u, this->annotation->num_labels());

    this->annotation->insert_rows({ 0, });
    ASSERT_EQ(0u, this->annotation->num_labels());

    this->annotation->set_labels(0, { "Label0", "Label2", "Label8" });

    ASSERT_EQ(3u, this->annotation->num_labels());
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
}

TYPED_TEST(AnnotatorDynamicTest, add_label_sequential) {
    size_t graph_half_size = 1000;
    TypeParam annotation(graph_half_size * 2);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_label(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TYPED_TEST(AnnotatorDynamicTest, add_label_random) {
    size_t graph_half_size = 1000;
    TypeParam annotation(graph_half_size * 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, labels[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TYPED_TEST(AnnotatorDynamicTest, set_labels) {
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_no_empty_rows) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({});

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_one_empty_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({ 4, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_first_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({ 0, });

    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_last_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({ 5, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_empty_rows) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({ 1, 2, 4, 5 });
    EXPECT_EQ(9u, this->annotation->num_objects());

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(6)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(7)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(8)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_empty_rows_many) {
    this->annotation.reset(new TypeParam(2'000'000));
    EXPECT_EQ(2'000'000u, this->annotation->num_objects());
    this->annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    this->annotation->set_labels(2, { "Label1", "Label2" });
    this->annotation->set_labels(4, { "Label8" });
    std::vector<uint64_t> indices(1'000'000);
    std::iota(indices.begin(), indices.end(), 5);
    this->annotation->add_labels(indices, { "Label9" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(4)));

    this->annotation->insert_rows({ 1, 2, 4, 5 });
    EXPECT_EQ(2'000'004u, this->annotation->num_objects());

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get_labels(6)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(7)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get_labels(8)));


    EXPECT_EQ(convert_to_set({ "Label9" }), convert_to_set(this->annotation->get_labels(9)));
    EXPECT_EQ(convert_to_set({ "Label9" }), convert_to_set(this->annotation->get_labels(1'000'000)));
}

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoad) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoadExtension1) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good + this->annotation->file_extension());
    this->annotation.reset(new TypeParam());
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad + this->annotation->file_extension()));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good + this->annotation->file_extension()));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoadExtension2) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good + this->annotation->file_extension());
    this->annotation.reset(new TypeParam());
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoadExtension3) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    this->annotation.reset(new TypeParam());
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad + this->annotation->file_extension()));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good + this->annotation->file_extension()));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SerializationAndLoad) {
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SerializationAndLoadExtension1) {
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SerializationAndLoadExtension2) {
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SerializationAndLoadExtension3) {
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadText) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, false));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    uint64_t size, num_set_bits, pos;
    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".text.annodbg");
        ASSERT_TRUE(fin.good());

        fin >> size >> num_set_bits;

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            fin >> pos;
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadTextParallel) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, false, 3));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    uint64_t size, num_set_bits, pos;
    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".text.annodbg");
        ASSERT_TRUE(fin.good());

        fin >> size >> num_set_bits;

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            fin >> pos;
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadBinary) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, true));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".raw.annodbg",
                          std::ios::binary);
        ASSERT_TRUE(fin.good());

        uint64_t size = load_number(fin);
        uint64_t num_set_bits = load_number(fin);

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            size_t pos = load_number(fin);
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadBinaryParallel) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, true, 3));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".raw.annodbg",
                          std::ios::binary);
        ASSERT_TRUE(fin.good());

        uint64_t size = load_number(fin);
        uint64_t num_set_bits = load_number(fin);

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            size_t pos = load_number(fin);
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, MergeLoadDisjoint) {
    std::filesystem::remove(test_dump_basename_vec_good + "_1");
    std::filesystem::remove(test_dump_basename_vec_good + "_1" + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good + "_2");
    std::filesystem::remove(test_dump_basename_vec_good + "_2" + this->annotation->file_extension());

    {
        TypeParam annotation(5);
        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label2" });
        annotation.set_labels(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_1");
    }
    {
        TypeParam annotation(5);
        annotation.set_labels(1, { "2_Label0", "2_Label2", "2_Label8" });
        annotation.set_labels(2, { "2_Label1", "2_Label9", "2_Label0" });
        annotation.set_labels(3, { "2_Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_2");
    }
    {
        TypeParam annotation(0);
        ASSERT_TRUE(annotation.merge_load({ test_dump_basename_vec_good + "_1",
                                            test_dump_basename_vec_good + "_2" }));

        EXPECT_EQ(5u, annotation.num_objects());
        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({ "2_Label0", "2_Label2", "2_Label8" }), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "2_Label1", "Label2", "2_Label9", "2_Label0" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({ "2_Label8" }), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, MergeLoad) {
    std::filesystem::remove(test_dump_basename_vec_good + "_1");
    std::filesystem::remove(test_dump_basename_vec_good + "_1" + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good + "_2");
    std::filesystem::remove(test_dump_basename_vec_good + "_2" + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good + "_1");

    this->annotation.reset(new TypeParam(5));
    this->annotation->set_labels(1, { "Label0", "Label2", "Label8" });
    this->annotation->set_labels(2, { "Label1", "Label9", "Label0" });
    this->annotation->set_labels(3, { "Label8" });
    this->annotation->serialize(test_dump_basename_vec_good + "_2");

    this->annotation.reset(new TypeParam(0));
    ASSERT_TRUE(this->annotation->merge_load({ test_dump_basename_vec_good + "_1",
                                               test_dump_basename_vec_good + "_2" }));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }),
              convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }),
              convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label9", "Label0" }),
              convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label8" }),
              convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),
              convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, NoRenameColumns) {
    this->annotation->rename_labels({});

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, RenameColumns) {
    this->annotation->rename_labels({ { "Label2", "Label2Renamed" },
                                      { "Label8", "Label8Renamed" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2Renamed", "Label8Renamed" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2Renamed" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8Renamed" }), convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SwapColumns) {
    this->annotation->rename_labels({ { "Label2", "Label8" },
                                      { "Label8", "Label2" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Label8", "Label2" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label8" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }), convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorStaticTest, RenameColumnsMerge) {
    annotate::ColumnCompressed<> column_annotator(5);
    column_annotator.set_labels(0, { "Label0", "Label2", "Label8" });
    column_annotator.set_labels(2, { "Label1", "Label2" });
    column_annotator.set_labels(4, { "Label8" });

    this->set(std::move(column_annotator));
    ASSERT_DEATH(
        this->annotation->rename_labels({ { "Label2", "Merged" },
                                          { "Label8", "Merged" } }),
        ""
    );
}

TYPED_TEST(AnnotatorStaticTest, RenameColumnsMergeAll) {
    annotate::ColumnCompressed<> column_annotator(5);
    column_annotator.set_labels(0, { "Label0", "Label2", "Label8" });
    column_annotator.set_labels(2, { "Label1", "Label2" });
    column_annotator.set_labels(4, { "Label8" });

    this->set(std::move(column_annotator));
    ASSERT_DEATH(
        this->annotation->rename_labels({ { "Label0", "Merged" },
                                          { "Label1", "Merged" },
                                          { "Label2", "Merged" },
                                          { "Label8", "Merged" }, }),
        ""
    );
}

TYPED_TEST(AnnotatorStaticLargeTest, CheckCache) {
    size_t num_rows = 20000;
    size_t num_columns = 200;
    BitVectorPtrArray columns, copy;
    annotate::LabelEncoder label_encoder;

    for (size_t j = 0; j < num_columns; ++j) {

        columns.emplace_back(new bit_vector_stat(num_rows));

        for (size_t i = 0; i < num_rows; ++i) {
            columns.back()->set(i, (i + 2 * j) % 1000);
        }
        copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));

        label_encoder.insert_and_encode(std::to_string(j));
    }

    auto annotator = TypeParam(
        std::make_unique<typename TypeParam::binary_matrix_type>(
            build_matrix_from_columns<typename TypeParam::binary_matrix_type>(
                std::move(copy), num_rows
            )
        ),
        label_encoder,
        1000000
    );

    std::vector<std::vector<std::string>> rows;
    for (size_t i = 0; i < num_rows; i += 1000) {
        rows.emplace_back(annotator.get_labels(i));
    }

    auto it = rows.begin();
    for (size_t i = 0; i < num_rows; i += 1000) {
        ASSERT_NE(rows.end(), it);
        EXPECT_EQ(*it++, annotator.get_labels(i));
    }
}

// This can be run with --gtest_also_run_disabled_tests
TYPED_TEST(AnnotatorStaticLargeTest, DISABLED_QueryRowsCached_LONG_TEST) {
    size_t num_rows = 2000000;
    size_t num_columns = 200;
    BitVectorPtrArray columns, copy;
    annotate::LabelEncoder label_encoder;

    for (size_t j = 0; j < num_columns; ++j) {

        columns.emplace_back(new bit_vector_stat(num_rows));

        for (size_t i = 0; i < num_rows; ++i) {
            columns.back()->set(i, (i + 2 * j) % 1000);
        }
        copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));

        label_encoder.insert_and_encode(std::to_string(j));
    }

    auto annotator = TypeParam(
        std::make_unique<typename TypeParam::binary_matrix_type>(
            build_matrix_from_columns<typename TypeParam::binary_matrix_type>(
                std::move(copy), num_rows
            )
        ),
        label_encoder
    );

    for (size_t cache_size : { 0, 1000, 10000, 100000, 1000000, 10000000 }) {
        annotator.reset_row_cache(cache_size);
        Timer timer;
        for (size_t i = 0; i < num_rows; ++i) {
            annotator.get_labels(i);
        }
        TEST_COUT << "Query all rows\t"
                  << "Cache size:\t" << cache_size << "\t\t"
                  << "Time:\t" << timer.elapsed();
    }

    for (size_t cache_size : { 0, 1000, 10000, 100000, 1000000, 10000000 }) {
        annotator.reset_row_cache(cache_size);
        Timer timer;
        for (size_t j = 0; j < 1000; ++j) {
            for (size_t i = 0; i < num_rows; i += 1000) {
                annotator.get_labels(i);
            }
        }
        TEST_COUT << "Query some rows repeatedly\t"
                  << "Cache size:\t" << cache_size << "\t\t"
                  << "Time:\t" << timer.elapsed();
    }
}

