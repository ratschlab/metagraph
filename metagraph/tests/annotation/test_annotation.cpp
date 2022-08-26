#include "gtest/gtest.h"

#include "test_annotation.hpp"
#include "../test_helpers.hpp"

#define private public
#include "test_matrix_helpers.hpp"
#include "common/unix_tools.hpp"
#include "common/data_generation.hpp"


namespace mtg {
namespace test {

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
              convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(std::vector<std::string>({}),
              this->annotation->get(1));
    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorPresetTest, CountLabels) {
    EXPECT_EQ(
        convert_to_set(std::vector<std::pair<uint64_t, size_t>>({
            {0, 1}, {3, 2}, {1, 4}, {2, 2}
        })),
        convert_to_set(this->annotation->get_matrix().sum_rows(
            std::vector<std::pair<uint64_t, size_t>>({
                {0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}
            })
        ))
    );

    EXPECT_EQ(
        convert_to_set(std::vector<std::pair<uint64_t, size_t>>({
            {0, 1}, {3, 2}, {1, 4}, {2, 2}
        })),
        convert_to_set(this->annotation->get_matrix().sum_rows(
            std::vector<std::pair<uint64_t, size_t>>({
                {0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}
            }), 0
        ))
    );

    EXPECT_EQ(
        convert_to_set(std::vector<std::pair<uint64_t, size_t>>({
            {0, 1}, {3, 2}, {1, 2}, {2, 2}
        })),
        convert_to_set(this->annotation->get_matrix().sum_rows(
            std::vector<std::pair<uint64_t, size_t>>({
                {0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}
            }), 0, 2
        ))
    );

    EXPECT_EQ(
        convert_to_set(std::vector<std::pair<uint64_t, size_t>>({
            {3, 2}, {1, 2}, {2, 2}
        })),
        convert_to_set(this->annotation->get_matrix().sum_rows(
            std::vector<std::pair<uint64_t, size_t>>({
                {0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}
            }), 2, 2
        ))
    );

    EXPECT_EQ(
        convert_to_set(std::vector<std::pair<uint64_t, size_t>>({})),
        convert_to_set(this->annotation->get_matrix().sum_rows(
            std::vector<std::pair<uint64_t, size_t>>({
                {0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}
            }), 0, 0
        ))
    );
}

TYPED_TEST(AnnotatorTest, GetLabelsOneRow) {
    annot::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels({ 0 }, {"Label 0", "Label 1", "Label 2", "Label 3", "Label 4"});
    this->set(std::move(column_annotator));
    EXPECT_EQ(convert_to_set({"Label 0", "Label 1", "Label 2", "Label 3", "Label 4"}),
              convert_to_set(this->annotation->get(0)));
}

TYPED_TEST(AnnotatorPresetTest, NumObjects) {
    EXPECT_EQ(5u, this->annotation->num_objects());
}

TYPED_TEST(AnnotatorStaticTest, NumLabelsOneEmptyRow) {
    annot::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels({ 0 }, {});
    this->set(std::move(column_annotator));
    EXPECT_EQ(0u, this->annotation->num_labels());
    EXPECT_EQ(0u, this->annotation->num_objects());
    EXPECT_EQ(0u, this->annotation->num_relations());
}

TYPED_TEST(AnnotatorTest, NumLabelsOneRow) {
    annot::ColumnCompressed<> column_annotator(5);
    column_annotator.add_labels({0}, {"Label 0", "Label 1", "Label 2",
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
        << this->annotation->get(4)[0];
    EXPECT_FALSE(this->annotation->has_label(4, "Label0"));
    EXPECT_FALSE(this->annotation->has_label(4, "Label1"))
        << this->annotation->get(4)[0];
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
        << this->annotation->get(4)[0];
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
    annotation.add_labels({ 0 }, { "0" });
    annotation.add_labels({ 1 }, { "0" });
    annotation.add_labels({ 2 }, { "1" });
    annotation.add_labels({ 1 }, { "2" });

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

    this->annotation->set(0, { "Label0", "Label2", "Label8" });

    ASSERT_EQ(3u, this->annotation->num_labels());
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
}

TYPED_TEST(AnnotatorDynamicTest, add_label_sequential) {
    size_t graph_half_size = 1000;
    TypeParam annotation(graph_half_size * 2);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_labels({ i }, { "Label1" });
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_labels({ i }, { "Label2" });
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TYPED_TEST(AnnotatorDynamicTest, add_label_random) {
    size_t graph_half_size = 1000;
    TypeParam annotation(graph_half_size * 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_labels({ i }, { labels[i % 2] });
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TYPED_TEST(AnnotatorDynamicTest, set) {
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_no_empty_rows) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({});

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_one_empty_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({ 4, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_first_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({ 0, });

    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_last_row) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({ 5, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_empty_rows) {
    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({ 1, 2, 4, 5 });
    EXPECT_EQ(9u, this->annotation->num_objects());

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(6)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(7)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(8)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, insert_empty_rows_many) {
    this->annotation.reset(new TypeParam(2'000'000));
    EXPECT_EQ(2'000'000u, this->annotation->num_objects());
    this->annotation->set(0, { "Label0", "Label2", "Label8" });
    this->annotation->set(2, { "Label1", "Label2" });
    this->annotation->set(4, { "Label8" });
    std::vector<uint64_t> indices(1'000'000);
    std::iota(indices.begin(), indices.end(), 5);
    this->annotation->add_labels(indices, { "Label9" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(4)));

    this->annotation->insert_rows({ 1, 2, 4, 5 });
    EXPECT_EQ(2'000'004u, this->annotation->num_objects());

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(5)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(this->annotation->get(6)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(this->annotation->get(7)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(this->annotation->get(8)));


    EXPECT_EQ(convert_to_set({ "Label9" }), convert_to_set(this->annotation->get(9)));
    EXPECT_EQ(convert_to_set({ "Label9" }), convert_to_set(this->annotation->get(1'000'000)));
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
    annot::ColumnCompressed<> column_annotator(5);
    column_annotator.set(0, { "Label0", "Label2", "Label8" });
    column_annotator.set(2, { "Label1", "Label2" });
    column_annotator.set(4, { "Label8" });

    this->set(std::move(column_annotator));
    ASSERT_DEATH_SILENT(
        this->annotation->rename_labels({ { "Label2", "Merged" },
                                          { "Label8", "Merged" } }),
        ""
    );
}

TYPED_TEST(AnnotatorStaticTest, RenameColumnsMergeAll) {
    annot::ColumnCompressed<> column_annotator(5);
    column_annotator.set(0, { "Label0", "Label2", "Label8" });
    column_annotator.set(2, { "Label1", "Label2" });
    column_annotator.set(4, { "Label8" });

    this->set(std::move(column_annotator));
    ASSERT_DEATH_SILENT(
        this->annotation->rename_labels({ { "Label0", "Merged" },
                                          { "Label1", "Merged" },
                                          { "Label2", "Merged" },
                                          { "Label8", "Merged" }, }),
        ""
    );
}

} // namespace test
} // namespace mtg
