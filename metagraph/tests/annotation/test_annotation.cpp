#include <random>

#include <filesystem>

#include "gtest/gtest.h"

#define protected public
#define private public

#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_annotation";


template <typename... Args>
class RowCompressedParallel : public annotate::RowCompressed<Args...> {
  public:
    template <typename... CArgs>
    RowCompressedParallel(CArgs&&... args)
          : annotate::RowCompressed<Args...>(std::forward<CArgs>(args)...) {}
};

template <typename... Args>
class RowCompressedDynamic : public annotate::RowCompressed<Args...> {
  public:
    template <typename... CArgs>
    RowCompressedDynamic(CArgs&&... args)
          : annotate::RowCompressed<Args...>(std::forward<CArgs>(args)...) {}
};

template <typename... Args>
class RowCompressedSparse : public annotate::RowCompressed<Args...> {
  public:
    RowCompressedSparse(uint64_t num_rows = 0)
          : annotate::RowCompressed<Args...>(num_rows, true) {}
};


template <typename Annotator>
class AnnotatorTest : public ::testing::Test {
  public:
    std::unique_ptr<Annotator> annotation;
    virtual void set(annotate::ColumnCompressed<>&& column_annotator) {
        annotation.reset(annotate::convert<Annotator>(
            std::move(column_annotator)
        ).release());
    }

    virtual void SetUp() { set(annotate::ColumnCompressed<>(0)); }
};

template<>
void AnnotatorTest<annotate::BRWTCompressed<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    annotation = annotate::convert_to_simple_BRWT<annotate::BRWTCompressed<>>(
        std::move(column_annotator)
    );
}

template<>
void AnnotatorTest<annotate::RowCompressed<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    annotation.reset(new annotate::RowCompressed<>(0));
    column_annotator.convert_to_row_annotator(annotation.get());
}

template<>
void AnnotatorTest<RowCompressedParallel<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    annotation.reset(new RowCompressedParallel<>(0));
    column_annotator.convert_to_row_annotator(annotation.get(), 10);
}

template<>
void AnnotatorTest<RowCompressedDynamic<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    annotation.reset(new RowCompressedDynamic<>(column_annotator.num_objects()));
    for (RowCompressedDynamic<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
        annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
    }
}

template<>
void AnnotatorTest<RowCompressedSparse<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    annotation.reset(new RowCompressedSparse<>(column_annotator.num_objects()));
    for (RowCompressedSparse<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
        annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
    }
}

template<>
void AnnotatorTest<annotate::ColumnCompressed<>>
::set(annotate::ColumnCompressed<>&& column_annotator) {
    // TODO: introduce move constructor for ColumnCompressed
    //annotation.reset(new annotate::ColumnCompressed<>(std::move(column_annotator)));
    annotation.reset(new annotate::ColumnCompressed<>(column_annotator.num_objects()));
    for (annotate::ColumnCompressed<>::Index i = 0; i < column_annotator.num_objects(); ++i) {
        annotation->add_labels(i, std::move(column_annotator.get_labels(i)));
    }
}

template <typename Annotator>
class AnnotatorStaticTest : public AnnotatorTest<Annotator> { };

template <typename Annotator>
class AnnotatorPresetTest : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(3, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels(4, {"Label2"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset2Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, { "Label0", "Label2", "Label8" });
        column_annotator.add_labels(2, { "Label1", "Label2" });
        column_annotator.add_labels(4, { "Label8" });
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorPreset3Test : public AnnotatorTest<Annotator> {
  public:
    virtual void SetUp() override {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(3, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels(4, {"Label2", "Label8"});
        this->set(std::move(column_annotator));
    }
};

template <typename Annotator>
class AnnotatorDynamicTest : public AnnotatorPreset2Test<Annotator> { };

template <typename Annotator>
class AnnotatorDynamicNoSparseTest : public AnnotatorPreset2Test<Annotator> { };


typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::BRWTCompressed<>,
                         annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorTypes;
typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::BRWTCompressed<>> AnnotatorStaticTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>,
                         RowCompressedSparse<>> AnnotatorDynamicTypes;
typedef ::testing::Types<annotate::ColumnCompressed<>,
                         annotate::RowCompressed<>,
                         RowCompressedParallel<>,
                         RowCompressedDynamic<>> AnnotatorDynamicNoSparseTypes;

TYPED_TEST_CASE(AnnotatorTest, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPresetTest, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPreset2Test, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorPreset3Test, AnnotatorTypes);
TYPED_TEST_CASE(AnnotatorStaticTest, AnnotatorStaticTypes);
TYPED_TEST_CASE(AnnotatorDynamicTest, AnnotatorDynamicTypes);
TYPED_TEST_CASE(AnnotatorDynamicNoSparseTest, AnnotatorDynamicNoSparseTypes);

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

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

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

TYPED_TEST(AnnotatorPresetTest, get_labels) {
    EXPECT_EQ(std::vector<std::string>({}),
              this->annotation->get_labels({}, 1));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(this->annotation->get_labels({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(this->annotation->get_labels({ 0, 1, 2, 3, 4 }, 1)));
}

TYPED_TEST(AnnotatorPreset3Test, get_top_labels) {
    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              this->annotation->get_top_labels({}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }));

    EXPECT_EQ(to_set(VectorCounts({ std::make_pair("Label1", 1),
                                    std::make_pair("Label2", 1) })),
              to_set(this->annotation->get_top_labels({ 2 })));

    EXPECT_EQ(VectorCounts({}),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              this->annotation->get_top_labels({ 0, 1, 2, 3, 4 }, 1000));
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
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);

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
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);

    this->annotation->serialize(test_dump_basename_vec_good + this->annotation->kExtension);
    this->annotation.reset(new TypeParam());
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad + this->annotation->kExtension));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good + this->annotation->kExtension));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoadExtension2) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);

    this->annotation->serialize(test_dump_basename_vec_good + this->annotation->kExtension);
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
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);

    this->annotation->serialize(test_dump_basename_vec_good);
    this->annotation.reset(new TypeParam());
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad + this->annotation->kExtension));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good + this->annotation->kExtension));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get_labels(4)));
}

TYPED_TEST(AnnotatorPreset2Test, SerializationAndLoad) {
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);

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
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);

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
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);

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
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->kExtension);

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(this->annotation->get(4)));
}

TYPED_TEST(AnnotatorDynamicNoSparseTest, MergeLoadDisjoint) {
    std::filesystem::remove(test_dump_basename_vec_good + "_1");
    std::filesystem::remove(test_dump_basename_vec_good + "_1" + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good + "_2");
    std::filesystem::remove(test_dump_basename_vec_good + "_2" + this->annotation->kExtension);

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
    std::filesystem::remove(test_dump_basename_vec_good + "_1" + this->annotation->kExtension);
    std::filesystem::remove(test_dump_basename_vec_good + "_2");
    std::filesystem::remove(test_dump_basename_vec_good + "_2" + this->annotation->kExtension);

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
