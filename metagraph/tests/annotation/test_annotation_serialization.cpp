#include <filesystem>
#include "gtest/gtest.h"

#include "test_annotation.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_matrix";


namespace {

using namespace mtg;
using namespace mtg::test;

TYPED_TEST(AnnotatorPresetTest, SerializationAndLoad) {
    std::filesystem::remove(test_dump_basename_vec_good);
    std::filesystem::remove(test_dump_basename_vec_good + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_bad);
    std::filesystem::remove(test_dump_basename_vec_bad + this->annotation->file_extension());

    this->annotation->serialize(test_dump_basename_vec_good);
    ASSERT_FALSE(this->annotation->load(test_dump_basename_vec_bad));
    ASSERT_TRUE(this->annotation->load(test_dump_basename_vec_good));

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get(4)));
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

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get(4)));
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

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get(4)));
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

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(this->annotation->get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(this->annotation->get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(this->annotation->get(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(this->annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }),                     convert_to_set(this->annotation->get(4)));
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

TYPED_TEST(AnnotatorDynamicNoSparseTest, MergeLoadDisjoint) {
    std::filesystem::remove(test_dump_basename_vec_good + "_1");
    std::filesystem::remove(test_dump_basename_vec_good + "_1" + this->annotation->file_extension());
    std::filesystem::remove(test_dump_basename_vec_good + "_2");
    std::filesystem::remove(test_dump_basename_vec_good + "_2" + this->annotation->file_extension());

    {
        TypeParam annotation(5);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_1");
    }
    {
        TypeParam annotation(5);
        annotation.set(1, { "2_Label0", "2_Label2", "2_Label8" });
        annotation.set(2, { "2_Label1", "2_Label9", "2_Label0" });
        annotation.set(3, { "2_Label8" });

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
    this->annotation->set(1, { "Label0", "Label2", "Label8" });
    this->annotation->set(2, { "Label1", "Label9", "Label0" });
    this->annotation->set(3, { "Label8" });
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

} // namespace
