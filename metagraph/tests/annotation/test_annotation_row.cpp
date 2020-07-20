#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#define private public
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_row_compressed";


TEST(RowCompressed, load_label_encoder) {
    {
        annot::RowCompressed<> annotation(5, false);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        auto label_encoder = annot::RowCompressed<>::load_label_encoder(test_dump_basename_vec_good);
        ASSERT_TRUE(label_encoder.get());
        EXPECT_EQ(4u, label_encoder->size());
    }
}

TEST(RowCompressed, load_shape) {
    {
        annot::RowCompressed<> annotation(5, false);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        uint64_t num_rows, num_relations;
        annot::RowCompressed<>::load_shape(test_dump_basename_vec_good,
                                              &num_rows, &num_relations);
        ASSERT_EQ(5u, num_rows);
        ASSERT_EQ(6u, num_relations);
    }
}

TEST(RowCompressed, load_label_encoder_and_load_shape) {
    {
        annot::RowCompressed<> annotation(5, false);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        uint64_t num_rows, num_relations;
        auto label_encoder = annot::RowCompressed<>::load_label_encoder(
            test_dump_basename_vec_good
        );
        annot::RowCompressed<>::load_shape(test_dump_basename_vec_good,
                                              &num_rows, &num_relations);
        ASSERT_EQ(5u, num_rows);
        ASSERT_EQ(6u, num_relations);
        ASSERT_TRUE(label_encoder.get());
        EXPECT_EQ(4u, label_encoder->size());
    }
}

TEST(RowCompressed, load_shape_and_load_label_encoder) {
    {
        annot::RowCompressed<> annotation(5, false);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        uint64_t num_rows, num_relations;
        annot::RowCompressed<>::load_shape(test_dump_basename_vec_good,
                                              &num_rows, &num_relations);
        auto label_encoder = annot::RowCompressed<>::load_label_encoder(
            test_dump_basename_vec_good
        );
        ASSERT_EQ(5u, num_rows);
        ASSERT_EQ(6u, num_relations);
        ASSERT_TRUE(label_encoder.get());
        EXPECT_EQ(4u, label_encoder->size());
    }
}

TEST(RowCompressed, SerializationExtension) {
    {
        annot::RowCompressed<> annotation(5, false);
        annotation.set(0, { "Label0", "Label2", "Label8" });
        annotation.set(2, { "Label1", "Label2" });
        annotation.set(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good
                                        + annotation.file_extension());
    }
    {
        annot::RowCompressed<> annotation(5, false);
        ASSERT_FALSE(annotation.load(test_dump_basename_vec_bad
                                        + annotation.file_extension()));
        ASSERT_TRUE(annotation.load(test_dump_basename_vec_good
                                        + annotation.file_extension()));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }),
                  convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),
                  convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }),
                  convert_to_set(annotation.get(4)));
    }
}

TEST(RowCompressed, RenameColumnsMerge) {
    annot::RowCompressed<> annotation(5);
    annotation.set(0, { "Label0", "Label2", "Label8" });
    annotation.set(2, { "Label1", "Label2" });
    annotation.set(4, { "Label8" });

    ASSERT_DEATH_SILENT(
        annotation.rename_labels({ { "Label2", "Merged" },
                                   { "Label8", "Merged" } }),
        ""
    );

    // EXPECT_EQ(convert_to_set({ "Label0", "Merged" }), convert_to_set(annotation.get(0)));
    // EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    // EXPECT_EQ(convert_to_set({ "Label1", "Merged" }), convert_to_set(annotation.get(2)));
    // EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    // EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(RowCompressed, RenameColumnsMergeAll) {
    annot::RowCompressed<> annotation(5);
    annotation.set(0, { "Label0", "Label2", "Label8" });
    annotation.set(2, { "Label1", "Label2" });
    annotation.set(4, { "Label8" });

    ASSERT_DEATH_SILENT(
        annotation.rename_labels({ { "Label0", "Merged" },
                                   { "Label1", "Merged" },
                                   { "Label2", "Merged" },
                                   { "Label8", "Merged" }, }),
        ""
    );

    // EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(0)));
    // EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    // EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(2)));
    // EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    // EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

} // namespace
