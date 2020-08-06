#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "common/threads/threading.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column_compressed";


TEST(ColumnCompressed, add_label_random_with_caching) {
    size_t graph_half_size = 1000;
    annot::ColumnCompressed<> annotation(graph_half_size * 2, 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_labels({ i }, { labels[i % 2] });
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TEST(ColumnCompressed, RenameColumnsMerge) {
    annot::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    annotation.rename_labels({ { "Label2", "Merged" },
                               { "Label8", "Merged" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, RenameColumnsMergeAll) {
    annot::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    annotation.rename_labels({ { "Label0", "Merged" },
                               { "Label1", "Merged" },
                               { "Label2", "Merged" },
                               { "Label8", "Merged" }, });

    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, GetColumn) {
    annot::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    EXPECT_EQ(1u, annotation.get_column("Label0").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label0")[0]);

    EXPECT_EQ(1u, annotation.get_column("Label1").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label1")[2]);

    EXPECT_EQ(2u, annotation.get_column("Label2").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label2")[0]);
    EXPECT_TRUE(annotation.get_column("Label2")[2]);

    EXPECT_EQ(2u, annotation.get_column("Label8").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label8")[0]);
    EXPECT_TRUE(annotation.get_column("Label8")[4]);
}

} // namespace
