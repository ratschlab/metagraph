#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "common/threads/threading.hpp"


namespace {

using namespace mtg;

const std::string test_dump_basename = test_dump_dir() + "/dump_test";
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
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TEST(ColumnCompressed, RenameColumnsMerge) {
    annot::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    annotation.rename_labels({ { "Label2", "Merged" },
                               { "Label8", "Merged" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Merged" }), convert_to_set(annotation.get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Merged" }), convert_to_set(annotation.get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get_labels(4)));
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

    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get_labels(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get_labels(4)));
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

TEST(ColumnCompressed, SparseUpdatesAfterCacheEviction) {
    const uint64_t rows = uint64_t(1) << 34;
    annot::ColumnCompressed<> annotation(rows, 1);
    annotation.add_labels({ 0, rows - 1 }, { "first" });
    annotation.add_labels({ rows / 2 }, { "second" });
    annotation.add_labels({ 0, 7, rows - 1 }, { "first" });
    annotation.add_labels({ 8, rows / 2 }, { "second" });
    const auto filename = test_dump_basename + "_sparse_revisit";
    annotation.serialize(filename);
    annot::ColumnCompressed<> restored;
    ASSERT_TRUE(restored.merge_load({ filename }));
    EXPECT_EQ(3u, restored.get_column("first").num_set_bits());
    EXPECT_EQ(2u, restored.get_column("second").num_set_bits());
    EXPECT_TRUE(restored.get_column("first")[rows - 1]);
    EXPECT_FALSE(restored.get_column("first")[8]);
    EXPECT_TRUE(restored.get_column("second")[8]);
    EXPECT_EQ(convert_to_set({ "first" }), convert_to_set(restored.get_labels(7)));
    EXPECT_EQ(convert_to_set({ "second" }), convert_to_set(restored.get_labels(8)));
    EXPECT_TRUE(restored.get_labels(9).empty());
}

TEST(ColumnCompressed, MergeSparseColumnsWithLargeRowCount) {
    const uint64_t rows = uint64_t(1) << 34;
    const auto first_file = test_dump_basename + "_sparse_merge_first";
    const auto second_file = test_dump_basename + "_sparse_merge_second";
    annot::ColumnCompressed<> first(rows, 1), second(rows, 1);
    first.add_labels({ 0, rows - 1 }, { "shared" });
    first.add_labels({ 5 }, { "first" });
    second.add_labels({ 7, rows - 1 }, { "shared" });
    second.add_labels({ 8 }, { "second" });
    first.serialize(first_file);
    second.serialize(second_file);
    for (const auto &files : std::vector<std::vector<std::string>>{
            {first_file, second_file}, {second_file, first_file},
            {first_file, second_file, first_file}}) {
        annot::ColumnCompressed<> merged(rows, 1);
        ASSERT_TRUE(merged.merge_load(files));
        EXPECT_EQ(3u, merged.get_column("shared").num_set_bits());
        EXPECT_EQ(convert_to_set({ "shared" }), convert_to_set(merged.get_labels(0)));
        EXPECT_EQ(convert_to_set({ "shared" }), convert_to_set(merged.get_labels(7)));
        EXPECT_EQ(convert_to_set({ "shared" }), convert_to_set(merged.get_labels(rows - 1)));
        EXPECT_EQ(convert_to_set({ "first" }), convert_to_set(merged.get_labels(5)));
        EXPECT_EQ(convert_to_set({ "second" }), convert_to_set(merged.get_labels(8)));
        EXPECT_TRUE(merged.get_labels(6).empty());
    }
}

} // namespace
