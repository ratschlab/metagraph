#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/file_utils.hpp"


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


TEST(ColumnCompressed, MergeLargeSparseColumnsWithoutDenseWorkspace) {
    const uint64_t rows = uint64_t(1) << 34;
    const uint64_t count = 4'000'001;
    std::vector<std::string> files;
    for (uint64_t offset : { 0u, 1u }) {
        files.push_back(test_dump_basename + "_large_sparse_" + std::to_string(offset));
        annot::ColumnCompressed<> part(rows);
        std::vector<uint64_t> positions(count);
        for (uint64_t i = 0; i < count; ++i)
            positions[i] = 2 * i + offset;
        part.add_labels(positions, { "shared" });
        part.serialize(files.back());
    }
    // Repeated inputs must not duplicate positions in the compressed output.
    files.push_back(files.front());
    annot::ColumnCompressed<> merged;
    ASSERT_TRUE(merged.merge_load(files));
    const auto &column = merged.get_column("shared");
    EXPECT_EQ(2 * count, column.num_set_bits());
    EXPECT_TRUE(column[0]);
    EXPECT_TRUE(column[2 * count - 1]);
    EXPECT_FALSE(column[2 * count]);
    EXPECT_FALSE(column[rows - 1]);

    const auto output = test_dump_basename + "_large_sparse_roundtrip";
    merged.serialize(output);
    annot::ColumnCompressed<> restored;
    ASSERT_TRUE(restored.load(output));
    for (const auto *result : { &column, &restored.get_column("shared") }) {
        uint64_t next = 0;
        bool ordered = true;
        result->call_ones([&](uint64_t position) { ordered &= position == next++; });
        EXPECT_TRUE(ordered);
        EXPECT_EQ(2 * count, next);
    }
}

TEST(ColumnCompressed, MergeManyFilesByLabel) {
    const uint64_t rows = uint64_t(1) << 34;
    std::vector<std::string> files;
    for (size_t f = 0; f < 48; ++f) {
        files.push_back(test_dump_basename + "_group_" + std::to_string(f));
        annot::ColumnCompressed<> part(rows, 1);
        // Vary the local column order; labels determine the destination.
        if (f % 2) part.add_labels({ rows - 1 }, { "tail" });
        part.add_labels({ 0, 100 + f, rows / 2 }, { "shared" });
        part.add_labels({ 1000 + f }, { "other" });
        for (size_t label = 0; label < 5; ++label)
            part.add_labels({ 2000 + 100 * label + f },
                            { "shared_" + std::to_string(label) });
        part.add_labels({ f }, { "unique_" + std::to_string(f) });
        if (!(f % 2)) part.add_labels({ rows - 1 }, { "tail" });
        part.serialize(files.back());
    }
    // Include repeats and reverse input order at different worker counts.
    files.push_back(files.front());
    const auto previous_threads = get_num_threads();
    for (size_t threads : { 1u, 4u, 8u }) {
        set_num_threads(threads);
        annot::ColumnCompressed<> merged;
        ASSERT_TRUE(merged.merge_load(files));
        EXPECT_EQ(56u, merged.num_labels());
        EXPECT_EQ(50u, merged.get_column("shared").num_set_bits());
        EXPECT_EQ(48u, merged.get_column("other").num_set_bits());
        EXPECT_EQ(1u, merged.get_column("tail").num_set_bits());
        for (size_t f = 0; f < 48; ++f) {
            EXPECT_TRUE(merged.get_column("shared")[100 + f]);
            EXPECT_TRUE(merged.get_column("other")[1000 + f]);
            const auto &unique = merged.get_column("unique_" + std::to_string(f));
            EXPECT_EQ(1u, unique.num_set_bits());
            EXPECT_TRUE(unique[f]);
            for (size_t label = 0; label < 5; ++label) {
                const auto &shared = merged.get_column("shared_" + std::to_string(label));
                EXPECT_EQ(48u, shared.num_set_bits());
                EXPECT_TRUE(shared[2000 + 100 * label + f]);
            }
        }
        const auto output = test_dump_basename + "_group_roundtrip";
        merged.serialize(output);
        annot::ColumnCompressed<> restored;
        ASSERT_TRUE(restored.merge_load({ output }));
        EXPECT_EQ(merged.num_relations(), restored.num_relations());
        for (const auto &label : merged.get_label_encoder().get_labels())
            EXPECT_EQ(merged.get_column(label), restored.get_column(label));
        std::reverse(files.begin(), files.end());
    }
    set_num_threads(previous_threads);
}

TEST(ColumnCompressed, MergeDenseAndEmptyColumns) {
    const uint64_t rows = 1024;
    std::vector<std::string> files;
    for (size_t f = 0; f < 16; ++f) {
        files.push_back(test_dump_basename + "_dense_group_" + std::to_string(f));
        annot::ColumnCompressed<> part(rows);
        std::vector<uint64_t> indexes;
        for (size_t i = f; i < rows; i += 16) indexes.push_back(i);
        part.add_labels(indexes, { "dense", "dense_second", "dense_third" });
        part.add_labels({ f }, { "becomes_dense" });
        part.add_labels({ 7 }, { "stays_sparse" });
        part.add_labels({}, { "empty" });
        part.serialize(files.back());
    }
    const auto previous_threads = get_num_threads();
    for (size_t threads : { 1u, 4u, 8u }) {
        set_num_threads(threads);
        annot::ColumnCompressed<> merged(0, 1);
        ASSERT_TRUE(merged.merge_load(files));
        for (const auto &label : { "dense", "dense_second", "dense_third" }) {
            EXPECT_EQ(rows, merged.get_column(label).num_set_bits());
            for (size_t i = 0; i < rows; ++i) EXPECT_TRUE(merged.get_column(label)[i]);
        }
        EXPECT_EQ(16u, merged.get_column("becomes_dense").num_set_bits());
        for (size_t i = 0; i < 16; ++i)
            EXPECT_TRUE(merged.get_column("becomes_dense")[i]);
        EXPECT_EQ(1u, merged.get_column("stays_sparse").num_set_bits());
        EXPECT_TRUE(merged.get_column("stays_sparse")[7]);
        EXPECT_EQ(0u, merged.get_column("empty").num_set_bits());
    }
    set_num_threads(previous_threads);
}

TEST(ColumnCompressed, MergeRejectsBadInputs) {
    const auto first = test_dump_basename + "_group_good";
    const auto other = test_dump_basename + "_group_wrong_rows";
    annot::ColumnCompressed<> a(100), b(101);
    a.add_labels({ 7 }, { "shared" });
    b.add_labels({ 8 }, { "shared" });
    a.serialize(first);
    b.serialize(other);
    annot::ColumnCompressed<> result;
    EXPECT_FALSE(result.merge_load({ first, other }));
    EXPECT_FALSE(result.merge_load({ first, first + "_missing" }));
    // A failed call must not prevent a subsequent valid merge.
    ASSERT_TRUE(result.merge_load({ first, first }));
    EXPECT_EQ(1u, result.get_column("shared").num_set_bits());
    EXPECT_TRUE(result.get_column("shared")[7]);
    // Preserve the existing treatment of files with no label columns.
    annot::ColumnCompressed<> empty(999);
    const auto empty_file = test_dump_basename + "_group_no_columns";
    empty.serialize(empty_file);
    ASSERT_TRUE(result.merge_load({ empty_file, first }));
    EXPECT_EQ(100u, result.num_objects());
}

TEST(ColumnCompressed, MergeEmptyAndTruncatedFiles) {
    const auto path = test_dump_basename + "_empty_group";
    annot::ColumnCompressed<> empty(0);
    empty.add_labels({}, { "empty" });
    empty.serialize(path);
    annot::ColumnCompressed<> restored;
    ASSERT_TRUE(restored.merge_load({ path, path }));
    EXPECT_EQ(0u, restored.num_objects());
    EXPECT_EQ(1u, restored.num_labels());
    EXPECT_EQ(0u, restored.num_relations());
    const auto truncated = path + "_truncated" + restored.file_extension();
    std::filesystem::copy_file(path + restored.file_extension(), truncated,
                              std::filesystem::copy_options::overwrite_existing);
    std::filesystem::resize_file(truncated, std::filesystem::file_size(truncated) - 1);
    EXPECT_FALSE(restored.merge_load({ path, truncated }));
}

TEST(ColumnCompressed, MergeMappedDenseColumns) {
    const uint64_t rows = 1024;
    const auto first_file = test_dump_basename + "_mapped_first";
    const auto second_file = test_dump_basename + "_mapped_second";
    annot::ColumnCompressed<> first(rows), second(rows);
    std::vector<uint64_t> even, odd;
    for (uint64_t i = 0; i < rows; ++i)
        (i % 2 ? odd : even).push_back(i);
    first.add_labels(even, { "shared" });
    second.add_labels(odd, { "shared" });
    first.serialize(first_file);
    second.serialize(second_file);

    const auto read_bytes = [](const std::string &path) {
        std::ifstream in(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
    };
    const auto first_bytes = read_bytes(first_file + first.file_extension());
    const auto second_bytes = read_bytes(second_file + second.file_extension());
    struct MmapSetting {
        bool previous = utils::with_mmap();
        MmapSetting() { utils::set_mmap(true); }
        ~MmapSetting() { utils::set_mmap(previous); }
    } mmap_setting;
    annot::ColumnCompressed<> merged;
    ASSERT_TRUE(merged.merge_load({ first_file, second_file, first_file }));
    EXPECT_EQ(rows, merged.get_column("shared").num_set_bits());
    for (uint64_t i = 0; i < rows; ++i)
        EXPECT_TRUE(merged.get_column("shared")[i]);
    const auto output = test_dump_basename + "_mapped_merged";
    merged.serialize(output);
    annot::ColumnCompressed<> restored;
    ASSERT_TRUE(restored.load(output));
    EXPECT_EQ(merged.get_column("shared"), restored.get_column("shared"));
    // Updating the merged bitmap must leave each mapped input unchanged.
    ASSERT_TRUE(restored.load(first_file));
    EXPECT_EQ(rows / 2, restored.get_column("shared").num_set_bits());
    for (uint64_t i = 0; i < rows; ++i)
        EXPECT_EQ(i % 2 == 0, restored.get_column("shared")[i]);
    EXPECT_EQ(first_bytes, read_bytes(first_file + first.file_extension()));
    EXPECT_EQ(second_bytes, read_bytes(second_file + second.file_extension()));
}

} // namespace
