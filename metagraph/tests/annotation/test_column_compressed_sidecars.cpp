#include <filesystem>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/string_utils.hpp"

namespace fs = std::filesystem;

namespace {

using namespace mtg;

using AC = annot::ColumnCompressed<>;

std::string sidecar_base(const char *suffix) {
    return test_dump_dir() + "/cc_sidecars_" + suffix;
}

void remove_cc_files(const std::string &base_no_ext) {
    std::error_code ec;
    const std::string main = utils::make_suffix(base_no_ext, AC::kExtension);
    const std::string stem = utils::remove_suffix(main, AC::kExtension);
    fs::remove(main, ec);
    fs::remove(stem + AC::kCoordExtension, ec);
    fs::remove(stem + AC::kCountExtension, ec);
}

std::string coord_path_for_main(const std::string &base_no_ext) {
    const std::string main = utils::make_suffix(base_no_ext, AC::kExtension);
    return utils::remove_suffix(main, AC::kExtension) + AC::kCoordExtension;
}

std::string counts_path_for_main(const std::string &base_no_ext) {
    const std::string main = utils::make_suffix(base_no_ext, AC::kExtension);
    return utils::remove_suffix(main, AC::kExtension) + AC::kCountExtension;
}

TEST(ColumnCompressedSidecars, EmptyAnnotatorCoordinateModeWritesCoordsOnly) {
    const std::string base = sidecar_base("empty_coord_only");
    remove_cc_files(base);

    // count_width 0 => no .counts; index_coordinates => empty .coords for Snakemake-style pipelines.
    AC anno(16, 1, "", uint64_t(10'000'000), /*count_width=*/0, /*index_coordinates=*/true, 2000);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_TRUE(fs::exists(coord_path_for_main(base)));
    EXPECT_FALSE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, EmptyAnnotatorCountsModeWritesCountsOnly) {
    const std::string base = sidecar_base("empty_counts_only");
    remove_cc_files(base);

    AC anno(16, 1, "", uint64_t(10'000'000), /*count_width=*/8, /*index_coordinates=*/false, 2000);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_FALSE(fs::exists(coord_path_for_main(base)));
    EXPECT_TRUE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, EmptyAnnotatorBothModesWriteBothSidecars) {
    const std::string base = sidecar_base("empty_both");
    remove_cc_files(base);

    AC anno(16, 1, "", uint64_t(10'000'000), /*count_width=*/8, /*index_coordinates=*/true, 2000);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_TRUE(fs::exists(coord_path_for_main(base)));
    EXPECT_TRUE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, EmptyAnnotatorNeitherModeWritesNoSidecars) {
    const std::string base = sidecar_base("empty_neither");
    remove_cc_files(base);

    AC anno(16, 1, "", uint64_t(10'000'000), /*count_width=*/0, /*index_coordinates=*/false, 2000);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_FALSE(fs::exists(coord_path_for_main(base)));
    EXPECT_FALSE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, ExplicitCountWidthWritesCountsSidecar) {
    const std::string base = sidecar_base("explicit_counts_ctor");
    remove_cc_files(base);

    AC anno(12, 1, "", uint64_t(10'000'000), /*count_width=*/8, false, 2000);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_FALSE(fs::exists(coord_path_for_main(base)));
    EXPECT_TRUE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, DefaultCtorWritesNoCountsSidecar) {
    const std::string base = sidecar_base("default_ctor_no_counts");
    remove_cc_files(base);

    AC anno(12);
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(utils::make_suffix(base, AC::kExtension)));
    EXPECT_FALSE(fs::exists(coord_path_for_main(base)));
    EXPECT_FALSE(fs::exists(counts_path_for_main(base)));
}

TEST(ColumnCompressedSidecars, CoordinateColumnSerializeAndBitmapRoundTrip) {
    const std::string base = sidecar_base("coord_one_col");
    remove_cc_files(base);

    AC anno(8, 1, "", uint64_t(10'000'000), 0, true, 2000);
    anno.add_labels({2}, {"colA"});
    anno.add_label_coords({{2, 42}}, {"colA"});
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(coord_path_for_main(base)));

    AC loaded(8, 1, "", uint64_t(10'000'000), 0, false, 2000);
    ASSERT_TRUE(loaded.load(utils::make_suffix(base, AC::kExtension)));

    EXPECT_EQ(1u, loaded.num_labels());
    EXPECT_EQ(convert_to_set({"colA"}), convert_to_set(loaded.get_labels(2)));
}

TEST(ColumnCompressedSidecars, CountColumnSerializeAndRoundTrip) {
    const std::string base = sidecar_base("counts_one_col");
    remove_cc_files(base);

    AC anno(8, 1, "", uint64_t(10'000'000), 8, false, 2000);
    anno.add_labels({4}, {"gene"});
    anno.add_label_counts({4}, {"gene"}, {5});
    anno.serialize(base);

    EXPECT_TRUE(fs::exists(counts_path_for_main(base)));

    AC loaded(8);
    ASSERT_TRUE(loaded.load(utils::make_suffix(base, AC::kExtension)));
    EXPECT_EQ(convert_to_set({"gene"}), convert_to_set(loaded.get_labels(4)));
}

TEST(ColumnCompressedSidecars, LoadThenSerializeAgainCoordinateMode) {
    const std::string base = sidecar_base("coord_reload_a");
    const std::string base2 = sidecar_base("coord_reload_b");
    remove_cc_files(base);
    remove_cc_files(base2);

    AC anno(8, 1, "", uint64_t(10'000'000), 0, true, 2000);
    anno.add_labels({1}, {"x"});
    anno.add_label_coords({{1, 7}}, {"x"});
    anno.serialize(base);

    ASSERT_TRUE(anno.load(utils::make_suffix(base, AC::kExtension)));
    ASSERT_NO_THROW(anno.serialize(base2));

    EXPECT_TRUE(fs::exists(utils::make_suffix(base2, AC::kExtension)));
    EXPECT_TRUE(fs::exists(coord_path_for_main(base2)));
}

TEST(ColumnCompressedSidecars, LoadClearsThenSerializeCountsStillWritten) {
    const std::string base = sidecar_base("counts_reload");
    remove_cc_files(base);

    AC anno(8, 1, "", uint64_t(10'000'000), 8, false, 2000);
    anno.add_labels({5}, {"y"});
    anno.add_label_counts({5}, {"y"}, {3});
    anno.serialize(base);

    ASSERT_TRUE(anno.load(utils::make_suffix(base, AC::kExtension)));
    ASSERT_NO_THROW(anno.serialize(base)); // overwrite same stem; relation_counts cleared by load

    EXPECT_TRUE(fs::exists(counts_path_for_main(base)));
}

} // namespace
