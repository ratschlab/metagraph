#include <gtest/gtest.h>

#include "annotation/coord_to_header.hpp"
#include "graph/alignment/alignment.hpp"
#include "annotation/representation/base/annotation.hpp"


namespace {

using namespace mtg;
using namespace mtg::annot;

TEST(CoordToHeader, MapSingleCoordFirstSequence) {
    // Column 0: two sequences, with 5 and 3 k-mers respectively
    //   seq0: coords 0..4, seq1: coords 5..7
    CoordToHeader cth(
        { { "seq0", "seq1" } },
        { { 5, 3 } }
    );

    EXPECT_EQ(cth.num_columns(), 1u);
    EXPECT_EQ(cth.num_headers(0), 2u);

    // coord 0 -> seq0, local 0
    auto [h0, lc0] = cth.map_single_coord(0, 0);
    EXPECT_EQ(h0, 0u);
    EXPECT_EQ(lc0, 0u);

    // coord 4 -> seq0, local 4
    auto [h4, lc4] = cth.map_single_coord(0, 4);
    EXPECT_EQ(h4, 0u);
    EXPECT_EQ(lc4, 4u);
}

TEST(CoordToHeader, MapSingleCoordSecondSequence) {
    CoordToHeader cth(
        { { "alpha", "beta" } },
        { { 5, 3 } }
    );

    // coord 5 -> beta, local 0
    auto [h5, lc5] = cth.map_single_coord(0, 5);
    EXPECT_EQ(h5, 1u);
    EXPECT_EQ(lc5, 0u);

    // coord 7 -> beta, local 2
    auto [h7, lc7] = cth.map_single_coord(0, 7);
    EXPECT_EQ(h7, 1u);
    EXPECT_EQ(lc7, 2u);
}

TEST(CoordToHeader, MapSingleCoordMultipleColumns) {
    // Column 0: "A" (4 kmers), "B" (6 kmers)
    // Column 1: "X" (3 kmers)
    CoordToHeader cth(
        { { "A", "B" }, { "X" } },
        { { 4, 6 }, { 3 } }
    );

    EXPECT_EQ(cth.num_columns(), 2u);
    EXPECT_EQ(cth.num_headers(0), 2u);
    EXPECT_EQ(cth.num_headers(1), 1u);

    // Column 0, coord 3 -> A, local 3
    auto [h0, lc0] = cth.map_single_coord(0, 3);
    EXPECT_EQ(h0, 0u);
    EXPECT_EQ(lc0, 3u);

    // Column 0, coord 4 -> B, local 0
    auto [h1, lc1] = cth.map_single_coord(0, 4);
    EXPECT_EQ(h1, 1u);
    EXPECT_EQ(lc1, 0u);

    // Column 1, coord 2 -> X, local 2
    auto [h2, lc2] = cth.map_single_coord(1, 2);
    EXPECT_EQ(h2, 0u);
    EXPECT_EQ(lc2, 2u);
}

TEST(CoordToHeader, MapSingleCoordOutOfRange) {
    CoordToHeader cth(
        { { "only" } },
        { { 5 } }
    );

    EXPECT_THROW(cth.map_single_coord(0, 5), std::runtime_error);
    EXPECT_THROW(cth.map_single_coord(0, 100), std::runtime_error);
}

TEST(CoordToHeader, MapSingleCoordConsistentWithBatch) {
    // Verify map_single_coord gives same results as map_to_local_coords
    CoordToHeader cth(
        { { "h0", "h1", "h2" } },
        { { 3, 4, 2 } }
    );

    // Batch mapping
    using RowTuples = CoordToHeader::RowTuples;
    std::vector<RowTuples> rows = {
        { { 0, { 0, 2, 3, 6, 8 } } }
    };
    cth.map_to_local_coords(&rows);

    // Single-coord mapping should produce same encoded values
    for (auto &[col, coords] : rows[0]) {
        size_t idx = 0;
        for (uint64_t global_coord : { 0u, 2u, 3u, 6u, 8u }) {
            auto [header_id, local_coord] = cth.map_single_coord(col, global_coord);
            uint64_t encoded = local_coord * cth.num_headers(col) + header_id;
            EXPECT_EQ(encoded, coords[idx])
                << "Mismatch at global coord " << global_coord;
            ++idx;
        }
    }
}

// Test Alignment::format_coords() with coord_to_header
TEST(AlignmentFormatCoords, WithCoordToHeader) {
    using namespace mtg::graph::align;

    CoordToHeader cth(
        { { "accession_A", "accession_B" } },
        { { 10, 20 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),  // 4-char sequence for coordinate range calculation
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0, 0 };
    aln.label_coordinates = {
        { 3 },   // column 0, coord 3 -> accession_A (local 3)
        { 12 },  // column 0, coord 12 -> accession_B (local 2)
    };

    // Without coord_to_header, needs label_encoder (not set), falls back to column index
    std::string without = aln.format_coords();
    EXPECT_NE(without.find("0:"), std::string::npos);

    // With coord_to_header
    aln.coord_to_header = &cth;
    std::string with_cth = aln.format_coords();

    EXPECT_NE(with_cth.find("accession_A"), std::string::npos);
    EXPECT_NE(with_cth.find("accession_B"), std::string::npos);
    EXPECT_EQ(with_cth.find("0:"), std::string::npos);  // no raw column index
}

TEST(AlignmentFormatCoords, CoordToHeaderGroupsByHeader) {
    using namespace mtg::graph::align;

    // Two sequences in column 0: "seqA" (5 kmers), "seqB" (5 kmers)
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 5, 5 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACG"),  // 3-char sequence
        0, {}, 0, false, 0
    );
    // One column entry with coords from different sequences
    aln.label_columns = { 0 };
    aln.label_coordinates = {
        { 2, 7 },  // coord 2 -> seqA (local 2), coord 7 -> seqB (local 2)
    };
    aln.coord_to_header = &cth;

    std::string result = aln.format_coords();

    // Both headers should appear, separated by ;
    EXPECT_NE(result.find("seqA"), std::string::npos);
    EXPECT_NE(result.find("seqB"), std::string::npos);
    EXPECT_NE(result.find(";"), std::string::npos);
}

TEST(AlignmentFormatCoords, WithoutCoordToHeaderUnchanged) {
    using namespace mtg::graph::align;

    LabelEncoder<> encoder;
    encoder.insert_and_encode("/path/to/file.fa");

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),
        0, {}, 0, false, 0
    );
    aln.label_encoder = &encoder;
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 42 } };

    std::string result = aln.format_coords();
    EXPECT_NE(result.find("/path/to/file.fa"), std::string::npos);
    EXPECT_NE(result.find("43-46"), std::string::npos);
}

} // namespace
