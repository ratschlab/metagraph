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

TEST(CoordToHeader, MapSingleCoordThreeSequences) {
    // Column 0: "a" (2 kmers), "b" (3 kmers), "c" (4 kmers) -> coords 0-1, 2-4, 5-8
    CoordToHeader cth(
        { { "a", "b", "c" } },
        { { 2, 3, 4 } }
    );

    EXPECT_EQ(cth.num_headers(0), 3u);

    // coord 0 -> a, local 0  (first seq, first coord)
    auto [h0, lc0] = cth.map_single_coord(0, 0);
    EXPECT_EQ(h0, 0u);
    EXPECT_EQ(lc0, 0u);

    // coord 1 -> a, local 1  (first seq, last coord)
    auto [h1, lc1] = cth.map_single_coord(0, 1);
    EXPECT_EQ(h1, 0u);
    EXPECT_EQ(lc1, 1u);

    // coord 2 -> b, local 0  (second seq boundary)
    auto [h2, lc2] = cth.map_single_coord(0, 2);
    EXPECT_EQ(h2, 1u);
    EXPECT_EQ(lc2, 0u);

    // coord 4 -> b, local 2  (second seq, last coord)
    auto [h4, lc4] = cth.map_single_coord(0, 4);
    EXPECT_EQ(h4, 1u);
    EXPECT_EQ(lc4, 2u);

    // coord 5 -> c, local 0  (third seq boundary)
    auto [h5, lc5] = cth.map_single_coord(0, 5);
    EXPECT_EQ(h5, 2u);
    EXPECT_EQ(lc5, 0u);

    // coord 8 -> c, local 3  (last valid coord in column)
    auto [h8, lc8] = cth.map_single_coord(0, 8);
    EXPECT_EQ(h8, 2u);
    EXPECT_EQ(lc8, 3u);
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

    // accession_A: 10 kmers (coords 0-9), accession_B: 20 kmers (coords 10-29)
    CoordToHeader cth(
        { { "accession_A", "accession_B" } },
        { { 10, 20 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),  // 4-char sequence -> range width = 4
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0, 0 };
    aln.label_coordinates = {
        { 3 },   // column 0, coord 3 -> accession_A, local 3 -> "4-7"
        { 12 },  // column 0, coord 12 -> accession_B, local 2 -> "3-6"
    };
    aln.coord_to_header = &cth;

    EXPECT_EQ(aln.format_coords(), "accession_A:4-7;accession_B:3-6");
}

TEST(AlignmentFormatCoords, CoordToHeaderGroupsByHeader) {
    using namespace mtg::graph::align;

    // seqA: 5 kmers (coords 0-4), seqB: 5 kmers (coords 5-9)
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 5, 5 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACG"),  // 3-char sequence -> range width = 3
        0, {}, 0, false, 0
    );
    // One column entry with coords from different sequences
    aln.label_columns = { 0 };
    aln.label_coordinates = {
        { 2, 7 },  // coord 2 -> seqA local 2, coord 7 -> seqB local 2
    };
    aln.coord_to_header = &cth;

    EXPECT_EQ(aln.format_coords(), "seqA:3-5;seqB:3-5");
}

TEST(AlignmentFormatCoords, CoordToHeaderSameHeaderMultipleCoords) {
    using namespace mtg::graph::align;

    // Single sequence with 20 kmers
    CoordToHeader cth(
        { { "ref" } },
        { { 20 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("AC"),  // 2-char sequence -> range width = 2
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = {
        { 3, 10 },  // two coords in the same header
    };
    aln.coord_to_header = &cth;

    // Same header gets both ranges appended
    EXPECT_EQ(aln.format_coords(), "ref:4-5:11-12");
}

TEST(AlignmentFormatCoords, EmptyCoordinates) {
    using namespace mtg::graph::align;

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),
        0, {}, 0, false, 0
    );
    // No coordinates -> empty string
    EXPECT_EQ(aln.format_coords(), "");
}

TEST(AlignmentFormatCoords, WithoutCoordToHeaderUsesLabelEncoder) {
    using namespace mtg::graph::align;

    LabelEncoder<> encoder;
    encoder.insert_and_encode("/path/to/file.fa");

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),  // 4-char sequence -> range width = 4
        0, {}, 0, false, 0
    );
    aln.label_encoder = &encoder;
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 42 } };  // coord 42 -> "43-46"

    EXPECT_EQ(aln.format_coords(), "/path/to/file.fa:43-46");
}

TEST(AlignmentFormatCoords, WithoutCoordToHeaderOrEncoder) {
    using namespace mtg::graph::align;

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),
        0, {}, 0, false, 0
    );
    // No label_encoder, no coord_to_header -> falls back to column index
    aln.label_columns = { 5 };
    aln.label_coordinates = { { 0 } };

    EXPECT_EQ(aln.format_coords(), "5:1-4");
}

// When coord_to_header_k is set and a range extends past the starting
// sequence's boundary, format_coords splits it across consecutive headers
// (joined by ';').  Without coord_to_header_k, the legacy unsplit behaviour
// is retained (nucleotide length is unknown).
TEST(AlignmentFormatCoords, CrossSequenceBoundaryWithK) {
    using namespace mtg::graph::align;

    // seqA: 10 kmers (coords 0-9), seqB: 10 kmers (coords 10-19).
    // With k=5, each sequence has 10+k-1 = 14 nt.
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 10, 10 } }
    );

    // 8-nt alignment starting at coord 7 in seqA (local 7).
    // seqA has 14 nt, so available from local 7 is 14-7 = 7 nt.
    // Remaining 1 nt spills into seqB starting at local 0.
    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGTACGT"),  // 8-char reference spelling
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 7 } };
    aln.coord_to_header = &cth;
    aln.coord_to_header_k = 5;

    EXPECT_EQ(aln.format_coords(), "seqA:8-14;seqB:1-1");
}

TEST(AlignmentFormatCoords, CrossSequenceBoundaryWithoutKUsesLegacy) {
    using namespace mtg::graph::align;

    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 10, 10 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGTACGT"),
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 7 } };
    aln.coord_to_header = &cth;
    // coord_to_header_k left as 0 -> legacy unsplit output

    EXPECT_EQ(aln.format_coords(), "seqA:8-15");
}

TEST(AlignmentFormatCoords, CrossThreeSequences) {
    using namespace mtg::graph::align;

    // Three sequences each with 5 kmers.  With k=3 each is 7 nt, so coords:
    //   seqA: global 0-4 (7 nt), seqB: global 5-9 (7 nt), seqC: global 10-14 (7 nt)
    CoordToHeader cth(
        { { "seqA", "seqB", "seqC" } },
        { { 5, 5, 5 } }
    );

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGTACGTACGTACGT"),  // 16 nt spelling
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 3 } };  // start in seqA at local 3
    aln.coord_to_header = &cth;
    aln.coord_to_header_k = 3;

    // seqA: local 3 to 6 (4 nt), seqB: local 0 to 6 (7 nt), seqC: local 0 to 4 (5 nt)
    EXPECT_EQ(aln.format_coords(), "seqA:4-7;seqB:1-7;seqC:1-5");
}

} // namespace
