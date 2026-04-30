#include <gtest/gtest.h>

#include "../test_helpers.hpp"
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
    EXPECT_EQ(cth.num_sequences(0), 2u);

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
    EXPECT_EQ(cth.num_sequences(0), 2u);
    EXPECT_EQ(cth.num_sequences(1), 1u);

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

    EXPECT_EQ(cth.num_sequences(0), 3u);

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

    // Out-of-range coord.
    EXPECT_THROW(cth.map_single_coord(0, 5), std::out_of_range);
    EXPECT_THROW(cth.map_single_coord(0, 100), std::out_of_range);
    // Out-of-range column.
    EXPECT_THROW(cth.map_single_coord(1, 0), std::out_of_range);
    EXPECT_THROW(cth.map_single_coord(99, 0), std::out_of_range);
}

TEST(CoordToHeader, NumKmersInSequence) {
    // Three columns with unequal sequence sizes, including singletons
    // at the start, middle, and end of a column.
    CoordToHeader cth(
        { { "only" }, { "a", "b", "c" }, { "head", "mid", "tail" } },
        { { 7 },      { 1, 4, 2 },       { 3, 1, 5 } }
    );

    // Column 0: single sequence, 7 k-mers.
    EXPECT_EQ(7u, cth.num_kmers_in_sequence(0, 0));

    // Column 1: singleton at start, then 4, then 2.
    EXPECT_EQ(1u, cth.num_kmers_in_sequence(1, 0));
    EXPECT_EQ(4u, cth.num_kmers_in_sequence(1, 1));
    EXPECT_EQ(2u, cth.num_kmers_in_sequence(1, 2));

    // Column 2: 3, singleton in the middle, then 5.
    EXPECT_EQ(3u, cth.num_kmers_in_sequence(2, 0));
    EXPECT_EQ(1u, cth.num_kmers_in_sequence(2, 1));
    EXPECT_EQ(5u, cth.num_kmers_in_sequence(2, 2));

    // Out-of-range inputs.
    EXPECT_THROW(cth.num_kmers_in_sequence(0, 1), std::out_of_range);
    EXPECT_THROW(cth.num_kmers_in_sequence(3, 0), std::out_of_range);
}

TEST(CoordToHeader, SerializeLoadRoundTrip) {
    // Build a non-trivial index, serialize, load into a fresh instance,
    // and verify every observable method returns the same value.
    CoordToHeader src(
        { { "only" },
          { "a", "b", "c" },
          { "head", "mid", "tail" } },
        { { 7 },
          { 1, 4, 2 },
          { 3, 1, 5 } }
    );

    std::string path_base = test_dump_dir() + "/cth_roundtrip";
    src.serialize(path_base);

    CoordToHeader dst;
    ASSERT_TRUE(dst.load(path_base));

    ASSERT_EQ(src.num_columns(), dst.num_columns());
    for (uint64_t col = 0; col < src.num_columns(); ++col) {
        EXPECT_EQ(src.num_sequences(col), dst.num_sequences(col))  << "col " << col;
        EXPECT_EQ(src.num_kmers(col),   dst.num_kmers(col))    << "col " << col;
        EXPECT_EQ(src.get_headers(col), dst.get_headers(col))  << "col " << col;
        for (size_t h = 0; h < src.num_sequences(col); ++h) {
            EXPECT_EQ(src.num_kmers_in_sequence(col, h),
                      dst.num_kmers_in_sequence(col, h))
                << "col " << col << " seq_id " << h;
        }
        for (uint64_t coord = 0; coord < src.num_kmers(col); ++coord) {
            EXPECT_EQ(src.map_single_coord(col, coord),
                      dst.map_single_coord(col, coord))
                << "col " << col << " coord " << coord;
        }
    }
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
            auto [seq_id, local_coord] = cth.map_single_coord(col, global_coord);
            uint64_t encoded = local_coord * cth.num_sequences(col) + seq_id;
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

    EXPECT_EQ(aln.format_coords(cth, 5), "accession_A:4-7;accession_B:3-6");
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

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA:3-5;seqB:3-5");
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

    // Same header gets both ranges appended
    EXPECT_EQ(aln.format_coords(cth, 5), "ref:4-5:11-12");
}

TEST(AlignmentFormatCoords, EmptyCoordinates) {
    using namespace mtg::graph::align;

    LabelEncoder<> encoder;
    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),
        0, {}, 0, false, 0
    );
    // No coordinates -> empty string
    EXPECT_EQ(aln.format_coords(encoder), "");
}

TEST(AlignmentFormatCoords, WithEncoder) {
    using namespace mtg::graph::align;

    LabelEncoder<> encoder;
    encoder.insert_and_encode("/path/to/file.fa");

    Alignment aln(
        std::string_view{},
        {},
        std::string("ACGT"),  // 4-char sequence -> range width = 4
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 42 } };  // coord 42 -> "43-46"

    EXPECT_EQ(aln.format_coords(encoder), "/path/to/file.fa:43-46");
}

// A range that extends past the starting sequence is split across
// consecutive headers (joined by ';').
TEST(AlignmentFormatCoords, CrossSequenceBoundaryWithK) {
    using namespace mtg::graph::align;

    // seqA/seqB: 10 kmers each -> 14 nt with k=5.
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 10, 10 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGT"),  // 8 nt
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 7 } };  // seqA local 7; 7 nt fit + 1 spills

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA:8-14;seqB:1-1");
}

TEST(AlignmentFormatCoords, CrossThreeSequences) {
    using namespace mtg::graph::align;

    // Three headers x 5 kmers = 7 nt each with k=3.
    CoordToHeader cth(
        { { "seqA", "seqB", "seqC" } },
        { { 5, 5, 5 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGTACGTACGT"),  // 16 nt, crosses both boundaries
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 3 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA:4-7;seqB:1-7;seqC:1-5");
}

TEST(AlignmentFormatCoords, TwoCoordsBothCrossIntoSameNextSequence) {
    // Two starting coords that both spill from seqA into seqB. Each coord
    // produces a range in seqA and a range in seqB; the per-header ranges
    // are aggregated into insertion-ordered lists.
    using namespace mtg::graph::align;

    // seqA/seqB: 8 k-mers each -> 12 nt each with k=5.
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 8, 8 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGT"),  // 8 nt
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    // Both coords are in seqA; both produce 8 nt ranges that cross into seqB.
    //   coord 6 -> seqA local 6, avail=6, spills 2 nt into seqB
    //   coord 7 -> seqA local 7, avail=5, spills 3 nt into seqB
    aln.label_coordinates = { { 6, 7 } };

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA:7-12:8-12;seqB:1-2:1-3");
}

TEST(AlignmentFormatCoords, OverflowPastLastHeader) {
    using namespace mtg::graph::align;

    // Total column capacity = 14 nt; alignment of 20 nt overruns by 6 nt,
    // exercising the warn+truncate branch.
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 5, 5 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGTACGTACGTACGT"),  // 20 nt
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 0 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA:1-7;seqB:1-7");
}

TEST(AlignmentFormatCoords, ExactlyFillsStartingHeader) {
    using namespace mtg::graph::align;

    // L == available in the starting header: emit a single range and stop.
    CoordToHeader cth(
        { { "seqA", "seqB" } },
        { { 5, 5 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTA"),  // 5 nt; seqA has 7 nt, local 2..6 fills to end
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 2 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA:3-7");
}

TEST(AlignmentFormatCoords, SpansEntireColumn) {
    using namespace mtg::graph::align;

    // L == total nt of the column; the walk must land on zero remaining
    // after emitting every header exactly once.
    CoordToHeader cth(
        { { "seqA", "seqB", "seqC" } },
        { { 3, 4, 2 } }  // 5 + 6 + 4 = 15 nt with k=3
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGTACGTACG"),  // 15 nt
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 0 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA:1-5;seqB:1-6;seqC:1-4");
}

} // namespace
