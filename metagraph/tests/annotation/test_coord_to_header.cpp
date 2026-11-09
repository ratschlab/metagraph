#include <gtest/gtest.h>

#include <atomic>
#include <fstream>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <thread>

#include "../test_helpers.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/alignment/alignment.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "coordinate_test_graph.hpp"


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

// sequence_range (one rank/select per sequence for sorted coordinates) agrees with
// map_single_coord at every coordinate and throws where it throws
TEST(CoordToHeader, SequenceRangeAgreesWithMapSingleCoord) {
    CoordToHeader cth(
        { { "only" }, { "a", "b", "c" }, { "head", "mid", "tail" } },
        { { 7 },      { 1, 4, 2 },       { 3, 1, 5 } }
    );
    for (CoordToHeader::Column col = 0; col < cth.num_columns(); ++col) {
        for (uint64_t coord = 0; coord < cth.num_kmers(col); ++coord) {
            const auto [seq, local] = cth.map_single_coord(col, coord);
            const auto range = cth.sequence_range(col, coord);
            EXPECT_EQ(seq, range.seq_id) << col << " " << coord;
            EXPECT_EQ(coord - local, range.first) << col << " " << coord;
            EXPECT_EQ(range.first + cth.num_kmers_in_sequence(col, seq) - 1, range.last);
            EXPECT_LE(range.first, coord);
            EXPECT_GE(range.last, coord);
        }
        EXPECT_THROW(cth.sequence_range(col, cth.num_kmers(col)), std::out_of_range);
    }
    EXPECT_THROW(cth.sequence_range(3, 0), std::out_of_range);
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

// the bytes of a file
std::string file_bytes(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

CoordToHeader three_columns() {
    return CoordToHeader({ { "only" }, { "a", "b", "c" }, { "head", "mid", "tail" } },
                         { { 7 }, { 1, 4, 2 }, { 3, 1, 5 } });
}

// The .seqs states its column labels after its columns (H4): the bytes before them are a file
// without labels, which loads as before (none stated), and a load reads them back
TEST(CoordToHeader, ColumnLabelsTrailerRoundTrip) {
    const std::vector<std::string> labels { "dir/x.fa", "y.fa.gz", "" };
    CoordToHeader plain = three_columns();
    CoordToHeader labelled = three_columns();
    EXPECT_FALSE(labelled.column_labels());
    labelled.set_column_labels(std::vector<std::string>(labels));
    ASSERT_TRUE(labelled.column_labels());
    EXPECT_EQ(labels, *labelled.column_labels());
    const std::string plain_base = test_dump_dir() + "/cth_plain";
    const std::string labelled_base = test_dump_dir() + "/cth_labelled";
    plain.serialize(plain_base);
    labelled.serialize(labelled_base);
    const std::string old_bytes = file_bytes(plain_base + CoordToHeader::kExtension);
    const std::string new_bytes = file_bytes(labelled_base + CoordToHeader::kExtension);
    ASSERT_LT(old_bytes.size(), new_bytes.size());
    EXPECT_EQ(old_bytes, new_bytes.substr(0, old_bytes.size()));
    // the trailer: the marker, then the labels as a string vector
    std::ostringstream trailer;
    serialize_number(trailer, CoordToHeader::kColumnLabelsMarker);
    serialize_string_vector(trailer, labels);
    EXPECT_EQ(trailer.str(), new_bytes.substr(old_bytes.size()));
    EXPECT_EQ("mtg-cols", new_bytes.substr(old_bytes.size(), 8));

    CoordToHeader old_file;
    ASSERT_TRUE(old_file.load(plain_base));
    EXPECT_FALSE(old_file.column_labels());
    CoordToHeader new_file;
    ASSERT_TRUE(new_file.load(labelled_base));
    ASSERT_TRUE(new_file.column_labels());
    EXPECT_EQ(labels, *new_file.column_labels());
    for (uint64_t col = 0; col < plain.num_columns(); ++col) {
        EXPECT_EQ(plain.get_headers(col), new_file.get_headers(col));
        EXPECT_EQ(plain.num_kmers(col), new_file.num_kmers(col));
    }
    // a load of a file without labels forgets the labels loaded before
    ASSERT_TRUE(new_file.load(plain_base));
    EXPECT_FALSE(new_file.column_labels());
    // labels of another number than the columns are refused
    EXPECT_THROW(plain.set_column_labels({ "a", "b" }), std::invalid_argument);
}

// check_columns: the count always; the order when the labels are stated, which it drops once
// they matched (column_order_verified); a mismatch keeps everything as it was
TEST(CoordToHeader, CheckColumnsAgainstTheAnnotationsLabels) {
    const std::vector<std::string> labels { "x.fa", "y.fa", "z.fa" };
    CoordToHeader unlabelled = three_columns();
    EXPECT_EQ("it has 3 columns, the annotation 2 labels",
              unlabelled.check_columns({ "x.fa", "y.fa" }));
    EXPECT_EQ("", unlabelled.check_columns({ "a", "b", "c" }));
    EXPECT_FALSE(unlabelled.column_order_verified());

    CoordToHeader labelled = three_columns();
    labelled.set_column_labels(std::vector<std::string>(labels));
    EXPECT_EQ("it has 3 columns, the annotation 4 labels",
              labelled.check_columns({ "x.fa", "y.fa", "z.fa", "w.fa" }));
    EXPECT_EQ("its column 1 is 'y.fa', the annotation's column 1 'z.fa' (the files were passed "
              "in another order or under other names)",
              labelled.check_columns({ "x.fa", "z.fa", "y.fa" }));
    EXPECT_EQ("its column 0 is 'x.fa', the annotation's column 0 '/abs/x.fa' (the files were "
              "passed in another order or under other names)",
              labelled.check_columns({ "/abs/x.fa", "y.fa", "z.fa" }));
    EXPECT_FALSE(labelled.column_order_verified());
    ASSERT_TRUE(labelled.column_labels());
    EXPECT_EQ("", labelled.check_columns(labels));
    EXPECT_TRUE(labelled.column_order_verified());
    EXPECT_FALSE(labelled.column_labels());
    // a copy keeps what was verified
    EXPECT_TRUE(CoordToHeader(labelled).column_order_verified());
}

// Data after the columns is the trailer or a corrupt file: anything else fails the load
TEST(CoordToHeader, LoadRefusesUnknownDataAfterTheColumns) {
    const std::string base = test_dump_dir() + "/cth_tail";
    const std::string path = base + CoordToHeader::kExtension;
    three_columns().serialize(base);
    const std::string columns = file_bytes(path);
    auto write = [&](const std::string &tail) {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out << columns << tail;
    };
    auto trailer = [](uint64_t marker, const std::vector<std::string> &labels) {
        std::ostringstream out;
        serialize_number(out, marker);
        serialize_string_vector(out, labels);
        return out.str();
    };
    CoordToHeader cth;
    write(trailer(CoordToHeader::kColumnLabelsMarker, { "a", "b", "c" }));
    EXPECT_TRUE(cth.load(base));
    // another marker, labels for another number of columns, data after the labels, a cut file
    write(trailer(CoordToHeader::kColumnLabelsMarker + 1, { "a", "b", "c" }));
    EXPECT_FALSE(cth.load(base));
    write(trailer(CoordToHeader::kColumnLabelsMarker, { "a", "b" }));
    EXPECT_FALSE(cth.load(base));
    write(trailer(CoordToHeader::kColumnLabelsMarker, { "a", "b", "c" }) + "x");
    EXPECT_FALSE(cth.load(base));
    write("tail");
    EXPECT_FALSE(cth.load(base));
    const std::string whole = trailer(CoordToHeader::kColumnLabelsMarker, { "a", "b", "c" });
    write(whole.substr(0, whole.size() - 1));
    EXPECT_FALSE(cth.load(base));
    write("");
    EXPECT_TRUE(cth.load(base));
    EXPECT_FALSE(cth.column_labels());
}

// A file cut anywhere is refused, also one written before the trailer and cut inside its last
// column: its load must not end in a failed read that looks like the end of an old file (H4
// review S1: cut 1, 9 or 40 bytes, it loaded as an old file and answered without the last
// column's records, or aborted on a coordinate it could not place)
TEST(CoordToHeader, LoadRefusesATruncatedFile) {
    const std::string base = test_dump_dir() + "/cth_cut";
    const std::string path = base + CoordToHeader::kExtension;
    for (bool labelled : { false, true }) {
        CoordToHeader cth = three_columns();
        if (labelled)
            cth.set_column_labels({ "x.fa", "y.fa", "z.fa" });
        cth.serialize(base);
        const std::string whole = file_bytes(path);
        CoordToHeader loaded;
        ASSERT_TRUE(loaded.load(base));
        for (size_t cut : std::initializer_list<size_t> { 1, 8, 9, 40 }) {
            ASSERT_LT(cut, whole.size());
            {
                std::ofstream out(path, std::ios::binary | std::ios::trunc);
                out << whole.substr(0, whole.size() - cut);
            }
            EXPECT_FALSE(loaded.load(base)) << cut << (labelled ? " labelled" : " old");
        }
        // every shorter prefix of an old file too
        if (!labelled) {
            for (size_t size = 0; size < whole.size(); ++size) {
                {
                    std::ofstream out(path, std::ios::binary | std::ios::trunc);
                    out << whole.substr(0, size);
                }
                EXPECT_FALSE(loaded.load(base)) << size << " of " << whole.size() << " bytes";
            }
        }
    }
}

// An empty file fails its load, also when the stream maps it (--mmap: mapping an empty file
// throws, which ended the process before its load could name the file)
TEST(CoordToHeader, AnEmptyFileFailsItsLoad) {
    const std::string base = test_dump_dir() + "/cth_empty";
    { std::ofstream out(base + CoordToHeader::kExtension, std::ios::binary | std::ios::trunc); }
    CoordToHeader cth;
    EXPECT_FALSE(cth.load(base));
    const bool mmap = utils::with_mmap();
    utils::set_mmap(true);
    EXPECT_FALSE(cth.load(base));
    utils::set_mmap(mmap);
}

// The load's check of a record mapping against its annotation (cli::check_record_mapping): the
// count, the stated labels, and the coordinates of a column against its k-mers where the
// annotation states them per column (column_coord, brwt_coord); the row-diff form only by count
TEST(CoordToHeader, RecordMappingCheckAtLoad) {
    // A: two records of 26 k-mers each, B: one of 8 (k = 5)
    const std::vector<test::CoordinateRecord> records {
        { "ACGTTGCATGACCTAGGATCCGATTACGGA", "A", "a1" },
        { "TTGACCGTAGGCATCAGTCAATGGCTTCAG", "A", "a2" },
        { "GGCATTACCAGT", "B", "b1" },
    };
    auto mapping = [](bool swapped) {
        return swapped ? CoordToHeader({ { "b1" }, { "a1", "a2" } }, { { 8 }, { 26, 26 } })
                       : CoordToHeader({ { "a1", "a2" }, { "b1" } }, { { 26, 26 }, { 8 } });
    };
    for (const auto &[row_diff, brwt] : std::vector<std::pair<bool, bool>> {
            { false, false }, { false, true }, { true, false }, { true, true } }) {
        SCOPED_TRACE(std::string(row_diff ? "row_diff " : "") + (brwt ? "brwt" : "column"));
        auto index = test::coordinate_test_graph(5, records, row_diff, brwt);
        const auto &annotation = index->get_annotator();
        ASSERT_EQ((std::vector<std::string> { "A", "B" }),
                  annotation.get_label_encoder().get_labels());
        // the right mapping, without and with its labels
        CoordToHeader right = mapping(false);
        EXPECT_EQ("", cli::check_record_mapping(annotation, &right));
        EXPECT_FALSE(right.column_order_verified());
        right.set_column_labels({ "A", "B" });
        EXPECT_EQ("", cli::check_record_mapping(annotation, &right));
        EXPECT_TRUE(right.column_order_verified());
        // another number of columns
        CoordToHeader one({ { "a1" } }, { { 26 } });
        EXPECT_EQ("it has 1 columns, the annotation 2 labels",
                  cli::check_record_mapping(annotation, &one));
        // the columns swapped, their labels stated
        CoordToHeader swapped_labelled = mapping(true);
        swapped_labelled.set_column_labels({ "B", "A" });
        EXPECT_EQ("its column 0 is 'B', the annotation's column 0 'A' (the files were passed in "
                  "another order or under other names)",
                  cli::check_record_mapping(annotation, &swapped_labelled));
        // the columns swapped, without labels: A's 52 coordinates do not fit B's 8 k-mers where
        // the annotation states them per column; the row-diff forms check the count only
        CoordToHeader swapped = mapping(true);
        if (row_diff) {
            EXPECT_EQ("", cli::check_record_mapping(annotation, &swapped));
        } else {
            EXPECT_EQ("column 0 ('A') of the annotation holds 52 coordinates, more than the 8 "
                      "k-mers the mapping states for it (the files were passed in another order)",
                      cli::check_record_mapping(annotation, &swapped));
        }
        EXPECT_FALSE(swapped.column_order_verified());
    }
}

// An index is never built with a mapping of another column count: no query has to check it
TEST(CoordToHeader, AnIndexRefusesAMappingOfAnotherColumnCount) {
    auto dbg = std::make_shared<graph::DBGSuccinct>(5);
    dbg->add_sequence("ACGTTGCATGACCTAGGA");
    auto annotation = [&]() {
        auto anno = std::make_unique<ColumnCompressed<>>(dbg->max_index());
        anno->add_labels({ 0 }, { "A", "B" });
        return anno;
    };
    EXPECT_THROW(graph::AnnotatedDBG(dbg, annotation(), false,
                                     std::make_shared<CoordToHeader>(
                                             std::vector<std::vector<std::string>> { { "a1" } },
                                             std::vector<std::vector<uint64_t>> { { 3 } })),
                 std::invalid_argument);
    EXPECT_NO_THROW(graph::AnnotatedDBG(dbg, annotation(), false,
                                        std::make_shared<CoordToHeader>(
                                                std::vector<std::vector<std::string>> {
                                                        { "a1" }, { "b1" } },
                                                std::vector<std::vector<uint64_t>> {
                                                        { 3 }, { 4 } })));
}

TEST(CoordToHeader, FindDuplicateHeadersNamesEveryRepeatAgainstTheFirst) {
    // a column with one header twice, one three times, a repeat that is not consecutive
    // and a header holding a byte that is not UTF-8 (a name, not a duplicate)
    const std::vector<std::string> headers
            = { "a", "b", "a", "c", "b", "b", std::string("d\xff", 2), "e" };
    std::vector<std::pair<size_t, size_t>> repeats;
    CoordToHeader::find_duplicate_headers(headers, [&](size_t first, size_t repeat) {
        repeats.emplace_back(first, repeat);
        return true;
    });
    EXPECT_EQ((std::vector<std::pair<size_t, size_t>>{ { 0, 2 }, { 1, 4 }, { 1, 5 } }), repeats);

    // the callback stops the scan
    repeats.clear();
    CoordToHeader::find_duplicate_headers(headers, [&](size_t first, size_t repeat) {
        repeats.emplace_back(first, repeat);
        return false;
    });
    EXPECT_EQ((std::vector<std::pair<size_t, size_t>>{ { 0, 2 } }), repeats);

    // distinct headers, an empty column: no repeat
    size_t calls = 0;
    CoordToHeader::find_duplicate_headers({ "x", "y", "" }, [&](size_t, size_t) {
        return ++calls, true;
    });
    CoordToHeader::find_duplicate_headers({}, [&](size_t, size_t) { return ++calls, true; });
    EXPECT_EQ(0u, calls);

    // find_header() answers a repeated header with its first record: what the scan reports
    CoordToHeader cth({ headers }, { { 1, 1, 1, 1, 1, 1, 1, 1 } });
    EXPECT_EQ(std::make_pair(CoordToHeader::Column(0), size_t(0)), *cth.find_header("a"));
    EXPECT_EQ(std::make_pair(CoordToHeader::Column(0), size_t(1)), *cth.find_header("b"));
    EXPECT_EQ(5u, cth.build_header_index());
}

// The header index (find_header): a sorted array of (hash, column, seq_id), searched by hash
// and settled by comparing the headers. Every test below holds for any hash function; the weak
// ones force the collisions a 64-bit hash makes rare.

// every header to one hash: one run of equal hashes, every lookup settled by the text
uint64_t constant_hash(std::string_view) { return 7; }
// the largest hash: the last run of the array, where a merge must still end
uint64_t max_hash(std::string_view) { return std::numeric_limits<uint64_t>::max(); }
// a few distinct hashes, each shared by many headers
uint64_t length_hash(std::string_view header) { return header.size(); }
// a hash that differs from the others, to merge with
uint64_t first_byte_hash(std::string_view header) {
    return header.empty() ? 0 : static_cast<unsigned char>(header[0]);
}
// 61 hashes over many headers: runs long enough to collide often, short enough to search
uint64_t weak_hash(std::string_view header) {
    return CoordToHeader::default_header_hash(header) % 61;
}

// what find_header must answer: the first occurrence in column order, each column in seq_id
// order
std::map<std::string, std::pair<CoordToHeader::Column, size_t>>
first_occurrences(const std::vector<std::vector<std::string>> &headers) {
    std::map<std::string, std::pair<CoordToHeader::Column, size_t>> first;
    for (size_t col = 0; col < headers.size(); ++col) {
        for (size_t s = 0; s < headers[col].size(); ++s) {
            first.emplace(headers[col][s], std::make_pair(CoordToHeader::Column(col), s));
        }
    }
    return first;
}

CoordToHeader make_cth(const std::vector<std::vector<std::string>> &headers) {
    std::vector<std::vector<uint64_t>> num_kmers;
    for (const auto &column : headers) {
        num_kmers.emplace_back(column.size(), 1);
    }
    return CoordToHeader(std::vector<std::vector<std::string>>(headers), std::move(num_kmers));
}

std::unique_ptr<CoordToHeader> new_cth(const std::vector<std::vector<std::string>> &headers) {
    std::vector<std::vector<uint64_t>> num_kmers;
    for (const auto &column : headers) {
        num_kmers.emplace_back(column.size(), 1);
    }
    return std::make_unique<CoordToHeader>(std::vector<std::vector<std::string>>(headers),
                                           std::move(num_kmers));
}

void expect_index_answers(const CoordToHeader &cth,
                          const std::vector<std::vector<std::string>> &headers,
                          const std::vector<std::string> &absent) {
    const auto first = first_occurrences(headers);
    for (const auto &[header, where] : first) {
        auto found = cth.find_header(header);
        ASSERT_TRUE(found) << header;
        EXPECT_EQ(where, *found) << header;
    }
    for (const std::string &header : absent) {
        ASSERT_EQ(0u, first.count(header)) << "not absent: " << header;
        EXPECT_FALSE(cth.find_header(header)) << header;
    }
    EXPECT_EQ(first.size(), cth.build_header_index());
}

TEST(CoordToHeaderIndex, FindsEveryHeaderAndNoOther) {
    const std::vector<std::vector<std::string>> headers {
        { "NC_000913.3", "NZ_CP009072.1", "plasmid pX" },
        {},
        { "contig_1", std::string("bin\0ary", 7), "", "contig_10" },
        { "x" },
    };
    auto cth = make_cth(headers);
    EXPECT_EQ(0u, cth.header_index_bytes());
    expect_index_answers(cth, headers, { "NC_000913", "NC_000913.30", "contig_", "contig_100",
                                         std::string("bin", 3), "X", " x", "plasmid" });
    // 16 bytes per header, the empty header included
    EXPECT_EQ(8u * 16, cth.header_index_bytes());
    EXPECT_EQ(1u, cth.num_header_index_builds());

    // an index without headers answers nothing
    auto empty = make_cth({ {}, {} });
    EXPECT_FALSE(empty.find_header(""));
    EXPECT_FALSE(empty.find_header("x"));
    EXPECT_EQ(0u, empty.build_header_index());
    EXPECT_EQ(0u, empty.header_index_bytes());
    CoordToHeader none;
    EXPECT_FALSE(none.find_header("x"));
    EXPECT_EQ(0u, none.build_header_index());
}

// A header in several columns, or twice in one column, answers its first occurrence in column
// order: the earliest column, there the earliest record, whatever the records' positions in
// the later columns
TEST(CoordToHeaderIndex, DuplicatesAnswerTheFirstOccurrenceInColumnOrder) {
    const std::vector<std::vector<std::string>> headers {
        { "b", "c", "z" },
        { "z", "a", "c" },
        { "a", "d", "a" },
        { "e", "z", "e" },
    };
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash, &constant_hash,
                                            &max_hash, &length_hash }) {
        auto cth = make_cth(headers);
        cth.set_header_hash(hash);
        using P = std::pair<CoordToHeader::Column, size_t>;
        EXPECT_EQ(P(0, 2), *cth.find_header("z"));
        EXPECT_EQ(P(1, 1), *cth.find_header("a"));
        EXPECT_EQ(P(0, 1), *cth.find_header("c"));
        EXPECT_EQ(P(2, 1), *cth.find_header("d"));
        EXPECT_EQ(P(3, 0), *cth.find_header("e"));
        EXPECT_FALSE(cth.find_header("f"));
        EXPECT_EQ(6u, cth.build_header_index());
    }
}

// Forced collisions: with every header on one hash, on the largest hash or on its length, every
// lookup is still exact (present, absent, first occurrence), and the number of distinct
// headers is counted on the headers, not on the hashes
TEST(CoordToHeaderIndex, HashCollisionsAreSettledByTheHeaders) {
    const std::vector<std::vector<std::string>> headers {
        { "aa", "bb", "cc", "a" },
        { "dd", "bb", "b", "ee" },
        { "aa", "c", "ff", "" },
    };
    for (CoordToHeader::HeaderHash hash : { &constant_hash, &max_hash, &length_hash,
                                            &first_byte_hash }) {
        auto cth = make_cth(headers);
        cth.set_header_hash(hash);
        expect_index_answers(cth, headers, { "ab", "ba", "d", "e", "aaa", "gg", "f" });
    }
}

// Setting the hash drops a built index (it is sorted by the old one); a copy hashes as its
// original and builds its own index
TEST(CoordToHeaderIndex, SettingTheHashRebuildsTheIndex) {
    auto cth = make_cth({ { "a", "b" }, { "c", "a" } });
    EXPECT_EQ(3u, cth.build_header_index());
    EXPECT_EQ(1u, cth.num_header_index_builds());
    cth.set_header_hash(&constant_hash);
    EXPECT_EQ(0u, cth.header_index_bytes());
    EXPECT_EQ(std::make_pair(CoordToHeader::Column(1), size_t(0)), *cth.find_header("c"));
    EXPECT_EQ(2u, cth.num_header_index_builds());

    CoordToHeader copy(cth);
    EXPECT_EQ(0u, copy.num_header_index_builds());
    EXPECT_EQ(std::make_pair(CoordToHeader::Column(0), size_t(0)), *copy.find_header("a"));
    // the copy's index merges with the original's: both are sorted by constant_hash
    size_t calls = 0;
    EXPECT_EQ(3u, CoordToHeader::find_repeated_headers({ &cth, &copy },
            [&](std::string_view, uint64_t, const auto &) { ++calls; }));
    EXPECT_EQ(3u, calls);
}

// A larger index against the reference (a std::map of first occurrences), built on several
// threads with the default hash and with weak ones: random headers, repeated within and
// across columns, columns of very different sizes, empty ones
TEST(CoordToHeaderIndex, AgreesWithTheReferenceOnRandomHeaders) {
    std::mt19937 rng(42);
    std::vector<std::vector<std::string>> headers(37);
    std::uniform_int_distribution<int> id(0, 40000);
    for (size_t col = 0; col < headers.size(); ++col) {
        const size_t n = col % 5 == 4 ? 0 : (col == 7 ? 100000 : rng() % 3000);
        for (size_t s = 0; s < n; ++s) {
            headers[col].push_back("seq_" + std::to_string(id(rng)));
        }
    }
    std::vector<std::string> absent;
    for (int i = 40001; i < 40100; ++i) {
        absent.push_back("seq_" + std::to_string(i));
    }
    absent.push_back("seq_");
    absent.push_back("");
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash,
                                            &weak_hash }) {
        for (size_t threads : { 1, 4 }) {
            auto cth = make_cth(headers);
            cth.set_header_hash(hash);
            const auto first = first_occurrences(headers);
            EXPECT_EQ(first.size(), cth.build_header_index(threads));
            expect_index_answers(cth, headers, absent);
            EXPECT_EQ(1u, cth.num_header_index_builds());
        }
    }
}

// Concurrent first lookups build the index once, and every one of them is answered from it
TEST(CoordToHeaderIndex, ConcurrentLookupsBuildItOnce) {
    std::vector<std::vector<std::string>> headers(4);
    for (size_t i = 0; i < 20000; ++i) {
        headers[i % 4].push_back("h" + std::to_string(i));
    }
    auto cth = make_cth(headers);
    std::atomic<size_t> wrong { 0 };
    std::vector<std::thread> threads;
    for (size_t t = 0; t < 8; ++t) {
        threads.emplace_back([&, t]() {
            for (size_t i = t; i < 20000; i += 97) {
                auto found = cth.find_header("h" + std::to_string(i));
                if (!found || found->first != i % 4 || found->second != i / 4)
                    ++wrong;
            }
        });
    }
    for (auto &thread : threads) {
        thread.join();
    }
    EXPECT_EQ(0u, wrong.load());
    EXPECT_EQ(1u, cth.num_header_index_builds());
}

// A reader holds the array it searches: while another thread swaps the hash again and again
// (set_header_hash drops the index, the next lookup builds one sorted by the new hash), every
// lookup is answered from one whole index, the one it took, and answers right
TEST(CoordToHeaderIndex, ReadersKeepTheirIndexWhileTheHashIsSwapped) {
    std::vector<std::vector<std::string>> headers(3);
    for (size_t i = 0; i < 3000; ++i) {
        headers[i % 3].push_back("h" + std::to_string(i));
    }
    auto cth = make_cth(headers);
    std::atomic<bool> done { false };
    std::atomic<size_t> wrong { 0 };
    std::atomic<size_t> lookups { 0 };
    std::vector<std::thread> readers;
    for (size_t t = 0; t < 4; ++t) {
        readers.emplace_back([&, t]() {
            for (size_t round = 0; !done || round == 0; ++round) {
                for (size_t i = t; i < 3000; i += 41) {
                    auto found = cth.find_header("h" + std::to_string(i));
                    if (!found || found->first != i % 3 || found->second != i / 3)
                        ++wrong;
                    if (cth.find_header("h" + std::to_string(3000 + i)))
                        ++wrong;
                    lookups += 2;
                }
            }
        });
    }
    const std::vector<CoordToHeader::HeaderHash> hashes {
        &weak_hash, &CoordToHeader::default_header_hash, &length_hash, &first_byte_hash };
    for (size_t swap = 0; swap < 200; ++swap) {
        cth.set_header_hash(hashes[swap % hashes.size()]);
        // a lookup between swaps, so that each hash's index is built and read
        auto found = cth.find_header("h" + std::to_string(swap));
        if (!found || found->first != swap % 3 || found->second != swap / 3)
            ++wrong;
    }
    done = true;
    for (auto &reader : readers) {
        reader.join();
    }
    EXPECT_EQ(0u, wrong.load());
    EXPECT_LT(0u, lookups.load());
    EXPECT_LE(200u, cth.num_header_index_builds());
    // the last hash's index answers as a fresh one
    expect_index_answers(cth, headers, { "h3000", "h", "" });
}

// The index stores a column and a seq_id in 32 bits: up to 2^32 columns, and 2^32 sequences
// in a column, fit; one more is refused, naming the limit
TEST(CoordToHeaderIndex, RefusesWhatExceedsItsThirtyTwoBitFields) {
    const uint64_t limit = uint64_t(1) << 32;
    EXPECT_EQ(limit, CoordToHeader::kMaxHeaderIndexColumns);
    EXPECT_EQ(limit, CoordToHeader::kMaxHeaderIndexRecords);
    EXPECT_NO_THROW(CoordToHeader::check_header_index_limits(0, 0));
    EXPECT_NO_THROW(CoordToHeader::check_header_index_limits(limit, limit));
    EXPECT_THROW(CoordToHeader::check_header_index_limits(limit + 1, 1), std::length_error);
    EXPECT_THROW(CoordToHeader::check_header_index_limits(1, limit + 1), std::length_error);
    try {
        CoordToHeader::check_header_index_limits(1, limit + 1);
        FAIL() << "not refused";
    } catch (const std::length_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("4294967297 sequences"))
                << e.what();
        EXPECT_NE(std::string::npos, std::string(e.what()).find("limit of 4294967296"))
                << e.what();
    }
    try {
        CoordToHeader::check_header_index_limits(limit + 1, 0);
        FAIL() << "not refused";
    } catch (const std::length_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("4294967297 columns"))
                << e.what();
    }
}

// find_repeated_headers: the headers that occur more than once over several objects, each once
// with every occurrence in (index, column, seq_id) order
using Occurrences = std::vector<CoordToHeader::HeaderOccurrence>;
std::map<std::string, Occurrences>
repeated_headers(const std::vector<const CoordToHeader*> &indexes, uint64_t *count = nullptr) {
    std::map<std::string, Occurrences> repeated;
    const uint64_t n = CoordToHeader::find_repeated_headers(indexes,
            [&](std::string_view header, uint64_t count, const Occurrences &occurrences) {
                EXPECT_EQ(count, occurrences.size()) << header;
                EXPECT_TRUE(repeated.emplace(std::string(header), occurrences).second)
                        << "reported twice: " << header;
            });
    EXPECT_EQ(repeated.size(), n);
    if (count)
        *count = n;
    return repeated;
}

TEST(CoordToHeaderRepeats, UniqueHeadersHaveNone) {
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash, &constant_hash,
                                            &max_hash, &length_hash }) {
        auto a = new_cth({ { "a1", "a2" }, { "a3" } });
        auto b = new_cth({ { "b1" }, {}, { "b2", "b3", "b44" } });
        auto c = new_cth({});
        a->set_header_hash(hash);
        b->set_header_hash(hash);
        c->set_header_hash(hash);
        EXPECT_TRUE(repeated_headers({ a.get(), b.get(), c.get() }).empty());
        EXPECT_TRUE(repeated_headers({ a.get() }).empty());
        EXPECT_TRUE(repeated_headers({}).empty());
    }
}

TEST(CoordToHeaderRepeats, NamesEveryRepeatWithAllItsOccurrences) {
    using O = CoordToHeader::HeaderOccurrence;
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash, &constant_hash,
                                            &max_hash, &length_hash, &first_byte_hash }) {
        // "x" across objects, "y" across the columns of one, "w" twice in one column, "v" in
        // all three; "u" only once (in a collision with others under the weak hashes)
        auto a = new_cth({ { "x", "y", "v" }, { "y", "w", "w" } });
        auto b = new_cth({ { "v" }, { "z" } });
        auto c = new_cth({ { "u", "x" }, { "v", "x" } });
        for (auto *cth : { a.get(), b.get(), c.get() }) {
            cth->set_header_hash(hash);
        }
        // a null entry holds no headers, and keeps the positions of the others
        const auto repeated = repeated_headers({ a.get(), nullptr, b.get(), c.get() });
        const std::map<std::string, Occurrences> expected {
            { "x", { O { 0, 0, 0 }, O { 3, 0, 1 }, O { 3, 1, 1 } } },
            { "y", { O { 0, 0, 1 }, O { 0, 1, 0 } } },
            { "w", { O { 0, 1, 1 }, O { 0, 1, 2 } } },
            { "v", { O { 0, 0, 2 }, O { 2, 0, 0 }, O { 3, 1, 0 } } },
        };
        EXPECT_EQ(expected.size(), repeated.size());
        for (const auto &[header, occurrences] : expected) {
            ASSERT_TRUE(repeated.count(header)) << header;
            EXPECT_EQ(occurrences, repeated.at(header)) << header;
        }
        // each object alone: its own repeats only
        EXPECT_EQ((std::vector<std::string>{ "w", "y" }),
                  [&]() {
                      std::vector<std::string> keys;
                      for (const auto &[header, _] : repeated_headers({ a.get() })) {
                          keys.push_back(header);
                      }
                      return keys;
                  }());
        EXPECT_EQ(1u, repeated_headers({ c.get() }).size());
        EXPECT_TRUE(repeated_headers({ b.get() }).empty());
    }
}

// The merge reads each object's index (building it once if needed) and refuses objects sorted
// by different hashes
TEST(CoordToHeaderRepeats, BuildsTheIndexesAndRefusesMixedHashes) {
    auto a = new_cth({ { "x" } });
    auto b = new_cth({ { "x" } });
    EXPECT_EQ(1u, repeated_headers({ a.get(), b.get() }).size());
    EXPECT_EQ(1u, a->num_header_index_builds());
    EXPECT_EQ(1u, repeated_headers({ a.get(), b.get() }).size());
    EXPECT_EQ(1u, b->num_header_index_builds());
    b->set_header_hash(&length_hash);
    EXPECT_THROW(CoordToHeader::find_repeated_headers({ a.get(), b.get() },
                    [](std::string_view, uint64_t, const auto &) {}),
                 std::invalid_argument);
}

// A header repeated in many records is counted in full and listed up to |max_listed|, its first
// occurrences in (index, column, seq_id) order, in one object and over several
TEST(CoordToHeaderRepeats, ListsTheFirstOccurrencesOfAFrequentHeader) {
    using O = CoordToHeader::HeaderOccurrence;
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash,
                                            &constant_hash }) {
        auto a = new_cth({ std::vector<std::string>(1000, "same"), { "same", "other" } });
        auto b = new_cth({ { "other" }, std::vector<std::string>(5, "same") });
        a->set_header_hash(hash);
        b->set_header_hash(hash);
        std::map<std::string, std::pair<uint64_t, Occurrences>> found;
        const auto collect = [&](std::string_view header, uint64_t count,
                                 const Occurrences &listed) {
            found[std::string(header)] = { count, listed };
        };
        EXPECT_EQ(1u, CoordToHeader::find_repeated_headers({ a.get() }, collect, 3));
        EXPECT_EQ(1001u, found.at("same").first);
        EXPECT_EQ((Occurrences { O { 0, 0, 0 }, O { 0, 0, 1 }, O { 0, 0, 2 } }),
                  found.at("same").second);
        found.clear();
        EXPECT_EQ(2u, CoordToHeader::find_repeated_headers({ b.get(), a.get() }, collect, 2));
        EXPECT_EQ(1006u, found.at("same").first);
        EXPECT_EQ((Occurrences { O { 0, 1, 0 }, O { 0, 1, 1 } }), found.at("same").second);
        EXPECT_EQ(2u, found.at("other").first);
        EXPECT_EQ((Occurrences { O { 0, 0, 0 }, O { 1, 1, 1 } }), found.at("other").second);
        found.clear();
        EXPECT_EQ(2u, CoordToHeader::find_repeated_headers({ b.get(), a.get() }, collect, 0));
        EXPECT_TRUE(found.at("same").second.empty());
    }
}

// Many objects with random headers against a reference count, with the default and a weak hash
TEST(CoordToHeaderRepeats, AgreesWithTheReferenceOnRandomHeaders) {
    std::mt19937 rng(7);
    std::vector<std::vector<std::vector<std::string>>> all(9);
    std::map<std::string, Occurrences> expected;
    std::map<std::string, Occurrences> seen;
    for (size_t i = 0; i < all.size(); ++i) {
        all[i].resize(1 + rng() % 4);
        for (size_t col = 0; col < all[i].size(); ++col) {
            const size_t n = rng() % 500;
            for (size_t s = 0; s < n; ++s) {
                all[i][col].push_back("r" + std::to_string(rng() % 6000));
                seen[all[i][col].back()].push_back({ i, static_cast<CoordToHeader::Column>(col), s });
            }
        }
    }
    for (const auto &[header, occurrences] : seen) {
        if (occurrences.size() > 1)
            expected.emplace(header, occurrences);
    }
    ASSERT_LT(10u, expected.size());
    for (CoordToHeader::HeaderHash hash : { &CoordToHeader::default_header_hash,
                                            &weak_hash, &length_hash }) {
        std::vector<std::unique_ptr<CoordToHeader>> objects;
        std::vector<const CoordToHeader*> pointers;
        for (const auto &headers : all) {
            objects.push_back(new_cth(headers));
            objects.back()->set_header_hash(hash);
            pointers.push_back(objects.back().get());
        }
        EXPECT_EQ(expected, repeated_headers(pointers));
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

// Test Alignment::format_coords(cth, k) renders per-header ranges
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

    EXPECT_EQ(aln.format_coords(cth, 5), "accession_A/14:4-7;accession_B/24:3-6");
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

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA/9:3-5;seqB/9:3-5");
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
    EXPECT_EQ(aln.format_coords(cth, 5), "ref/24:4-5:11-12");
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
    aln.label_coordinates = { { 7 } };  // Three k-mers in seqA, one in seqB.

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA/14:8-14;seqB/14:1-5");
}

TEST(AlignmentFormatCoords, PartialLastKmer) {
    using namespace mtg::graph::align;
    CoordToHeader cth({ { "seqA", "seqB" } }, { { 5, 5 } });
    for (size_t length = 0; length <= 5; ++length) {
        Alignment aln(std::string_view{}, {}, std::string(length, 'A'), 0, {}, 0, false, 0);
        aln.label_columns = { 0 };
        aln.label_coordinates = { { 4 } };
        EXPECT_EQ(aln.format_coords(cth, 5),
                  length ? "seqA/9:5-" + std::to_string(4 + length) : "");
    }
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
        std::string("ACGTACGTACGT"),  // 12 nt, crosses both boundaries
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 3 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA/7:4-7;seqB/7:1-7;seqC/7:1-5");
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
    //   coord 6 -> two k-mers in seqA, two in seqB
    //   coord 7 -> one k-mer in seqA, three in seqB
    aln.label_coordinates = { { 6, 7 } };

    EXPECT_EQ(aln.format_coords(cth, 5), "seqA/12:7-12:8-12;seqB/12:1-6:1-7");
}

TEST(AlignmentFormatCoords, OverflowPastLastHeader) {
    using namespace mtg::graph::align;

    // The ten k-mers spell 12 nt; alignment of 20 nt overruns by 8 nt,
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

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA/7:1-7;seqB/7:1-7");
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

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA/7:3-7");
}

TEST(AlignmentFormatCoords, SpansEntireColumn) {
    using namespace mtg::graph::align;

    // The nine k-mers spell 11 nt with k=3. Each header's covered range
    // includes the shared boundary bases.
    CoordToHeader cth(
        { { "seqA", "seqB", "seqC" } },
        { { 3, 4, 2 } }
    );

    Alignment aln(
        std::string_view{}, {},
        std::string("ACGTACGTACG"),  // 11 nt
        0, {}, 0, false, 0
    );
    aln.label_columns = { 0 };
    aln.label_coordinates = { { 0 } };

    EXPECT_EQ(aln.format_coords(cth, 3), "seqA/5:1-5;seqB/6:1-6;seqC/4:1-4");
}

} // namespace
