// A coordinate annotation cut short is refused, whether the load maps the file or not.
//
// Until master's #727, TupleCSCMatrix kept its column values in sdsl::int_vector<>, and an owned
// int_vector's load from a short stream only sets failbit. So a load without mmap, the path of the
// server's per-request in_ram copies, accepted a row_diff_brwt_coord annotation cut inside its
// coordinate values: its last 8 bytes or 4 KiB cut off, `stats -a` answered rc 0 with the tail
// unread. The values are now sdsl::packed_vector_view, whose load throws on a short read, and
// TupleCSCMatrix::load_tuples turns the throw into a refusal. This test keeps that so.
#include "gtest/gtest.h"
#include "../test_helpers.hpp"
#include "test_matrix_helpers.hpp"

#include <fstream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

const std::string test_dump_basename = test_dump_dir() + "/annotation_cut";

// Sets the global mmap mode and restores it on scope exit, also after a failed assertion.
class MmapMode {
  public:
    explicit MmapMode(bool mapped) : old_(utils::with_mmap()) { utils::set_mmap(mapped); }
    ~MmapMode() { utils::set_mmap(old_); }

  private:
    bool old_;
};

// The deployed coordinate form (row_diff_brwt_coord): every cut inside the coordinates is refused,
// by a load without mmap (as the in_ram copies load) and by a mapped load.
TEST(AnnotationCut, RowDiffBRWTCoordCutInItsCoordinatesIsRefused) {
    using Annotator = annot::RowDiffBRWTCoordAnnotator;
    using Matrix = annot::matrix::TupleRowDiff<annot::matrix::TupleCSCMatrix<annot::matrix::BRWT>>;
    const uint64_t num_rows = 3000;
    const size_t num_columns = 5;
    std::mt19937_64 rng(5);
    BitVectorPtrArray columns;
    std::vector<bit_vector_smart> delimiters;
    std::vector<sdsl::int_vector<>> values;
    for (size_t j = 0; j < num_columns; ++j) {
        sdsl::bit_vector bits(num_rows, 0);
        // per set row a one, then a zero per coordinate; a one at the end
        std::vector<bool> delims;
        std::vector<uint64_t> coords;
        for (uint64_t i = 0; i < num_rows; ++i) {
            // the last row is set in every column, so the file ends in coordinate values
            if (rng() % 3 && i + 1 < num_rows)
                continue;
            bits[i] = 1;
            delims.push_back(true);
            for (size_t t = 1 + rng() % 3; t > 0; --t) {
                coords.push_back(rng() % (1u << 20));
                delims.push_back(false);
            }
        }
        delims.push_back(true);
        columns.emplace_back(new bit_vector_stat(std::move(bits)));
        sdsl::bit_vector delim_bits(delims.size(), 0);
        for (size_t i = 0; i < delims.size(); ++i)
            delim_bits[i] = delims[i];
        delimiters.emplace_back(delim_bits);
        sdsl::int_vector<> column_values(coords.size(), 0, 20);
        for (size_t i = 0; i < coords.size(); ++i)
            column_values[i] = coords[i];
        values.push_back(std::move(column_values));
    }
    auto matrix = std::make_unique<Matrix>(
            nullptr, build_matrix_from_columns<annot::matrix::BRWT>(columns, num_rows),
            std::move(delimiters), std::move(values));

    std::vector<uint64_t> rows(num_rows);
    std::iota(rows.begin(), rows.end(), 0);
    // the coordinates as stored (no graph, so no row-diff decoding)
    const auto expected = matrix->diffs().get_row_tuples(rows);
    std::stringstream tuples;
    matrix->diffs().serialize_tuples(tuples);
    const size_t tuple_bytes = tuples.str().size();

    annot::LabelEncoder<std::string> labels;
    for (size_t j = 0; j < num_columns; ++j) {
        labels.insert_and_encode("column " + std::to_string(j));
    }
    const std::string filename = test_dump_basename + Annotator::kExtension;
    Annotator(std::move(matrix), labels).serialize(filename);
    std::ifstream in(filename, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    // the coordinates are the file's tail, long enough for every cut below to fall inside them
    ASSERT_LT(tuple_bytes, bytes.size());
    ASSERT_GT(tuple_bytes, 3 * 4096u);

    for (bool mapped : { false, true }) {
        MmapMode mode(mapped);
        Annotator loaded;
        ASSERT_TRUE(loaded.load(filename)) << (mapped ? "mapped" : "not mapped");
        EXPECT_EQ(expected, loaded.get_matrix().diffs().get_row_tuples(rows));
    }

    const std::string cut_filename = test_dump_basename + "_cut" + Annotator::kExtension;
    for (size_t cut : { size_t(1), size_t(8), size_t(4096), tuple_bytes / 2, tuple_bytes }) {
        {
            std::ofstream out(cut_filename, std::ios::binary);
            out.write(bytes.data(), bytes.size() - cut);
        }
        for (bool mapped : { false, true }) {
            MmapMode mode(mapped);
            Annotator loaded;
            EXPECT_FALSE(loaded.load(cut_filename))
                    << "cut " << cut << " of " << bytes.size() << " bytes, "
                    << (mapped ? "mapped" : "not mapped");
        }
    }
}

} // namespace
