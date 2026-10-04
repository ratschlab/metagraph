#ifndef __COORD_TO_HEADER_HPP__
#define __COORD_TO_HEADER_HPP__

#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vectors/bit_vector_sd.hpp"


namespace mtg {
namespace annot {

/**
 * An index mapping k-mer coordinates in annotation columns to sequence headers.
 * A column is typically the output of `annotate --anno-filename` (all sequences
 * in a FASTA file are collapsed into one label, with consecutive global
 * coordinates). This index remembers, for each column, the FASTA header of
 * each sequence and how many k-mers it contributed, so downstream tools can
 * report per-sequence labels instead of whole-file offsets:
 *
 *   metagraph query (coords): transforms
 *       "<file.fa>:0-0-3:0-23-27"   (pos-first[-last], tab-separated labels)
 *     into
 *       "<seq1>:0-0-3\t<seq3>:0-1-5"
 *
 *   metagraph align (coord-aware): transforms
 *       "file.fa:4-30"              (label:1-based-inclusive nt range)
 *     into
 *       "seq1:4-10;seq3:1-5"        (';'-separated per-sequence ranges)
 *
 * Terminology: a "header" is the sequence's label (FASTA header line), and
 * `seq_id` is the 0-based index of a sequence within its column.
 *
 * Built by `annotate --index-header-coords` and stored alongside the column
 * annotation as a ".seqs" file.
 */
class CoordToHeader {
  public:
    using Column = mtg::annot::matrix::BinaryMatrix::Column;
    using Tuple = mtg::annot::matrix::MultiIntMatrix::Tuple;
    using RowTuples = mtg::annot::matrix::MultiIntMatrix::RowTuples;

    // out of line, as the destructor: the header index type is private to the .cpp
    CoordToHeader();
    // Constructor from sequence names and number of k-mers in them, per annotation column.
    CoordToHeader(std::vector<std::vector<std::string>> &&headers,
                  std::vector<std::vector<uint64_t>> &&num_kmers);
    ~CoordToHeader();
    // a copy or a moved-to object starts without a header index of its own (find_header)
    CoordToHeader(const CoordToHeader &other);
    CoordToHeader(CoordToHeader &&other);
    CoordToHeader& operator=(const CoordToHeader &other);
    CoordToHeader& operator=(CoordToHeader &&other);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    // global_coord -> (seq_id, local_coord)
    std::pair<size_t, uint64_t> map_single_coord(Column col, uint64_t coord) const;

    // The sequence a global coordinate of a column falls in, with the first and last global
    // coordinates of that sequence: every coordinate c in [first, last] maps to (seq_id,
    // c - first), so a caller walking sorted coordinates maps a whole sequence's run for
    // the cost of one map_single_coord. Throws as map_single_coord does.
    struct SequenceRange {
        size_t seq_id;
        uint64_t first;
        uint64_t last;
    };
    SequenceRange sequence_range(Column col, uint64_t coord) const;
    // The last global coordinate of a sequence (one select; its first is the previous
    // sequence's last + 1). Throws as num_kmers_in_sequence does.
    uint64_t last_coord(Column col, size_t seq_id) const;

    /**
     * Batch variant of `map_single_coord` that replaces each global coord
     * in `rows_tuples` with a packed `(seq_id, local_coord)` pair
     * encoded as a single `uint64_t`. Used by the query pipeline to reuse
     * the existing `uint64_t`-valued row/tuple data structures instead of
     * introducing a parallel structure for pairs.
     *
     * Packing:   `packed = local_coord * num_sequences[col] + seq_id`
     * Unpacking: `seq_id   = packed % num_sequences[col]`
     *            `local_coord = packed / num_sequences[col]`
     *
     * For single-coord lookups prefer `map_single_coord`, which returns
     * the pair directly without packing.
     */
    void map_to_local_coords(std::vector<RowTuples> *rows_tuples) const;

    uint64_t num_columns() const { return coord_offsets_.size(); }
    size_t num_sequences(Column column) const { return headers_[column].size(); }
    // Number of k-mers in a sequence.
    uint64_t num_kmers_in_sequence(Column col, size_t seq_id) const;
    // Total number of k-mers / coordinates in a column (across all sequences).
    uint64_t num_kmers(Column column) const { return coord_offsets_[column].size(); }
    // FASTA headers of the sequences in a column, indexed by seq_id.
    const std::vector<std::string>& get_headers(Column column) const { return headers_[column]; }

    /**
     * The (column, seq_id) of a FASTA header: its first occurrence in column order (a
     * header is unique within a column, not within the index), or nothing. The reverse
     * index is built on the first call and is owned by this object, so it lives and dies
     * with the headers it indexes: a server builds it once for its loaded index, and an
     * object constructed later at a freed address can never be answered from another
     * object's index. That is what a process-wide cache keyed by the object's address
     * did, and a fingerprint of a few headers could not tell two such objects apart (GPT
     * review of stage 2, finding 7: [A,X,Y,Z] cached, [A,Y,X,Z] rebuilt at the same
     * address, and X resolved to the sequence now named Y). load() drops the index.
     * Thread-safe.
     */
    std::optional<std::pair<Column, size_t>> find_header(std::string_view header) const;
    // how many times this object built its reverse index (to test that it is built once)
    size_t num_header_index_builds() const;

    static constexpr auto kExtension = ".seqs";

  private:
    struct HeaderIndex;

    std::vector<std::vector<std::string>> headers_;
    std::vector<bit_vector_sd> coord_offsets_;

    // find_header(): keys are views into |headers_|, which only load() changes, and load()
    // drops the index first
    mutable std::mutex header_index_mutex_;
    mutable std::unique_ptr<const HeaderIndex> header_index_;
    mutable size_t header_index_builds_ = 0;
};

} // namespace annot
} // namespace mtg

#endif // __COORD_TO_HEADER_HPP__
