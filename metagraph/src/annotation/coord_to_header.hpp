#ifndef __COORD_TO_HEADER_HPP__
#define __COORD_TO_HEADER_HPP__

#include <functional>
#include <limits>
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
 * annotation as a ".seqs" file: the number of columns, then per column its headers
 * and coordinate offsets, then (since H4, train 3) a trailer stating the label of
 * each column (kColumnLabelsMarker, then the labels), which binaries before it never
 * read: they stop after the columns.
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
    // a copy starts without a header index of its own (find_header); not assignable
    CoordToHeader(const CoordToHeader &other);

    // load() reads the trailer of column labels when the file has one (column_labels()); data
    // after the columns that is not that trailer, or a trailer of another number of labels, fails
    // the load. serialize() writes the trailer when the labels are set (set_column_labels)
    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    // The labels of the annotation columns this mapping was built for, in column order, as
    // `annotate --index-header-coords` states them: for its i-th file the label --anno-filename
    // gives it (the file name as passed). Null when none are stated (a file written before
    // H4, or a mapping constructed without them).
    const std::vector<std::string>* column_labels() const {
        return column_labels_ ? &*column_labels_ : nullptr;
    }
    // |labels|: one per column (std::invalid_argument otherwise)
    void set_column_labels(std::vector<std::string> &&labels);
    /**
     * Why this mapping is not the record mapping of an annotation whose labels are |labels|, in
     * column order, or "" when it is (as far as the mapping can tell): another number of
     * columns ("it has n columns, the annotation m labels"), or, when the mapping states its
     * column labels, the first column whose label differs ("its column i is 'x', the
     * annotation's column i 'y' ..."). On success the stated labels are dropped (the annotation
     * holds the same bytes) and column_order_verified() turns true if there were any: a loader
     * calls it once, before the mapping is shared.
     */
    std::string check_columns(const std::vector<std::string> &labels);
    // Whether check_columns() matched every column's stated label with the annotation's: the
    // column order of this mapping was verified at load (SPEC-multigraph-envelope.md section
    // 9.1, seqs_order "verified"). False for a mapping that states no labels: its column count
    // alone was checked ("count_only")
    bool column_order_verified() const { return column_order_verified_; }

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
     * object's index. A process-wide cache keyed by the object's address could do that,
     * and a fingerprint of a few headers cannot tell two such objects apart. load()
     * drops the index. Thread-safe (with each other, not with load()).
     *
     * The reverse index is one array of 16 bytes per header, (64-bit hash of the header,
     * column, seq_id), sorted by (hash, column, seq_id): a lookup is a binary search for the
     * hash, then a comparison of the stored headers with that hash, in column order, against
     * |header|. A hash collision therefore costs one more comparison, never a wrong answer,
     * and the first header that compares equal is the first occurrence in column order. A
     * hash table of views into the headers took about 40 bytes per header (21.5 GB for the
     * 447.6M headers of AllTheBacteria) and a single thread to build.
     */
    std::optional<std::pair<Column, size_t>> find_header(std::string_view header) const;
    // Builds the reverse index now, if not built yet, and returns its number of distinct
    // headers: the server builds it when the index is loaded, so that no request pays it (on
    // refseq33m, 33M headers, a first header name took seconds inside a 1 s budget). The
    // hashing and the sort run on |num_threads| threads (0: get_num_threads()). Throws
    // std::length_error when the columns or the records of a column do not fit the index's
    // 32-bit fields (check_header_index_limits).
    size_t build_header_index(size_t num_threads = 0) const;
    // how many times this object built its reverse index (to test that it is built once)
    size_t num_header_index_builds() const;
    // the bytes the reverse index holds (16 per header), 0 while it is not built
    size_t header_index_bytes() const;

    // The reverse index stores a column and a seq_id in 32 bits each: an index with more
    // columns than kMaxHeaderIndexColumns, or a column with more sequences than
    // kMaxHeaderIndexRecords, is refused with std::length_error, naming the limit, rather
    // than indexed with truncated numbers
    static constexpr uint64_t kMaxHeaderIndexColumns = uint64_t(1) << 32;
    static constexpr uint64_t kMaxHeaderIndexRecords = uint64_t(1) << 32;
    static void check_header_index_limits(uint64_t num_columns, uint64_t max_sequences);

    // The hash of a header the reverse index sorts by. Two headers with equal hashes are
    // told apart by their text, so the hash decides speed, never an answer: a test sets a
    // weak one to force collisions. Setting it drops a built index.
    using HeaderHash = uint64_t (*)(std::string_view header);
    static uint64_t default_header_hash(std::string_view header);
    void set_header_hash(HeaderHash hash);

    // One occurrence of a header among several CoordToHeader objects (find_repeated_headers):
    // the position of its object in the list passed, its column there and its seq_id
    struct HeaderOccurrence {
        size_t index;
        Column column;
        size_t seq_id;

        bool operator==(const HeaderOccurrence &other) const {
            return index == other.index && column == other.column && seq_id == other.seq_id;
        }
    };
    /**
     * The headers that occur more than once over |indexes|: in two of them, in two columns
     * of one, or twice in one column. |callback|(header, count, listed) is called once per such
     * header with its number of occurrences and the first |max_listed| of them, ordered by
     * (index, column, seq_id); the headers come in the order of their hashes. Returns their
     * number. Builds each object's reverse index if not built yet (on |num_threads| threads,
     * as build_header_index), then merges the sorted arrays (a k-way merge over the objects'
     * hashes): equal hashes are compared on the headers, so a hash collision is never
     * reported as a repeat. A header repeated many times costs its count, not a list. Null
     * entries hold no headers. |indexes| lists each object once: an object listed twice
     * repeats every one of its headers (a caller with one object under several names, as a
     * multi-graph server's pair under two names, passes it once). Every object must hash with
     * the same function (std::invalid_argument).
     */
    static uint64_t find_repeated_headers(
            const std::vector<const CoordToHeader*> &indexes,
            const std::function<void(std::string_view header, uint64_t count,
                                     const std::vector<HeaderOccurrence> &listed)> &callback,
            size_t max_listed = std::numeric_limits<size_t>::max(),
            size_t num_threads = 0);

    /**
     * The records of |headers| whose header repeats an earlier record's, in record order:
     * |callback|(first, repeat) with the 0-based indices of the first occurrence and of the
     * repeat, once per repeat, until the callback returns false. One pass with a hash set over
     * the record indices (no header is copied), so a column of 33M headers is checked in O(n).
     * find_header() answers a header with its first record only, so a column whose headers
     * are not distinct hides every repeat from it: `annotate --index-header-coords` refuses
     * such a column and `stats --check-names` reports it.
     */
    static void find_duplicate_headers(const std::vector<std::string> &headers,
                                       const std::function<bool(size_t, size_t)> &callback);

    static constexpr auto kExtension = ".seqs";
    // the number that opens the trailer of column labels (the bytes "mtg-cols", written big-endian
    // as every number of the file)
    static constexpr uint64_t kColumnLabelsMarker = 0x6d74672d636f6c73;

  private:
    struct HeaderIndex;
    // the reverse index, built now if not built yet: shared, so that a reader keeps the array
    // it searches although set_header_hash() drops it meanwhile (load() is not thread-safe
    // with readers: it replaces the headers the entries name, find_header)
    std::shared_ptr<const HeaderIndex> header_index(size_t num_threads = 0) const;

    std::vector<std::vector<std::string>> headers_;
    std::vector<bit_vector_sd> coord_offsets_;
    // the trailer's labels while they are kept (until check_columns), and what it verified
    std::optional<std::vector<std::string>> column_labels_;
    bool column_order_verified_ = false;

    // find_header(): the entries name records of |headers_|, which only load() changes, and
    // load() drops the index first
    mutable std::mutex header_index_mutex_;
    mutable std::shared_ptr<const HeaderIndex> header_index_;
    mutable size_t header_index_builds_ = 0;
    HeaderHash header_hash_ = &default_header_hash;
};

} // namespace annot
} // namespace mtg

#endif // __COORD_TO_HEADER_HPP__
