#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <atomic>
#include <cstdlib>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
#include <variant>
#include <optional>

#include <json/json.h>
#include <sdsl/int_vector.hpp>

#include "common/vector.hpp"
#include "graph/alignment/alignment.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {

namespace seq_io {
    class FastaParser;
}

namespace graph {
    class AnnotatedDBG;
    namespace align {
        struct DBGAlignerConfig;
    }
}


namespace cli {

class Config;

using StringGenerator = std::function<void(std::function<void(const std::string &)>)>;


/**
 * Given a vector of kmer label matched coordinates, collapse continuous ranges of coordinates
 * to start-end tuples.
 *
 * @param coords    the vector of tuples as stored in a SeqSearchResult
 *                  originally returned from AnnotatedDBG::get_kmer_coordinates
 * @return vector of 'begin-end' range string representations
 */
std::vector<std::string>
collapse_coord_ranges(const std::vector<SmallVector<uint64_t>> &tuples);

/**
 * collapse_coord_ranges of the coordinates per query position |coords| stands for (a position
 * without coordinates ends every range open), computed from the positions with coordinates.
 */
std::vector<std::string>
collapse_coord_ranges(const graph::PositionCoords &coords);

// The text collapse_coord_ranges(coords) makes, without making it: |bytes|, its strings joined
// by ':' (as the JSON and text outputs write them); |ranges|, the number of strings;
// |long_ranges_bytes|, the bytes of the strings longer than a std::string holds inline (15),
// each with its terminator (a /search's answer account prices them)
struct CoordRangesText {
    uint64_t bytes = 0;
    uint64_t ranges = 0;
    uint64_t long_ranges_bytes = 0;
};
CoordRangesText collapsed_coords_text(const graph::PositionCoords &coords);
// The most the ranges collapse_coord_ranges and collapsed_coords_text keep open hold at once,
// for coordinates whose longest list (of a position) has |longest| coordinates: two vectors
// of (position, first, last), each grown by doubling to at most a list's coordinates (to the
// power of two at least |longest|), and at a growth the old buffer beside the new one, each
// buffer at its size class (alloc_bytes)
uint64_t coord_ranges_bytes(uint64_t longest);

/**
 * Encode a k-mer presence mask as alternating x<len>o<len>... runs (x = ones, o = zeros),
 * the trailing zeros included: 11100110 -> x3o2x2o1. The mask is given as a bit vector or as
 * its runs of ones (PresenceRuns), which give the same string.
 */
std::string encode_presence_mask(const sdsl::bit_vector &mask);
std::string encode_presence_mask(const graph::PresenceRuns &mask);

/**
 * The runs of one nonzero value of k-mer counts per query position, "<i>=<v>" for a run of
 * one position and "<i>-<j>=<v>" otherwise, in order (the runs of zeros left out).
 */
std::vector<std::string> count_run_strings(const graph::CountRuns &counts);


// Simple struct to wrap a query sequence
struct QuerySequence {
    size_t id;            // Sequence ID
    std::string name;     // Sequence name
    std::string sequence; // Sequence string representation
};

// Simple struct to wrap alignment results for a query sequence
struct Alignment {
    graph::align::DBGAlignerConfig::score_t score;
    graph::align::DBGAlignerConfig::score_t max_score;
    std::string cigar;
    bool orientation;
};


/**
 * Search result for an individual sequence. Stores the search result in one of several possible
 * vector types (kmer counts, label counts, quantiles, coordinates, signatures) represented
 * as a typesafe union std::variant.
 *
 * Note that this is an intermediate storage class and does not (1) store the actual sequence
 * nor (2) perform any of the search itself.
 */
class SeqSearchResult {
  public:
    typedef std::string Label;
    typedef std::vector<Label> LabelVec;
    typedef std::vector<std::pair<Label, size_t>> LabelCountVec;
    // the values per query position held where a label matches (PresenceRuns, CountRuns,
    // PositionCoords), written out as the vectors they stand for
    typedef std::vector<std::tuple<Label, size_t, graph::PresenceRuns>> LabelSigVec;
    typedef std::vector<std::tuple<Label, size_t, graph::CountRuns>> LabelCountAbundancesVec;
    // (label, num_kmer_matches, kmer_coordinates, num_kmers_in_target). The last
    // element is 0 when no CoordToHeader mapping is loaded (the label refers to
    // a whole annotation column rather than a single indexed sequence).
    typedef std::vector<std::tuple<Label, size_t, graph::PositionCoords, size_t>> LabelCountCoordsVec;

    typedef std::variant<LabelVec,
                         LabelCountVec,
                         LabelSigVec,
                         LabelCountAbundancesVec,
                         LabelCountCoordsVec> result_type;

    // JSON Field Keys
    static constexpr auto SEQ_DESCRIPTION_JSON_FIELD = "seq_description";
    static constexpr auto LABEL_SAMPLE_FIELD = "sample";
    static constexpr auto KMER_COUNT_FIELD = "kmer_count";
    static constexpr auto KMER_COORDINATE_FIELD = "kmer_coords";
    static constexpr auto KMERS_IN_TARGET_FIELD = "kmers_in_target";
    static constexpr auto SIGNATURE_FIELD = "signature";
    static constexpr auto KMER_ABUNDANCE_FIELD = "kmer_abundances";
    static constexpr auto SCORE_JSON_FIELD = "score";
    static constexpr auto MAX_SCORE_JSON_FIELD = "max_score";
    static constexpr auto SEQUENCE_JSON_FIELD = "sequence";
    static constexpr auto ALIGNMENT_JSON_FIELD = "alignments";
    static constexpr auto CIGAR_JSON_FIELD = "cigar";
    static constexpr auto ORIENTATION_JSON_FIELD = "orientation";

    /**
     * Construct an instance of SeqSearchResult by providing QuerySequence instance describing
     * sequence and result from that query on an annotated DBG. Will move both into this instance.
     *
     * @param sequence rvalue reference to QuerySequence
     * @param result rvalue reference to result_type
     * @param alignment alignment to graph (how the query was transformed)
     * @param cut the top-N cut the result was taken with
     */
    SeqSearchResult(QuerySequence&& sequence,
                    result_type&& result,
                    std::optional<Alignment>&& alignment = {},
                    const graph::SearchCut &cut = {})
          : sequence_(std::move(sequence)),
            alignment_(std::move(alignment)),
            result_(std::move(result)),
            cut_(cut) {}

    /** Const reference getters */
    const QuerySequence& get_sequence() const { return sequence_; }
    const std::optional<Alignment>& get_alignment() const { return alignment_; }
    std::optional<Alignment>& get_alignment() { return alignment_; }
    const result_type& get_result() const { return result_; }
    const graph::SearchCut& get_cut() const { return cut_; }

    /**
     * Returns a Json object representing the individual query result for the
     * represented sequence.
     *
     * @param verbose_output    do not collapse continuous ranges of coords (or counts)
     * @param k                 k-mer length for kmer presence mask scoring
     * @return  Json::Value instance representing sequence result
     */
    Json::Value to_json(bool verbose_output, size_t k) const;

    /**
     * Returns the cut of this sequence's search as a /search answer states it per pair
     * (SPEC-multigraph-envelope.md section 3), without the pair's tags: returned, complete,
     * nth_score, tied, tied_returned, and on an annotation with k-mer counts nth_median and
     * tied_same_median. A complete cut states nth_score (and nth_median) null and its counts
     * 0; a cut whose N-th score has one label (tied 1) states nth_median null, no median
     * having been read for it. A unit refused alone (refuse()) states no row (returned 0,
     * complete false, nth_score and nth_median null, the counts 0, the member set of its
     * computed cut) and `refused`: {account, code, max_memory_mb, and with the answer account
     * priced_mb}.
     */
    Json::Value cut_to_json() const;

    /**
     * Refuses this result's unit (its sequence on its pair) alone, with a /search's field
     * `partial_refusal` (lane T7-P1): its rows are dropped (their buffers freed; the variant
     * keeps its type), its alignment and sequence are kept, and its cut is marked with the
     * account that refused it, the bound and, for the answer account, the unit's own object's
     * price (both in MiB).
     */
    void refuse(graph::SearchCut::Refused account, uint64_t max_mb, uint64_t priced_mb = 0);

    /**
     * Returns a string representing the individual query result for the represented sequence.
     * Follows the format '<seq_name>\t<label>:<info>\t<label>:<info>\t ...'
     *
     * Done this way to allow for some output related config (as opposed to ostream override)
     *
     * @param delimiter             the delimiter between labels for that sequence
     * @param suppress_unlabeled    do not print seq_name if sequence is unlabeled
     * @param verbose_output        do not collapse continuous ranges of coords or counts
     * @param k                     k-mer length for kmer presence mask scoring
     */
    std::string to_string(std::string delimiter,
                          bool suppress_unlabeled,
                          bool verbose_output,
                          size_t k) const;

  private:
    QuerySequence sequence_;             // query sequence this result represents
    std::optional<Alignment> alignment_; // optional wrapper for alignment struct
    result_type result_;                 // result vector of labels and additional info
    graph::SearchCut cut_;               // the top-N cut the result was taken with
};


/**
 * The stop of a /search (lane T7-P2): its client gone or the server stopping (the server's
 * predicate, which answers the reason the work stops, or null while it goes on; with the result
 * cache, a shutdown only). Asked before each record's work starts (query_fasta's tasks), before
 * each batch is read, each of its alignments, its query graph and each of its records
 * (batched_query_fasta) and before each pair's search (the server's fan-out). Once the predicate
 * answered a reason it answers true without asking again: no further record, batch or pair
 * starts; those already running finish. Thread-safe. It keeps the first reason, and how far the
 * request got, for the log: the most records any one pair started, and the pairs started.
 */
class SearchStop {
  public:
    explicit SearchStop(std::function<const char *()> predicate)
          : predicate_(std::move(predicate)) {}
    SearchStop(const SearchStop &) = delete;
    SearchStop& operator=(const SearchStop &) = delete;

    // whether the work stops now: asks the predicate until it answered a reason once
    bool operator()() {
        if (stopped())
            return true;
        const char *reason = predicate_ ? predicate_() : nullptr;
        if (!reason)
            return false;
        // the first reason answered (two threads may ask at once)
        const char *none = nullptr;
        reason_.compare_exchange_strong(none, reason, std::memory_order_acq_rel);
        return true;
    }
    // whether it answered true (no call of the predicate)
    bool stopped() const { return reason_.load(std::memory_order_acquire) != nullptr; }
    // the reason the predicate answered first (null while not stopped)
    const char* reason() const { return reason_.load(std::memory_order_acquire); }

    // a pair's search started |records| of its records (the most of any pair is kept)
    void note_records_started(size_t records) {
        size_t most = records_started_.load(std::memory_order_relaxed);
        while (records > most && !records_started_.compare_exchange_weak(
                    most, records, std::memory_order_relaxed)) {}
    }
    size_t records_started() const { return records_started_.load(std::memory_order_relaxed); }
    // a pair's search began
    void note_pair_started() { pairs_started_.fetch_add(1, std::memory_order_relaxed); }
    size_t pairs_started() const { return pairs_started_.load(std::memory_order_relaxed); }

  private:
    std::function<const char *()> predicate_;
    std::atomic<const char *> reason_ { nullptr };
    std::atomic<size_t> records_started_ { 0 };
    std::atomic<size_t> pairs_started_ { 0 };
};

// A /search stopped by its SearchStop (lane T7-P2): thrown once the records already running
// have finished, whatever else they raised; the server writes nothing (its client is gone, or
// the server stops)
class SearchStopped : public std::runtime_error {
  public:
    SearchStopped() : std::runtime_error("the /search stopped: its client is gone") {}
};

// The records of the FASTA text |fasta| as query_fasta reads them (with |forward_and_reverse|
// each record and its reverse complement): a stopped /search's log line states them
size_t count_query_records(const std::string &fasta, bool forward_and_reverse);


enum QueryMode {
    LABELS = 0,
    MATCHES,
    COUNTS_SUM,
    COUNTS,
    COORDS,
    SIGNATURE
};

/**
 * Query sequences from a FASTA file.
 *
 * Provide a callback which receives each SeqSearchResult instance as an rvalue, which the
 * callback function should use as desired (typically generate a string/JSON output, or move
 * it into a container: a callback that keeps the result moves it, so that its per-position
 * values are never held twice; one taking a const reference reads it in place).
 *
 * @param file_path     path to FASTA file
 * @param callback      callback function
 * @param config        config object
 * @param anno_graph    annotated de Bruijn graph
 * @param aligner_config alignment config
 * @param budget        a /search's answer budget (null: none, the CLI's;
 *                      SPEC-multigraph-envelope.md section 3): each record's task makes the
 *                      values account of its unit (the record on this pair, named |pair| in
 *                      the log) on its own stack; the per-position modes charge their rows'
 *                      values to it as they grow, and each record's rows are charged to it as
 *                      its result is handed over. Past the bound the query throws
 *                      graph::SearchAnswerTooLarge (the request refused); with the request's
 *                      `partial_refusal` the unit is refused alone and its callback receives
 *                      the refused stub (SeqSearchResult::refuse). The callback consumes the
 *                      result before the task returns
 * @param pair          the pair, as the log names it ("pair (<name>, <graph>, <annotation>)"
 *                      or "the graph")
 * @param pair_index    the pair's place in the answer (the order of a record's `cuts`
 *                      entries; the order of the units refused alone in the log)
 * @param stop          a /search's stop (null: none, the CLI's; lane T7-P2): asked before each
 *                      record's work starts (and, on the batched path, before each batch, each
 *                      alignment and each query graph); once it answered true no further record
 *                      starts, the records running finish, the records this call started are
 *                      noted in it (on the batched path, those whose query on their batch's
 *                      query graph began), and SearchStopped is thrown, whatever else the
 *                      records raised
 *
 * @return the number of base pairs (characters) in query file
 */
size_t query_fasta(const std::string &file_path,
                   const std::function<void(SeqSearchResult &&)> &callback,
                   const Config &config,
                   const graph::AnnotatedDBG &anno_graph,
                   const graph::align::DBGAlignerConfig *aligner_config,
                   graph::SearchBudget *budget = nullptr,
                   const std::string &pair = "",
                   size_t pair_index = 0,
                   SearchStop *stop = nullptr);

int query_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __QUERY_GRAPH_HPP__
