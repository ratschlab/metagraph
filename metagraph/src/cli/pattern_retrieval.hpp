#ifndef __METAGRAPH_CLI_PATTERN_RETRIEVAL_HPP__
#define __METAGRAPH_CLI_PATTERN_RETRIEVAL_HPP__

/**
 * The labelled retrieval of POST /pattern (docs/DESIGN-pattern-search.md §4.3, §5.2-§5.5,
 * §7.2; increment 3): output.labels "all" for patterns of L <= k. The engine
 * (graph::pattern::PatternSearch::enumerate) releases the graph contexts; this module reads
 * their labels and places them, in two steps on the traversal's budget-aware, paced classes:
 *
 *  1. discovery: a LabelRecorder (column labels, at most max_labels_per_anchor per row, the
 *     row's true total beside the list) over the contexts' annotation keys;
 *  2. placement: a LabelQuery on the discovered labels, with coordinates, over the same keys
 *     (BASIC indexes with coordinates only): a context's first k-mer coordinate c is mapped
 *     to (seq_id, local) by the CoordToHeader FIRST, the offset added AFTER (start = local +
 *     offset + 1, 1-based); without a .seqs the occurrence is (kmer_coord, offset), placed
 *     nowhere.
 *
 * Both steps read through LabelOracle::decode_charged() annotations with a DecodeBudget (the
 * account's remainder) and ReadPacing (the request's deadline); an unbudgeted backend is
 * read only when the request allows it (allow_unbudgeted_annotation), and the answer says so.
 * They read one row at a time, the time and the work checked before each row, every row whose
 * read began charged its units (refused and interrupted reads too), so that the reads pass the
 * work budget by one row at most (review GPT-2, findings 1 and 2).
 * One memory account per request, over all its patterns, holds what the reads return, the
 * label dictionary, the contexts' descriptors (charged as the engine releases them, before
 * their result objects are built), the statements of refused and truncated rows (reserved
 * before the read that can produce them), the deduplication state and the labels built for
 * the answer (a deterministic model, as the walker's, never a measurement; an unbudgeted read
 * may take it past its maximum by the label names it returned, after which nothing more fits
 * and the reads stop, stated); the work
 * account counts the oracle's units (8 per row, 1 per entry and coordinate, and the rows'
 * row-diff dependencies; a refused row what its read decoded, at least 8). Every refusal,
 * truncation, cut and stop is stated (the owner's
 * guarantee rule); a count is exact only when everything behind it was read.
 */

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <json/json.h>

#include "graph/alignment/pattern_search.hpp"


namespace mtg {

namespace annot {
class CoordToHeader;
}

namespace graph {
class AnnotatedDBG;
namespace traversal {
class LabelOracle;
}
}

namespace cli {

/**
 * What the answer of one /pattern request has built so far, and the time writing it is
 * expected to take (review of 2026-10-07, X-EFFICIENCY-04). The finalisation reserve
 * (--pattern-finalize-ms) is the floor of the time kept back from the work; the work stops
 * earlier by this estimate, so that a request whose patterns buffered many results is still
 * answered within its budget (a time stop, its counts kept) rather than 503 "deadline" with
 * all of its work lost. As /traverse's delivery reserve (traverse_attempts.hpp,
 * Attempt::reserve_ms), but from configured rates only:
 *
 *   finalize_ms = margin x scale x ((T + P) / (build x 1000) + (T + P) / (compress x 1000)
 *                                   + P / (build x 1000))
 *
 * T the bytes of compact JSON text of the objects built (each result as the route builds it,
 * counted exactly by compact_json_bytes), P those of the label objects about to be built
 * (estimated from above before they are built; counted once more, at the build rate, for
 * building them), rates in MB/s (bytes per microsecond), scale the writer's text per compact
 * byte (1 on the server; 2 for the CLI's indented text). Used by one request's thread only.
 */
class AnswerVolume {
  public:
    AnswerVolume() = default;
    AnswerVolume(double build_mbps, double compress_mbps, double text_scale = 1)
          : build_mbps_(build_mbps), compress_mbps_(compress_mbps), scale_(text_scale) {}

    // |bytes| of compact JSON text built (to be written and compressed)
    void add(uint64_t bytes) { text_ += bytes; }
    // |bytes| of compact JSON text whose objects are about to be built
    void add_pending(uint64_t bytes) { pending_ += bytes; }
    void drop_pending(uint64_t bytes) { pending_ -= std::min(bytes, pending_); }
    // the pending objects were built: from then on written only
    void settle(uint64_t bytes) {
        drop_pending(bytes);
        add(bytes);
    }

    uint64_t text_bytes() const { return text_; }
    uint64_t pending_bytes() const { return pending_; }
    // the estimated time to write (and build the pending part of) what is buffered (ms); 0
    // for a rate of 0 or infinity (no model: the reserve alone, as before the review)
    double finalize_ms() const;

  private:
    double build_mbps_ = 0;
    double compress_mbps_ = 0;
    double scale_ = 1;
    uint64_t text_ = 0;
    uint64_t pending_ = 0;
};

// the margin of AnswerVolume::finalize_ms over its model (as /traverse's kReserveMargin)
constexpr double kAnswerVolumeMargin = 1.25;

// The length of |value| written as compact JSON by jsoncpp's StreamWriter, from above: a
// string's every control, quote, backslash or non-ASCII byte counted as an escape of 6 bytes,
// a double as 24 characters
uint64_t compact_json_bytes(const Json::Value &value);

/**
 * The annotation limits of one request, effective (after the server's caps, §5.3).
 */
struct RetrievalLimits {
    // the labels kept per row (LabelRecorder's cap); a row with more is truncated and stated
    uint64_t max_labels_per_anchor = 64;
    // the oracle's work units over all reads of the request (8 per row, 1 per entry and
    // coordinate, plus the rows' row-diff dependency units; checked before every row)
    uint64_t max_annotation_work = 100'000'000;
    // the request's memory account (bytes)
    uint64_t max_memory_bytes = uint64_t(256) << 20;
    // partial only: the labels listed per pattern, in label order (contexts desc, column asc)
    uint64_t max_labels = 1'000;
    // partial only: the placed occurrences listed per label, over its deduplicated union
    uint64_t max_occurrences_per_label = 16;
    // read an annotation without the budget-aware decode (no memory bound on the reads
    // themselves; the deadline checked between chunks of keys)
    bool allow_unbudgeted = false;
    // output.occurrences: place the labels (step 2) where the index can
    bool occurrences = true;
    // the expected duration of one uninterruptible piece of a read under the deadline (the
    // reads are decoded in chunks sized from the observed per-row time; 0: one piece)
    double chunk_target_ms = 50;
};

/**
 * What an index's annotation gives a labelled retrieval (§4.3), as the capabilities state it.
 *  placement  record (BASIC, coordinates, .seqs) | global (BASIC, coordinates, no .seqs) |
 *             none (no coordinates) | none_canonical (CANONICAL or PRIMARY)
 *  support    record_verified | label_intersection (the best per-label support, for L > k)
 *  budgeted   LabelOracle::decode_charged()
 */
struct AnnotationDescription {
    const char *placement = "none";
    const char *support = "label_intersection";
    bool budgeted = false;
};

AnnotationDescription describe_annotation(const graph::traversal::LabelOracle &oracle,
                                          graph::pattern::GraphMode mode);

/**
 * One released context of a pattern with L <= k, as the route collected it from the engine
 * (in answer order): its orientation and offset, its k-mer as the graph spells it, and the
 * annotation key of the k-mer (row + 1; the canonical k-mer's on a native CANONICAL graph,
 * the stored k-mer's on a wrapped PRIMARY one), npos when it has no row.
 */
struct RetrievalContext {
    graph::pattern::Orientation orientation = graph::pattern::Orientation::FORWARD;
    uint32_t offset = 0;
    std::string kmer;
    uint64_t key = 0;
};

// Tests: a record mapping instead of the index's; a hook called before every annotation read
// with its row count (to slow the reads down on a virtual clock); the memory account in bytes
// instead of the request's MiB (0: the request's); a hook asked at every charge of a read's
// DecodeBudget, refusing it when true (DecodeBudget::deny); a hook called before the work time
// is read for the labels of each context built for the answer, with the context's index (to
// move a virtual clock past the work time in the middle of the output); and a hook called
// once the route's work is done, before the answer is assembled (to move it into the
// finalisation reserve or past the deadline)
struct RetrievalHooks {
    const annot::CoordToHeader *coord_to_header = nullptr;
    std::function<void(size_t rows)> read_hook;
    uint64_t max_memory_bytes = 0;
    std::function<bool(uint64_t ordinal)> deny_decode;
    std::function<void(size_t context)> output_hook;
    std::function<void()> work_done_hook;
};

/**
 * The labels of one pattern, read and built during the work phase, merged into the
 * pattern's entry during the finalisation (apply()).
 */
struct LabelsAnswer {
    // the entry's fields written as they are: placement, annotation, by_label, rows_refused,
    // anchors_truncated, labels_cut, occurrences_cut
    Json::Value fields = Json::Value(Json::objectValue);
    Json::Value labels_count;
    Json::Value occurrences_count;
    // added to `work` and `timing`
    Json::Value work = Json::Value(Json::objectValue);
    Json::Value timing = Json::Value(Json::objectValue);
    // per admitted context, in order: support, labels_status, labels_total, labels (merged
    // into the result objects); empty when nothing was read
    std::vector<Json::Value> result_fields;
    // all_or_count: the results are withheld for this reason (rows refused, a truncated
    // anchor, a budget or deadline stop of the reads)
    std::optional<std::string> withheld;
    // partial: the contexts list was cut for this reason (max_memory)
    std::optional<std::string> cut;
    // a stop of the annotation work: {phase, reason}
    std::optional<std::pair<std::string, std::string>> stop;
    bool time_limited = false;
    // every released context came with all its labels and placements, nothing cut
    bool complete = false;
    // added after the engine's notes
    std::vector<std::string> notes;
};

/**
 * One request's labelled retrieval: the oracle, the label dictionary (LabelRecorder), the
 * memory and work accounts, shared by the request's patterns in request order. Used by one
 * thread. The deadline is the request's (|budget|): a time stop of the reads is recorded in
 * |budget| too, so that the patterns after it answer as after any time stop. |volume| (may be
 * null) is the request's answer volume: the labels built for the answer are added to it, and
 * the work time is read before each context's labels are built (with |budget|'s deadline,
 * which reads the volume: a time stop there is stop {output, time}).
 */
class PatternRetrieval {
  public:
    PatternRetrieval(const graph::AnnotatedDBG &anno_graph, graph::pattern::GraphMode mode,
                     const RetrievalLimits &limits, graph::pattern::Budget &budget,
                     const RetrievalHooks *hooks = nullptr, AnswerVolume *volume = nullptr);
    ~PatternRetrieval();

    const AnnotationDescription& description() const { return description_; }
    // the placement this request's answers give: the index's, or "not_requested" when the
    // request set output.occurrences false
    const char* placement() const;
    const RetrievalLimits& limits() const { return limits_; }

    /**
     * Before the engine releases one pattern's contexts (mode ALL_OR_COUNT or PARTIAL): opens
     * the pattern's allowance for their descriptors, what the account has left in
     * all_or_count, half of it in partial (the other half kept for the pattern's reads and
     * labels, so that a memory cut of the list still returns labelled contexts).
     */
    void begin_release(graph::pattern::Mode mode);
    /**
     * One released context, before the route builds its result object: charges its
     * descriptor (512 + 2k bytes of the model). False when it does not fit the allowance, and
     * for every context of the pattern after the first that did not: the route builds none of
     * them (partial: the list is cut, max_memory; all_or_count: withheld, output_budget).
     */
    bool admit_context();

    /**
     * Reads the labels of the admitted contexts of one pattern (L <= k, or a long pattern's
     * empty release) and builds their JSON. |released| is the number of contexts the engine
     * released (|contexts| are the first of them, admitted); |extraction| is the engine's
     * (withheld: nothing is read); |length| is L; |graph_name| names the graph in by_label
     * (null when unnamed). Mode ALL_OR_COUNT or PARTIAL.
     */
    LabelsAnswer retrieve(const std::vector<RetrievalContext> &contexts, uint64_t released,
                          size_t length, graph::pattern::Mode mode,
                          const graph::pattern::Extraction &extraction,
                          const Json::Value &graph_name);

    // the memory account's peak so far (bytes of the model) and what it holds now
    uint64_t memory_peak() const;
    uint64_t memory_held() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    RetrievalLimits limits_;
    AnnotationDescription description_;
};

/**
 * Merges |answer| into the pattern's |entry| (written by the route for the engine's result,
 * with every released context in `results`): the label fields, the counts, the results'
 * labels, and — when the reads did not complete — withheld (all_or_count), the cut
 * (partial), retrieval_complete, stop and determinism.
 */
void apply_labels(Json::Value *entry, LabelsAnswer &&answer, graph::pattern::Mode mode);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_RETRIEVAL_HPP__
