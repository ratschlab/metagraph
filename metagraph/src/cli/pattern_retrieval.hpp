#ifndef __METAGRAPH_CLI_PATTERN_RETRIEVAL_HPP__
#define __METAGRAPH_CLI_PATTERN_RETRIEVAL_HPP__

/**
 * The labelled retrieval of POST /pattern (docs/DESIGN-pattern-search.md §4.3, §5.2-§5.5,
 * §7.2): output.labels "all" for patterns of L <= k, and (long_search "paths") for the paths of
 * patterns longer than k (retrieve_paths: the labels on every k-mer of a path, each with its
 * support). The engine (graph::pattern::PatternSearch::enumerate) releases the graph contexts;
 * this module reads their labels and places them, in two steps on the traversal's budget-aware,
 * paced classes:
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
 * account's remainder) and ReadPacing (the request's deadline); an unbudgeted backend is read
 * only when the request allows it (allow_unbudgeted_annotation), and the answer says so. They
 * read one row at a time, the time and the work checked before each row, every row whose read
 * began charged its units (refused and interrupted reads too), so that the reads pass the work
 * budget by one row at most.
 * One memory account per request, over all its patterns, holds what the reads return, the label
 * dictionary, the contexts' descriptors (charged as the engine releases them, before their
 * result objects are built), the statements of refused and truncated rows (reserved before the
 * read that can produce them), the deduplication state, the runs of occurrences a path's
 * verification keeps for its output and the labels built for the answer (a deterministic model,
 * as the walker's, never a measurement; an unbudgeted read may take it past its maximum by the
 * label names it returned, after which nothing more fits and the reads stop, stated); the work
 * account counts the oracle's units (8 per row, 1 per entry and coordinate, and the rows'
 * row-diff dependencies; a refused row what its read decoded, at least 8). The work between the
 * reads (the occurrences, the paths' label lists and their verification) reads the clock at
 * least every Budget::kClockStride units, before the work. Each placed occurrence is counted in
 * its label's union, but a label's list holds, and the answer's volume and account are charged
 * for, only what it can list (partial: max_occurrences_per_label of the context's or path's
 * first), cut to the union's first once it is complete. Every refusal, truncation, cut and stop
 * is stated; a count is exact only when everything behind it was read.
 *
 * Predicates (pattern_selection.cpp; the internals shared through pattern_retrieval_impl.hpp):
 * the selection of a request's predicate for patterns of L <= k — bind() once per request, then
 * per pattern begin_selection(), admit_tested() and tested_context() as the engine releases the
 * raw contexts, select() (the pass: rows and reverse-complement lookups, one restricted read
 * per row under max_predicate_work, decisions in answer order, the relations of SPEC §19.7, the
 * selection admission) and the projection of the chosen contexts: none, predicate_only
 * (retrieve_given: the pass's rows instead of a discovery, placed) or all (retrieve()).
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
#include "graph/traversal/traversal_types.hpp"


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

namespace predicate {
struct Predicate;
class Bound;
struct Binding;
}

/**
 * What the answer of one /pattern request has built so far, and the time writing it is expected
 * to take. The finalisation reserve (--pattern-finalize-ms) is the floor of the time kept back
 * from the work; the work stops earlier by this estimate, so that a request whose patterns
 * buffered many results is still answered within its budget (a time stop, its counts kept)
 * rather than 503 "deadline" with all of its work lost. As /traverse's delivery reserve
 * (traverse_attempts.hpp, Attempt::reserve_ms), but from configured rates only:
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
    // the estimated time to write (and build the pending part of) what is buffered (ms); 0 for
    // a rate of 0 or infinity (no model: the reserve alone)
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

/**
 * One released path of a pattern longer than k (long_search "paths"; §4.2), as the route
 * collected it from the engine (in answer order): its orientation, the L bases it spells, and
 * the annotation key of each of its n = L - k + 1 k-mers in reading order (row + 1, as
 * RetrievalContext::key; npos when one has no row).
 */
struct RetrievalPath {
    graph::pattern::Orientation orientation = graph::pattern::Orientation::FORWARD;
    std::string sequence;
    std::vector<uint64_t> keys;
};

/**
 * The memory model's price of one released path (§5.3, item 3: the retained paths are charged
 * before they are built): its result object — sequence, instance and anchor_kmer, its nodes
 * and rows arrays (n = L - k + 1 entries each, an array element of the JSON library and its
 * text) — its descriptor, and the retrieval's copy of its sequence and keys. A deterministic
 * model, as context_bytes: 512 + 2k + 3L + 192n.
 */
uint64_t path_descriptor_bytes(size_t k, size_t length);

// Tests: a record mapping instead of the index's; a hook called before every annotation read
// with its row count (to slow the reads down on a virtual clock); the memory account in bytes
// instead of the request's MiB (0: the request's); a hook asked at every charge of a read's
// DecodeBudget, refusing it when true (DecodeBudget::deny); a hook called before the work time
// is read for the labels of each context built for the answer, with the context's index (to
// move a virtual clock past the work time in the middle of the output); a hook called once
// the route's work is done, before the answer is assembled (to move it into the finalisation
// reserve or past the deadline); and a hook called before the occurrences of each context read
// are made, and before each path's labels are verified (once the clock was read there), with
// its index (to move the clock past the work time inside that work, which reads it itself)
struct RetrievalHooks {
    const annot::CoordToHeader *coord_to_header = nullptr;
    std::function<void(size_t rows)> read_hook;
    uint64_t max_memory_bytes = 0;
    std::function<bool(uint64_t ordinal)> deny_decode;
    std::function<void(size_t context)> output_hook;
    std::function<void()> work_done_hook;
    std::function<void(size_t item)> occurrences_hook;
};

/**
 * Counters of the labelled retrieval, additive: one pattern's in its LabelsAnswer, the
 * request's so far in PatternRetrieval::counters(). Not in the answer unless the route states
 * them.
 */
struct RetrievalCounters {
    // the distinct annotation rows whose labels were read (LabelRecorder, complete or
    // truncated): each row once per pattern, however many of its contexts or of its paths'
    // k-mers share it; the placement's second read of a row is not counted again (the
    // request's sum counts a row once per pattern that read it)
    uint64_t rows_distinct = 0;
    // paths: the time spent intersecting the label lists of each path's rows (ms)
    double label_intersection_ms = 0;
    // paths: the time spent verifying the labels carrying the paths, the chains of every
    // (path, label) in its rows' coordinates with their record placement (ms)
    double verification_ms = 0;
    // paths: the verification's units of work, one per k-mer row looked up for a label
    // carrying a path, per list ordered, per galloping seek, per run of chains extended in a
    // list and per record a run crosses — a homopolymer's run of chains is a few units per
    // k-mer of the path, not one per coordinate
    uint64_t verification_steps = 0;

    RetrievalCounters& operator+=(const RetrievalCounters &other) {
        rows_distinct += other.rows_distinct;
        label_intersection_ms += other.label_intersection_ms;
        verification_ms += other.verification_ms;
        verification_steps += other.verification_steps;
        return *this;
    }
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
    // this pattern's (not merged by apply_labels: the route states them)
    RetrievalCounters counters;
};

/**
 * The selection of a predicate (SPEC-pattern-search.md §19): which contexts of a pattern of
 * L <= k satisfy the request's predicate, read from their annotation rows by the selection pass
 * (PatternRetrieval::select, pattern_selection.cpp).
 */

// The request-wide options of the selection, effective (after the server's caps, §19.2)
struct SelectionLimits {
    // the work units of the selection's annotation reads, lookups and decisions over the
    // request (§19.9), a budget of its own beside max_annotation_work
    uint64_t max_predicate_work = 100'000'000;
    // predicate_strands "either" (true): on a BASIC graph a label is present for a context
    // when it annotates the context's k-mer or its reverse complement; "context" (false): the
    // context's own k-mer only. On CANONICAL and PRIMARY graphs one row serves both
    // orientations: "either" whatever is asked
    bool either = true;
};

// output.labels of a predicate request (§19.2): what the selected contexts are returned with
enum class Projection { NONE, PREDICATE_ONLY, ALL };

// selection.pass (§19.7)
enum class SelectionPass { COMPLETED, STOPPED, NOT_ADMITTED, NOT_STARTED, CONSTANT };
const char* to_string(SelectionPass pass);

/**
 * One raw context of a pattern released into the selection pass, in answer order: what the
 * engine released (orientation, offset, node, base_node) and the annotation key of its k-mer
 * (row + 1, as RetrievalContext::key; PatternRetrieval::tested_context computes it as the
 * route names a context's row), then what the pass made of it. 64 bytes in the memory model
 * (admit_tested), whatever its layout.
 */
struct TestedContext {
    graph::pattern::Orientation orientation = graph::pattern::Orientation::FORWARD;
    uint32_t offset = 0;
    uint64_t node = 0;
    uint64_t base_node = 0;
    uint64_t key = 0;
    // the pass's: its row in the pass's table, whether every row it needs was read and it
    // was decided (tested), and the decision
    uint32_t row = UINT32_MAX;
    bool decided = false;
    bool selected = false;
};

// What one pattern's selection is asked for, beside the request-wide SelectionLimits
struct SelectionRequest {
    // the request's mode: COUNT (the counts only), ALL_OR_COUNT or PARTIAL
    graph::pattern::Mode mode = graph::pattern::Mode::ALL_OR_COUNT;
    // the selection admission's threshold (all_or_count), the list's cap (partial), and
    // stop_at_threshold's threshold on the selected count
    uint64_t max_contexts = 10'000;
    bool stop_at_threshold = false;
    // what the selected contexts are returned with: NONE keeps no row and no
    // selection_labels; PREDICATE_ONLY keeps the chosen contexts' own rows for
    // retrieve_given and builds their selection_labels; ALL their selection_labels (the
    // projection reads the rows again, retrieve())
    Projection projection = Projection::NONE;
};

/**
 * What the selection pass made of one pattern (§19.6 steps 3 and 4, §19.7, §19.8). The
 * route writes counts.tested and counts.selected, selection.pass, work.predicate_rows and
 * work.predicate_units, timing.selection_ms, and composes withheld, cut and stop with the
 * engine's.
 */
struct SelectionAnswer {
    SelectionPass pass = SelectionPass::NOT_STARTED;
    // unit graph_contexts: tested exact (the decisions made; a constant: the raw count),
    // selected with the relation of §19.7
    graph::pattern::Count tested
            = graph::pattern::Count::unknown(graph::pattern::Unit::GRAPH_CONTEXTS);
    graph::pattern::Count selected
            = graph::pattern::Count::unknown(graph::pattern::Unit::GRAPH_CONTEXTS);
    // the selected contexts the answer lists, indices into the pass's tested contexts in
    // answer order: all_or_count all of them or none, partial the first max_contexts, count
    // none
    std::vector<uint32_t> chosen;
    // with a projection that reads labels: per chosen context its selection_labels, the
    // predicate's labels in the set it was evaluated on (its row's, and with "either" its
    // reverse complement's), as ids into Bound::labels(), in label order (contexts desc over
    // the chosen, column asc); charged when the context was decided
    std::vector<std::vector<graph::traversal::LabelId>> selection_labels;
    // beside each of them, which row it was found on (per selected result and label the
    // orientation that supported it): kOnContext (the context's k-mer x as spelled),
    // kOnReverseComplement (rc(x), "either" on a BASIC graph), or both (a palindromic x under
    // "either": x is rc(x)). selection_strands_json names them
    static constexpr uint8_t kOnContext = 1;
    static constexpr uint8_t kOnReverseComplement = 2;
    std::vector<std::vector<uint8_t>> selection_label_rows;
    // all_or_count: why nothing is listed (predicate_budget, deadline, threshold_crossed,
    // selected_above_threshold, output_budget)
    std::optional<std::string> withheld;
    // partial: why the list may be shorter than the selected contexts, by the pass itself
    // (max_predicate_work, max_memory, time, max_contexts for a stop_at_threshold stop, or
    // max_memory for refused rows and descriptors the account could not hold)
    std::optional<std::string> cut;
    // partial: more contexts were selected than the list holds (cut max_contexts, after the
    // engine's and the pass's cuts, §19.8)
    bool list_cut = false;
    // the pass's stop, {selection, max_predicate_work | max_memory | time | max_contexts},
    // or {output, max_memory | time} while its selection_labels were built
    std::optional<std::pair<std::string, std::string>> stop;
    bool time_limited = false;
    // the statements of the rows the account could not hold (phase "selection")
    Json::Value rows_refused = Json::Value(Json::arrayValue);
    // work.predicate_rows: the rows the pass read; work.predicate_units: its units (refused
    // and interrupted reads included); the reverse-complement lookups made
    uint64_t rows = 0;
    uint64_t units = 0;
    uint64_t lookups = 0;
    double ms = 0;
};

// The selection of a pattern without a pass (§19.7): its normal form a constant (pass
// "constant": false selects nothing, tested the raw count; true selects the raw count, with
// its relation), or a pass that did not run (not_admitted, not_started: selected unknown,
// exact 0 when the raw count is exact 0)
SelectionAnswer constant_selection(bool value, const graph::pattern::Count &raw);
SelectionAnswer selection_without_pass(SelectionPass pass, const graph::pattern::Count &raw);

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
    // the same for one released path of a pattern of |length| bases (path_descriptor_bytes)
    bool admit_path(size_t length);

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

    /**
     * The labels of the admitted paths of one pattern longer than k (long_search "paths";
     * DESIGN §4.3 "Label consistency for long"), as retrieve() for contexts, with the same
     * budgets, statements and modes:
     *  1. discovery: the rows of every k-mer of every admitted path (each distinct row read
     *     once, one row per read, in answer order of first appearance); a path's labels are
     *     the labels present on EVERY one of its k-mers (support label_intersection);
     *  2. verification (placement record, occurrences requested): the coordinates of the rows
     *     of the paths with at least one such label; a label is record_verified on a path when
     *     its coordinates show one contiguous occurrence of the whole path in ONE record: a
     *     column coordinate c of the first k-mer with c + i a coordinate of the i-th k-mer for
     *     every i, c mapped to (seq_id, local) first, and local + n - 1 inside the record's
     *     k-mers (a chain crossing into the next record of the column is not one). Each such
     *     occurrence is placed: (seq_id, 1-based local + 1, the path's strand), nt_coords over
     *     the L bases. Placement global: the chains (kmer_coord, offset 0), record bounds
     *     unknown, nothing verified. Elsewhere label_intersection only. The chains are the
     *     intersection of the k-mers' coordinate lists shifted by -i (a leapfrog join from the
     *     smallest, consecutive chains kept as runs), made once per (path, label) and kept for
     *     the output as runs of occurrences.
     * |require_verified| (require_support "record_verified"; placement record, checked by the
     * route): only the verified labels are listed, the others counted per path and per entry
     * (labels_excluded_unverified). Mode ALL_OR_COUNT or PARTIAL.
     */
    LabelsAnswer retrieve_paths(const std::vector<RetrievalPath> &paths, uint64_t released,
                                size_t length, graph::pattern::Mode mode,
                                const graph::pattern::Extraction &extraction,
                                const Json::Value &graph_name, bool require_verified);

    // the memory account's peak so far (bytes of the model)
    uint64_t memory_peak() const;
    // the request's counters so far (the sum of its patterns')
    const RetrievalCounters& counters() const;

    // ---- the selection of a predicate (pattern_selection.cpp; SPEC §19)

    /**
     * Binds the request's |predicate| to this index's columns (predicate::Bound::bind, the
     * work time read every 4,096 names), its bytes (Bound::bytes) charged to the memory
     * account for the whole request, once, before the first pattern; |limits| are the
     * request's. The binding's stop ("time", "max_memory") stops every pattern's pass after
     * it (not_started, stop {selection, <stop>}). Called once per request, before any select().
     */
    const predicate::Binding& bind(const predicate::Predicate &predicate,
                                   const SelectionLimits &limits);
    // the bound predicate (null before bind() and when bind() stopped)
    const predicate::Bound* bound() const;
    // predicate.strands as evaluated: "either" on CANONICAL and PRIMARY graphs whatever was
    // asked, else as asked ("either" | "context")
    const char* selection_strands() const;
    // selection.access: "rows" (budget-aware, or unbudgeted rows) or "columns" (an
    // unbudgeted annotation with direct access and at most 16 known labels)
    const char* selection_access() const;
    // the request's selection units so far (max_predicate_work)
    uint64_t predicate_units() const;

    /**
     * Before the engine releases one pattern's raw contexts into the pass (the route runs
     * enumerate() with max_contexts = max_predicate_contexts, mode ALL_OR_COUNT for a count
     * request): frees what the previous pattern's pass still held (end_selection) and opens
     * the allowance of the descriptors, as begin_release: what the account has left
     * (all_or_count, count), half of it (partial).
     */
    void begin_selection(graph::pattern::Mode mode);
    // One released raw context, before the route stores it: charges its descriptor (64 bytes
    // of the model, §19.9). False when it does not fit the allowance, and for every one after
    // the first that did not: the route stores none of them (the pass's set ends there)
    bool admit_tested();
    // the descriptor of the released context |c|: its key named as the route names a
    // context's row (§7.10: the stored k-mer's on BASIC and wrapped PRIMARY graphs, the
    // canonical k-mer's on a native CANONICAL graph, spelled; npos when it has no row)
    TestedContext tested_context(const graph::pattern::Context &c) const;

    /**
     * The selection pass of one pattern (§19.6 steps 3 and 4). |tested| are the admitted raw
     * contexts in answer order (admit_tested), the first of |released| the engine released
     * (more when the account could not hold their descriptors: the pass's set ended there, stop
     * {selection, max_memory}); |raw| is the raw count (counts.contexts.total) and |x| the
     * engine's extraction (withheld: no pass, not_admitted for count_above_threshold, else
     * not_started). In order:
     *  1. the rows: each context's row and, with "either" on a BASIC graph, its reverse
     *     complement's — the k-mer spelled, reverse-complemented and looked up
     *     (LabelOracle::keys_of_sequence), k units charged and the gate checked before each;
     *     one lookup per distinct row, its result kept (and the mirror's mirror known: rc is an
     *     involution); every distinct row charged kSelectionRowBytes before it is added;
     *  2. the reads: one row per read, in the order of first appearance (a context's row before
     *     its mirror's), the time and max_predicate_work checked before each, a statement
     *     reserved, a LabelQuery over Bound::labels() (Access::ROWS for the budget-aware fetch,
     *     never AUTO; unbudgeted: rows, or single cells for at most 16 labels when the
     *     annotation has direct access) with the account's remainder; 8 + KeyCost::entries +
     *     dependency units charged on every outcome (a refused read what its decode reached,
     *     at least 8; unbudgeted: 8 + the hits); a refused row stated (rows_refused, phase
     *     "selection"; all_or_count ends the pass there, the other modes go on);
     *  3. the decisions, in answer order as their rows are read: Bound::eval on the labels of
     *     the context's row and its mirror's, its units charged (never refused), the clock
     *     every 64 lookups and decisions; stop_at_threshold ends the pass once more than
     *     max_contexts are selected; after a stop (but time) every context whose rows were
     *     read is decided all the same. A selected context among the first max_contexts gets
     *     its selection_labels (with a projection that reads labels) and, for predicate_only,
     *     its row kept, both charged when it is decided;
     *  4. the counts and their relations (§19.7), the selection admission (all_or_count: all
     *     selected or none; partial: the first max_contexts), the label order of the
     *     selection_labels; every row's hits freed but the chosen contexts' own rows
     *     (predicate_only, for retrieve_given).
     * Throws std::logic_error before bind(), for a constant normal form (no pass: see
     * constant_selection), and std::runtime_error for a context without a row.
     */
    SelectionAnswer select(std::vector<TestedContext> &tested, uint64_t released,
                           const graph::pattern::Count &raw,
                           const graph::pattern::Extraction &x,
                           const SelectionRequest &request);
    // the j-th chosen context's selection_labels (the names; their bytes were charged by the
    // pass)
    Json::Value selection_labels_json(const SelectionAnswer &answer, size_t j) const;
    // beside them, per label the orientation whose row carries it (§19.10): "context" (the
    // context's k-mer as spelled), "reverse_complement" (its reverse complement's row only,
    // "either" on a BASIC graph), "both", or "either" on CANONICAL and PRIMARY graphs (one row
    // serves both orientations: no strand is known)
    Json::Value selection_strands_json(const SelectionAnswer &answer, size_t j) const;
    // frees what the last pattern's pass still holds: its descriptors (after the route built
    // the chosen contexts' results) and its kept rows (when retrieve_given did not take them)
    void end_selection();

    /**
     * The projection predicate_only (§19.10) of the chosen contexts of the last select()
     * (|contexts|, admitted by admit_context after begin_release as for retrieve()): retrieve()
     * with the discovery replaced by the pass's rows — each context's labels are the
     * predicate's labels on its own row, as the pass read them (never truncated; labels_total
     * their number) — and the placement of those labels (LabelQuery with coordinates over the
     * labels found on the rows, under max_annotation_work), with retrieve()'s statements,
     * cuts, modes and answer. Takes the kept rows (end_selection is then needed for the
     * descriptors only).
     */
    LabelsAnswer retrieve_given(const std::vector<RetrievalContext> &contexts, uint64_t released,
                                size_t length, graph::pattern::Mode mode,
                                const graph::pattern::Extraction &extraction,
                                const Json::Value &graph_name);

  private:
    struct Impl;
    struct Given;
    std::unique_ptr<Impl> impl_;
    RetrievalLimits limits_;
    AnnotationDescription description_;

    bool admit(uint64_t bytes);
    // retrieve() and retrieve_given(): |given| null reads the rows' labels (discovery)
    LabelsAnswer retrieve_rows(const std::vector<RetrievalContext> &contexts, uint64_t released,
                               size_t length, graph::pattern::Mode mode,
                               const graph::pattern::Extraction &extraction,
                               const Json::Value &graph_name, Given *given);
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
