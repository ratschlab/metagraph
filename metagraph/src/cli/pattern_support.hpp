#ifndef __METAGRAPH_CLI_PATTERN_SUPPORT_HPP__
#define __METAGRAPH_CLI_PATTERN_SUPPORT_HPP__

/**
 * The supported-path search of POST /pattern (long_search "supported_paths", SPEC
 * §20): the label and trace trackers the engine's extension asks per level
 * (graph::pattern::SupportTracker) and the list of supported paths it hands its complete walks
 * to (graph::pattern::PathSink), over the request's labelled retrieval (its oracle, memory
 * account, annotation work and deadline: PathSupportEnv).
 *
 * A supported path is a walk spelling the pattern that some label supports along its whole
 * length, in ONE orientation as a whole: the label annotates every k-mer of the walk as the walk
 * spells it (label_intersection), and at the record level one of its records holds the walk
 * whole (record_verified: a column coordinate c of the first k-mer, c + i of the i-th, all in
 * one record). Each oriented pattern is its own walk (FORWARD spells P, REVERSE rc(P)), so a
 * label holding some k-mers of a walk on one strand and the rest on the other supports neither
 * (on CANONICAL and PRIMARY graphs a k-mer and its reverse complement share one row: no strand
 * is known there, and only the label level exists).
 *
 * Per pattern:
 *  1. before the first anchor is extended (SupportTracker::prepare), the anchors' rows are read
 *     WHOLE (a LabelRecorder over every column, never truncated); the union of their labels is
 *     the pattern's permitted set — only an anchor's labels can support a walk from it — and
 *     every later row is read for those labels only (a LabelQuery, with coordinates at the
 *     record level);
 *  2. each anchor opens frame 0 from its row, each step of the DFS narrows the top frame by the
 *     next k-mer's row into the next frame (graph::traversal::support_step: label merges, and at
 *     the record level the chains of consecutive coordinates as runs, clipped at their record's
 *     end); a frame whose support is empty is DEAD and the branch pruned;
 *  3. a complete walk is a supported path: the sink keeps it (its retention rule, its memory
 *     admitted first) with its labels from the frames — every label carrying the walk, each
 *     record_verified when a chain survived, with the chains' occurrences —, nothing read again.
 *
 * Rows are read once per pattern while the row cache holds them: a fixed allotment of the
 * account (a quarter of what is left when the pattern starts, at most 64 MiB) holds the rows as
 * the steps use them (labels and coordinate runs), evicted wholesale when the next one does not
 * fit, and the row-diff path cache in what the rows leave of it; a row read again after an
 * eviction is charged again.
 *
 * Use, per request: a PathTracker and a SupportedPathSink over the request's PathSupportEnv
 * (both destroyed before what the env refers to); per pattern longer than k:
 *   tracker.begin_pattern(L); sink.begin_pattern(mode, max_paths);
 *   request.extend_paths = true; request.support = &tracker; request.sink = &sink;
 *   result = search.count(pattern, request, budget);      // enumerate() refuses a sink
 *   x = supported_release(result, request, sink, &named);  // nullopt in mode count
 *   results: sink.paths()[0, x->returned), as the paths of long_search "paths" are written;
 *   labels (PathSupportOptions::labels): sink.labels_answer(*x, named, L, graph name);
 *   sink.end_pattern(the results listed); tracker.end_pattern();
 * A stop {extension, external} of the engine is the tracker's (stop_reason(): time,
 * max_annotation_work, max_memory) or the sink's (time).
 *
 * Budgets: before every row read the work time is read (Budget::check_time, so that the engine
 * states a time stop as {extension, time}), the request's annotation units are compared with
 * max_annotation_work and room for the row's statement is checked; a read charges 8 units, the
 * decoded row's entries and its row-diff dependency units on every outcome, and turning its
 * hits into runs 1 per hit and coordinate; a step charges the labels and runs its merges
 * compare (support_step's units), its frame's bytes admitted before it is built. A row the
 * account cannot hold stops the search (max_memory, stated in rows_refused); it is never taken
 * as empty. Long merges read the clock every 4,096 rounds.
 */

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <json/json.h>

#include "pattern_retrieval.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/traversal/support_step.hpp"
#include "graph/traversal/traversal_types.hpp"


namespace mtg {

namespace graph {
namespace traversal {
class LabelOracle;
class LabelQuery;
class LabelRecorder;
}
}

namespace cli {

namespace retrieval {
class Account;
}

/**
 * What the supported-path search reads with: the request's labelled retrieval's oracle, memory
 * account, annotation units and deadline (PatternRetrieval's, so that one request holds one
 * account and one work budget over all its patterns, whatever each is answered by), its limits
 * and how its annotation is read.
 */
struct PathSupportEnv {
    const graph::traversal::LabelOracle &oracle;
    graph::pattern::Budget &budget;
    retrieval::Account &account;
    // the request's annotation units so far, compared with limits.max_annotation_work
    uint64_t &units;
    const RetrievalLimits &limits;
    graph::pattern::GraphMode mode = graph::pattern::GraphMode::BASIC;
    // LabelOracle::decode_charged(): the budget-aware reads; else the unbudgeted ones (the route
    // requires allow_unbudgeted_annotation for them), their results charged after the fact
    bool budgeted = false;
    // the request's answer volume (null: none): the text of the labels built for the answer
    AnswerVolume *volume = nullptr;
    // tests: DecodeBudget::deny of every budget-aware read
    std::function<bool(uint64_t)> deny_decode;
};

// What a request asks of the supported-path search
struct PathSupportOptions {
    /**
     * The level the search prunes on: TRACE (record_verified: one record of the label holds
     * the walk whole) on a BASIC index with coordinates and the record mapping, else KMER
     * (label_intersection: the label annotates every k-mer of the walk).
     */
    graph::traversal::Support level = graph::traversal::Support::KMER;
    // output.labels "all": every kept path with its labels (and their occurrences where the
    // index places them and output.occurrences asks for them); false: the paths only
    bool labels = false;
    // require_support "record_verified" (TRACE only): a path lists its verified labels only,
    // the others counted (labels_excluded_unverified)
    bool require_verified = false;
    // the labels a kept path lists, by column (output.labels "predicate_only": the predicate's
    // columns); a label it does not accept is neither listed nor counted (labels_total). Empty:
    // every label carrying the path
    std::function<bool(graph::traversal::Column)> listed;
};

/**
 * What one pattern's supported-path search did (work beside the engine's steps; exact as
 * work, whether or not the search completed).
 */
struct PathSupportWork {
    // annotation rows read: the anchors' whole rows, then the walk's rows (a row read again
    // after the cache evicted it counted again)
    uint64_t rows = 0;
    // of them, the anchors' whole rows (the permitted set)
    uint64_t anchor_rows = 0;
    // the steps (frames opened or narrowed) whose row came from the cache
    uint64_t cache_hits = 0;
    // the times the row cache was emptied to make room (the rows after it may be read again),
    // and the rows larger than its allotment (read for the one step that needed them)
    uint64_t evictions = 0;
    uint64_t uncached_rows = 0;
    // the labels of the permitted set (the anchors' rows' union)
    uint64_t permitted_labels = 0;
    // this pattern's annotation units (reads, conversions, merges)
    uint64_t units = 0;
    // the frames' bytes at their most at once (the memory model)
    uint64_t frames_peak_bytes = 0;
    // the time spent in the tracker: reads, merges, the occurrences of the kept paths (ms)
    double ms = 0;
};

/**
 * The label and trace trackers (SupportTracker), one per request, used by its patterns in
 * request order (begin_pattern before each pattern's count(), end_pattern after it). Used by
 * one thread.
 */
class PathTracker final : public graph::pattern::SupportTracker {
  public:
    /**
     * Throws std::invalid_argument for the TRACE level on an index without coordinates or the
     * record mapping, or on a graph other than BASIC (coordinates carry no strand elsewhere).
     */
    PathTracker(const PathSupportEnv &env, const PathSupportOptions &options);
    ~PathTracker() override;

    PathTracker(const PathTracker&) = delete;
    PathTracker& operator=(const PathTracker&) = delete;

    // before one pattern's search: the pattern of |length| bases (> k)
    void begin_pattern(size_t length);
    // after it: frees what the pattern held (its dictionary and the row cache's allotment)
    // and restores the oracle's path cache
    void end_pattern();
    /**
     * After the pattern's search, before end_pattern, when other reads follow it (the mirror
     * walks of a predicate's selection, pattern_supported.hpp): empties the row cache, returns
     * its allotment to the account and restores the oracle's path cache, so that the next
     * reader's allotment replaces this one instead of adding to it. The dictionary stays (the
     * labels' answer names its labels); no row of the pattern is read after it.
     */
    void release_rows();

    // SupportTracker
    bool prepare(const std::vector<graph::pattern::Context> &anchors) override;
    Verdict open(const graph::pattern::SearchState &anchor) override;
    Verdict push(const graph::pattern::SearchState &child) override;
    void pop() override;
    Verdict complete(const graph::pattern::PathView &walk) override;
    const char* stop_reason() const override { return stop_; }

    /**
     * The support of the walk the frames hold (complete() answered ALIVE, before its pop): per
     * label carrying every k-mer of it (ascending dictionary ids), whether a chain survived
     * (record_verified) and, with |occurrences| and a chain level, the walk's occurrences as
     * runs of consecutive starts. Each start is mapped to its record (seq_id, 1-based local
     * start) when the index has the record mapping, else stated as its column coordinate
     * (record bounds unknown). Reads the clock every 4,096 runs mapped: false when the work
     * time passed first (stop_reason() "time"; |out| incomplete).
     */
    struct Occurrences {
        // (seq_id, 1-based start) with the record mapping, else (column coordinate, 0)
        uint64_t a = 0;
        uint64_t b = 0;
        // consecutive starts from (a, b) on: b + i (a + i without the record mapping)
        uint64_t count = 0;
    };
    struct WalkLabel {
        graph::traversal::LabelId label = 0;
        bool verified = false;
        uint64_t occurrences = 0;
        std::vector<Occurrences> runs;
    };
    bool walk_labels(bool occurrences, std::vector<WalkLabel> *out);
    // before walk_labels(): how many labels it gives, and the runs of label i (bytes before)
    size_t walk_num_labels() const;
    size_t walk_num_runs(size_t i) const;

    /**
     * The labels the top frame's support holds, at any depth (ascending dictionary ids): at
     * the record level those with a chain (a record holding the branch so far), else every
     * label of the frame (carrying every k-mer of it). Empty when no frame is open. What a
     * predicate is evaluated on (pattern_supported.hpp): a branch's support only shrinks as
     * it is extended, so a monotone predicate false on it is false on every walk below it.
     */
    void frame_support(std::vector<graph::traversal::LabelId> *out) const;

    // the pattern's dictionary: the labels of the anchors' rows (LabelId == index)
    const std::vector<graph::traversal::LabelRef>& labels() const;
    // whether the frames carry chains (TRACE frames), and whether the search prunes on them
    bool chains() const { return chains_; }
    bool prunes_on_chains() const { return prune_on_chains_; }
    graph::traversal::Support level() const { return options_.level; }
    const PathSupportOptions& options() const { return options_; }
    const PathSupportEnv& env() const { return env_; }
    // the statements of the rows the account could not hold (phase "extension")
    const Json::Value& rows_refused() const { return rows_refused_; }
    const PathSupportWork& work() const { return work_; }

    // the row cache's allotment for a pattern starting when the account has |left| bytes left
    static uint64_t cache_allotment(uint64_t left);

  private:
    struct Impl;
    const PathSupportEnv env_;
    const PathSupportOptions options_;
    std::unique_ptr<Impl> impl_;
    bool chains_ = false;
    bool prune_on_chains_ = false;
    const char *stop_ = nullptr;
    Json::Value rows_refused_ = Json::Value(Json::arrayValue);
    PathSupportWork work_;

    // the row of |key| (the k-mer |kmer|) as the steps read it; null: stopped (stop_)
    const graph::traversal::support_step::RowRuns* row(uint64_t key, std::string_view kmer);
};

/**
 * The supported paths of one pattern (PathSink): every complete supported walk the extension
 * hands over, kept under the request's mode — count: none; partial: the first max_paths;
 * all_or_count: all while they are at most max_paths, none once they are more — with its labels
 * from the tracker's frames (PathSupportOptions::labels). What a kept path holds (its descriptor,
 * path_descriptor_bytes, and its labels and occurrence runs) is admitted before it is copied,
 * within an allowance of half of what the account has left when the pattern starts (the other
 * half kept for the frames and rows); the first path that does not fit ends the list
 * (output_cut(): partial lists the paths before it, all_or_count none). Never stops the search
 * itself, but for the work time while it maps a path's occurrences (stop_reason() "time").
 */
class SupportedPathSink final : public graph::pattern::PathSink {
  public:
    explicit SupportedPathSink(PathTracker &tracker);
    ~SupportedPathSink() override;

    SupportedPathSink(const SupportedPathSink&) = delete;
    SupportedPathSink& operator=(const SupportedPathSink&) = delete;

    // before one pattern's search, after PathTracker::begin_pattern; |mode| and |max_paths|
    // as the request's
    void begin_pattern(graph::pattern::Mode mode, uint64_t max_paths);
    // after the pattern's answer was built (labels_answer before it, PathTracker::end_pattern
    // after it): frees what its paths still hold, but the descriptors of the first |returned|
    // (the result objects of the paths the answer lists stay in the account, as a released
    // path's do)
    void end_pattern(size_t returned = 0);

    /**
     * After the search, before labels_answer: keeps the kept paths |indices| (ascending
     * indices into paths()) and frees what the others hold, so that paths() is those, in that
     * order. A predicate's selection among held paths (pattern_supported.hpp) lists the
     * selected ones only.
     */
    void retain(const std::vector<size_t> &indices);

    bool accept(const graph::pattern::PathView &path,
                const graph::pattern::SupportTracker *support) override;
    const char* stop_reason() const override { return stop_; }

    // the kept paths, in answer order
    const std::vector<graph::pattern::Context>& paths() const { return paths_; }
    // the supported paths handed over (kept or not)
    uint64_t accepted() const { return accepted_; }
    // a path did not fit the allowance: the list ends before it (partial), or holds none
    // (all_or_count)
    bool output_cut() const { return output_cut_; }
    // all_or_count: more than max_paths were handed over, none is kept
    bool dropped() const { return dropped_; }

    /**
     * The labels of the paths |x| releases (the first |x|.returned kept paths), as
     * PatternRetrieval::retrieve_paths states a released path's labels: the label order,
     * by_label, the counts (labels by support, occurrences), partial's max_labels and
     * max_occurrences_per_label cuts, each path's labels with its support and occurrences, the
     * statements and the work, withheld (all_or_count) or cut (partial) when they are not
     * complete. Their objects are charged to the account and their text to the answer's volume
     * before they are built, the work time read before each path's (stop {output, time |
     * max_memory}). |released|: the paths the release names (more than |x|.returned when the
     * list was cut by the allowance: stop {output, max_memory}). Only with
     * PathSupportOptions::labels.
     */
    LabelsAnswer labels_answer(const graph::pattern::Extraction &x, uint64_t released,
                               size_t length, const Json::Value &graph_name);

  private:
    struct Kept;
    PathTracker &tracker_;
    graph::pattern::Mode mode_ = graph::pattern::Mode::ALL_OR_COUNT;
    uint64_t max_paths_ = 0;
    uint64_t allowance_ = 0;
    uint64_t held_ = 0;
    uint64_t accepted_ = 0;
    bool output_cut_ = false;
    bool dropped_ = false;
    bool keeping_ = true;
    const char *stop_ = nullptr;
    std::vector<graph::pattern::Context> paths_;
    std::vector<Kept> kept_;

    void drop_all();
};

/**
 * What the release of a pattern's supported paths is (SPEC §20.5), from the engine's
 * count() with the tracker and the sink, and the sink: mode all_or_count every supported path
 * or none — withheld ANCHORS_ABOVE_THRESHOLD (the extension not admitted), DISCOVERY_BUDGET
 * (max_steps), DEADLINE (time), THRESHOLD_CROSSED (stop_at_threshold), EXTERNAL (the
 * tracker's or the sink's stop: annotation_budget), COUNT_ABOVE_THRESHOLD (the supported
 * paths exact and more than max_paths) —; partial the supported paths completed before a stop
 * of the extension, a prefix in answer order, cut by the stop's reason (EXTERNAL: the
 * tracker's), else MAX_PATHS when there are more than max_paths; none after a stop before the
 * extension. Count: no release (nullopt). |returned| is what the release names (the sink's
 * kept paths, or as many as it would have kept: the sink's own memory cut is the labels'
 * answer's, labels_answer).
 */
std::optional<graph::pattern::Extraction>
supported_release(const graph::pattern::Result &result, const graph::pattern::Request &request,
                  const SupportedPathSink &sink, uint64_t *returned);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_SUPPORT_HPP__
