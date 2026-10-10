#ifndef __METAGRAPH_CLI_PATTERN_SUPPORTED_HPP__
#define __METAGRAPH_CLI_PATTERN_SUPPORTED_HPP__

/**
 * A request's predicate on the supported paths of a pattern longer than k (long_search
 * "supported_paths" with a predicate; SPEC-pattern-search.md §20.9). The supported-path search
 * (pattern_support.hpp: PathTracker, SupportedPathSink) finds the walks some label supports;
 * PathSelection sits between it and the engine, as the engine's SupportTracker and PathSink at
 * once, and decides which supported walks the predicate selects:
 *
 *  - a walk's support is the labels supporting it at the search's level (record level: the
 *    labels a record of which holds it whole; label level: the labels carrying every k-mer),
 *    strand-consistent: the walk as spelled, never mixed with its reverse complement;
 *  - predicate_strands "context", and every walk on CANONICAL and PRIMARY graphs (a k-mer and
 *    its reverse complement share one row, so a walk and its reverse-complement walk have one
 *    support): the walk is decided when it completes, on its own support; a selected one is
 *    handed to the SupportedPathSink, which keeps it under max_paths (the threshold is on the
 *    selected paths: the selection owns stop_at_threshold's max_paths, the engine's is off).
 *    A MONOTONE normal form (no none, no not) also prunes: the support only shrinks along a
 *    walk, so a branch on whose support the normal form is already false holds no selected
 *    walk, and it is not followed (counted apart from the support's prunes);
 *  - predicate_strands "either" on a BASIC graph: a walk is decided on its support and its
 *    reverse-complement walk's (the mirror), so the decision waits for the search to end: the
 *    supported walks are held (at most max_predicate_contexts, each charged before it is held;
 *    the SupportedPathSink keeps their paths and labels), and decided after it — the mirror's
 *    support taken from the held walks when the search covered the mirrors (strands "both", or
 *    a palindromic pattern: a mirror not held has no support there), else read: the mirror
 *    walk's k-mers looked up and its rows read by a PathTracker of its own, as the search
 *    reads (n rows, under max_annotation_work, the account and the deadline).
 *
 * Use per pattern (the request's PathTracker and SupportedPathSink, begun by the caller):
 *   selection.begin_pattern(options, L); request.support = request.sink = &selection;
 *   request.max_paths = no limit (the selection's own threshold); result = search.count(...);
 *   answer = selection.finish(result); results from sink.paths()[0, answer.listed), labels
 *   from sink.labels_answer(...) as without a predicate; selection.end_pattern().
 */

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <json/json.h>

#include "graph/alignment/pattern_search.hpp"
#include "pattern_retrieval.hpp"
#include "pattern_support.hpp"


namespace mtg {
namespace cli {

namespace predicate {
class Bound;
}

// What one pattern's selection of supported paths is asked
struct PathSelectionOptions {
    graph::pattern::Mode mode = graph::pattern::Mode::ALL_OR_COUNT;
    // the selected paths listed (all_or_count's threshold, partial's cap) and, with
    // stop_at_threshold, the selected count that stops the search
    uint64_t max_paths = 1'000;
    bool stop_at_threshold = false;
    // "either" on a BASIC graph: the supported walks held for their decisions, at most this
    // many (all_or_count and count: more is not admitted; partial: the first ones)
    uint64_t max_predicate_contexts = 100'000;
    // predicate_strands "either" on a BASIC graph: a walk's mirror is part of its decision
    bool either = false;
    // the search reaches every mirror walk (strands "both", or a palindromic pattern): a mirror
    // it did not hand over has no support; false: the mirrors are read
    bool mirrors_searched = false;
    // a projection that reads labels: each listed path's selection_labels are kept
    bool selection_labels = false;
};

/**
 * What one pattern's selection of supported paths did (SPEC §20.9), beside the engine's result:
 * the route writes counts.tested and counts.selected (unit paths), selection.pass, the
 * release of the selected paths (withheld, cut, the paths listed), the selection's own stop and
 * its work.
 */
struct PathSelectionAnswer {
    SelectionPass pass = SelectionPass::NOT_STARTED;
    graph::pattern::Count tested = graph::pattern::Count::unknown(graph::pattern::Unit::PATHS);
    graph::pattern::Count selected = graph::pattern::Count::unknown(graph::pattern::Unit::PATHS);
    // the retrieval modes: the selected paths listed (the sink's paths()[0, listed) once the
    // selection retained them), the paths the release names (more than listed when the memory
    // cut the list: labels_answer's |released|), whether every selected path is listed, and
    // why not (all_or_count: withheld; partial: cut), by the selection itself — the engine's
    // own stop is the route's to compose before these
    uint64_t listed = 0;
    uint64_t named = 0;
    bool complete = false;
    std::optional<std::string> withheld;
    std::optional<std::string> cut;
    // the account ended the list (cut max_memory, stop {output, max_memory}): a cut of the
    // results themselves, which replaces the engine's
    bool output_cut = false;
    // the rows of the mirror walks the account could not hold (phase "selection")
    Json::Value rows_refused = Json::Value(Json::arrayValue);
    // the selection's own stop after the search ({selection, time | max_annotation_work |
    // max_memory}, or {output, max_memory} for the list), and whether the clock touched it
    std::optional<std::pair<std::string, std::string>> stop;
    bool time_limited = false;
    // per listed path, with selection_labels: the predicate's labels in the set it was decided
    // on (ids into Bound::labels(), label order over the listed paths) and beside each the
    // walks whose support carries it (kOnWalk, kOnMirror)
    static constexpr uint8_t kOnWalk = 1;
    static constexpr uint8_t kOnMirror = 2;
    std::vector<std::vector<graph::traversal::LabelId>> selection_labels;
    std::vector<std::vector<uint8_t>> selection_rows;
    // the branches the predicate pruned (a monotone normal form already false on their
    // support), per orientation (the supported paths of an orientation pruned in are at_least)
    uint64_t pruned = 0;
    std::map<graph::pattern::Orientation, uint64_t> pruned_by;
    // work: the decisions' and the pruning's units (max_predicate_work, charged by the caller
    // to the request's selection work), the mirror walks looked up and their rows read (part
    // of the search's annotation work, max_annotation_work)
    uint64_t units = 0;
    uint64_t lookups = 0;
    uint64_t rows = 0;
    uint64_t mirror_units = 0;
    double ms = 0;
};

/**
 * The selection of one request's supported paths by its bound predicate (see the top of this
 * file). One per request, used by its patterns in request order; one thread.
 */
class PathSelection final : public graph::pattern::SupportTracker,
                            public graph::pattern::PathSink {
  public:
    // |tracker| and |sink| are the request's supported-path search; |bound| its predicate,
    // not a constant
    PathSelection(PathTracker &tracker, SupportedPathSink &sink, const predicate::Bound &bound);
    ~PathSelection() override;

    PathSelection(const PathSelection&) = delete;
    PathSelection& operator=(const PathSelection&) = delete;

    /**
     * Before one pattern's search, after the tracker's and the sink's begin_pattern (the sink
     * begun with the mode and the threshold of sink_mode() and sink_threshold()): the pattern
     * of |length| bases.
     */
    void begin_pattern(const PathSelectionOptions &options, size_t length);
    // what the sink is begun with: the request's mode, and its max_paths (decided at
    // completion: the sink keeps the selected paths) or ("either" on BASIC: the sink keeps the
    // held walks) max_predicate_contexts
    static graph::pattern::Mode sink_mode(const PathSelectionOptions &options);
    static uint64_t sink_threshold(const PathSelectionOptions &options);

    /**
     * After the search: the held walks decided (the mirrors read where needed: the work time
     * read before each, the work checked), the selected paths chosen for the release and the
     * sink told to keep those only (SupportedPathSink::retain). |result| is the engine's.
     */
    PathSelectionAnswer finish(const graph::pattern::Result &result);
    // after the pattern's answer was built: frees what the selection still holds
    void end_pattern();

    // SupportTracker (the request's PathTracker's, with the predicate's pruning)
    bool prepare(const std::vector<graph::pattern::Context> &anchors) override;
    Verdict open(const graph::pattern::SearchState &anchor) override;
    Verdict push(const graph::pattern::SearchState &child) override;
    void pop() override;
    Verdict complete(const graph::pattern::PathView &walk) override;
    // PathSink: a complete supported walk, decided or held
    bool accept(const graph::pattern::PathView &path,
                const graph::pattern::SupportTracker *support) override;
    // the stop of the selection (max_paths, max_predicate_contexts), else the tracker's or
    // the sink's
    const char* stop_reason() const override;
    // the selection's own stops are stop_at_threshold's thresholds (the engine diagnoses
    // kNoteLowComplexity beside them, as beside its own); the tracker's and the sink's are
    // budgets
    bool stopped_at_threshold() const override;

  private:
    struct Held;
    PathTracker &tracker_;
    SupportedPathSink &sink_;
    const predicate::Bound &bound_;
    // the predicate's label of each column (Bound::labels())
    std::unordered_map<graph::traversal::Column, graph::traversal::LabelId> label_of_;
    // the mirrors' own tracker (one strand, "either"), made at its first use
    std::unique_ptr<PathTracker> mirror_;

    PathSelectionOptions options_;
    size_t length_ = 0;
    bool prune_ = false;
    bool hold_ = false;
    bool keeping_ = true;
    const char *stop_ = nullptr;
    // the supported walks handed over, decided, selected
    uint64_t supported_ = 0;
    uint64_t tested_ = 0;
    uint64_t selected_ = 0;
    // decided at completion: the selection labels of the paths the sink keeps (aligned with its
    // paths()), their bytes held (kept_bytes_; of them, the text of their strings in the answer,
    // kept_text_) and the text of the listed paths' (answer_bytes_), which stays with the answer
    // when the pattern ends while the rest is given back; the list ended by the account
    // (output_cut_)
    std::vector<std::vector<graph::traversal::LabelId>> kept_labels_;
    uint64_t kept_bytes_ = 0;
    uint64_t kept_text_ = 0;
    uint64_t answer_bytes_ = 0;
    bool output_cut_ = false;
    // "either": the walks held, their bytes, and why holding ended (over_: more than
    // max_predicate_contexts in all_or_count or count; hold_cut_: partial's reason;
    // hold_memory_: the account refused one)
    std::vector<Held> held_;
    uint64_t held_bytes_ = 0;
    bool holding_ = true;
    bool over_ = false;
    const char *hold_cut_ = nullptr;
    bool hold_memory_ = false;
    PathSelectionAnswer answer_;
    // scratch
    std::vector<graph::traversal::LabelId> frame_;
    std::vector<graph::traversal::LabelId> present_;

    // the predicate's labels among the top frame's support (ascending ids into Bound::labels())
    void support_of_frame(std::vector<graph::traversal::LabelId> *out);
    // whether a monotone normal form can still hold below the top frame (its units counted)
    bool may_hold();
    bool hold(const graph::pattern::PathView &path);
    bool decide_now(const graph::pattern::PathView &path);
    void drop_held();
    // the mirror walk of |sequence| read by mirror_: its support (true), or a stop (false,
    // answer_.stop set)
    bool read_mirror(const std::string &sequence, std::vector<graph::traversal::LabelId> *out);
};

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_SUPPORTED_HPP__
