#ifndef __METAGRAPH_CLI_PATTERN_RETRIEVAL_IMPL_HPP__
#define __METAGRAPH_CLI_PATTERN_RETRIEVAL_IMPL_HPP__

/**
 * The labelled retrieval's internals (pattern_retrieval.hpp), shared by the files that make up
 * PatternRetrieval: pattern_retrieval.cpp (the labels of contexts and paths) and
 * pattern_selection.cpp (the selection pass of a predicate); and by the supported-path search,
 * which reads through the retrieval's oracle and account (pattern_support.cpp,
 * pattern_supported.cpp), and its tests. Not an interface of the route: pattern.cpp reaches
 * these through PatternRetrieval and the supported-path search only. It holds the memory and
 * text models, the request's account, the rows' states and PatternRetrieval::Impl.
 */

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <json/json.h>

#include "json_helpers.hpp"
#include "pattern_retrieval.hpp"
#include "pattern_predicate.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"


namespace mtg {
namespace cli {
namespace retrieval {

/**
 * The names of a pattern's label dictionary as the answer may carry them (SPEC §8.9,
 * unrepresentable_label_name): a name is written only when it is valid UTF-8. jsoncpp
 * replaces the bytes of any other by U+FFFD, and the replaced name can be another column's,
 * so a client that fed it into a predicate would select the wrong column. Each dictionary
 * label is validated once, where its name is first copied into the answer (at the cost of
 * that copy; an answer that lists no name validates nothing), and the verdict kept per
 * label. The first label found unrepresentable is recorded: the lists are still built and
 * their charges released as for an answered entry, and the route refuses the entry in its
 * slot, naming the column and never its bytes.
 */
class NameCheck {
  public:
    explicit NameCheck(const std::vector<graph::traversal::LabelRef> &dict)
          : dict_(dict), state_(dict.size(), kUnknown) {}

    // the name of dictionary label |id| when the answer may carry it; else an empty string,
    // the label's column recorded if none was before
    const std::string& name(graph::traversal::LabelId id) {
        uint8_t &state = state_.at(id);
        if (state == kUnknown)
            state = utils::valid_utf8(dict_[id].name) ? kValid : kInvalid;
        if (state == kValid)
            return dict_[id].name;
        if (!column_)
            column_ = dict_[id].column;
        return empty_;
    }
    const std::optional<graph::traversal::Column>& unrepresentable() const { return column_; }

  private:
    static constexpr uint8_t kUnknown = 0;
    static constexpr uint8_t kValid = 1;
    static constexpr uint8_t kInvalid = 2;
    const std::vector<graph::traversal::LabelRef> &dict_;
    std::vector<uint8_t> state_;
    std::optional<graph::traversal::Column> column_;
    std::string empty_;
};

/**
 * The memory model of the account (§5.3): deterministic prices of what the answer holds,
 * never a measurement, so that where a request stops does not depend on the allocator. What
 * the reads return is charged by the reads themselves (DecodeBudget, the model of
 * decode_budget.hpp); the rest is priced here.
 */
// a released context: its result object (k-mer, instance, offset, strand, node, row) and its
// descriptor, beside its k-mer twice (the result's strings); charged before the object is built
inline uint64_t context_bytes(size_t k) { return 512 + 2 * k; }
// Every copy of a label's name the request holds is priced where it is made, at its length (a
// flat price would let one long name in many contexts build an answer far past max_memory_mb):
// a dictionary label: its LabelRef and name in the recorder, its counters here, and the copy of
// its name in each pattern's placement (the LabelQuery of step 2 copies the dictionary; one at
// a time); charged inside the read that names it
inline uint64_t label_name_bytes(std::string_view name) { return 192 + 2 * name.size(); }
// a label object of a result (column, support, the occurrences count, the list) and its copy
// of the label's name
inline uint64_t label_entry_bytes(std::string_view name) { return 256 + name.size(); }
// a by_label entry of a pattern (graph, column, contexts, contexts_suffix, occurrences) and its
// copy of the label's name
inline uint64_t by_label_bytes(std::string_view name) { return 512 + name.size(); }
// a placed occurrence object (seq_id, record, strand, nt_coords, nt_length), its record name
// beside it
inline uint64_t occurrence_bytes(std::string_view record) { return 256 + record.size(); }
// a (kmer_coord, offset, strand) object of the global placement
constexpr uint64_t kGlobalOccurrenceBytes = 192;
// a placed occurrence in a label's deduplication set (§5.4)
constexpr uint64_t kDedupBytes = 64;
// a statement of a row: a rows_refused entry (k-mer, row, phase, reason, needed_bytes,
// available_bytes) or an anchors_truncated entry (k-mer, row, cap, total); reserved before
// the read that can produce it, so that it always fits, and held with the answer
inline uint64_t statement_bytes(size_t k) { return 384 + k; }
// a released path (long_search "paths"): see path_descriptor_bytes in the header
constexpr uint64_t kPathKmerBytes = 192;

/**
 * The compact JSON text of the labels built for the answer, from above (AnswerVolume):
 * estimated before the objects are built, so that the work time is read with them counted.
 * Integers at their widest (20 digits), strings as string_text_bytes.
 */
// a count {"relation":"at_least","unit":"placed_occurrences","value":<20 digits>}: 79
constexpr uint64_t kCountText = 96;
// a result's label fields: ,"labels":[],"labels_status":"output_budget","labels_total":<20
// digits>,"support":"kmer" (about 100)
constexpr uint64_t kResultLabelsText = 112;
// a label object without its name and list: {"column":,"occurrence_list":[],"occurrences":
// <count>,"support":"kmer"},
constexpr uint64_t kLabelText = 64 + kCountText;
// a placed occurrence without its record name: {"nt_coords":"<20>-<20>","nt_length":<20>,
// "record":,"seq_id":<20>,"strand":"+"}, (142)
constexpr uint64_t kOccurrenceText = 144;
// a global one: {"kmer_coord":<20>,"offset":<20>,"strand":"+"}, (79)
constexpr uint64_t kGlobalOccurrenceText = 80;
// a by_label entry without its names: {"column":,"contexts":<count>,"contexts_suffix":
// <count>,"graph":,"occurrences":<count>}, (68 and three counts)
constexpr uint64_t kByLabelText = 72 + 3 * kCountText;
// paths (long_search "paths"): a result's label fields, ,"labels":[],"labels_excluded_
// unverified":<20 digits>,"labels_status":"output_budget","labels_total":<20 digits>,
// "support":"label_intersection" (about 160)
constexpr uint64_t kPathResultLabelsText = 176;
// a path's label object without its name and list: {"column":,"occurrence_list":[],
// "occurrences":<count>,"support":"label_intersection"},
constexpr uint64_t kPathLabelText = 80 + kCountText;
// a path's by_label entry without its names: {"column":,"graph":,"occurrences":<count>,
// "paths":<count>,"paths_record_verified":<count>}, (73 and three counts)
constexpr uint64_t kPathByLabelText = 80 + 3 * kCountText;

// a string's text: quoted, every byte jsoncpp may escape (a control character, a quote, a
// backslash, a byte of a non-ASCII character: \uXXXX) counted as 6
inline uint64_t string_text_bytes(const char *begin, const char *end) {
    uint64_t bytes = 2;
    for (const char *c = begin; c != end; ++c) {
        const unsigned char u = static_cast<unsigned char>(*c);
        bytes += u < 0x20 || u >= 0x80 || u == '"' || u == '\\' ? 6 : 1;
    }
    return bytes;
}

inline uint64_t string_text_bytes(std::string_view s) {
    return string_text_bytes(s.data(), s.data() + s.size());
}

constexpr uint64_t kNoKey = graph::traversal::npos;

inline Json::Value count_json(graph::pattern::Relation relation, uint64_t value,
                              graph::pattern::Unit unit) {
    Json::Value v;
    v["value"] = relation == graph::pattern::Relation::UNKNOWN ? Json::Value() : uint_json(value);
    v["relation"] = graph::pattern::to_string(relation);
    v["unit"] = graph::pattern::to_string(unit);
    return v;
}

inline const char* strand_of(graph::pattern::Orientation o) {
    return graph::pattern::strand_symbol(o);
}

inline uint8_t strand_rank(graph::pattern::Orientation o) {
    switch (o) {
        case graph::pattern::Orientation::FORWARD: return 0;
        case graph::pattern::Orientation::REVERSE: return 1;
        case graph::pattern::Orientation::PALINDROMIC: return 2;
    }
    return 3;
}

// the account of §5.3: one per request, every item charged before it is held (what an
// unbudgeted read returned excepted: force())
class Account {
  public:
    explicit Account(uint64_t max) : max_(max) {}
    // refused when it does not fit what is left; once an unbudgeted read's names were forced
    // past the maximum nothing fits any more (left() is 0: no unsigned wrap of max_ - held_)
    bool charge(uint64_t bytes) {
        if (bytes > left())
            return false;
        held_ += bytes;
        peak_ = std::max(peak_, held_);
        return true;
    }
    // what an unbudgeted read returned, held whether or not it fits (the read is done); the
    // answer states the unbudgeted access, and the peak shows any excess
    void force(uint64_t bytes) {
        held_ += bytes;
        peak_ = std::max(peak_, held_);
    }
    void release(uint64_t bytes) {
        assert(bytes <= held_);
        held_ -= std::min(bytes, held_);
    }
    uint64_t left() const { return held_ < max_ ? max_ - held_ : 0; }
    uint64_t held() const { return held_; }
    uint64_t peak() const { return peak_; }

  private:
    uint64_t max_;
    uint64_t held_ = 0;
    uint64_t peak_ = 0;
};

// What the reads made of a row (the annotation key of one or more contexts)
enum class RowStatus {
    PENDING,
    // read, every label of the row listed
    COMPLETE,
    // read, more labels than max_labels_per_anchor: the first that many listed, the total
    // stated
    TRUNCATED,
    // the account could not hold the row's read (stated in rows_refused)
    REFUSED,
    // not read: the reads stopped before it (time, work), or all_or_count stopped reading
    // once its results could no longer be published
    NOT_READ,
};

inline const char* to_string(RowStatus s) {
    switch (s) {
        case RowStatus::COMPLETE: return "complete";
        case RowStatus::TRUNCATED: return "truncated";
        case RowStatus::REFUSED: return "refused";
        case RowStatus::PENDING:
        case RowStatus::NOT_READ: return "not_read";
    }
    return "not_read";
}

struct RowState {
    uint64_t key = kNoKey;
    // the first context of the pattern with this key: its k-mer names the row
    size_t first = 0;
    RowStatus status = RowStatus::PENDING;
    graph::traversal::LabelRecorder::NodeLabels labels;
    // placement (step 2): read, refused by the account, or not reached
    bool placed = false;
    bool place_refused = false;
    graph::traversal::LabelQuery::NodeHits hits;
};

/**
 * The memory model of the selection pass (SPEC §19.9), beside the items above.
 */
// a raw context kept for the pass (TestedContext: orientation, offset, node, row key, the
// reverse complement's key, row indices, its decision); charged as the engine releases it
constexpr uint64_t kTestedBytes = 64;
// a distinct row of the pass: its entry in the pass's row table (key, first context, mirror,
// status, the hits' list) and in the index of its key; charged before the row is added
constexpr uint64_t kSelectionRowBytes = 96;
// a selected context's selection_labels (the list and its strings), the names' lengths beside
inline uint64_t selection_labels_bytes(uint64_t names_length) { return 32 + names_length; }
// beside them its selection_strands (the list; per label the row it was found on, a byte kept,
// and its string in the answer, at most 18 characters): 32 + 24 per label
inline uint64_t selection_strands_bytes(uint64_t labels) { return 32 + 24 * labels; }
// a label of a pattern's motif union (SelectionRequest::motif): its entry in the pass's index
// of the union (label, contexts, first, the context that last counted it, rows), its
// MotifLabel in the answer and its object in the motif block (column, contexts, strands) with
// its copy of the label's name; charged before the entry is made, held with the answer
inline uint64_t motif_label_bytes(std::string_view name) { return 64 + label_entry_bytes(name); }
// an id of the motif's list of undecided labels (Bound::eval3's maybes), held while it runs
constexpr uint64_t kMotifUndecidedBytes = sizeof(graph::traversal::LabelId);

/**
 * A row of the selection pass held for the projection predicate_only (retrieve_given): the
 * predicate's labels on the row (hits restricted to Bound::labels(), no coordinates) and what
 * the account holds for them (the read's hits and the row's entry, kSelectionRowBytes).
 */
struct KeptRow {
    graph::traversal::LabelQuery::NodeHits hits;
    uint64_t bytes = 0;
};

/**
 * The text of the labels about to be built for one pattern's answer (AnswerVolume's pending
 * part, from above): added before they are built, so that the work time is read with them
 * counted, dropped for what is not built, settled once built.
 */
class PendingText {
  public:
    explicit PendingText(AnswerVolume *volume) : volume_(volume) {}

    void add(uint64_t bytes) {
        pending_ += bytes;
        if (volume_)
            volume_->add_pending(bytes);
    }
    void drop(uint64_t bytes) {
        pending_ -= std::min(bytes, pending_);
        if (volume_)
            volume_->drop_pending(bytes);
    }
    // none of it is built
    void drop_all() { drop(pending_); }
    // all of it was built: from now on it is written only
    void settle() const {
        if (volume_)
            volume_->settle(pending_);
    }

  private:
    AnswerVolume *volume_;
    uint64_t pending_ = 0;
};

} // namespace retrieval


struct PatternRetrieval::Impl {
    Impl(const graph::AnnotatedDBG &anno_graph, const RetrievalLimits &limits,
         graph::pattern::Budget &budget, const RetrievalHooks *hooks, AnswerVolume *volume)
          : oracle(anno_graph, hooks ? hooks->coord_to_header : nullptr),
            limits(limits), budget(budget),
            account(hooks && hooks->max_memory_bytes ? hooks->max_memory_bytes
                                                     : limits.max_memory_bytes),
            volume(volume) {
        if (hooks) {
            deny_decode = hooks->deny_decode;
            output_hook = hooks->output_hook;
            occurrences_hook = hooks->occurrences_hook;
        }
        // the reads under the deadline are decoded in paced chunks (as /traverse's)
        oracle.pacer().target_ms = limits.chunk_target_ms;
        if (hooks && hooks->read_hook)
            oracle.test_read_hook = hooks->read_hook;
        recorder = std::make_unique<graph::traversal::LabelRecorder>(
                oracle, graph::traversal::LabelKind::COLUMN,
                std::max<uint64_t>(1, limits.max_labels_per_anchor));
        // no cache: every row is read once per step, its keys deduplicated, and nothing a
        // cache held would be in the account (a fixed allotment of zero)
        recorder->set_max_cache_bytes(0);
    }

    graph::traversal::LabelOracle oracle;
    const RetrievalLimits limits;
    graph::pattern::Budget &budget;
    retrieval::Account account;
    std::unique_ptr<graph::traversal::LabelRecorder> recorder;
    bool budgeted = false;
    bool place = false;          // placement record or global, and requested
    bool records = false;        // placement record
    // the request's annotation work (units) so far
    uint64_t units = 0;
    // tests: DecodeBudget::deny of every read
    std::function<bool(uint64_t)> deny_decode;
    // the request's answer volume (null: none), and the tests' hook before the work time is
    // read for a context's labels
    AnswerVolume *volume = nullptr;
    std::function<void(size_t)> output_hook;
    // tests: called before each context's occurrences are made and each path's verified
    std::function<void(size_t)> occurrences_hook;
    // the request's counters (RetrievalCounters)
    RetrievalCounters counters;
    // the pattern being released (admit_context): its descriptors' bytes, held, and what they
    // may hold (all_or_count: what the account had left; partial: half of it); false once a
    // descriptor did not fit
    uint64_t descriptors = 0;
    uint64_t allowance = 0;
    bool admitting = false;

    uint64_t statement() const { return retrieval::statement_bytes(oracle.get_k()); }

    // what a read may hold: the account's remainder, less the statements reserved for its rows
    annot::matrix::DecodeBudget decode_budget(uint64_t reserve) const {
        annot::matrix::DecodeBudget decode(account.left() - std::min(reserve, account.left()));
        if (deny_decode)
            decode.deny = deny_decode;
        return decode;
    }

    // the stop of the current pattern's annotation work (the first one), and whether the
    // reads themselves stopped (time, work, an unbudgeted read the account could not hold)
    std::optional<std::pair<std::string, std::string>> stop;
    bool time_stop = false;
    bool read_stop = false;

    graph::traversal::ReadPacing pacing() {
        graph::traversal::ReadPacing p;
        const graph::pattern::Deadline *d = &budget.deadline();
        p.ms_left = [d]() {
            return d->time_budget_ms() - d->finalize_reserve_ms() - d->elapsed_ms();
        };
        p.stop = [d]() { return d->work_expired(); };
        return p;
    }

    void set_stop(const char *phase, const char *reason) {
        if (!stop)
            stop = std::make_pair(std::string(phase), std::string(reason));
        read_stop |= std::string_view(phase) != "output";
    }

    // the time and work checks before a read; false: stopped (recorded)
    bool may_read(const char *phase) {
        if (!budget.check_time()) {
            set_stop(phase, "time");
            time_stop = true;
            return false;
        }
        if (units >= limits.max_annotation_work) {
            set_stop(phase, "max_annotation_work");
            return false;
        }
        return true;
    }

    // the work of a read, charged to the request and to the pattern
    void charge_work(uint64_t u, uint64_t *pattern_units) {
        units += u;
        *pattern_units += u;
    }

    /**
     * The clock of the retrieval's own work between the reads: the occurrences made for the
     * counts and the output, the paths' label lists and their verification. Asked before |u|
     * more units are done: the clock is read when they would take the units since its last
     * reading past Budget::kClockStride, as the engine's steps read it; false when the work
     * time passed (the caller states the stop and does not do the work).
     */
    uint64_t unclocked = 0;
    bool may_work(uint64_t u) {
        if (unclocked + u > graph::pattern::Budget::kClockStride) {
            unclocked = 0;
            if (!budget.check_time())
                return false;
        }
        unclocked += u;
        return true;
    }

    // What a refused read is charged: the units of what it decoded (FetchRefusal::units: the
    // row was read, and refused for its demand or its names), or a key's units when its read
    // itself did not fit (its units are not known): a refused row is work like a row read, so
    // that refusals cannot go on past the work budget
    static uint64_t refused_units(const graph::traversal::FetchRefusal &r) {
        return r.units ? r.units : graph::traversal::kKeyUnits;
    }

    Json::Value refusal_json(const retrieval::RowState &row,
                             const std::vector<RetrievalContext> &contexts,
                             const char *phase, const graph::traversal::FetchRefusal *r) const {
        return refusal_json(contexts[row.first].kmer, row.key, phase, r);
    }

    // the statement of a refused row: its k-mer, its row, the phase, what reading the row
    // alone needed (at least) against what the account had left
    Json::Value refusal_json(const std::string &kmer, uint64_t key, const char *phase,
                             const graph::traversal::FetchRefusal *r) const {
        using graph::traversal::FetchRefusal;
        Json::Value v;
        v["kmer"] = kmer;
        v["row"] = uint_json(graph::AnnotatedDBG::graph_to_anno_index(key));
        v["phase"] = phase;
        v["reason"] = "max_memory";
        uint64_t need = 0;
        if (r) {
            need = r->cause == FetchRefusal::DECODE ? r->need
                 : r->cause == FetchRefusal::NAMES ? r->demand + r->names_bytes : r->demand;
        }
        v["needed_bytes"] = r ? Json::Value(uint_json(need)) : Json::Value();
        v["available_bytes"] = uint_json(r ? r->left : account.left());
        return v;
    }

    // step 1: the labels of every row, at most max_labels_per_anchor each
    void discover(std::vector<retrieval::RowState> &rows,
                  const std::vector<RetrievalContext> &contexts,
                  graph::pattern::Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                  uint64_t *held_lists, Json::Value *refused);

    // step 2: the coordinates of the rows' labels; |needed| (paths): only the rows it marks
    // (those of a path with at least one label on every k-mer), every row read with a label
    // when null (contexts); |dict| the labels the rows' ids name (the recorder's dictionary
    // when null; retrieve_given: the predicate's labels found on the selected rows)
    void place_rows(std::vector<retrieval::RowState> &rows,
                    const std::vector<RetrievalContext> &contexts,
                    graph::pattern::Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                    uint64_t *held_hits, Json::Value *refused,
                    const std::vector<bool> *needed = nullptr,
                    const std::vector<graph::traversal::LabelRef> *dict = nullptr);

    // ---- the selection of a predicate (pattern_selection.cpp)

    graph::pattern::GraphMode mode = graph::pattern::GraphMode::BASIC;
    // bind(): the request's predicate bound to the index (its bytes held for the request), and
    // the request's selection limits; null until bind() was called
    std::unique_ptr<predicate::Binding> binding;
    SelectionLimits selection_limits;
    // the selection's reads, one LabelQuery over Bound::labels() for the request (built at its
    // first read: the bound predicate's model prices this copy of its labels), and how it reads
    std::unique_ptr<graph::traversal::LabelQuery> selection_query;
    const char *selection_access = "rows";
    // the request's selection work (units, max_predicate_work) so far
    uint64_t predicate_units = 0;
    // the pattern being released into the pass (admit_tested): its descriptors' bytes, held
    // until end_selection(), what they may hold, and false once one did not fit
    uint64_t tested_bytes = 0;
    uint64_t tested_allowance = 0;
    bool admitting_tested = false;
    // the selected rows held for the projection predicate_only (retrieve_given), by key
    std::unordered_map<uint64_t, retrieval::KeptRow> kept_rows;
    uint64_t kept_bytes = 0;

    // the predicate's labels found on the kept rows (ids into Bound::labels(), ascending) and
    // their names' bytes held (label_name_bytes each: the labels' dictionary for
    // retrieve_given, as a recorder's dictionary label is priced)
    std::vector<graph::traversal::LabelId> kept_labels;

    // frees what the pass of the last pattern still held (its descriptors, its kept rows)
    void end_selection();

    // ---- the steps retrieve_rows() and retrieve_paths() share (pattern_retrieval.cpp)

    // the pattern's annotation work into |a|, on every return: its rows and units, the
    // account's peak, the steps' times, the statements' text (they stay in the answer whatever
    // else does: AnswerVolume), and the pattern's counters into the request's
    void finish_work(LabelsAnswer &a, uint64_t rows_read, uint64_t pattern_units,
                     double discovery_ms, double placement_ms);
    // Before the reads, over the |keep| of the |released| items the route built (their
    // descriptors charged as the engine released them, before their result objects were
    // built: the first that did not fit ended the list). True when nothing is read: the engine
    // withheld the release, or all_or_count's list was cut (withheld output_budget, stop
    // {output, max_memory}); |a| is then the answer, finished. Partial's cut is stated in |a|.
    bool skip_reads(LabelsAnswer &a, const graph::pattern::Extraction &x, size_t keep,
                    uint64_t released, uint64_t descriptors, bool partial);
    // Every row's statement after the reads: the rows read (counters.rows_distinct) and each
    // truncated row's anchors_truncated entry (|named|[row.first] names it by its k-mer).
    // Returns whether every row was read completely; |truncated|: whether one was truncated
    bool state_rows(LabelsAnswer &a, const std::vector<retrieval::RowState> &rows,
                    const std::vector<RetrievalContext> &named, bool *truncated);
    // The output of item |i|'s labels (a context's or a path's): their |bytes| charged first
    // (where the output stops does not depend on the machine), then the work time read with
    // their |text| counted in the answer's volume. False when the item was |late| or |refused|
    // already, or the account or the time refuses it now: stop {output, time | max_memory}
    // (a time stop is recorded in the budget too: the later patterns answer as after any time
    // stop), nothing of it charged or pending
    bool commit_output(size_t i, uint64_t bytes, uint64_t text, bool late, bool refused,
                       retrieval::PendingText *pending);
    // all_or_count without every label of every item: withheld, the reason named (the reads'
    // own budgets — a refused row, the work, an unbudgeted read the account could not hold —
    // before the answer's memory, before a truncated anchor), what was |held| for the output
    // released, nothing pending; what explains it (rows_refused, anchors_truncated) stays
    void withhold_all(LabelsAnswer &a, bool rows_refused, bool output_stopped, bool truncated,
                      uint64_t held, retrieval::PendingText *pending);
};

/**
 * What retrieve_given() gives the shared body of retrieve() in place of the discovery: the
 * dictionary the rows' labels are renamed into (the predicate's labels found on the kept rows:
 * dict[i] is Bound::labels()[(*ids)[i]], |ids| ascending) and the kept rows by key (their hits
 * name Bound::labels()).
 */
struct PatternRetrieval::Given {
    std::vector<graph::traversal::LabelRef> dict;
    const std::vector<graph::traversal::LabelId> *ids = nullptr;
    const std::unordered_map<uint64_t, retrieval::KeptRow> *rows = nullptr;
};

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_RETRIEVAL_IMPL_HPP__
