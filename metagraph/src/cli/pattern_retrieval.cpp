#include "pattern_retrieval.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>

#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/coord_to_header.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph::pattern;
using graph::AnnotatedDBG;
using graph::traversal::LabelOracle;
using graph::traversal::LabelRecorder;
using graph::traversal::LabelQuery;
using graph::traversal::LabelRef;
using graph::traversal::LabelKind;
using graph::traversal::LabelId;
using graph::traversal::KeyCost;
using graph::traversal::FetchRefusal;
using graph::traversal::ReadPacing;
using graph::traversal::Column;
using graph::traversal::Coord;
using annot::matrix::DecodeBudget;

namespace {

/**
 * The memory model of the account (§5.3): deterministic prices of what the answer holds,
 * never a measurement, so that where a request stops does not depend on the allocator. What
 * the reads return is charged by the reads themselves (DecodeBudget, the model of
 * decode_budget.hpp); the rest is priced here.
 */
// a released context: its result object (k-mer, instance, offset, strand, node, row) and its
// descriptor, beside its k-mer twice (the result's strings); charged before the object is built
uint64_t context_bytes(size_t k) { return 512 + 2 * k; }
// Every copy of a label's name the request holds is priced where it is made, at its length
// (review GPT-2, finding 3: a result's label was priced 256 whatever its name, and one label
// of 512 KiB in 64 contexts built a 34 MB answer under max_memory_mb 2):
// a dictionary label: its LabelRef and name in the recorder, its counters here, and the copy
// of its name in each pattern's placement (the LabelQuery of step 2 copies the dictionary;
// one at a time); charged inside the read that names it
uint64_t label_name_bytes(std::string_view name) { return 192 + 2 * name.size(); }
// a label object of a result (column, support, the occurrences count, the list) and its copy
// of the label's name
uint64_t label_entry_bytes(std::string_view name) { return 256 + name.size(); }
// a by_label entry of a pattern (graph, column, contexts, contexts_suffix, occurrences) and its
// copy of the label's name
uint64_t by_label_bytes(std::string_view name) { return 512 + name.size(); }
// a placed occurrence object (seq_id, record, strand, nt_coords, nt_length), its record name
// beside it
uint64_t occurrence_bytes(std::string_view record) { return 256 + record.size(); }
// a (kmer_coord, offset, strand) object of the global placement
constexpr uint64_t kGlobalOccurrenceBytes = 192;
// a placed occurrence in a label's deduplication set (§5.4)
constexpr uint64_t kDedupBytes = 64;
// a statement of a row: a rows_refused entry (k-mer, row, phase, reason, needed_bytes,
// available_bytes) or an anchors_truncated entry (k-mer, row, cap, total); reserved before
// the read that can produce it, so that it always fits, and held with the answer
uint64_t statement_bytes(size_t k) { return 384 + k; }

/**
 * The compact JSON text of the labels built for the answer, from above (AnswerVolume, review
 * of 2026-10-07, X-EFFICIENCY-04): estimated before the objects are built, so that the work
 * time is read with them counted. Integers at their widest (20 digits), strings as
 * string_text_bytes.
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

// a string's text: quoted, every byte jsoncpp may escape (a control character, a quote, a
// backslash, a byte of a non-ASCII character: \uXXXX) counted as 6
uint64_t string_text_bytes(const char *begin, const char *end) {
    uint64_t bytes = 2;
    for (const char *c = begin; c != end; ++c) {
        const unsigned char u = static_cast<unsigned char>(*c);
        bytes += u < 0x20 || u >= 0x80 || u == '"' || u == '\\' ? 6 : 1;
    }
    return bytes;
}

uint64_t string_text_bytes(std::string_view s) {
    return string_text_bytes(s.data(), s.data() + s.size());
}

uint64_t decimal_digits(uint64_t x) {
    uint64_t digits = 1;
    while (x >= 10) {
        x /= 10;
        ++digits;
    }
    return digits;
}

constexpr uint64_t kNoKey = graph::traversal::npos;

Json::Value uint_json(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

Json::Value count_json(Relation relation, uint64_t value, Unit unit) {
    Json::Value v;
    v["value"] = relation == Relation::UNKNOWN ? Json::Value() : uint_json(value);
    v["relation"] = to_string(relation);
    v["unit"] = to_string(unit);
    return v;
}

Json::Value reason_json(const std::string &reason) {
    Json::Value r;
    r["reason"] = reason;
    return r;
}

const char* strand_of(Orientation o) { return strand_symbol(o); }

uint8_t strand_rank(Orientation o) {
    switch (o) {
        case Orientation::FORWARD: return 0;
        case Orientation::REVERSE: return 1;
        case Orientation::PALINDROMIC: return 2;
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

const char* to_string(RowStatus s) {
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
    LabelRecorder::NodeLabels labels;
    // placement (step 2): read, refused by the account, or not reached
    bool placed = false;
    bool place_refused = false;
    LabelQuery::NodeHits hits;
};

// one occurrence of a label in a context: (seq_id, 1-based start) with a record mapping,
// (kmer_coord, offset) without one; the strand of the context
struct Occurrence {
    uint64_t a = 0;
    uint64_t b = 0;
    uint8_t strand = 0;
    bool operator<(const Occurrence &o) const {
        return std::tie(a, b, strand) < std::tie(o.a, o.b, o.strand);
    }
    bool operator==(const Occurrence &o) const {
        return a == o.a && b == o.b && strand == o.strand;
    }
};

struct ContextLabel {
    LabelId label = 0;
    // placed: the occurrences of the label in this context (sorted, distinct); not placed
    // (placement not applicable, refused, not reached): none
    bool placed = false;
    std::vector<Occurrence> occurrences;
    // the occurrences of this context-label before the partial cut
    uint64_t total = 0;
};

} // namespace


double AnswerVolume::finalize_ms() const {
    // MB/s are bytes per microsecond: x 1000 bytes per ms. A rate of 0 or infinity: no model
    auto ms = [](double bytes, double mbps) {
        return mbps > 0 && std::isfinite(mbps) ? bytes / (mbps * 1000) : 0.0;
    };
    const double all = static_cast<double>(text_ + pending_) * scale_;
    const double pending = static_cast<double>(pending_) * scale_;
    return kAnswerVolumeMargin * (ms(all, build_mbps_) + ms(all, compress_mbps_)
                                  + ms(pending, build_mbps_));
}

uint64_t compact_json_bytes(const Json::Value &v) {
    switch (v.type()) {
        case Json::nullValue:
            return 4;
        case Json::booleanValue:
            return v.asBool() ? 4 : 5;
        case Json::intValue: {
            const int64_t x = v.asInt64();
            return x < 0 ? 1 + decimal_digits(static_cast<uint64_t>(-(x + 1)) + 1)
                         : decimal_digits(static_cast<uint64_t>(x));
        }
        case Json::uintValue:
            return decimal_digits(v.asUInt64());
        case Json::realValue:
            // at most 17 significant digits, a sign, a point and an exponent
            return 24;
        case Json::stringValue: {
            const char *begin = nullptr, *end = nullptr;
            v.getString(&begin, &end);
            return string_text_bytes(begin, end);
        }
        case Json::arrayValue: {
            uint64_t bytes = 2 + (v.size() ? v.size() - 1 : 0);
            for (const Json::Value &x : v) {
                bytes += compact_json_bytes(x);
            }
            return bytes;
        }
        case Json::objectValue: {
            // {"name":value,...}
            uint64_t bytes = 2 + (v.size() ? v.size() - 1 : 0);
            for (auto it = v.begin(); it != v.end(); ++it) {
                const char *end = nullptr;
                const char *begin = it.memberName(&end);
                bytes += string_text_bytes(begin, end) + 1 + compact_json_bytes(*it);
            }
            return bytes;
        }
    }
    return 0;
}


AnnotationDescription describe_annotation(const LabelOracle &oracle, GraphMode mode) {
    AnnotationDescription d;
    const bool basic = mode == GraphMode::BASIC;
    const bool coords = oracle.has_coordinates();
    const bool records = coords && oracle.coord_to_header();
    d.placement = !basic ? "none_canonical" : records ? "record" : coords ? "global" : "none";
    d.support = basic && records ? "record_verified" : "label_intersection";
    d.budgeted = oracle.decode_charged();
    return d;
}


struct PatternRetrieval::Impl {
    Impl(const AnnotatedDBG &anno_graph, const RetrievalLimits &limits, Budget &budget,
         const RetrievalHooks *hooks, AnswerVolume *volume)
          : oracle(anno_graph, hooks ? hooks->coord_to_header : nullptr),
            limits(limits), budget(budget),
            account(hooks && hooks->max_memory_bytes ? hooks->max_memory_bytes
                                                     : limits.max_memory_bytes),
            volume(volume) {
        if (hooks) {
            deny_decode = hooks->deny_decode;
            output_hook = hooks->output_hook;
        }
        // the reads under the deadline are decoded in paced chunks (as /traverse's)
        oracle.pacer().target_ms = limits.chunk_target_ms;
        if (hooks && hooks->read_hook)
            oracle.test_read_hook = hooks->read_hook;
        recorder = std::make_unique<LabelRecorder>(oracle, LabelKind::COLUMN,
                                                   std::max<uint64_t>(1, limits.max_labels_per_anchor));
        // no cache: every row is read once per step, its keys deduplicated, and nothing a
        // cache held would be in the account (a fixed allotment of zero)
        recorder->set_max_cache_bytes(0);
    }

    LabelOracle oracle;
    const RetrievalLimits limits;
    Budget &budget;
    Account account;
    std::unique_ptr<LabelRecorder> recorder;
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
    // the pattern being released (admit_context): its descriptors' bytes, held, and what they
    // may hold (all_or_count: what the account had left; partial: half of it); false once a
    // descriptor did not fit
    uint64_t descriptors = 0;
    uint64_t allowance = 0;
    bool admitting = false;

    uint64_t statement() const { return statement_bytes(oracle.get_k()); }

    // what a read may hold: the account's remainder, less the statements reserved for its rows
    DecodeBudget decode_budget(uint64_t reserve) const {
        DecodeBudget decode(account.left() - std::min(reserve, account.left()));
        if (deny_decode)
            decode.deny = deny_decode;
        return decode;
    }

    // the stop of the current pattern's annotation work (the first one), and whether the
    // reads themselves stopped (time, work, an unbudgeted read the account could not hold)
    std::optional<std::pair<std::string, std::string>> stop;
    bool time_stop = false;
    bool read_stop = false;

    ReadPacing pacing() {
        ReadPacing p;
        const Deadline *d = &budget.deadline();
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

    // What a refused read is charged (review GPT-2, finding 1): the units of what it decoded
    // (FetchRefusal::units: the row was read, and refused for its demand or its names), or 8
    // when its read itself did not fit (its units are not known): a refused row is work like
    // a row read, so that refusals cannot go on past the work budget
    static uint64_t refused_units(const FetchRefusal &r) { return r.units ? r.units : 8; }

    Json::Value refusal_json(const RowState &row, const std::vector<RetrievalContext> &contexts,
                             const char *phase, const FetchRefusal *r) const {
        Json::Value v;
        v["kmer"] = contexts[row.first].kmer;
        v["row"] = uint_json(AnnotatedDBG::graph_to_anno_index(row.key));
        v["phase"] = phase;
        v["reason"] = "max_memory";
        // what reading the row alone needed (at least), against what the account had left
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
    void discover(std::vector<RowState> &rows, const std::vector<RetrievalContext> &contexts,
                  Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                  uint64_t *held_lists, Json::Value *refused);

    // step 2: the coordinates of the rows' labels
    void place_rows(std::vector<RowState> &rows, const std::vector<RetrievalContext> &contexts,
                    Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                    uint64_t *held_hits, Json::Value *refused);
};

// Both steps read one row at a time (review GPT-2, finding 2): the time and the work are
// checked before every row, and every row whose read began is charged its units (also when it
// is refused or interrupted), so that a read passes the work budget by its one row at most —
// a batch sized by the rows read before it overshot by every wider row in it (7 rows of one
// label, then 8 of 100: 927 units under a budget of 150). The units of a row do not depend on
// how the rows are cut into reads (KeyCost: its whole row-diff path), so where the reads stop
// depends only on the index and the request.
void PatternRetrieval::Impl::discover(std::vector<RowState> &rows,
                                      const std::vector<RetrievalContext> &contexts,
                                      Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                                      uint64_t *held_lists, Json::Value *refused) {
    auto name_bytes = [](std::string_view name) { return label_name_bytes(name); };
    for (RowState &row : rows) {
        if (!may_read("label_discovery"))
            break;
        // room for the row's statement (refused or truncated), reserved before its read
        if (account.left() < statement()) {
            set_stop("label_discovery", "max_memory");
            break;
        }
        ReadPacing pace = pacing();
        if (budgeted) {
            DecodeBudget decode = decode_budget(statement());
            std::vector<LabelRecorder::NodeLabels> out;
            std::vector<KeyCost> costs;
            out.reserve(1);
            costs.reserve(1);
            size_t refused_at = 0;
            if (!recorder->fetch(&row.key, 1, decode, &out, &costs, &refused_at, name_bytes,
                                 &pace)) {
                if (pace.interrupted) {
                    // decoded work, though nothing was returned (none: the deadline came
                    // before the row's read)
                    charge_work(pace.units, pattern_units);
                    budget.check_time();
                    set_stop("label_discovery", "time");
                    time_stop = true;
                    break;
                }
                // the row does not fit what the account has left: refused, stated (the
                // statement in its reserve), and its read charged as work
                charge_work(refused_units(recorder->refusal()), pattern_units);
                row.status = RowStatus::REFUSED;
                const bool stated = account.charge(statement());
                assert(stated);
                (void)stated;
                refused->append(refusal_json(row, contexts, "label_discovery",
                                             &recorder->refusal()));
                if (mode == Mode::ALL_OR_COUNT)
                    break;      // the results cannot be published: no further reads
                continue;
            }
            // the list and the names the call gave, held within what was left
            const bool held = account.charge(decode.held());
            assert(held);
            (void)held;
            // the names stay with the dictionary (the request's); the list and the
            // provisional naming charges are this pattern's
            *held_lists += decode.held() - recorder->last_names_bytes();
            row.labels = std::move(out[0]);
            row.status = row.labels.truncated() ? RowStatus::TRUNCATED : RowStatus::COMPLETE;
            if (row.status == RowStatus::TRUNCATED) {
                // its anchors_truncated entry, in the read's reserve
                const bool stated = account.charge(statement());
                assert(stated);
                (void)stated;
            }
            charge_work(8 + row.labels.total + costs[0].dependency_units, pattern_units);
        } else {
            const size_t named = recorder->labels().size();
            std::vector<LabelRecorder::NodeLabels> out = recorder->fetch({ row.key }, &pace);
            if (pace.interrupted) {
                charge_work(pace.units, pattern_units);
                budget.check_time();
                set_stop("label_discovery", "time");
                time_stop = true;
                break;
            }
            // read without a budget: what it returned is held (charged after the fact), the
            // statement of a truncated row with it; the names stay with the dictionary whether
            // or not the rest fits (forced: the account may go past its maximum by them, and
            // then nothing more fits and the reads stop)
            uint64_t list = LabelRecorder::held_bytes(out[0]), names = 0;
            if (out[0].truncated())
                list += statement();
            for (size_t id = named; id < recorder->labels().size(); ++id) {
                names += label_name_bytes(recorder->labels()[id].name);
            }
            account.force(names);
            charge_work(8 + out[0].total, pattern_units);
            if (!account.charge(list)) {
                set_stop("label_discovery", "max_memory");
                break;
            }
            row.labels = std::move(out[0]);
            row.status = row.labels.truncated() ? RowStatus::TRUNCATED : RowStatus::COMPLETE;
            // the statement stays with the answer; the list is the pattern's
            *held_lists += list - (row.status == RowStatus::TRUNCATED ? statement() : 0);
        }
        ++*rows_read;
    }
    for (RowState &row : rows) {
        if (row.status == RowStatus::PENDING)
            row.status = RowStatus::NOT_READ;
    }
}

void PatternRetrieval::Impl::place_rows(std::vector<RowState> &rows,
                                        const std::vector<RetrievalContext> &contexts,
                                        Mode mode, uint64_t *rows_read, uint64_t *pattern_units,
                                        uint64_t *held_hits, Json::Value *refused) {
    // the rows read with at least one label
    std::vector<size_t> todo;
    for (size_t i = 0; i < rows.size(); ++i) {
        if ((rows[i].status == RowStatus::COMPLETE || rows[i].status == RowStatus::TRUNCATED)
                && !rows[i].labels.labels.empty()) {
            todo.push_back(i);
        }
    }
    if (todo.empty())
        return;
    // the permitted set: every label discovered so far (LabelQuery's ids are the
    // dictionary's, so that a hit names its label by the recorder's id)
    LabelQuery query(oracle, recorder->labels(), true);
    query.set_max_cache_bytes(0);
    // the units of a row's hits: 8, 1 per label and 1 per coordinate, and its dependency rows'
    auto hits_units = [](const LabelQuery::NodeHits &hits, const KeyCost &cost) {
        uint64_t u = 8 + hits.size() + cost.dependency_units;
        for (const auto &hit : hits) {
            u += hit.coords.size();
        }
        return u;
    };
    for (size_t t : todo) {
        RowState &row = rows[t];
        if (!may_read("placement"))
            break;
        // room for the row's statement (refused), reserved
        if (account.left() < statement()) {
            set_stop("placement", "max_memory");
            break;
        }
        ReadPacing pace = pacing();
        std::vector<LabelQuery::NodeHits> out;
        std::vector<KeyCost> costs;
        if (budgeted) {
            DecodeBudget decode = decode_budget(statement());
            out.reserve(1);
            costs.reserve(1);
            size_t refused_at = 0;
            if (!query.fetch(&row.key, 1, decode, &out, &costs, &refused_at, &pace)) {
                if (pace.interrupted) {
                    charge_work(pace.units, pattern_units);
                    budget.check_time();
                    set_stop("placement", "time");
                    time_stop = true;
                    break;
                }
                charge_work(refused_units(query.refusal()), pattern_units);
                row.place_refused = true;
                const bool stated = account.charge(statement());
                assert(stated);
                (void)stated;
                refused->append(refusal_json(row, contexts, "placement", &query.refusal()));
                if (mode == Mode::ALL_OR_COUNT)
                    break;
                continue;
            }
            const bool held = account.charge(decode.held());
            assert(held);
            (void)held;
            *held_hits += decode.held();
        } else {
            out = query.fetch({ row.key }, &pace);
            if (pace.interrupted) {
                charge_work(pace.units, pattern_units);
                budget.check_time();
                set_stop("placement", "time");
                time_stop = true;
                break;
            }
            costs.resize(1);
            const uint64_t bytes = LabelQuery::held_bytes(out[0]);
            if (!account.charge(bytes)) {
                charge_work(hits_units(out[0], costs[0]), pattern_units);
                set_stop("placement", "max_memory");
                break;
            }
            *held_hits += bytes;
        }
        charge_work(hits_units(out[0], costs[0]), pattern_units);
        row.hits = std::move(out[0]);
        row.placed = true;
        ++*rows_read;
    }
}


PatternRetrieval::PatternRetrieval(const AnnotatedDBG &anno_graph, GraphMode mode,
                                   const RetrievalLimits &limits, Budget &budget,
                                   const RetrievalHooks *hooks, AnswerVolume *volume)
      : impl_(std::make_unique<Impl>(anno_graph, limits, budget, hooks, volume)),
        limits_(limits) {
    description_ = describe_annotation(impl_->oracle, mode);
    impl_->budgeted = description_.budgeted;
    const std::string placement = description_.placement;
    impl_->place = limits.occurrences && (placement == "record" || placement == "global");
    impl_->records = impl_->place && placement == "record";
}

PatternRetrieval::~PatternRetrieval() {}

const char* PatternRetrieval::placement() const {
    return limits_.occurrences ? description_.placement : "not_requested";
}

void PatternRetrieval::begin_release(Mode mode) {
    Impl &m = *impl_;
    // a previous release whose descriptors no retrieve() took (none: the engine releases
    // nothing to a refused pattern) is not held further
    m.account.release(m.descriptors);
    m.descriptors = 0;
    // partial keeps half of what is left for the pattern's reads and labels, so that a memory
    // cut of the list still returns labelled contexts; all_or_count needs all of them anyway
    m.allowance = mode == Mode::PARTIAL ? m.account.left() / 2 : m.account.left();
    m.admitting = true;
}

bool PatternRetrieval::admit_context() {
    Impl &m = *impl_;
    if (!m.admitting)
        return false;
    const uint64_t bytes = context_bytes(m.oracle.get_k());
    if (m.descriptors + bytes > m.allowance || !m.account.charge(bytes)) {
        m.admitting = false;
        return false;
    }
    m.descriptors += bytes;
    return true;
}

uint64_t PatternRetrieval::memory_peak() const { return impl_->account.peak(); }
uint64_t PatternRetrieval::memory_held() const { return impl_->account.held(); }

LabelsAnswer PatternRetrieval::retrieve(const std::vector<RetrievalContext> &contexts,
                                        uint64_t released, size_t length, Mode mode,
                                        const Extraction &x, const Json::Value &graph_name) {
    Impl &m = *impl_;
    // the descriptors admit_context() charged for |contexts|, the pattern's from here on
    const uint64_t descriptors = m.descriptors;
    m.descriptors = 0;
    m.admitting = false;
    assert(contexts.size() <= released);
    assert(descriptors == contexts.size() * context_bytes(m.oracle.get_k()));
    m.stop.reset();
    m.time_stop = false;
    m.read_stop = false;
    const size_t k = m.oracle.get_k();
    const bool partial = mode == Mode::PARTIAL;

    LabelsAnswer a;
    a.fields["placement"] = placement();
    a.fields["annotation"] = description_.budgeted ? "budgeted" : "unbudgeted";
    a.fields["rows_refused"] = Json::Value(Json::arrayValue);
    a.fields["anchors_truncated"] = Json::Value(Json::arrayValue);
    a.fields["labels_cut"] = Json::Value();
    a.fields["occurrences_cut"] = Json::Value();
    a.fields["by_label"] = Json::Value();
    a.labels_count = count_json(Relation::UNKNOWN, 0, Unit::LABELS);
    a.occurrences_count = count_json(Relation::UNKNOWN, 0, Unit::PLACED_OCCURRENCES);
    if (!description_.budgeted)
        a.notes.push_back("annotation_unbudgeted");
    if (m.place && !m.records)
        a.notes.push_back("record_bounds_unknown");

    uint64_t rows_read = 0, pattern_units = 0;
    double discovery_ms = 0, placement_ms = 0;
    auto finish_work = [&]() {
        a.work["annotation_rows"] = uint_json(rows_read);
        a.work["annotation_units"] = uint_json(pattern_units);
        a.work["memory_bytes"] = uint_json(m.account.peak());
        a.timing["label_discovery_ms"] = discovery_ms;
        a.timing["placement_ms"] = placement_ms;
        // the statements stay in the answer, whatever else does (AnswerVolume)
        if (m.volume) {
            m.volume->add(compact_json_bytes(a.fields["rows_refused"])
                          + compact_json_bytes(a.fields["anchors_truncated"]));
        }
    };

    if (x.withheld) {
        // nothing was released: nothing is read (§5.2: no annotation work on a pattern
        // whose count was not admitted)
        m.account.release(descriptors);
        finish_work();
        return a;
    }

    // the descriptors of the released contexts were charged as the engine released them,
    // before their result objects were built (admit_context, §5.3 items 3 and 5): the first
    // that did not fit ended the list, and no result object was built after it
    const size_t keep = contexts.size();
    if (keep < released) {
        m.set_stop("output", "max_memory");
        if (!partial) {
            m.account.release(descriptors);
            a.withheld = "output_budget";
            a.stop = m.stop;
            finish_work();
            return a;
        }
        a.cut = "max_memory";
    }

    // the rows: the contexts' distinct keys, in answer order of their first context
    std::vector<RowState> rows;
    std::vector<size_t> row_of(keep);
    {
        std::unordered_map<uint64_t, size_t> index;
        for (size_t i = 0; i < keep; ++i) {
            if (contexts[i].key == kNoKey) {
                throw std::runtime_error("pattern: the k-mer " + contexts[i].kmer + " has no "
                                         "annotation row: the annotation does not describe "
                                         "this graph");
            }
            auto [it, inserted] = index.emplace(contexts[i].key, rows.size());
            if (inserted) {
                RowState row;
                row.key = contexts[i].key;
                row.first = i;
                rows.push_back(std::move(row));
            }
            row_of[i] = it->second;
        }
    }

    Json::Value &refused = a.fields["rows_refused"];
    uint64_t held_lists = 0, held_hits = 0;
    const auto t0 = std::chrono::steady_clock::now();
    m.discover(rows, contexts, mode, &rows_read, &pattern_units, &held_lists, &refused);
    const auto t1 = std::chrono::steady_clock::now();
    discovery_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    bool truncated = false;
    bool rows_complete = true;
    for (const RowState &row : rows) {
        if (row.status == RowStatus::TRUNCATED) {
            truncated = true;
            Json::Value t;
            t["kmer"] = contexts[row.first].kmer;
            t["row"] = uint_json(AnnotatedDBG::graph_to_anno_index(row.key));
            t["cap"] = uint_json(m.limits.max_labels_per_anchor);
            t["total"] = uint_json(row.labels.total);
            a.fields["anchors_truncated"].append(std::move(t));
        }
        rows_complete &= row.status == RowStatus::COMPLETE;
    }

    // all_or_count publishes nothing once a row is cut, refused or not read: no placement
    // (partial: on the rows read, unless a stop ended the reads)
    if (m.place && !m.read_stop && (partial || (rows_complete && !m.stop))) {
        const auto t2 = std::chrono::steady_clock::now();
        m.place_rows(rows, contexts, mode, &rows_read, &pattern_units, &held_hits, &refused);
        placement_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t2).count();
    }

    // the label order of §5.5 over the contexts read: contexts desc, column asc
    const std::vector<LabelRef> &dict = m.recorder->labels();
    std::vector<uint64_t> label_contexts(dict.size(), 0), label_suffix(dict.size(), 0);
    for (size_t i = 0; i < keep; ++i) {
        const RowState &row = rows[row_of[i]];
        if (row.status != RowStatus::COMPLETE && row.status != RowStatus::TRUNCATED)
            continue;
        for (LabelId id : row.labels.labels) {
            label_contexts[id]++;
            if (contexts[i].offset + length == k)
                label_suffix[id]++;
        }
    }
    std::vector<LabelId> order;
    for (LabelId id = 0; id < dict.size(); ++id) {
        if (label_contexts[id])
            order.push_back(id);
    }
    std::sort(order.begin(), order.end(), [&](LabelId x, LabelId y) {
        if (label_contexts[x] != label_contexts[y])
            return label_contexts[x] > label_contexts[y];
        return dict[x].name < dict[y].name;
    });
    const uint64_t num_labels = order.size();
    std::vector<uint64_t> rank(dict.size(), std::numeric_limits<uint64_t>::max());
    size_t kept_labels = order.size();
    if (partial && kept_labels > m.limits.max_labels) {
        kept_labels = m.limits.max_labels;
        Json::Value cut = reason_json("max_labels");
        cut["returned"] = uint_json(kept_labels);
        a.fields["labels_cut"] = std::move(cut);
    }
    for (size_t r = 0; r < kept_labels; ++r) {
        rank[order[r]] = r;
    }

    // the text the labels built for the answer will write (AnswerVolume, review of
    // 2026-10-07, X-EFFICIENCY-04), from above and before they are built: every result's label
    // fields and by_label first, then each context's labels as the loop below takes them; the
    // work time is read with them counted before each context's labels are built
    uint64_t pending = 0;
    auto pend = [&](uint64_t bytes) {
        pending += bytes;
        if (m.volume)
            m.volume->add_pending(bytes);
    };
    auto unpend = [&](uint64_t bytes) {
        pending -= std::min(bytes, pending);
        if (m.volume)
            m.volume->drop_pending(bytes);
    };
    std::vector<uint64_t> name_text(dict.size(), 0);
    auto label_text = [&](LabelId id) {
        if (!name_text[id])
            name_text[id] = string_text_bytes(dict[id].name);
        return name_text[id];
    };
    uint64_t summary_text = 0;
    {
        const uint64_t graph_text = compact_json_bytes(graph_name);
        for (size_t r = 0; r < kept_labels; ++r) {
            summary_text += kByLabelText + label_text(order[r]) + graph_text;
        }
        pend(keep * kResultLabelsText + summary_text);
    }

    // by_label (its entries with their copies of the names) is charged before any context's
    // labels: the answer holds it whatever is cut after it. When it does not fit, no label of
    // the pattern is built: stop {output, max_memory} (unless an earlier stop is stated),
    // all_or_count withholds (output_budget), partial returns the contexts read with
    // labels_status output_budget and by_label null
    uint64_t summary = 0;
    for (size_t r = 0; r < kept_labels; ++r) {
        summary += by_label_bytes(dict[order[r]].name);
    }
    const bool summary_held = m.account.charge(summary);
    if (!summary_held)
        unpend(summary_text);

    // the contexts' label lists with their occurrences; each label's deduplicated union
    // (§5.4: (column, seq_id, start, strand); the label is the column)
    std::vector<std::vector<ContextLabel>> lists(keep);
    std::vector<std::set<Occurrence>> unions(dict.size());
    std::vector<bool> output_cut(keep, !summary_held);
    uint64_t output = 0, dedup = 0;
    bool placement_complete = true;
    bool output_stopped = !summary_held;
    if (!summary_held)
        m.set_stop("output", "max_memory");
    for (size_t i = 0; i < keep && !output_stopped; ++i) {
        const RowState &row = rows[row_of[i]];
        if (row.status != RowStatus::COMPLETE && row.status != RowStatus::TRUNCATED)
            continue;
        const RetrievalContext &c = contexts[i];
        std::vector<ContextLabel> list;
        uint64_t bytes = 0, new_dedup = 0, text = 0;
        std::vector<std::pair<LabelId, Occurrence>> inserted;
        for (LabelId id : row.labels.labels) {
            // a label partial's max_labels cut is not listed, but its occurrences are
            // counted all the same (counts.occurrences is over every label)
            const bool listed = rank[id] != std::numeric_limits<uint64_t>::max();
            ContextLabel cl;
            cl.label = id;
            if (listed) {
                bytes += label_entry_bytes(dict[id].name);
                text += kLabelText + label_text(id);
            }
            if (m.place) {
                auto hit = std::lower_bound(row.hits.begin(), row.hits.end(), id,
                                            [](const LabelQuery::Hit &h, LabelId l) {
                                                return h.label < l;
                                            });
                if (row.placed && hit != row.hits.end() && hit->label == id) {
                    cl.placed = true;
                    const Column column = dict[id].column;
                    const uint8_t strand = strand_rank(c.orientation);
                    if (m.records) {
                        // the record first, the offset after (§4.3)
                        m.oracle.map_coords(column, hit->coords.data(), hit->coords.size(),
                                            [&](Coord, uint64_t seq_id, Coord local) {
                            const uint64_t nt_length
                                    = m.oracle.num_kmers_in_sequence(column, seq_id) + k - 1;
                            if (local + c.offset + length > nt_length) {
                                throw std::runtime_error(
                                        "pattern: a coordinate of " + c.kmer + " maps past the "
                                        "end of its record: the record mapping does not "
                                        "describe this annotation");
                            }
                            cl.occurrences.push_back(
                                    Occurrence { seq_id, local + c.offset + 1, strand });
                        });
                    } else {
                        for (Coord coord : hit->coords) {
                            cl.occurrences.push_back(Occurrence { coord, c.offset, strand });
                        }
                    }
                    std::sort(cl.occurrences.begin(), cl.occurrences.end());
                    cl.occurrences.erase(std::unique(cl.occurrences.begin(),
                                                     cl.occurrences.end()),
                                         cl.occurrences.end());
                    cl.total = cl.occurrences.size();
                    for (const Occurrence &o : cl.occurrences) {
                        if (listed && m.records) {
                            const std::string_view record = m.oracle.header_name(column, o.a);
                            bytes += occurrence_bytes(record);
                            text += kOccurrenceText + string_text_bytes(record);
                        } else if (listed) {
                            bytes += kGlobalOccurrenceBytes;
                            text += kGlobalOccurrenceText;
                        }
                        if (unions[id].insert(o).second) {
                            inserted.emplace_back(id, o);
                            new_dedup += kDedupBytes;
                        }
                    }
                } else {
                    // the row's placement refused or not reached, or no coordinates for a
                    // label the row carries: not placed, stated as unknown
                    placement_complete = false;
                }
            }
            if (listed)
                list.push_back(std::move(cl));
        }
        // the memory first (where it stops does not depend on the machine), then the time:
        // the work time, read with this context's labels counted in the answer's volume
        const bool held = m.account.charge(bytes + new_dedup);
        bool late = false;
        if (held) {
            if (m.output_hook)
                m.output_hook(i);
            pend(text);
            late = !m.budget.check_time();
            if (late) {
                unpend(text);
                m.account.release(bytes + new_dedup);
            }
        }
        if (!held || late) {
            for (const auto &[id, o] : inserted) {
                unions[id].erase(o);
            }
            output_stopped = true;
            for (size_t j = i; j < keep; ++j) {
                output_cut[j] = true;
            }
            // a time stop of the output is recorded in the budget too (Budget::check_time):
            // the later patterns answer as after any time stop
            m.set_stop("output", late ? "time" : "max_memory");
            m.time_stop |= late;
            break;
        }
        output += bytes;
        dedup += new_dedup;
        std::sort(list.begin(), list.end(), [&](const ContextLabel &x, const ContextLabel &y) {
            return rank[x.label] < rank[y.label];
        });
        lists[i] = std::move(list);
    }

    // partial: each label lists the first max_occurrences_per_label occurrences of its union
    uint64_t occurrences_total = 0;
    for (LabelId id : order) {
        occurrences_total += unions[id].size();
    }
    if (partial && m.place) {
        uint64_t cut_labels = 0;
        for (size_t r = 0; r < kept_labels; ++r) {
            const LabelId id = order[r];
            if (unions[id].size() <= m.limits.max_occurrences_per_label)
                continue;
            cut_labels++;
            auto last = unions[id].begin();
            std::advance(last, m.limits.max_occurrences_per_label);
            const std::set<Occurrence> kept(unions[id].begin(), last);
            for (auto &list : lists) {
                for (ContextLabel &cl : list) {
                    if (cl.label != id)
                        continue;
                    uint64_t dropped = 0;
                    const Column column = dict[id].column;
                    cl.occurrences.erase(std::remove_if(cl.occurrences.begin(),
                                                        cl.occurrences.end(),
                                                        [&](const Occurrence &o) {
                        if (kept.count(o))
                            return false;
                        dropped += m.records
                                ? occurrence_bytes(m.oracle.header_name(column, o.a))
                                : kGlobalOccurrenceBytes;
                        return true;
                    }), cl.occurrences.end());
                    m.account.release(dropped);
                    output -= std::min(output, dropped);
                }
            }
        }
        if (cut_labels) {
            Json::Value cut = reason_json("max_occurrences_per_label");
            cut["labels"] = uint_json(cut_labels);
            a.fields["occurrences_cut"] = std::move(cut);
        }
    }

    // what the pattern's reads held is freed; the dictionary, the descriptors and the labels
    // built for the answer stay (the buffered answer)
    m.account.release(held_lists + held_hits + dedup);

    // ---- the statements
    const bool all_returned = x.complete && keep == released;
    const bool anything_read = rows_read > 0;
    const bool contexts_exact = all_returned && rows_complete;
    bool any_placed = false;
    for (const RowState &row : rows) {
        any_placed |= row.placed;
    }
    for (const RowState &row : rows) {
        if ((row.status == RowStatus::COMPLETE || row.status == RowStatus::TRUNCATED)
                && !row.labels.labels.empty() && m.place && !row.placed) {
            placement_complete = false;
        }
    }
    const bool occurrences_exact = contexts_exact && (!m.place || placement_complete)
                                    && !output_stopped;
    const Relation labels_relation = contexts_exact ? Relation::EXACT
                                   : anything_read || keep ? Relation::AT_LEAST
                                                           : Relation::UNKNOWN;
    const Relation occ_relation = !m.records ? Relation::UNKNOWN
                                : occurrences_exact ? Relation::EXACT
                                : any_placed ? Relation::AT_LEAST : Relation::UNKNOWN;
    // an empty, complete release: nothing to read, the counts exact zeros
    const bool empty_complete = released == 0 && x.complete;

    a.complete = all_returned && rows_complete && (!m.place || placement_complete)
                    && !output_stopped && !m.stop && a.fields["labels_cut"].isNull()
                    && a.fields["occurrences_cut"].isNull();
    a.stop = m.stop;
    a.time_limited = m.time_stop;

    if (!partial && !a.complete) {
        // all_or_count: all or nothing (§5.2), the reason named; what explains it
        // (rows_refused, anchors_truncated) stays
        // the reads' own budgets (a refused row, the work, an unbudgeted read the account
        // could not hold) before the answer's memory, before a truncated anchor
        a.withheld = m.time_stop ? "deadline"
                   : !refused.empty() || m.read_stop ? "annotation_budget"
                   : output_stopped || m.stop ? "output_budget"
                   : truncated ? "anchor_labels_truncated"
                               : "annotation_budget";
        m.account.release(output + descriptors + (summary_held ? summary : 0));
        // no label of this pattern is built for the answer
        unpend(pending);
        finish_work();
        return a;
    }

    a.labels_count = count_json(empty_complete ? Relation::EXACT : labels_relation, num_labels,
                                Unit::LABELS);
    a.occurrences_count = count_json(empty_complete && m.records ? Relation::EXACT : occ_relation,
                                     occurrences_total, Unit::PLACED_OCCURRENCES);

    // by_label: the per-label summary over the returned contexts (§7.2), in label order (null
    // when the account could not hold it, above)
    Json::Value by_label(Json::arrayValue);
    for (size_t r = 0; r < kept_labels && summary_held; ++r) {
        const LabelId id = order[r];
        Json::Value b;
        b["graph"] = graph_name;
        b["column"] = dict[id].name;
        const Relation rel = contexts_exact ? Relation::EXACT : Relation::AT_LEAST;
        b["contexts"] = count_json(rel, label_contexts[id], Unit::GRAPH_CONTEXTS);
        b["contexts_suffix"] = count_json(rel, label_suffix[id], Unit::GRAPH_CONTEXTS);
        b["occurrences"] = count_json(!m.records ? Relation::UNKNOWN
                                      : occurrences_exact ? Relation::EXACT
                                      : any_placed ? Relation::AT_LEAST : Relation::UNKNOWN,
                                      unions[id].size(), Unit::PLACED_OCCURRENCES);
        by_label.append(std::move(b));
    }
    if (summary_held)
        a.fields["by_label"] = std::move(by_label);

    // the results' labels
    a.result_fields.reserve(keep);
    for (size_t i = 0; i < keep; ++i) {
        const RetrievalContext &c = contexts[i];
        const RowState &row = rows[row_of[i]];
        Json::Value f;
        f["support"] = "kmer";
        const bool read = row.status == RowStatus::COMPLETE || row.status == RowStatus::TRUNCATED;
        if (output_cut[i] && read) {
            f["labels_status"] = "output_budget";
        } else {
            f["labels_status"] = to_string(row.status);
        }
        f["labels_total"] = read ? uint_json(row.labels.total) : Json::Value();
        if (!read || output_cut[i]) {
            f["labels"] = Json::Value();
            a.result_fields.push_back(std::move(f));
            continue;
        }
        Json::Value labels(Json::arrayValue);
        for (const ContextLabel &cl : lists[i]) {
            Json::Value l;
            const Column column = dict[cl.label].column;
            l["column"] = dict[cl.label].name;
            l["support"] = "kmer";
            if (m.place) {
                if (!cl.placed) {
                    if (m.records)
                        l["occurrences"] = count_json(Relation::UNKNOWN, 0,
                                                      Unit::PLACED_OCCURRENCES);
                    l["occurrence_list"] = Json::Value();
                } else {
                    if (m.records)
                        l["occurrences"] = count_json(Relation::EXACT, cl.total,
                                                      Unit::PLACED_OCCURRENCES);
                    Json::Value occ(Json::arrayValue);
                    for (const Occurrence &o : cl.occurrences) {
                        Json::Value e;
                        if (m.records) {
                            e["seq_id"] = uint_json(o.a);
                            e["record"] = m.oracle.header_name(column, o.a);
                            e["strand"] = strand_of(c.orientation);
                            e["nt_coords"] = std::to_string(o.b) + "-"
                                                + std::to_string(o.b + length - 1);
                            e["nt_length"] = uint_json(
                                    m.oracle.num_kmers_in_sequence(column, o.a) + k - 1);
                        } else {
                            e["kmer_coord"] = uint_json(o.a);
                            e["offset"] = uint_json(o.b);
                            e["strand"] = strand_of(c.orientation);
                        }
                        occ.append(std::move(e));
                    }
                    l["occurrence_list"] = std::move(occ);
                }
            }
            labels.append(std::move(l));
        }
        f["labels"] = std::move(labels);
        a.result_fields.push_back(std::move(f));
    }
    // the labels are built: from now on they are written only
    if (m.volume)
        m.volume->settle(pending);
    finish_work();
    return a;
}


void apply_labels(Json::Value *entry, LabelsAnswer &&answer, Mode mode) {
    Json::Value &e = *entry;
    if (e.isMember("error"))
        return;
    for (const std::string &name : answer.fields.getMemberNames()) {
        e[name] = std::move(answer.fields[name]);
    }
    e["counts"]["labels"] = std::move(answer.labels_count);
    e["counts"]["occurrences"] = std::move(answer.occurrences_count);
    for (const std::string &name : answer.work.getMemberNames()) {
        e["work"][name] = std::move(answer.work[name]);
    }
    for (const std::string &name : answer.timing.getMemberNames()) {
        e["timing"][name] = std::move(answer.timing[name]);
    }
    for (const std::string &note : answer.notes) {
        e["notes"].append(note);
    }
    if (answer.time_limited)
        e["determinism"] = "time_limited";
    // the engine's stop first: the annotation's comes after a completed (or cut) release
    if (answer.stop && e["stop"].isNull()) {
        Json::Value s;
        s["phase"] = answer.stop->first;
        s["reason"] = answer.stop->second;
        e["stop"] = std::move(s);
    }
    if (e["withheld"].isObject()) {
        // the engine withheld the release: nothing was read
        return;
    }

    Json::Value &results = e["results"];
    if (answer.withheld) {
        // all_or_count: all or nothing (§5.2)
        e["withheld"] = reason_json(*answer.withheld);
        e["retrieval_complete"] = false;
        e["returned"] = 0;
        e["cut"] = Json::Value();
        results = Json::Value(Json::arrayValue);
        e["by_label"] = Json::Value();
        return;
    }
    // the memory cut is the one that set the list's length (the engine's max_contexts cut,
    // when there was one, released more); the route built the admitted contexts only
    if (answer.cut && mode == Mode::PARTIAL)
        e["cut"] = reason_json(*answer.cut);
    for (Json::ArrayIndex i = 0; i < results.size() && i < answer.result_fields.size(); ++i) {
        for (const std::string &name : answer.result_fields[i].getMemberNames()) {
            results[i][name] = std::move(answer.result_fields[i][name]);
        }
    }
    // complete only when the release was and every label of every context came with it
    e["retrieval_complete"] = e["retrieval_complete"].asBool() && answer.complete;
}

} // namespace cli
} // namespace mtg
