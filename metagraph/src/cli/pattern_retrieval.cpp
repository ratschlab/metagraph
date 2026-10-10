#include "pattern_retrieval.hpp"
#include "pattern_retrieval_impl.hpp"

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
using graph::traversal::kKeyUnits;
using graph::traversal::hits_units;
using graph::traversal::Column;
using graph::traversal::Coord;
using annot::matrix::DecodeBudget;
// the memory and text models, the account and the rows (pattern_retrieval_impl.hpp)
using namespace retrieval;

namespace {

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
    // a path's label (long_search "paths"): record_verified, one record holds the whole path
    bool verified = false;
};

// What the verification of a path keeps for its output, so that the chains of every (path,
// label) are not made again for it: its occurrences as runs, |count| consecutive starts in one
// record from |first| on (global placement: consecutive chains, the coordinate |first|.a on) —
// a homopolymer's are one run. Charged before they are held: the entry of a label with
// occurrences (PathRuns) and each run, both below the deduplication state the output then
// charges for the occurrences (kDedupBytes each, a run holding one at least), so that keeping
// them does not raise the peak of the account
struct OccurrenceRun {
    Occurrence first;
    uint64_t count = 0;
};
struct PathRuns {
    LabelId label = 0;
    uint64_t total = 0;
    std::vector<OccurrenceRun> runs;
};
constexpr uint64_t kRunsBytes = 32;
constexpr uint64_t kRunBytes = 32;

/**
 * |o| into |set|, the deduplicated union of a label (§5.4), the occurrences coming in sorted
 * order: hinted at the position after the previous one's (|hint|, advanced), which makes a
 * sorted stream O(1) an occurrence where the union already holds its neighbours (a
 * homopolymer's contexts: every context's occurrences but its last are in the union). The
 * iterator of the new element, or set.end() when the union held it.
 */
std::set<Occurrence>::iterator insert_sorted(std::set<Occurrence> &set,
                                             std::set<Occurrence>::iterator *hint,
                                             const Occurrence &o) {
    const size_t before = set.size();
    auto it = set.insert(*hint, o);
    *hint = std::next(it);
    return set.size() > before ? it : set.end();
}

/**
 * One k-mer's coordinates of a label for the chains of a path (§4.3): the column coordinates
 * of the j-th k-mer (sorted, distinct), shifted by -j, and where the join stands in them.
 */
struct ShiftedList {
    const Coord *coords = nullptr;
    size_t size = 0;
    Coord shift = 0;
    size_t pos = 0;
};

// the first position at or after l.pos whose coordinate is at least |target|: galloping from
// l.pos (the join's targets only grow), then a binary search in the last stride
size_t seek(const ShiftedList &l, Coord target) {
    size_t lo = l.pos;
    if (lo >= l.size || l.coords[lo] >= target)
        return lo;
    size_t step = 1;
    while (lo + step < l.size && l.coords[lo + step] < target) {
        lo += step;
        step *= 2;
    }
    const size_t hi = std::min(lo + step, l.size);
    return std::lower_bound(l.coords + lo + 1, l.coords + hi, target) - l.coords;
}

// the largest d <= |limit| with coords[pos + d] == coords[pos] + d: how far the coordinates
// from l.pos on are consecutive (sorted and distinct, coords[pos + d] - coords[pos] >= d, equal
// exactly as long as they are: a binary search). Distinct: each is the position of one k-mer
// in its column (a repeated one would make the test unsound, a missing chain hidden by it)
uint64_t consecutive(const ShiftedList &l, uint64_t limit) {
    const Coord *c = l.coords + l.pos;
    uint64_t hi = std::min<uint64_t>(limit, l.size - 1 - l.pos);
    if (!hi || c[1] != c[0] + 1)
        return 0;
    if (c[hi] - c[0] == hi)
        return hi;
    // c[lo] consecutive, c[hi] not
    uint64_t lo = 1;
    while (hi - lo > 1) {
        const uint64_t mid = lo + (hi - lo) / 2;
        (c[mid] - c[0] == mid ? lo : hi) = mid;
    }
    return lo;
}

/**
 * The chains of a path in one label (§4.3): the column coordinates c with c + j a coordinate of
 * the j-th k-mer of the path for every j, the intersection of |lists| shifted. A leapfrog join,
 * driven by the smallest list (binary-searching each coordinate of the first list in every
 * other one would take 44 million searches for a 30,000-base homopolymer and a 1,500-base
 * path): its next coordinate is the candidate, the other lists (smallest first) seek it from
 * where they stood (galloping), the first holding a larger coordinate moves the candidate
 * there, and a candidate every list holds is extended into the run of consecutive chains as far
 * as every list is consecutive from it (a homopolymer's chains: one run, one check per list).
 * emit(first, last) takes the runs in order and may end the join (false). |work|(units) is
 * asked before every seek and every extension, and false ends the join (returned: false).
 */
template <class Emit, class Work>
bool join_chains(std::vector<ShiftedList> &lists, const Emit &emit, const Work &work) {
    if (lists.empty())
        return true;
    if (!work(lists.size()))
        return false;
    std::sort(lists.begin(), lists.end(), [](const ShiftedList &x, const ShiftedList &y) {
        return x.size < y.size;
    });
    ShiftedList &driver = lists[0];
    // a chain is not negative: the driver's coordinates from its shift on
    driver.pos = std::lower_bound(driver.coords, driver.coords + driver.size, driver.shift)
                    - driver.coords;
    while (driver.pos < driver.size) {
        const Coord c = driver.coords[driver.pos] - driver.shift;
        bool held = true;
        for (size_t t = 1; t < lists.size(); ++t) {
            ShiftedList &l = lists[t];
            if (!work(1))
                return false;
            l.pos = seek(l, c + l.shift);
            if (l.pos == l.size)
                return true;
            if (l.coords[l.pos] != c + l.shift) {
                // the next candidate: no chain before this list's next coordinate
                if (!work(1))
                    return false;
                driver.pos = seek(driver, l.coords[l.pos] - l.shift + driver.shift);
                held = false;
                break;
            }
        }
        if (!held)
            continue;
        uint64_t run = consecutive(driver, std::numeric_limits<uint64_t>::max());
        for (size_t t = 1; t < lists.size() && run; ++t) {
            if (!work(1))
                return false;
            run = consecutive(lists[t], run);
        }
        if (!emit(c, c + run))
            return true;
        // past the run (and a coordinate repeated in it)
        driver.pos = seek(driver, c + run + 1 + driver.shift);
    }
    return true;
}

/**
 * Partial (§5.4): each label lists the first |cap| occurrences of its union, those up to its
 * cap-th. |lists| (each context's or path's listed labels, each with at most its first |cap|
 * occurrences) are cut to them, and what the cut takes is given back by refund(bytes of the
 * memory model, bytes of answer text), for the labels whose union holds more than |cap| among
 * the first |kept_labels| of |order| (the labels listed), which |cut| states. Not clocked: it
 * runs after the work, over occurrences the account admitted (kGlobalOccurrenceBytes or more
 * each) and whose text the answer's volume still holds pending, so the finalisation reserve
 * already counts the time to write them, more than dropping them takes.
 */
template <class Refund>
void trim_to_unions(const std::vector<std::set<Occurrence>> &unions,
                    const std::vector<LabelId> &order, size_t kept_labels, uint64_t cap,
                    bool records, const LabelOracle &oracle,
                    const std::vector<LabelRef> &dict,
                    std::vector<std::vector<ContextLabel>> *lists, Json::Value *cut,
                    const Refund &refund) {
    // the last occurrence each label lists (null: all of them)
    std::vector<const Occurrence*> bound(dict.size(), nullptr);
    uint64_t cut_labels = 0;
    for (size_t r = 0; r < kept_labels; ++r) {
        const std::set<Occurrence> &u = unions[order[r]];
        if (u.size() <= cap)
            continue;
        cut_labels++;
        // (a cap of 0: the lists hold none)
        if (cap)
            bound[order[r]] = &*std::next(u.begin(), cap - 1);
    }
    if (!cut_labels)
        return;
    uint64_t bytes = 0, text = 0;
    for (auto &list : *lists) {
        for (ContextLabel &cl : list) {
            if (!bound[cl.label])
                continue;
            auto end = std::upper_bound(cl.occurrences.begin(), cl.occurrences.end(),
                                        *bound[cl.label]);
            for (auto it = end; it != cl.occurrences.end(); ++it) {
                if (records) {
                    const std::string_view record = oracle.header_name(dict[cl.label].column,
                                                                       it->a);
                    bytes += occurrence_bytes(record);
                    text += kOccurrenceText + string_text_bytes(record);
                } else {
                    bytes += kGlobalOccurrenceBytes;
                    text += kGlobalOccurrenceText;
                }
            }
            cl.occurrences.erase(end, cl.occurrences.end());
        }
    }
    refund(bytes, text);
    *cut = reason_json("max_occurrences_per_label");
    (*cut)["labels"] = uint_json(cut_labels);
}

// What the verification made of a label carrying a path (every k-mer of it annotated)
enum class PathSupport : uint8_t {
    // not verified: no record placement, or a row of the path not placed (refused, not
    // reached), or no coordinates of the label in a placed row
    NOT_PLACED,
    // the coordinates were read: no contiguous occurrence of the whole path in one record
    // (record placement), or record bounds unknown (global: the chains only)
    PLACED,
    // record placement: one record of the label holds the whole path (record_verified)
    VERIFIED,
};

// a label carrying a path, with what the verification made of it (8 bytes in the model)
struct PathLabel {
    LabelId label = 0;
    PathSupport support = PathSupport::NOT_PLACED;
};
// the memory model of a path's labels (the labels on every k-mer of it, held until the
// pattern's labels are built): its list and 8 bytes per label
constexpr uint64_t kPathLabelsBytes = 32;
constexpr uint64_t kPathLabelBytes = 8;

// an item's occurrences inserted into the unions, undone when its output is refused
using Inserted = std::vector<std::pair<LabelId, std::set<Occurrence>::iterator>>;

// The label fields of a pattern's answer before anything is read: the counts unknown, no
// statement and no cut, and the notes on the annotation's access and placement
LabelsAnswer begin_labels(const char *placement, bool budgeted, bool place, bool records) {
    LabelsAnswer a;
    a.fields["placement"] = placement;
    a.fields["annotation"] = budgeted ? "budgeted" : "unbudgeted";
    a.fields["rows_refused"] = Json::Value(Json::arrayValue);
    a.fields["anchors_truncated"] = Json::Value(Json::arrayValue);
    a.fields["labels_cut"] = Json::Value();
    a.fields["occurrences_cut"] = Json::Value();
    a.fields["by_label"] = Json::Value();
    a.labels_count = count_json(Relation::UNKNOWN, 0, Unit::LABELS);
    a.occurrences_count = count_json(Relation::UNKNOWN, 0, Unit::PLACED_OCCURRENCES);
    if (!budgeted)
        a.notes.push_back("annotation_unbudgeted");
    if (place && !records)
        a.notes.push_back("record_bounds_unknown");
    return a;
}

/**
 * The label order of §5.5 over the labels the items carry (|counts|: per label, the contexts
 * or the paths): items desc, column asc; and the labels the results list: all of them, or
 * partial's first max_labels (the cut stated in |labels_cut|). rank[id] is the label's place
 * among the listed ones, uint64_t's maximum when it is not listed.
 */
struct LabelOrder {
    std::vector<LabelId> order;
    std::vector<uint64_t> rank;
    size_t listed = 0;
};

LabelOrder order_labels(const std::vector<uint64_t> &counts, const std::vector<LabelRef> &dict,
                        bool partial, uint64_t max_labels, Json::Value *labels_cut) {
    LabelOrder o;
    for (LabelId id = 0; id < dict.size(); ++id) {
        if (counts[id])
            o.order.push_back(id);
    }
    std::sort(o.order.begin(), o.order.end(), [&](LabelId x, LabelId y) {
        if (counts[x] != counts[y])
            return counts[x] > counts[y];
        return dict[x].name < dict[y].name;
    });
    o.rank.assign(dict.size(), std::numeric_limits<uint64_t>::max());
    o.listed = o.order.size();
    if (partial && o.listed > max_labels) {
        o.listed = max_labels;
        Json::Value cut = reason_json("max_labels");
        cut["returned"] = uint_json(o.listed);
        *labels_cut = std::move(cut);
    }
    for (size_t r = 0; r < o.listed; ++r) {
        o.rank[o.order[r]] = r;
    }
    return o;
}

// The text of each label's name in the answer (string_text_bytes), computed once per label
class NameText {
  public:
    explicit NameText(const std::vector<LabelRef> &dict) : dict_(dict), bytes_(dict.size(), 0) {}

    uint64_t operator()(LabelId id) {
        if (!bytes_[id])
            bytes_[id] = string_text_bytes(dict_[id].name);
        return bytes_[id];
    }

  private:
    const std::vector<LabelRef> &dict_;
    std::vector<uint64_t> bytes_;
};

// by_label in the account: an entry with its copy of the name per listed label
uint64_t summary_bytes(const LabelOrder &labels, const std::vector<LabelRef> &dict) {
    uint64_t bytes = 0;
    for (size_t r = 0; r < labels.listed; ++r) {
        bytes += by_label_bytes(dict[labels.order[r]].name);
    }
    return bytes;
}

// by_label's text: per listed label |entry_text| (an entry without its names), the label's
// name and the graph's
uint64_t summary_text(const LabelOrder &labels, NameText &names, const Json::Value &graph_name,
                      uint64_t entry_text) {
    const uint64_t graph_text = compact_json_bytes(graph_name);
    uint64_t text = 0;
    for (size_t r = 0; r < labels.listed; ++r) {
        text += entry_text + names(labels.order[r]) + graph_text;
    }
    return text;
}

// a refused item's occurrences out of the unions again, it and every later item cut
void undo_output(std::vector<std::set<Occurrence>> *unions, const Inserted &inserted,
                 std::vector<bool> *output_cut, size_t i) {
    for (const auto &[id, it] : inserted) {
        (*unions)[id].erase(it);
    }
    for (size_t j = i; j < output_cut->size(); ++j) {
        (*output_cut)[j] = true;
    }
}

// One label object of a result (§7.2): its column (|name|, |label|'s as NameCheck gives it),
// its |support| and, with placement, its occurrences (their count with record placement, and
// their list; null when not placed)
Json::Value label_json(const ContextLabel &cl, const LabelRef &label, const std::string &name,
                       const char *support, const char *strand, size_t length,
                       const LabelOracle &oracle, bool place, bool records) {
    Json::Value l;
    l["column"] = name;
    l["support"] = support;
    if (!place)
        return l;
    if (!cl.placed) {
        if (records)
            l["occurrences"] = count_json(Relation::UNKNOWN, 0, Unit::PLACED_OCCURRENCES);
        l["occurrence_list"] = Json::Value();
        return l;
    }
    if (records)
        l["occurrences"] = count_json(Relation::EXACT, cl.total, Unit::PLACED_OCCURRENCES);
    const size_t k = oracle.get_k();
    Json::Value occ(Json::arrayValue);
    for (const Occurrence &o : cl.occurrences) {
        Json::Value e;
        if (records) {
            e["seq_id"] = uint_json(o.a);
            e["record"] = oracle.header_name(label.column, o.a);
            e["strand"] = strand;
            e["nt_coords"] = std::to_string(o.b) + "-" + std::to_string(o.b + length - 1);
            e["nt_length"] = uint_json(oracle.num_kmers_in_sequence(label.column, o.a) + k - 1);
        } else {
            e["kmer_coord"] = uint_json(o.a);
            e["offset"] = uint_json(o.b);
            e["strand"] = strand;
        }
        occ.append(std::move(e));
    }
    l["occurrence_list"] = std::move(occ);
    return l;
}

} // namespace


uint64_t path_descriptor_bytes(size_t k, size_t length) {
    const uint64_t n = length >= k ? length - k + 1 : 1;
    return context_bytes(k) + 3 * static_cast<uint64_t>(length) + kPathKmerBytes * n;
}

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


void PatternRetrieval::Impl::finish_work(LabelsAnswer &a, uint64_t rows_read,
                                         uint64_t pattern_units, double discovery_ms,
                                         double placement_ms) {
    a.work["annotation_rows"] = uint_json(rows_read);
    a.work["annotation_units"] = uint_json(pattern_units);
    a.work["memory_bytes"] = uint_json(account.peak());
    a.timing["label_discovery_ms"] = discovery_ms;
    a.timing["placement_ms"] = placement_ms;
    if (volume) {
        volume->add(compact_json_bytes(a.fields["rows_refused"])
                    + compact_json_bytes(a.fields["anchors_truncated"]));
    }
    counters += a.counters;
}

bool PatternRetrieval::Impl::skip_reads(LabelsAnswer &a, const Extraction &x, size_t keep,
                                        uint64_t released, uint64_t descriptors, bool partial) {
    if (x.withheld) {
        // nothing was released: nothing is read (§5.2: no annotation work on a pattern whose
        // count was not admitted)
        account.release(descriptors);
        finish_work(a, 0, 0, 0, 0);
        return true;
    }
    if (keep < released) {
        set_stop("output", "max_memory");
        if (!partial) {
            account.release(descriptors);
            a.withheld = "output_budget";
            a.stop = stop;
            finish_work(a, 0, 0, 0, 0);
            return true;
        }
        a.cut = "max_memory";
    }
    return false;
}

bool PatternRetrieval::Impl::state_rows(LabelsAnswer &a, const std::vector<RowState> &rows,
                                        const std::vector<RetrievalContext> &named,
                                        bool *truncated) {
    bool complete = true;
    for (const RowState &row : rows) {
        a.counters.rows_distinct += row.status == RowStatus::COMPLETE
                                    || row.status == RowStatus::TRUNCATED;
        if (row.status == RowStatus::TRUNCATED) {
            *truncated = true;
            Json::Value t;
            t["kmer"] = named[row.first].kmer;
            t["row"] = uint_json(AnnotatedDBG::graph_to_anno_index(row.key));
            t["cap"] = uint_json(limits.max_labels_per_anchor);
            t["total"] = uint_json(row.labels.total);
            a.fields["anchors_truncated"].append(std::move(t));
        }
        complete &= row.status == RowStatus::COMPLETE;
    }
    return complete;
}

bool PatternRetrieval::Impl::commit_output(size_t i, uint64_t bytes, uint64_t text, bool late,
                                           bool refused, PendingText *pending) {
    const bool held = !late && !refused && account.charge(bytes);
    if (held) {
        if (output_hook)
            output_hook(i);
        pending->add(text);
        late = !budget.check_time();
        if (late) {
            pending->drop(text);
            account.release(bytes);
        }
    }
    if (held && !late)
        return true;
    set_stop("output", late ? "time" : "max_memory");
    time_stop |= late;
    return false;
}

void PatternRetrieval::Impl::withhold_all(LabelsAnswer &a, bool rows_refused,
                                          bool output_stopped, bool truncated, uint64_t held,
                                          PendingText *pending) {
    a.withheld = time_stop ? "deadline"
               : rows_refused || read_stop ? "annotation_budget"
               : output_stopped || stop ? "output_budget"
               : truncated ? "anchor_labels_truncated"
                           : "annotation_budget";
    account.release(held);
    pending->drop_all();
}


// Both steps read one row at a time: the time and the work are checked before every row, and
// every row whose read began is charged its units (also when it is refused or interrupted), so
// that a read passes the work budget by its one row at most — a batch sized by the rows read
// before it would overshoot by every wider row in it. The units of a row do not depend on how
// the rows are cut into reads (KeyCost: its whole row-diff path), so where the reads stop
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
            charge_work(kKeyUnits + row.labels.total + costs[0].dependency_units, pattern_units);
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
            charge_work(kKeyUnits + out[0].total, pattern_units);
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
                                        uint64_t *held_hits, Json::Value *refused,
                                        const std::vector<bool> *needed,
                                        const std::vector<LabelRef> *dict) {
    // the rows read with at least one label (and, for paths, needed)
    std::vector<size_t> todo;
    for (size_t i = 0; i < rows.size(); ++i) {
        if ((rows[i].status == RowStatus::COMPLETE || rows[i].status == RowStatus::TRUNCATED)
                && !rows[i].labels.labels.empty() && (!needed || (*needed)[i])) {
            todo.push_back(i);
        }
    }
    if (todo.empty())
        return;
    // the permitted set: every label discovered so far (LabelQuery's ids are the
    // dictionary's, so that a hit names its label by the recorder's id), or the labels
    // retrieve_given's rows name
    LabelQuery query(oracle, dict ? *dict : recorder->labels(), true);
    query.set_max_cache_bytes(0);
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
                charge_work(hits_units(out[0]) + costs[0].dependency_units, pattern_units);
                set_stop("placement", "max_memory");
                break;
            }
            *held_hits += bytes;
        }
        charge_work(hits_units(out[0]) + costs[0].dependency_units, pattern_units);
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
    impl_->mode = mode;
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

bool PatternRetrieval::admit(uint64_t bytes) {
    Impl &m = *impl_;
    if (!m.admitting)
        return false;
    if (m.descriptors + bytes > m.allowance || !m.account.charge(bytes)) {
        m.admitting = false;
        return false;
    }
    m.descriptors += bytes;
    return true;
}

bool PatternRetrieval::admit_context() {
    return admit(context_bytes(impl_->oracle.get_k()));
}

bool PatternRetrieval::admit_path(size_t length) {
    return admit(path_descriptor_bytes(impl_->oracle.get_k(), length));
}

uint64_t PatternRetrieval::memory_peak() const { return impl_->account.peak(); }
const RetrievalCounters& PatternRetrieval::counters() const { return impl_->counters; }

LabelsAnswer PatternRetrieval::retrieve(const std::vector<RetrievalContext> &contexts,
                                        uint64_t released, size_t length, Mode mode,
                                        const Extraction &x, const Json::Value &graph_name) {
    return retrieve_rows(contexts, released, length, mode, x, graph_name, nullptr);
}

LabelsAnswer PatternRetrieval::retrieve_rows(const std::vector<RetrievalContext> &contexts,
                                             uint64_t released, size_t length, Mode mode,
                                             const Extraction &x, const Json::Value &graph_name,
                                             Given *given) {
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

    LabelsAnswer a = begin_labels(placement(), description_.budgeted, m.place, m.records);
    uint64_t rows_read = 0, pattern_units = 0;
    double discovery_ms = 0, placement_ms = 0;

    // the contexts' descriptors were charged by admit_context (§5.3 items 3 and 5)
    const size_t keep = contexts.size();
    if (m.skip_reads(a, x, keep, released, descriptors, partial))
        return a;

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
    if (!given) {
        m.discover(rows, contexts, mode, &rows_read, &pattern_units, &held_lists, &refused);
    } else {
        // retrieve_given: the labels the selection pass read on each row (the predicate's,
        // never truncated), renamed into the given dictionary; the lists are smaller than the
        // kept rows the account holds for them until retrieve_given ends. Light work, clocked
        // as the retrieval's own: a time stop leaves the later rows not read, as discovery's
        for (RowState &row : rows) {
            auto it = given->rows->find(row.key);
            if (it == given->rows->end())
                throw std::logic_error("pattern: retrieve_given for a row the selection did "
                                       "not keep");
            const LabelQuery::NodeHits &hits = it->second.hits;
            if (!m.may_work(1 + hits.size())) {
                m.set_stop("label_discovery", "time");
                m.time_stop = true;
                break;
            }
            row.labels.labels.reserve(hits.size());
            for (const LabelQuery::Hit &h : hits) {
                auto id = std::lower_bound(given->ids->begin(), given->ids->end(), h.label);
                if (id == given->ids->end() || *id != h.label)
                    throw std::logic_error("pattern: a kept row's label is not in the given "
                                           "dictionary");
                row.labels.labels.push_back(static_cast<LabelId>(id - given->ids->begin()));
            }
            row.labels.total = row.labels.labels.size();
            row.status = RowStatus::COMPLETE;
        }
        for (RowState &row : rows) {
            if (row.status == RowStatus::PENDING)
                row.status = RowStatus::NOT_READ;
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    discovery_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    bool truncated = false;
    const bool rows_complete = m.state_rows(a, rows, contexts, &truncated);

    // all_or_count publishes nothing once a row is cut, refused or not read: no placement
    // (partial: on the rows read, unless a stop ended the reads)
    if (m.place && !m.read_stop && (partial || (rows_complete && !m.stop))) {
        const auto t2 = std::chrono::steady_clock::now();
        m.place_rows(rows, contexts, mode, &rows_read, &pattern_units, &held_hits, &refused,
                     nullptr, given ? &given->dict : nullptr);
        placement_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t2).count();
    }

    // the label order of §5.5 over the contexts read: contexts desc, column asc
    const std::vector<LabelRef> &dict = given ? given->dict : m.recorder->labels();
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
    const LabelOrder labels = order_labels(label_contexts, dict, partial, m.limits.max_labels,
                                           &a.fields["labels_cut"]);
    const std::vector<LabelId> &order = labels.order;
    const std::vector<uint64_t> &rank = labels.rank;
    const size_t kept_labels = labels.listed;
    const uint64_t num_labels = order.size();

    // the text the labels built for the answer will write (AnswerVolume), from above and before
    // they are built: every result's label fields and by_label first, then each context's
    // labels as the loop below takes them; the work time is read with them counted before each
    // context's labels are built
    PendingText pending(m.volume);
    NameText label_text(dict);
    const uint64_t by_label_text = summary_text(labels, label_text, graph_name, kByLabelText);
    pending.add(keep * kResultLabelsText + by_label_text);

    // by_label (its entries with their copies of the names) is charged before any context's
    // labels: the answer holds it whatever is cut after it. When it does not fit, no label of
    // the pattern is built: stop {output, max_memory} (unless an earlier stop is stated),
    // all_or_count withholds (output_budget), partial returns the contexts read with
    // labels_status output_budget and by_label null
    const uint64_t summary = summary_bytes(labels, dict);
    const bool summary_held = m.account.charge(summary);
    if (!summary_held)
        pending.drop(by_label_text);

    // the contexts' label lists with their occurrences; each label's deduplicated union (§5.4:
    // (column, seq_id, start, strand); the label is the column). Every occurrence is counted in
    // its label's union, but a listed label's list holds only the first
    // max_occurrences_per_label of its context's occurrences in partial, so that the ones past
    // the cap are never built, charged or estimated in the answer's volume: the label lists the
    // first that many of its union, and those of one context are a prefix of its sorted
    // occurrences. Once the unions are complete, the lists are cut to them, and what they no
    // longer hold refunded to the account and the volume
    std::vector<std::vector<ContextLabel>> lists(keep);
    std::vector<std::set<Occurrence>> unions(dict.size());
    std::vector<bool> output_cut(keep, !summary_held);
    uint64_t output = 0, dedup = 0;
    bool placement_complete = true;
    bool output_stopped = !summary_held;
    if (!summary_held)
        m.set_stop("output", "max_memory");
    const uint64_t cap = partial ? m.limits.max_occurrences_per_label
                                 : std::numeric_limits<uint64_t>::max();
    for (size_t i = 0; i < keep && !output_stopped; ++i) {
        const RowState &row = rows[row_of[i]];
        if (row.status != RowStatus::COMPLETE && row.status != RowStatus::TRUNCATED)
            continue;
        const RetrievalContext &c = contexts[i];
        const uint8_t strand = strand_rank(c.orientation);
        if (m.occurrences_hook)
            m.occurrences_hook(i);
        std::vector<ContextLabel> list;
        uint64_t bytes = 0, new_dedup = 0, text = 0;
        Inserted inserted;
        bool late = false, refused = false;
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
                    // one occurrence per coordinate (distinct), in their order: the record
                    // mapping keeps the order of a column's coordinates and the offset is the
                    // context's (without a mapping the coordinate is the occurrence)
                    const auto &coords = hit->coords;
                    cl.total = coords.size();
                    const uint64_t first = listed ? std::min<uint64_t>(cap, cl.total) : 0;
                    cl.occurrences.reserve(first);
                    auto hint = unions[id].begin();
                    Occurrence previous;
                    bool any = false;
                    // false: what the context holds no longer fits what the account has left
                    // (the context's charge below would be refused: it stops here, before the
                    // union grows past the account)
                    auto add = [&](const Occurrence &o) {
                        // (the list keeps the first: their order is checked; a coordinate
                        // repeated is one occurrence)
                        if (any && !(previous < o)) {
                            if (o == previous) {
                                cl.total--;
                                return true;
                            }
                            throw std::runtime_error("pattern: the coordinates of " + c.kmer
                                                     + " are not in their record order: the "
                                                     "record mapping does not describe this "
                                                     "annotation");
                        }
                        previous = o;
                        any = true;
                        if (cl.occurrences.size() < first) {
                            cl.occurrences.push_back(o);
                            if (m.records) {
                                const std::string_view record = m.oracle.header_name(column,
                                                                                      o.a);
                                bytes += occurrence_bytes(record);
                                text += kOccurrenceText + string_text_bytes(record);
                            } else {
                                bytes += kGlobalOccurrenceBytes;
                                text += kGlobalOccurrenceText;
                            }
                        }
                        auto it = insert_sorted(unions[id], &hint, o);
                        if (it != unions[id].end()) {
                            inserted.emplace_back(id, it);
                            new_dedup += kDedupBytes;
                        }
                        return bytes + new_dedup <= m.account.left();
                    };
                    // in pieces between the clock's readings
                    uint64_t last_seq = std::numeric_limits<uint64_t>::max(), nt_length = 0;
                    for (size_t b = 0; b < coords.size() && !refused; b += Budget::kClockStride) {
                        const size_t e = std::min<size_t>(b + Budget::kClockStride,
                                                          coords.size());
                        if (!m.may_work(e - b)) {
                            late = true;
                            break;
                        }
                        if (m.records) {
                            // the record first, the offset after (§4.3)
                            m.oracle.map_coords(column, coords.data() + b, e - b,
                                                [&](Coord, uint64_t seq_id, Coord local) {
                                if (refused)
                                    return;
                                if (seq_id != last_seq) {
                                    nt_length = m.oracle.num_kmers_in_sequence(column, seq_id)
                                                    + k - 1;
                                    last_seq = seq_id;
                                }
                                if (local + c.offset + length > nt_length) {
                                    throw std::runtime_error(
                                            "pattern: a coordinate of " + c.kmer + " maps past "
                                            "the end of its record: the record mapping does "
                                            "not describe this annotation");
                                }
                                refused = !add(Occurrence { seq_id, local + c.offset + 1,
                                                            strand });
                            });
                        } else {
                            for (size_t t = b; t < e && !refused; ++t) {
                                refused = !add(Occurrence { coords[t], c.offset, strand });
                            }
                        }
                    }
                } else {
                    // the row's placement refused or not reached, or no coordinates for a
                    // label the row carries: not placed, stated as unknown
                    placement_complete = false;
                }
            }
            if (late || refused)
                break;
            if (listed)
                list.push_back(std::move(cl));
        }
        if (!m.commit_output(i, bytes + new_dedup, text, late, refused, &pending)) {
            undo_output(&unions, inserted, &output_cut, i);
            output_stopped = true;
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
        trim_to_unions(unions, order, kept_labels, cap, m.records, m.oracle, dict, &lists,
                       &a.fields["occurrences_cut"], [&](uint64_t bytes, uint64_t text) {
            m.account.release(bytes);
            output -= std::min(output, bytes);
            pending.drop(text);
        });
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
        // all_or_count: all or nothing (§5.2); no label of this pattern is built
        m.withhold_all(a, !refused.empty(), output_stopped, truncated,
                       output + descriptors + (summary_held ? summary : 0), &pending);
        m.finish_work(a, rows_read, pattern_units, discovery_ms, placement_ms);
        return a;
    }

    a.labels_count = count_json(empty_complete ? Relation::EXACT : labels_relation, num_labels,
                                Unit::LABELS);
    a.occurrences_count = count_json(empty_complete && m.records ? Relation::EXACT : occ_relation,
                                     occurrences_total, Unit::PLACED_OCCURRENCES);

    // by_label: the per-label summary over the returned contexts (§7.2), in label order (null
    // when the account could not hold it, above). The names as the answer may carry them
    retrieval::NameCheck names(dict);
    Json::Value by_label(Json::arrayValue);
    for (size_t r = 0; r < kept_labels && summary_held; ++r) {
        const LabelId id = order[r];
        Json::Value b;
        b["graph"] = graph_name;
        b["column"] = names.name(id);
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
        const char *strand = strand_of(c.orientation);
        for (const ContextLabel &cl : lists[i]) {
            labels.append(label_json(cl, dict[cl.label], names.name(cl.label), "kmer", strand,
                                     length, m.oracle, m.place, m.records));
        }
        f["labels"] = std::move(labels);
        a.result_fields.push_back(std::move(f));
    }
    a.unrepresentable = names.unrepresentable();
    // the labels are built: from now on they are written only
    pending.settle();
    m.finish_work(a, rows_read, pattern_units, discovery_ms, placement_ms);
    return a;
}


LabelsAnswer PatternRetrieval::retrieve_paths(const std::vector<RetrievalPath> &paths,
                                              uint64_t released, size_t length, Mode mode,
                                              const Extraction &x, const Json::Value &graph_name,
                                              bool require_verified) {
    Impl &m = *impl_;
    const size_t k = m.oracle.get_k();
    // the descriptors admit_path() charged for |paths|, the pattern's from here on
    const uint64_t descriptors = m.descriptors;
    m.descriptors = 0;
    m.admitting = false;
    if (length <= k)
        throw std::logic_error("pattern: paths of a pattern not longer than k");
    assert(paths.size() <= released);
    assert(descriptors == paths.size() * path_descriptor_bytes(k, length));
    // a label is verified only with record placement (BASIC, coordinates, the record mapping,
    // occurrences requested): the route refuses require_support "record_verified" elsewhere
    if (require_verified && !m.records) {
        throw std::logic_error("pattern: require_support record_verified without record "
                               "placement");
    }
    m.stop.reset();
    m.time_stop = false;
    m.read_stop = false;
    const size_t n = length - k + 1;
    const bool partial = mode == Mode::PARTIAL;

    // counts.labels.by_support: the labels by the support this answer gives them (the
    // strongest over the returned paths); record_verified is unknown where nothing can be
    // verified (no record placement)
    auto support_split = [&](Relation verified_relation, uint64_t verified,
                             Relation intersection_relation, uint64_t intersection) {
        Json::Value v;
        v["record_verified"] = count_json(m.records ? verified_relation : Relation::UNKNOWN,
                                          verified, Unit::LABELS);
        v["label_intersection"] = count_json(intersection_relation, intersection, Unit::LABELS);
        return v;
    };

    LabelsAnswer a = begin_labels(placement(), description_.budgeted, m.place, m.records);
    if (require_verified)
        a.fields["labels_excluded_unverified"] = count_json(Relation::UNKNOWN, 0, Unit::LABELS);
    a.labels_count["by_support"] = support_split(Relation::UNKNOWN, 0, Relation::UNKNOWN, 0);
    // no coordinates read (placement none, none_canonical, not_requested): every label of a
    // path is supported by the intersection of its k-mers' labels only (DESIGN §4.3)
    if (!m.place)
        a.notes.push_back("label_intersection_only");
    uint64_t rows_read = 0, pattern_units = 0;
    double discovery_ms = 0, placement_ms = 0;

    // the paths' descriptors were charged by admit_path
    const size_t keep = paths.size();
    if (m.skip_reads(a, x, keep, released, descriptors, partial))
        return a;

    // the rows: the distinct keys of the paths' k-mers, in answer order of their first
    // appearance, each named by its k-mer in the statements (rows_refused, anchors_truncated)
    std::vector<RowState> rows;
    std::vector<RetrievalContext> namer;
    std::vector<std::vector<size_t>> path_rows(keep);
    {
        std::unordered_map<uint64_t, size_t> index;
        for (size_t i = 0; i < keep; ++i) {
            const RetrievalPath &p = paths[i];
            if (p.sequence.size() != length || p.keys.size() != n) {
                throw std::logic_error("pattern: a released path of " + std::to_string(
                                       p.sequence.size()) + " bases and "
                                       + std::to_string(p.keys.size()) + " k-mers for a "
                                       "pattern of " + std::to_string(length) + " bases");
            }
            path_rows[i].reserve(n);
            for (size_t j = 0; j < n; ++j) {
                if (p.keys[j] == kNoKey) {
                    throw std::runtime_error("pattern: the k-mer " + p.sequence.substr(j, k)
                                             + " has no annotation row: the annotation does "
                                             "not describe this graph");
                }
                auto [it, inserted] = index.emplace(p.keys[j], rows.size());
                if (inserted) {
                    RowState row;
                    row.key = p.keys[j];
                    row.first = namer.size();
                    rows.push_back(std::move(row));
                    RetrievalContext c;
                    c.orientation = p.orientation;
                    c.kmer = p.sequence.substr(j, k);
                    c.key = p.keys[j];
                    namer.push_back(std::move(c));
                }
                path_rows[i].push_back(it->second);
            }
        }
    }

    Json::Value &refused = a.fields["rows_refused"];
    uint64_t held_lists = 0, held_hits = 0;
    const auto t0 = std::chrono::steady_clock::now();
    m.discover(rows, namer, mode, &rows_read, &pattern_units, &held_lists, &refused);
    discovery_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();

    bool truncated = false;
    const bool rows_complete = m.state_rows(a, rows, namer, &truncated);
    // the paths intersect the rows' label lists: sorted once per row (LabelRecorder gives them
    // in ascending ids), not copied and sorted for every path through the row
    for (RowState &row : rows) {
        const bool read = row.status == RowStatus::COMPLETE || row.status == RowStatus::TRUNCATED;
        if (read && !std::is_sorted(row.labels.labels.begin(), row.labels.labels.end()))
            std::sort(row.labels.labels.begin(), row.labels.labels.end());
    }

    // ---- each path's labels: those on EVERY one of its k-mers (the intersection of its
    // rows' lists; a truncated row's kept labels are true ones, so the intersection of the
    // kept lists is a true, possibly incomplete, list). A path's status is its worst row's:
    // refused, then not read, then truncated
    auto severity = [](RowStatus s) {
        switch (s) {
            case RowStatus::COMPLETE: return 0;
            case RowStatus::TRUNCATED: return 1;
            case RowStatus::PENDING:
            case RowStatus::NOT_READ: return 2;
            case RowStatus::REFUSED: return 3;
        }
        return 2;
    };
    std::vector<RowStatus> path_status(keep, RowStatus::COMPLETE);
    std::vector<std::vector<PathLabel>> carried(keep);
    std::vector<bool> needed(rows.size(), false);
    // the lists are held until the pattern's labels are built: charged as each is made; the
    // first that does not fit ends the labels of the pattern there (stop {output, max_memory},
    // as a context's labels that do not fit)
    uint64_t carried_bytes = 0;
    size_t lists_end = keep;
    auto is_read = [](RowStatus s) {
        return s == RowStatus::COMPLETE || s == RowStatus::TRUNCATED;
    };
    const auto t_lists = std::chrono::steady_clock::now();
    std::vector<LabelId> common;
    for (size_t i = 0; i < keep; ++i) {
        // the clock before each path's list (O(n) row lists to intersect, read again as they
        // are): a time stop ends the lists there, as the output of the labels (stop {output,
        // time})
        if (!m.budget.check_time() || !m.may_work(n)) {
            m.set_stop("output", "time");
            m.time_stop = true;
            lists_end = i;
            break;
        }
        RowStatus status = RowStatus::COMPLETE;
        for (size_t r : path_rows[i]) {
            if (severity(rows[r].status) > severity(status))
                status = rows[r].status;
        }
        path_status[i] = status == RowStatus::PENDING ? RowStatus::NOT_READ : status;
        if (!is_read(path_status[i]))
            continue;
        // from the smallest list, filtered in place by each other row's (sorted once per row,
        // after the reads; a row repeated at consecutive k-mers taken once)
        size_t smallest = path_rows[i][0];
        for (size_t r : path_rows[i]) {
            if (rows[r].labels.labels.size() < rows[smallest].labels.labels.size())
                smallest = r;
        }
        common.assign(rows[smallest].labels.labels.begin(), rows[smallest].labels.labels.end());
        size_t previous = smallest;
        bool late = false;
        for (size_t r : path_rows[i]) {
            if (common.empty())
                break;
            if (r == previous)
                continue;
            previous = r;
            const std::vector<LabelId> &next = rows[r].labels.labels;
            if (!m.may_work(common.size())) {
                late = true;
                break;
            }
            size_t w = 0;
            auto it = next.begin();
            for (size_t t = 0; t < common.size(); ++t) {
                it = std::lower_bound(it, next.end(), common[t]);
                if (it == next.end())
                    break;
                if (*it == common[t])
                    common[w++] = common[t];
            }
            common.resize(w);
        }
        if (late) {
            m.set_stop("output", "time");
            m.time_stop = true;
            lists_end = i;
            break;
        }
        const uint64_t bytes = kPathLabelsBytes + kPathLabelBytes * common.size();
        if (!m.account.charge(bytes)) {
            lists_end = i;
            break;
        }
        carried_bytes += bytes;
        carried[i].reserve(common.size());
        for (LabelId id : common) {
            carried[i].push_back(PathLabel { id, PathSupport::NOT_PLACED });
        }
        if (!common.empty()) {
            for (size_t r : path_rows[i]) {
                needed[r] = true;
            }
        }
    }
    const bool lists_complete = lists_end == keep;
    a.counters.label_intersection_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t_lists).count();

    // ---- step 2, the placement: the coordinates of the rows of the paths that carry a label
    // (all_or_count: only when every row was read completely, as for contexts)
    if (m.place && !m.read_stop && lists_complete && (partial || (rows_complete && !m.stop))) {
        const auto t2 = std::chrono::steady_clock::now();
        m.place_rows(rows, namer, mode, &rows_read, &pattern_units, &held_hits, &refused,
                     &needed);
        placement_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t2).count();
    }

    const std::vector<LabelRef> &dict = m.recorder->labels();

    // ---- step 3, the verification of every label carrying a path, once per (path, label), its
    // chains kept for the output and its work clocked. The chains of a path in a label:
    // join_chains over the coordinates of its k-mers (NOT_PLACED when a row of the path was not
    // placed or holds no coordinates of the label). Their occurrences, as runs (OccurrenceRun):
    // with record placement the chains whose whole path lies in one record, (seq_id, local +
    // 1), the record mapping first (§4.3; a chain whose first k-mer is in one record and whose
    // last is past that record's k-mers crosses into the next record of the column: not an
    // occurrence), and the label is verified when there is one; global, every chain
    // (kmer_coord, offset 0), nothing verified. The runs are kept for the output (PathRuns,
    // charged before they are held): when the account cannot hold a path's, its labels and the
    // later paths' are verified all the same (a label's first occurrence is enough), and their
    // labels are not output (stop {output, max_memory}: holding no more than the runs, the
    // account could not hold their occurrences' deduplication). The work is clocked inside
    // (stop {placement, time}: the path and the later ones not verified)
    const uint64_t cap = partial ? m.limits.max_occurrences_per_label
                                 : std::numeric_limits<uint64_t>::max();
    bool late = false;
    auto work = [&](uint64_t u) {
        a.counters.verification_steps += u;
        if (m.may_work(u))
            return true;
        late = true;
        return false;
    };
    enum class Chains { NOT_PLACED, DONE, LATE };
    std::vector<ShiftedList> shifted;
    // the occurrences of path |i| in label |id|, run by run in order, to on_run(run) (false
    // ends them)
    auto chains_of = [&](size_t i, LabelId id, const auto &on_run) {
        if (!work(n))
            return Chains::LATE;
        shifted.clear();
        for (size_t j = 0; j < n; ++j) {
            const RowState &row = rows[path_rows[i][j]];
            auto hit = std::lower_bound(row.hits.begin(), row.hits.end(), id,
                                        [](const LabelQuery::Hit &h, LabelId l) {
                                            return h.label < l;
                                        });
            if (!row.placed || hit == row.hits.end() || hit->label != id)
                return Chains::NOT_PLACED;
            shifted.push_back(ShiftedList { hit->coords.data(), hit->coords.size(), j });
        }
        const Column column = dict[id].column;
        const uint8_t strand = strand_rank(paths[i].orientation);
        auto emit = [&](Coord first, Coord last) {
            if (!m.records)
                return on_run(OccurrenceRun { Occurrence { first, 0, strand }, last - first + 1 });
            // record by record: the run's chains in one record are consecutive locals
            for (Coord c = first; ; ) {
                if (!work(1))
                    return false;
                const auto [seq_id, local] = m.oracle.map_coord(column, c);
                const uint64_t kmers = m.oracle.num_kmers_in_sequence(column, seq_id);
                if (local >= kmers) {
                    throw std::runtime_error("pattern: a coordinate of " + paths[i].sequence
                                             + " maps past the end of its record: the record "
                                             "mapping does not describe this annotation");
                }
                const Coord record_last = c - local + kmers - 1;
                const Coord end = std::min(last, record_last);
                if (local + n - 1 < kmers) {
                    const Coord whole = std::min(end, record_last - (n - 1));
                    if (!on_run(OccurrenceRun { Occurrence { seq_id, local + 1, strand },
                                                whole - c + 1 }))
                        return false;
                }
                if (end == last)
                    return true;
                c = end + 1;
            }
        };
        join_chains(shifted, emit, work);
        return late ? Chains::LATE : Chains::DONE;
    };

    std::vector<std::vector<PathRuns>> found(keep);
    std::vector<uint64_t> kept_bytes(keep, 0);
    // the paths whose runs were kept (found: every label of theirs with occurrences)
    size_t kept_end = keep;
    bool placement_complete = true;
    bool any_placed = false;
    const auto t_verify = std::chrono::steady_clock::now();
    for (size_t i = 0; i < lists_end; ++i) {
        // the clock before each path's verification: a time stop leaves the later paths'
        // labels not verified (stop {placement, time})
        if (m.place && !carried[i].empty()) {
            if (!m.budget.check_time()) {
                m.set_stop("placement", "time");
                m.time_stop = true;
                placement_complete = false;
                break;
            }
            if (m.occurrences_hook)
                m.occurrences_hook(i);
        }
        uint64_t kept = 0;
        for (PathLabel &pl : carried[i]) {
            if (!m.place)
                continue;
            PathRuns pr;
            pr.label = pl.label;
            uint64_t bytes = 0;
            const Chains r = chains_of(i, pl.label, [&](const OccurrenceRun &run) {
                pr.total += run.count;
                // not kept: the first occurrence verifies the label
                if (i >= kept_end)
                    return false;
                const uint64_t more = (pr.runs.empty() ? kRunsBytes : 0) + kRunBytes;
                if (!m.account.charge(more)) {
                    // nothing of this path kept: its labels and the later ones' not output
                    m.account.release(bytes + kept);
                    bytes = kept = 0;
                    found[i].clear();
                    kept_end = i;
                    m.set_stop("output", "max_memory");
                    return false;
                }
                bytes += more;
                pr.runs.push_back(run);
                return true;
            });
            if (r == Chains::LATE) {
                m.account.release(bytes);
                break;
            }
            if (r == Chains::NOT_PLACED) {
                placement_complete = false;
                continue;
            }
            any_placed = true;
            pl.support = m.records && pr.total ? PathSupport::VERIFIED : PathSupport::PLACED;
            if (i < kept_end && pr.total) {
                kept += bytes;
                found[i].push_back(std::move(pr));
            }
        }
        if (late) {
            m.account.release(kept);
            found[i].clear();
            for (PathLabel &pl : carried[i]) {
                pl.support = PathSupport::NOT_PLACED;
            }
            m.set_stop("placement", "time");
            m.time_stop = true;
            placement_complete = false;
            break;
        }
        kept_bytes[i] = kept;
    }
    a.counters.verification_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t_verify).count();
    // what a path lists: every label carrying it, or (require_support "record_verified") the
    // verified ones only
    auto listed_label = [&](const PathLabel &pl) {
        return !require_verified || pl.support == PathSupport::VERIFIED;
    };

    // the label order of §5.5 over the paths read: paths desc, column asc
    std::vector<uint64_t> label_paths(dict.size(), 0), label_verified(dict.size(), 0);
    std::vector<bool> label_carries(dict.size(), false);
    for (size_t i = 0; i < lists_end; ++i) {
        for (const PathLabel &pl : carried[i]) {
            label_carries[pl.label] = true;
            if (pl.support == PathSupport::VERIFIED)
                label_verified[pl.label]++;
            if (listed_label(pl))
                label_paths[pl.label]++;
        }
    }
    const LabelOrder labels = order_labels(label_paths, dict, partial, m.limits.max_labels,
                                           &a.fields["labels_cut"]);
    const std::vector<LabelId> &order = labels.order;
    const std::vector<uint64_t> &rank = labels.rank;
    const size_t kept_labels = labels.listed;
    const uint64_t num_labels = order.size();
    uint64_t excluded_labels = 0;
    for (LabelId id = 0; id < dict.size(); ++id) {
        if (require_verified && label_carries[id] && !label_verified[id])
            excluded_labels++;
    }

    // the text the labels built for the answer will write (AnswerVolume), from above and
    // before they are built
    PendingText pending(m.volume);
    NameText label_text(dict);
    const uint64_t by_label_text = summary_text(labels, label_text, graph_name,
                                                kPathByLabelText);
    pending.add(keep * kPathResultLabelsText + by_label_text);

    // by_label before any path's labels (as for contexts)
    const uint64_t summary = summary_bytes(labels, dict);
    const bool summary_held = m.account.charge(summary);
    if (!summary_held)
        pending.drop(by_label_text);

    // the paths' label lists with their occurrences; each label's deduplicated union (§5.4),
    // and, as for contexts, each listed label's list with its first max_occurrences_per_label
    // occurrences only in partial, cut to the unions once they are complete
    std::vector<std::vector<ContextLabel>> lists(keep);
    std::vector<std::set<Occurrence>> unions(dict.size());
    std::vector<bool> output_cut(keep, !summary_held);
    for (size_t i = lists_end; i < keep; ++i) {
        output_cut[i] = true;
    }
    uint64_t output = 0, dedup = 0;
    bool output_stopped = !summary_held || !lists_complete;
    if (output_stopped) {
        m.set_stop("output", "max_memory");
        // no path's labels are built: every path read is returned without them
        // (output_budget), not with an empty list
        output_cut.assign(keep, true);
    }
    for (size_t i = 0; i < keep && !output_stopped; ++i) {
        if (!is_read(path_status[i]))
            continue;
        if (i >= kept_end) {
            // the account could not hold this path's runs (stop {output, max_memory}, above)
            for (size_t j = i; j < keep; ++j) {
                output_cut[j] = true;
            }
            output_stopped = true;
            break;
        }
        std::vector<ContextLabel> list;
        uint64_t bytes = 0, new_dedup = 0, text = 0;
        Inserted inserted;
        late = false;
        bool refused = false;
        size_t f = 0;
        for (const PathLabel &pl : carried[i]) {
            // its runs from the verification (found[i] is in carried[i]'s order)
            const PathRuns *pr = f < found[i].size() && found[i][f].label == pl.label
                    ? &found[i][f++] : nullptr;
            if (!listed_label(pl))
                continue;
            // a label partial's max_labels cut does not list: its occurrences are counted
            const bool listed = rank[pl.label] != std::numeric_limits<uint64_t>::max();
            ContextLabel cl;
            cl.label = pl.label;
            cl.verified = pl.support == PathSupport::VERIFIED;
            if (listed) {
                bytes += label_entry_bytes(dict[pl.label].name);
                text += kPathLabelText + label_text(pl.label);
            }
            if (pl.support != PathSupport::NOT_PLACED) {
                cl.placed = true;
                const Column column = dict[pl.label].column;
                const uint64_t first = listed ? cap : 0;
                auto hint = unions[pl.label].begin();
                // a run's occurrences in pieces between the clock's readings
                auto add_run = [&](const OccurrenceRun &run) {
                    for (uint64_t t = 0; t < run.count; ) {
                        const uint64_t piece = std::min<uint64_t>(run.count - t,
                                                                  Budget::kClockStride);
                        if (!m.may_work(piece)) {
                            late = true;
                            return false;
                        }
                        for (const uint64_t e = t + piece; t < e; ++t) {
                            Occurrence o = run.first;
                            (m.records ? o.b : o.a) += t;
                            cl.total++;
                            if (cl.occurrences.size() < first) {
                                cl.occurrences.push_back(o);
                                if (m.records) {
                                    const std::string_view record
                                            = m.oracle.header_name(column, o.a);
                                    bytes += occurrence_bytes(record);
                                    text += kOccurrenceText + string_text_bytes(record);
                                } else {
                                    bytes += kGlobalOccurrenceBytes;
                                    text += kGlobalOccurrenceText;
                                }
                            }
                            auto it = insert_sorted(unions[pl.label], &hint, o);
                            if (it != unions[pl.label].end()) {
                                inserted.emplace_back(pl.label, it);
                                new_dedup += kDedupBytes;
                            }
                            // (the path's charge below, the runs it holds released first,
                            // would be refused: stopped before the union grows past it)
                            if (bytes + new_dedup > m.account.left() + kept_bytes[i]) {
                                refused = true;
                                return false;
                            }
                        }
                    }
                    return true;
                };
                // (none kept: no occurrence)
                for (size_t r = 0; pr && r < pr->runs.size() && add_run(pr->runs[r]); ++r) {}
            }
            if (late || refused)
                break;
            if (listed)
                list.push_back(std::move(cl));
        }
        // the runs kept for this path are no longer needed
        m.account.release(kept_bytes[i]);
        kept_bytes[i] = 0;
        if (!m.commit_output(i, bytes + new_dedup, text, late, refused, &pending)) {
            undo_output(&unions, inserted, &output_cut, i);
            output_stopped = true;
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
        trim_to_unions(unions, order, kept_labels, cap, m.records, m.oracle, dict, &lists,
                       &a.fields["occurrences_cut"], [&](uint64_t bytes, uint64_t text) {
            m.account.release(bytes);
            output -= std::min(output, bytes);
            pending.drop(text);
        });
    }

    // what the pattern's reads, its paths' label lists and the runs kept for paths not output
    // held is freed; the dictionary, the descriptors and the labels built for the answer stay
    uint64_t kept_left = 0;
    for (uint64_t bytes : kept_bytes) {
        kept_left += bytes;
    }
    m.account.release(held_lists + held_hits + dedup + carried_bytes + kept_left);

    // ---- the statements
    const bool all_returned = x.complete && keep == released;
    const bool anything_read = rows_read > 0;
    // every path returned, every row of every path read completely, every path's list made
    const bool paths_read = all_returned && rows_complete && lists_complete;
    // every label carrying a path verified or refuted (where verification applies)
    const bool verification_complete = !m.place || placement_complete;
    const bool labels_exact = paths_read && (!require_verified || verification_complete);
    const bool verified_exact = paths_read && verification_complete;
    const bool occurrences_exact = paths_read && verification_complete && !output_stopped;
    const Relation fallback = anything_read || keep ? Relation::AT_LEAST : Relation::UNKNOWN;
    const bool empty_complete = released == 0 && x.complete;

    a.complete = all_returned && rows_complete && lists_complete && verification_complete
                    && !output_stopped && !m.stop && a.fields["labels_cut"].isNull()
                    && a.fields["occurrences_cut"].isNull();
    a.stop = m.stop;
    a.time_limited = m.time_stop;

    if (!partial && !a.complete) {
        // all_or_count: all or nothing (§5.2), as for contexts
        m.withhold_all(a, !refused.empty(), output_stopped, truncated,
                       output + descriptors + (summary_held ? summary : 0), &pending);
        m.finish_work(a, rows_read, pattern_units, discovery_ms, placement_ms);
        return a;
    }

    uint64_t verified_labels = 0;
    for (LabelId id : order) {
        verified_labels += label_verified[id] > 0;
    }
    const Relation labels_relation = empty_complete || labels_exact ? Relation::EXACT
                                                                    : fallback;
    a.labels_count = count_json(labels_relation, num_labels, Unit::LABELS);
    if (m.records) {
        const bool exact = empty_complete || verified_exact;
        a.labels_count["by_support"] = support_split(
                exact ? Relation::EXACT : fallback, verified_labels,
                exact ? Relation::EXACT : Relation::UNKNOWN, num_labels - verified_labels);
    } else {
        a.labels_count["by_support"] = support_split(Relation::UNKNOWN, 0, labels_relation,
                                                     num_labels);
    }
    if (require_verified) {
        a.fields["labels_excluded_unverified"] = count_json(
                empty_complete || verified_exact ? Relation::EXACT : Relation::UNKNOWN,
                excluded_labels, Unit::LABELS);
    }
    a.occurrences_count = count_json(!m.records ? Relation::UNKNOWN
                                     : empty_complete || occurrences_exact ? Relation::EXACT
                                     : any_placed ? Relation::AT_LEAST : Relation::UNKNOWN,
                                     occurrences_total, Unit::PLACED_OCCURRENCES);

    // by_label: the per-label summary over the returned paths, in label order; the names as
    // the answer may carry them
    retrieval::NameCheck names(dict);
    Json::Value by_label(Json::arrayValue);
    for (size_t r = 0; r < kept_labels && summary_held; ++r) {
        const LabelId id = order[r];
        Json::Value b;
        b["graph"] = graph_name;
        b["column"] = names.name(id);
        b["paths"] = count_json(labels_exact ? Relation::EXACT : Relation::AT_LEAST,
                                label_paths[id], Unit::PATHS);
        b["paths_record_verified"] = count_json(
                !m.records ? Relation::UNKNOWN
                : verified_exact ? Relation::EXACT : Relation::AT_LEAST,
                label_verified[id], Unit::PATHS);
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
        const RetrievalPath &p = paths[i];
        Json::Value f;
        const bool read = is_read(path_status[i]) && i < lists_end;
        f["labels_status"] = output_cut[i] && is_read(path_status[i])
                ? "output_budget" : to_string(path_status[i]);
        // the labels on every k-mer of the path, their true number: known when every row of
        // the path was read completely (a truncated row's labels are only partly known)
        f["labels_total"] = read && path_status[i] == RowStatus::COMPLETE
                ? uint_json(carried[i].size()) : Json::Value();
        if (require_verified) {
            // the labels of the path left out for not being verified, their true number: known
            // when every row of the path was read completely (a truncated row's labels are
            // only partly known, as for labels_total) and every label carrying it was verified
            // or refuted (a placement stopped or refused before it leaves a label neither:
            // excluded, but not shown unverified); null otherwise, never a smaller integer
            uint64_t excluded = 0;
            bool decided = read && path_status[i] == RowStatus::COMPLETE;
            for (const PathLabel &pl : carried[i]) {
                excluded += pl.support != PathSupport::VERIFIED;
                decided &= pl.support != PathSupport::NOT_PLACED;
            }
            f["labels_excluded_unverified"] = decided ? uint_json(excluded) : Json::Value();
        }
        if (!read || output_cut[i]) {
            f["support"] = Json::Value();
            f["labels"] = Json::Value();
            a.result_fields.push_back(std::move(f));
            continue;
        }
        Json::Value labels(Json::arrayValue);
        bool any_verified = false, any_unverified = false;
        const char *strand = strand_of(p.orientation);
        for (const ContextLabel &cl : lists[i]) {
            (cl.verified ? any_verified : any_unverified) = true;
            labels.append(label_json(cl, dict[cl.label], names.name(cl.label),
                                     cl.verified ? "record_verified" : "label_intersection",
                                     strand, length, m.oracle, m.place, m.records));
        }
        // the path's support: its listed labels' (null when it lists none)
        f["support"] = !any_verified && !any_unverified ? Json::Value()
                     : any_verified && any_unverified ? Json::Value("mixed")
                     : Json::Value(any_verified ? "record_verified" : "label_intersection");
        f["labels"] = std::move(labels);
        a.result_fields.push_back(std::move(f));
    }
    a.unrepresentable = names.unrepresentable();
    pending.settle();
    m.finish_work(a, rows_read, pattern_units, discovery_ms, placement_ms);
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
