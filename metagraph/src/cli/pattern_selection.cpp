#include "pattern_retrieval.hpp"
#include "pattern_retrieval_impl.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "pattern_predicate.hpp"

/**
 * The selection pass of a predicate for patterns of L <= k (SPEC-pattern-search.md
 * §19.5-§19.9): which raw contexts of a pattern satisfy the request's bound predicate, read
 * from their annotation rows restricted to the predicate's labels, under max_predicate_work,
 * the request's memory account and its deadline; and the projection predicate_only's rows
 * (retrieve_given, pattern_retrieval.cpp). Every refusal, cut and stop is stated; a count is
 * exact only when every context behind it was decided.
 */

namespace mtg {
namespace cli {

using namespace mtg::graph::pattern;
using namespace retrieval;
using graph::AnnotatedDBG;
using graph::DeBruijnGraph;
using graph::traversal::LabelOracle;
using graph::traversal::LabelQuery;
using graph::traversal::LabelRef;
using graph::traversal::LabelId;
using graph::traversal::KeyCost;
using graph::traversal::ReadPacing;
using annot::matrix::DecodeBudget;

const char* to_string(SelectionPass pass) {
    switch (pass) {
        case SelectionPass::COMPLETED: return "completed";
        case SelectionPass::STOPPED: return "stopped";
        case SelectionPass::NOT_ADMITTED: return "not_admitted";
        case SelectionPass::NOT_STARTED: return "not_started";
        case SelectionPass::CONSTANT: return "constant";
    }
    return "not_started";
}

namespace {

constexpr uint32_t kNoRow = std::numeric_limits<uint32_t>::max();

// a pass that did not run (§19.7): tested unknown, selected unknown, exact 0 when the raw
// count is exact 0 (no context can pass)
void without_pass(SelectionAnswer *a, SelectionPass pass, const Count &raw) {
    a->pass = pass;
    a->tested = Count::unknown(Unit::GRAPH_CONTEXTS);
    a->selected = raw.relation == Relation::EXACT && raw.value == 0
            ? Count::exact(Unit::GRAPH_CONTEXTS, 0) : Count::unknown(Unit::GRAPH_CONTEXTS);
}

// the selection_labels of a listed context: the list the pass keeps (ids, 4 bytes each) and,
// for the answer, its strings (32 + the names' lengths), and its selection_strands (32 + 24
// per label, §19.9)
uint64_t listed_labels_bytes(size_t ids, uint64_t names_length) {
    return selection_labels_bytes(names_length) + 4 * static_cast<uint64_t>(ids)
            + selection_strands_bytes(ids);
}
// a label of the listed contexts' label order: its count and rank (an entry of a hash map)
constexpr uint64_t kOrderEntryBytes = 48;

// the pass's row of an annotation key
struct PassRow {
    uint64_t key = kNoKey;
    // the first context needing it: its k-mer (reverse-complemented for a mirror row) names
    // the row in a statement
    uint32_t first = 0;
    // "either" on BASIC: the reverse complement's row (kNoRow: none, or not looked up), once
    // looked_up (by a lookup of its k-mer, or because the row is another row's mirror)
    uint32_t mirror = kNoRow;
    bool looked_up = false;
    bool as_mirror = false;
    // the own row of a listed context (predicate_only): held for retrieve_given
    bool keep = false;
    RowStatus status = RowStatus::PENDING;
    // the contexts needing it not yet decided or given up: its hits are freed at 0
    uint32_t needs = 0;
    LabelQuery::NodeHits hits;
    // what the account holds for the hits
    uint64_t bytes = 0;
};

} // namespace


SelectionAnswer constant_selection(bool value, const Count &raw) {
    SelectionAnswer a;
    a.pass = SelectionPass::CONSTANT;
    a.tested = raw;
    a.tested.unit = Unit::GRAPH_CONTEXTS;
    a.selected = value ? a.tested : Count::exact(Unit::GRAPH_CONTEXTS, 0);
    return a;
}

SelectionAnswer selection_without_pass(SelectionPass pass, const Count &raw) {
    if (pass == SelectionPass::COMPLETED || pass == SelectionPass::STOPPED
            || pass == SelectionPass::CONSTANT) {
        throw std::logic_error("pattern: a selection that ran (or a constant) has a pass");
    }
    SelectionAnswer a;
    without_pass(&a, pass, raw);
    return a;
}


const predicate::Binding& PatternRetrieval::bind(const predicate::Predicate &predicate,
                                                 const SelectionLimits &limits) {
    Impl &m = *impl_;
    if (m.binding)
        throw std::logic_error("pattern: the request's predicate bound twice");
    m.selection_limits = limits;
    // the bound predicate is held for the whole request: its bytes charged here, once
    m.binding = std::make_unique<predicate::Binding>(predicate::Bound::bind(
            predicate, m.oracle, m.budget,
            [&m](uint64_t bytes) { return m.account.charge(bytes); }));
    return *m.binding;
}

const predicate::Bound* PatternRetrieval::bound() const {
    const Impl &m = *impl_;
    return m.binding && m.binding->bound ? &*m.binding->bound : nullptr;
}

const char* PatternRetrieval::selection_strands() const {
    const Impl &m = *impl_;
    return m.mode != GraphMode::BASIC || m.selection_limits.either ? "either" : "context";
}

// how the selection reads the predicate's labels: the budget-aware fetch decodes whole rows
// (ROWS, never AUTO: AUTO picks DIRECT for at most 16 labels, which the budget-aware fetch
// does not read); an unbudgeted annotation with direct access is read cell by cell for at most
// 16 labels (as AUTO would), else by rows
static LabelOracle::Access selection_access_of(const LabelOracle &oracle, bool budgeted,
                                               const predicate::Bound *bound) {
    return !budgeted && oracle.supports_direct() && bound && bound->labels().size() <= 16
            ? LabelOracle::Access::DIRECT : LabelOracle::Access::ROWS;
}

const char* PatternRetrieval::selection_access() const {
    const Impl &m = *impl_;
    return selection_access_of(m.oracle, m.budgeted, bound()) == LabelOracle::Access::DIRECT
            ? "columns" : "rows";
}

uint64_t PatternRetrieval::predicate_units() const { return impl_->predicate_units; }

void PatternRetrieval::Impl::end_selection() {
    account.release(tested_bytes);
    tested_bytes = 0;
    tested_allowance = 0;
    admitting_tested = false;
    account.release(kept_bytes);
    kept_bytes = 0;
    kept_rows.clear();
    kept_labels.clear();
}

void PatternRetrieval::end_selection() { impl_->end_selection(); }

void PatternRetrieval::begin_selection(Mode mode) {
    Impl &m = *impl_;
    // what a previous pattern's pass held is not held further
    m.end_selection();
    // partial keeps half of what is left for the pass's rows and the projection, as
    // begin_release; all_or_count and count need every descriptor anyway
    m.tested_allowance = mode == Mode::PARTIAL ? m.account.left() / 2 : m.account.left();
    m.admitting_tested = true;
}

bool PatternRetrieval::admit_tested() {
    Impl &m = *impl_;
    if (!m.admitting_tested)
        return false;
    if (m.tested_bytes + kTestedBytes > m.tested_allowance || !m.account.charge(kTestedBytes)) {
        m.admitting_tested = false;
        return false;
    }
    m.tested_bytes += kTestedBytes;
    return true;
}

TestedContext PatternRetrieval::tested_context(const Context &c) const {
    const Impl &m = *impl_;
    TestedContext t;
    t.orientation = c.orientation;
    t.offset = c.offset;
    t.node = c.node;
    t.base_node = c.base_node;
    const DeBruijnGraph &graph = m.oracle.anno_graph().get_graph();
    DeBruijnGraph::node_index key = c.base_node;
    if (m.mode == GraphMode::CANONICAL) {
        // the canonical k-mer's row (the route's context_json names it alike)
        key = DeBruijnGraph::npos;
        graph.map_to_nodes(graph.get_node_sequence(c.node),
                           [&](DeBruijnGraph::node_index n) { key = n; });
    }
    const uint64_t num_rows = m.oracle.anno_graph().get_annotator().num_objects();
    t.key = key != DeBruijnGraph::npos && AnnotatedDBG::graph_to_anno_index(key) < num_rows
            ? key : kNoKey;
    return t;
}

Json::Value PatternRetrieval::selection_labels_json(const SelectionAnswer &answer,
                                                    size_t j) const {
    const predicate::Bound *b = bound();
    if (!b || j >= answer.selection_labels.size())
        throw std::logic_error("pattern: no selection_labels for this context");
    Json::Value v(Json::arrayValue);
    for (LabelId id : answer.selection_labels[j]) {
        v.append(b->labels().at(id).name);
    }
    return v;
}

Json::Value PatternRetrieval::selection_strands_json(const SelectionAnswer &answer,
                                                     size_t j) const {
    const Impl &m = *impl_;
    if (j >= answer.selection_label_rows.size()
            || answer.selection_label_rows[j].size() != answer.selection_labels.at(j).size()) {
        throw std::logic_error("pattern: no selection_strands for this context");
    }
    Json::Value v(Json::arrayValue);
    for (uint8_t on : answer.selection_label_rows[j]) {
        if (m.mode != GraphMode::BASIC) {
            // one row serves x and rc(x): no strand is known
            v.append("either");
        } else if (on == (SelectionAnswer::kOnContext | SelectionAnswer::kOnReverseComplement)) {
            v.append("both");
        } else if (on == SelectionAnswer::kOnReverseComplement) {
            v.append("reverse_complement");
        } else {
            v.append("context");
        }
    }
    return v;
}


SelectionAnswer PatternRetrieval::select(std::vector<TestedContext> &tested, uint64_t released,
                                         const Count &raw, const Extraction &x,
                                         const SelectionRequest &request) {
    Impl &m = *impl_;
    const auto t0 = std::chrono::steady_clock::now();
    if (!m.binding)
        throw std::logic_error("pattern: a selection before the request's predicate was bound");
    if (tested.size() > released)
        throw std::logic_error("pattern: more contexts tested than the engine released");
    // the release into the pass is over
    m.admitting_tested = false;
    const predicate::Bound *bound = this->bound();
    if (bound && bound->constant()) {
        throw std::logic_error("pattern: a constant predicate has no selection pass "
                               "(constant_selection)");
    }
    const size_t k = m.oracle.get_k();
    const bool partial = request.mode == Mode::PARTIAL;
    const bool all_or_count = request.mode == Mode::ALL_OR_COUNT;
    const uint64_t max_work = m.selection_limits.max_predicate_work;

    SelectionAnswer a;
    uint64_t pattern_units = 0;
    auto charge = [&](uint64_t u) {
        m.predicate_units += u;
        pattern_units += u;
    };
    auto finish = [&]() {
        a.units = pattern_units;
        a.ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
        // the statements stay in the answer, whatever else does (AnswerVolume)
        if (m.volume)
            m.volume->add(compact_json_bytes(a.rows_refused));
        return std::move(a);
    };
    // the first stop is the answer's; |ended|: the pass goes no further (every stop but the
    // descriptors' admission cut, after which partial and count test what was admitted)
    bool ended = false;
    auto set_stop = [&](const char *phase, const char *reason, bool ends = true) {
        if (!a.stop)
            a.stop = std::make_pair(std::string(phase), std::string(reason));
        if (std::string_view(reason) == "time")
            a.time_limited = true;
        ended |= ends;
    };

    if (x.withheld) {
        // nothing was released: nothing is read (§19.6 step 2)
        without_pass(&a, *x.withheld == Withheld::COUNT_ABOVE_THRESHOLD
                            ? SelectionPass::NOT_ADMITTED : SelectionPass::NOT_STARTED, raw);
        return finish();
    }
    if (tested.empty() && released == 0 && x.complete) {
        // an empty, complete release: every context (none) tested
        a.pass = SelectionPass::COMPLETED;
        a.tested = Count::exact(Unit::GRAPH_CONTEXTS, 0);
        a.selected = Count::exact(Unit::GRAPH_CONTEXTS, 0);
        return finish();
    }
    // a stop before the pass: the predicate not bound (its own stop), or the request's
    // selection work spent by an earlier pattern (sticky, as max_annotation_work, §19.8)
    const char *before = !bound ? (m.binding->stop && std::string_view(m.binding->stop) == "time"
                                        ? "time" : "max_memory")
                       : m.predicate_units >= max_work ? "max_predicate_work" : nullptr;
    if (before) {
        without_pass(&a, SelectionPass::NOT_STARTED, raw);
        set_stop("selection", before);
        if (all_or_count) {
            a.withheld = std::string_view(before) == "time" ? "deadline" : "predicate_budget";
        } else if (partial) {
            a.cut = before;
        }
        return finish();
    }
    if (tested.empty() && released == 0 && !x.complete
            && (!x.cut || *x.cut != StopReason::MAX_CONTEXTS)) {
        // partial: the engine stopped (max_steps, time) before it released a context. Nothing
        // reached the pass, as when all_or_count's and count's release is withheld for the same
        // stop: not_started in every mode (§19.7). A release cut by max_predicate_contexts
        // alone is the pass's (stopped, its bounds)
        without_pass(&a, SelectionPass::NOT_STARTED, raw);
        return finish();
    }

    const bool either = m.mode == GraphMode::BASIC && m.selection_limits.either;
    const DeBruijnGraph &graph = m.oracle.anno_graph().get_graph();
    const uint64_t num_rows = m.oracle.anno_graph().get_annotator().num_objects();
    const Projection projection = request.projection;
    const std::vector<LabelRef> &labels = bound->labels();

    // the descriptors the account could not hold ended the pass's set (§19.9): all_or_count
    // can publish nothing (no read), the other modes test what was admitted
    const bool admission_cut = tested.size() < released;
    if (admission_cut)
        set_stop("selection", "max_memory", all_or_count);

    // ---- 1. the rows: each context's own row and, with "either" on BASIC, its reverse
    // complement's (a lookup per distinct row, kept), in the order of first appearance
    std::vector<PassRow> rows;
    std::unordered_map<uint64_t, uint32_t> index;
    uint64_t row_bytes = 0;
    // the contexts whose rows are all known (a stop in the lookups ends the set)
    size_t known = admission_cut && all_or_count ? 0 : tested.size();
    // the clock of the lookups and decisions (heavy: a spelling, or an evaluation), read before
    // the first of them and then every kReleaseClockStride, and before every read
    uint64_t since_clock = Budget::kReleaseClockStride - 1;
    auto clocked = [&]() {
        since_clock = 0;
        return m.budget.check_time();
    };
    auto heavy = [&]() {
        if (++since_clock < Budget::kReleaseClockStride)
            return true;
        return clocked();
    };
    // the row of |key|, added (charged before) when new; kNoRow when the account cannot hold it
    auto add_row = [&](uint64_t key, uint32_t first, bool as_mirror) -> uint32_t {
        auto it = index.find(key);
        if (it != index.end())
            return it->second;
        if (!m.account.charge(kSelectionRowBytes))
            return kNoRow;
        row_bytes += kSelectionRowBytes;
        PassRow row;
        row.key = key;
        row.first = first;
        row.as_mirror = as_mirror;
        rows.push_back(std::move(row));
        const uint32_t r = static_cast<uint32_t>(rows.size() - 1);
        index.emplace(key, r);
        return r;
    };
    auto valid_key = [&](DeBruijnGraph::node_index key) {
        return key != DeBruijnGraph::npos && AnnotatedDBG::graph_to_anno_index(key) < num_rows
                ? key : kNoKey;
    };
    for (size_t i = 0; i < known; ++i) {
        TestedContext &c = tested[i];
        c.row = kNoRow;
        c.decided = false;
        c.selected = false;
        if (c.key == kNoKey) {
            throw std::runtime_error("pattern: the k-mer of node " + std::to_string(c.node)
                                     + " has no annotation row: the annotation does not "
                                     "describe this graph");
        }
        // light work (a hash lookup), the clock every kClockStride contexts
        if (i && i % Budget::kClockStride == 0 && !clocked()) {
            set_stop("selection", "time");
            known = i;
            break;
        }
        const uint32_t r = add_row(c.key, i, false);
        if (r == kNoRow) {
            set_stop("selection", "max_memory");
            known = i;
            break;
        }
        c.row = r;
        if (either && !rows[r].looked_up) {
            // the gate before a lookup: the clock, then the work (k units: the spelling,
            // k - 1 BOSS steps, and the reverse complement's lookup)
            if (!heavy()) {
                set_stop("selection", "time");
                known = i;
                break;
            }
            if (m.predicate_units >= max_work) {
                set_stop("selection", "max_predicate_work");
                known = i;
                break;
            }
            charge(k);
            ++a.lookups;
            std::string kmer = graph.get_node_sequence(c.node);
            reverse_complement(kmer.begin(), kmer.end());
            const std::vector<DeBruijnGraph::node_index> keys = m.oracle.keys_of_sequence(kmer);
            const uint64_t mirror = keys.size() == 1 ? valid_key(keys[0]) : kNoKey;
            rows[r].looked_up = true;
            if (mirror == kNoKey) {
                rows[r].mirror = kNoRow;
            } else if (mirror == rows[r].key) {
                // a palindromic k-mer: its own row
                rows[r].mirror = r;
            } else {
                const uint32_t mr = add_row(mirror, i, true);
                if (mr == kNoRow) {
                    set_stop("selection", "max_memory");
                    known = i;
                    break;
                }
                rows[r].mirror = mr;
                // the reverse complement is an involution: the mirror's mirror is this row
                if (!rows[mr].looked_up) {
                    rows[mr].looked_up = true;
                    rows[mr].mirror = r;
                }
            }
        }
        rows[r].needs++;
        if (either && rows[r].mirror != kNoRow && rows[r].mirror != r)
            rows[rows[r].mirror].needs++;
    }
    // a stop in the lookups ends the pass (|ended|): no row is read, as its gate would refuse
    // it (the work, the time, or an account that could not hold a row's entry, let alone its
    // statement)

    // the mirror row of a context's row (kNoRow: none, or the row itself)
    auto mirror_of = [&](uint32_t r) {
        return either && rows[r].mirror != r ? rows[r].mirror : kNoRow;
    };
    enum class State { READY, PENDING, NEVER };
    auto state_of = [&](const TestedContext &c) {
        auto one = [](const PassRow &row) {
            return row.status == RowStatus::COMPLETE ? State::READY
                 : row.status == RowStatus::REFUSED ? State::NEVER : State::PENDING;
        };
        State s = one(rows[c.row]);
        const uint32_t mr = mirror_of(c.row);
        if (mr != kNoRow) {
            const State t = one(rows[mr]);
            if (t == State::NEVER || (t == State::PENDING && s == State::READY))
                s = t;
        }
        return s;
    };
    // frees a row's hits once no context needs them (unless kept for retrieve_given)
    auto drop = [&](uint32_t r) {
        PassRow &row = rows[r];
        if (row.keep || !row.bytes)
            return;
        m.account.release(row.bytes);
        row.bytes = 0;
        row.hits = LabelQuery::NodeHits();
    };
    auto resolve = [&](const TestedContext &c) {
        for (uint32_t r : { c.row, mirror_of(c.row) }) {
            if (r == kNoRow)
                continue;
            if (rows[r].needs)
                rows[r].needs--;
            if (!rows[r].needs)
                drop(r);
        }
    };

    // ---- 3. the decisions (in answer order, as the rows are read), and the list: the first
    // max_contexts selected with their selection_labels (a projection that reads labels) and,
    // for predicate_only, their own rows kept
    uint64_t selected = 0, decided = 0;
    std::vector<uint32_t> listed;
    // a listed label and the rows it was found on (SelectionAnswer::kOnContext, ...)
    using Tagged = std::pair<LabelId, uint8_t>;
    std::vector<std::vector<Tagged>> listed_labels;
    // what the list holds: the selection_labels, the label order's entries, the kept labels'
    // names (the dictionary of retrieve_given)
    uint64_t listed_bytes = 0, kept_names_bytes = 0;
    std::unordered_map<LabelId, uint64_t> label_count;
    std::unordered_set<LabelId> kept_names;
    bool listing = request.mode != Mode::COUNT;
    std::vector<LabelId> present;
    auto decide = [&](size_t i) {
        // the clock every kReleaseClockStride lookups and decisions; not decided after it
        if (!heavy()) {
            set_stop("selection", "time");
            return false;
        }
        TestedContext &c = tested[i];
        const PassRow &own = rows[c.row];
        present.clear();
        for (const LabelQuery::Hit &h : own.hits) {
            present.push_back(h.label);
        }
        const uint32_t mr = mirror_of(c.row);
        if (mr != kNoRow) {
            for (const LabelQuery::Hit &h : rows[mr].hits) {
                present.push_back(h.label);
            }
        }
        uint64_t units = 0;
        c.selected = bound->eval(present.data(), present.size(), &units);
        c.decided = true;
        charge(units);
        ++decided;
        if (c.selected)
            ++selected;
        if (c.selected && listing && selected <= request.max_contexts) {
            // listed: its selection_labels (the set it was evaluated on, ascending ids) with
            // the row each was found on (per label the orientation that supported it), its own
            // row kept for predicate_only with the names of its labels, charged before
            std::vector<Tagged> ids;
            uint64_t bytes = 0;
            if (projection != Projection::NONE) {
                // a palindromic k-mer under "either" is its own reverse complement
                const uint8_t on_own = either && rows[c.row].mirror == c.row
                        ? SelectionAnswer::kOnContext | SelectionAnswer::kOnReverseComplement
                        : SelectionAnswer::kOnContext;
                for (const LabelQuery::Hit &h : own.hits) {
                    ids.emplace_back(h.label, on_own);
                }
                if (mr != kNoRow) {
                    for (const LabelQuery::Hit &h : rows[mr].hits) {
                        ids.emplace_back(h.label, SelectionAnswer::kOnReverseComplement);
                    }
                }
                std::sort(ids.begin(), ids.end());
                // one entry per label, the rows it was found on merged
                size_t n = 0;
                for (size_t t = 0; t < ids.size(); ++t) {
                    if (n && ids[n - 1].first == ids[t].first) {
                        ids[n - 1].second |= ids[t].second;
                    } else {
                        ids[n++] = ids[t];
                    }
                }
                ids.resize(n);
                uint64_t names = 0;
                for (const auto &[id, on] : ids) {
                    names += labels[id].name.size();
                    if (!label_count.count(id))
                        bytes += kOrderEntryBytes;
                }
                bytes += listed_labels_bytes(ids.size(), names);
            }
            uint64_t name_bytes = 0;
            if (projection == Projection::PREDICATE_ONLY) {
                for (const LabelQuery::Hit &h : own.hits) {
                    if (!kept_names.count(h.label))
                        name_bytes += label_name_bytes(labels[h.label].name);
                }
            }
            if (!m.account.charge(bytes + name_bytes)) {
                // the list ends here (stated): the pass cannot hold what it would list
                set_stop("selection", "max_memory");
                listing = false;
            } else {
                listed_bytes += bytes;
                kept_names_bytes += name_bytes;
                for (const auto &[id, on] : ids) {
                    label_count[id]++;
                }
                if (projection == Projection::PREDICATE_ONLY) {
                    rows[c.row].keep = true;
                    for (const LabelQuery::Hit &h : own.hits) {
                        kept_names.insert(h.label);
                    }
                }
                listed.push_back(static_cast<uint32_t>(i));
                listed_labels.push_back(std::move(ids));
            }
        }
        resolve(c);
        // stop_at_threshold: an early exit once more than max_contexts are selected
        if (c.selected && request.stop_at_threshold && selected > request.max_contexts)
            set_stop("selection", "max_contexts");
        return true;
    };
    size_t cursor = 0;
    auto advance = [&]() {
        while (cursor < known && !ended) {
            const TestedContext &c = tested[cursor];
            const State s = state_of(c);
            if (s == State::PENDING)
                return;
            if (s == State::NEVER) {
                resolve(c);
            } else if (!decide(cursor)) {
                return;
            }
            ++cursor;
        }
    };

    // ---- 2. the reads: one row per read, in the order of first appearance
    bool refused_any = false;
    if (known && !ended) {
        if (!m.selection_query) {
            // the request's selection query, built once (its copy of the predicate's labels
            // priced by the bound predicate's model)
            const LabelOracle::Access access = selection_access_of(m.oracle, m.budgeted, bound);
            m.selection_query = std::make_unique<LabelQuery>(m.oracle, labels, false, access);
            m.selection_query->set_max_cache_bytes(0);
            m.selection_access = access == LabelOracle::Access::DIRECT ? "columns" : "rows";
        }
        LabelQuery &query = *m.selection_query;
        auto kmer_of = [&](const PassRow &row) {
            std::string kmer = graph.get_node_sequence(tested[row.first].node);
            if (row.as_mirror)
                reverse_complement(kmer.begin(), kmer.end());
            return kmer;
        };
        for (size_t r = 0; r < rows.size() && !ended; ++r) {
            PassRow &row = rows[r];
            // every context needing it was given up (a refused row of theirs): not read
            if (!row.needs)
                continue;
            // the gate before a read: the time, the work, room for the row's statement
            if (!clocked()) {
                set_stop("selection", "time");
                break;
            }
            if (m.predicate_units >= max_work) {
                set_stop("selection", "max_predicate_work");
                break;
            }
            if (m.account.left() < m.statement()) {
                set_stop("selection", "max_memory");
                break;
            }
            ReadPacing pace = m.pacing();
            if (m.budgeted) {
                DecodeBudget decode = m.decode_budget(m.statement());
                std::vector<LabelQuery::NodeHits> out;
                std::vector<KeyCost> costs;
                out.reserve(1);
                costs.reserve(1);
                size_t refused_at = 0;
                if (!query.fetch(&row.key, 1, decode, &out, &costs, &refused_at, &pace)) {
                    if (pace.interrupted) {
                        // the units its decode reached, at least 8
                        charge(std::max<uint64_t>(8, pace.units));
                        m.budget.check_time();
                        set_stop("selection", "time");
                        break;
                    }
                    // the row does not fit what the account has left: refused, stated (the
                    // statement in its reserve), its read charged as work (what its decode
                    // reached: FetchRefusal::units, at least 8)
                    charge(Impl::refused_units(query.refusal()));
                    row.status = RowStatus::REFUSED;
                    refused_any = true;
                    const bool stated = m.account.charge(m.statement());
                    assert(stated);
                    (void)stated;
                    a.rows_refused.append(m.refusal_json(kmer_of(row), row.key, "selection",
                                                         &query.refusal()));
                    if (all_or_count)
                        break;      // the results cannot be published: no further reads
                    advance();
                    continue;
                }
                const bool held = m.account.charge(decode.held());
                assert(held);
                (void)held;
                row.bytes = decode.held();
                row.hits = std::move(out[0]);
                // the decoded row's whole size: every entry, not only the predicate's
                charge(8 + static_cast<uint64_t>(costs[0].entries) + costs[0].dependency_units);
            } else {
                std::vector<LabelQuery::NodeHits> out = query.fetch({ row.key }, &pace);
                if (pace.interrupted) {
                    charge(std::max<uint64_t>(8, pace.units));
                    m.budget.check_time();
                    set_stop("selection", "time");
                    break;
                }
                // read without a budget: what it returned is held, charged after the fact
                // (the row's size is not known: 8 and its hits)
                const uint64_t bytes = LabelQuery::held_bytes(out[0]);
                charge(8 + out[0].size());
                if (!m.account.charge(bytes)) {
                    set_stop("selection", "max_memory");
                    break;
                }
                row.bytes = bytes;
                row.hits = std::move(out[0]);
            }
            row.status = RowStatus::COMPLETE;
            ++a.rows;
            advance();
        }
        advance();
    }
    // after a stop (but time): every context whose rows were read is decided all the same
    if (!a.time_limited) {
        for (size_t i = cursor; i < known; ++i) {
            if (tested[i].decided || state_of(tested[i]) != State::READY)
                continue;
            if (!decide(i))
                break;
        }
    }
    // what the rows still hold is freed, but the kept ones (predicate_only)
    for (uint32_t r = 0; r < rows.size(); ++r) {
        rows[r].needs = 0;
        drop(r);
    }

    // ---- 4. the counts (§19.7), the selection admission, the label order
    const bool completed = !a.stop && !refused_any && x.complete
            && raw.relation == Relation::EXACT && decided == raw.value
            && tested.size() == released;
    a.pass = completed ? SelectionPass::COMPLETED : SelectionPass::STOPPED;
    a.tested = Count::exact(Unit::GRAPH_CONTEXTS, decided);
    if (completed) {
        a.selected = Count::exact(Unit::GRAPH_CONTEXTS, selected);
    } else if (raw.relation == Relation::EXACT || raw.relation == Relation::BOUNDS) {
        // the untested contexts may all pass: S <= true <= S + R_upper - T
        const uint64_t upper = raw.relation == Relation::BOUNDS ? raw.upper : raw.value;
        a.selected = Count::bounds(Unit::GRAPH_CONTEXTS, selected,
                                   selected + (upper > decided ? upper - decided : 0));
    } else {
        a.selected = Count::at_least(Unit::GRAPH_CONTEXTS, selected);
    }

    const bool time = a.stop && a.stop->second == "time";
    bool keep_list = false;
    if (all_or_count) {
        if (!completed) {
            a.withheld = time ? "deadline"
                       : a.stop && a.stop->second == "max_contexts" ? "threshold_crossed"
                                                                    : "predicate_budget";
        } else if (selected > request.max_contexts) {
            a.withheld = "selected_above_threshold";
        } else {
            keep_list = true;
        }
    } else if (partial) {
        if (a.stop && a.stop->first == "selection") {
            a.cut = a.stop->second;
        } else if (refused_any || admission_cut) {
            a.cut = "max_memory";
        }
        a.list_cut = selected > request.max_contexts;
        keep_list = true;
    }
    if (keep_list && projection != Projection::NONE && !listed.empty()) {
        // the label order of §5.5 over the listed contexts' selection_labels: contexts desc,
        // column asc. The clock is read before it (§19.9), its sorts then counted as the
        // retrieval's light work (n log n comparisons each, Impl::may_work: the next light work
        // reads the clock once they pass its stride)
        std::vector<LabelId> order;
        order.reserve(label_count.size());
        for (const auto &[id, count] : label_count) {
            order.push_back(id);
        }
        auto sort_units = [](uint64_t n) {
            uint64_t u = n;
            for (uint64_t h = n; h > 1; h >>= 1) {
                u += n;
            }
            return u;
        };
        uint64_t sorting = sort_units(order.size());
        for (const auto &list : listed_labels) {
            sorting += sort_units(list.size());
        }
        const bool in_time = m.budget.check_time();
        m.unclocked = std::min<uint64_t>(sorting, Budget::kClockStride);
        if (!in_time) {
            set_stop("output", "time");
            keep_list = false;
            if (all_or_count) {
                a.withheld = "deadline";
            } else {
                a.cut = "time";
            }
        } else {
            std::sort(order.begin(), order.end(), [&](LabelId p, LabelId q) {
                const uint64_t cp = label_count[p], cq = label_count[q];
                if (cp != cq)
                    return cp > cq;
                return labels[p].name < labels[q].name;
            });
            // the rank of each label replaces its count
            for (size_t rank = 0; rank < order.size(); ++rank) {
                label_count[order[rank]] = rank;
            }
            for (auto &list : listed_labels) {
                std::sort(list.begin(), list.end(), [&](const Tagged &p, const Tagged &q) {
                    return label_count[p.first] < label_count[q.first];
                });
            }
        }
    }
    if (keep_list) {
        a.chosen = std::move(listed);
        if (projection != Projection::NONE) {
            a.selection_labels.reserve(listed_labels.size());
            a.selection_label_rows.reserve(listed_labels.size());
            for (const auto &list : listed_labels) {
                a.selection_labels.emplace_back();
                a.selection_label_rows.emplace_back();
                for (const auto &[id, on] : list) {
                    a.selection_labels.back().push_back(id);
                    a.selection_label_rows.back().push_back(on);
                }
            }
        }
    }

    // what the answer holds: the chosen contexts' selection_labels (and the label order's
    // entries); for predicate_only their own rows and the names of their labels, taken by
    // retrieve_given. The rest of the pass is freed (its descriptors at end_selection)
    if (a.chosen.empty() || projection == Projection::NONE) {
        m.account.release(listed_bytes);
        m.account.release(kept_names_bytes);
        kept_names_bytes = 0;
    }
    uint64_t kept_slots = 0;
    if (projection == Projection::PREDICATE_ONLY && !a.chosen.empty()) {
        for (PassRow &row : rows) {
            if (!row.keep)
                continue;
            KeptRow kept;
            kept.hits = std::move(row.hits);
            kept.bytes = row.bytes + kSelectionRowBytes;
            m.kept_bytes += kept.bytes;
            kept_slots += kSelectionRowBytes;
            m.kept_rows.emplace(row.key, std::move(kept));
            row.bytes = 0;
        }
        m.kept_labels.assign(kept_names.begin(), kept_names.end());
        std::sort(m.kept_labels.begin(), m.kept_labels.end());
        m.kept_bytes += kept_names_bytes;
    } else {
        for (PassRow &row : rows) {
            if (row.bytes) {
                m.account.release(row.bytes);
                row.bytes = 0;
            }
        }
    }
    m.account.release(row_bytes - kept_slots);
    return finish();
}


LabelsAnswer PatternRetrieval::retrieve_given(const std::vector<RetrievalContext> &contexts,
                                              uint64_t released, size_t length, Mode mode,
                                              const Extraction &x,
                                              const Json::Value &graph_name) {
    Impl &m = *impl_;
    const predicate::Bound *b = bound();
    if (!b || b->constant())
        throw std::logic_error("pattern: retrieve_given without a bound, non-constant predicate");
    // the dictionary: the predicate's labels found on the kept rows (their names held since
    // the pass kept them), and each row's labels renamed into it
    Given given;
    // (a copy of their LabelRefs, priced with the names when the pass kept them)
    given.dict.reserve(m.kept_labels.size());
    for (LabelId id : m.kept_labels) {
        given.dict.push_back(b->labels()[id]);
    }
    given.ids = &m.kept_labels;
    given.rows = &m.kept_rows;
    LabelsAnswer a = retrieve_rows(contexts, released, length, mode, x, graph_name, &given);
    // the kept rows and the dictionary's names are freed: the labels built for the answer hold
    // their own copies (label_entry_bytes, by_label_bytes)
    m.account.release(m.kept_bytes);
    m.kept_bytes = 0;
    m.kept_rows.clear();
    m.kept_labels.clear();
    return a;
}

} // namespace cli
} // namespace mtg
