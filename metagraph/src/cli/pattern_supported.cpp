#include "pattern_supported.hpp"
#include "pattern_retrieval_impl.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string_view>

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "pattern_predicate.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph::pattern;
using graph::DeBruijnGraph;
using graph::traversal::Column;
using graph::traversal::LabelId;
using graph::traversal::LabelRef;

namespace {

using node_index = DeBruijnGraph::node_index;

// The memory model of the selection (SPEC §20.9), twice the elements as everywhere in the
// account: a held walk (its descriptor, PathSelection::Held, which the vector of the held
// walks grows by; its sequence; its support's ids); a listed path's selection labels: the ids
// and, beside each, the walks carrying it (the selection's list, freed with the pattern), and
// the text the answer keeps of them, its selection_labels strings (32 + the names' lengths)
// and its selection_strands (32 + 24 per label), the model of a listed context's (§19.9),
// which stay with the result object until the request ends as every label built for the
// answer does (§14.4); a held walk's entry in the index of the held walks by sequence (the
// mirrors' lookup)
constexpr uint64_t kHeldBytes = 192;
constexpr uint64_t kHeldLabelBytes = 8;
constexpr uint64_t kSelectionLabelsBytes = 64;
constexpr uint64_t kSelectionLabelBytes = 10;
constexpr uint64_t kSequenceIndexBytes = 8;

uint64_t held_bytes(size_t length, size_t labels) {
    return kHeldBytes + 2 * length + kHeldLabelBytes * labels;
}

// the names' lengths of the labels |ids| of |dict|: what their copies in the answer hold
uint64_t names_length(const std::vector<LabelRef> &dict, const std::vector<LabelId> &ids) {
    uint64_t length = 0;
    for (LabelId id : ids) {
        length += dict[id].name.size();
    }
    return length;
}

// the part of a listed path's selection labels the answer keeps
uint64_t selection_text_bytes(size_t labels, uint64_t names_length) {
    return retrieval::selection_labels_bytes(names_length)
            + retrieval::selection_strands_bytes(labels);
}

uint64_t selection_bytes(size_t labels, uint64_t names_length) {
    return kSelectionLabelsBytes + kSelectionLabelBytes * labels
            + selection_text_bytes(labels, names_length);
}

double ms_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t)
            .count();
}

std::string reverse_complement_of(std::string_view s) {
    std::string out(s);
    reverse_complement(out.begin(), out.end());
    return out;
}

} // namespace


// ------------------------------------------------------------------ PatternRetrieval's part

PathSupportEnv PatternRetrieval::path_support_env() {
    Impl &m = *impl_;
    return PathSupportEnv { m.oracle, m.budget, m.account, m.units, m.limits, m.mode,
                            m.budgeted, m.volume, m.deny_decode };
}

void PatternRetrieval::charge_predicate_units(uint64_t units) {
    impl_->predicate_units += units;
}


// ------------------------------------------------------------------ the selection

struct PathSelection::Held {
    std::string sequence;
    Orientation orientation = Orientation::FORWARD;
    // the predicate's labels among the walk's support (ascending ids into Bound::labels())
    std::vector<LabelId> support;
    uint64_t bytes = 0;
    // its index among the sink's paths, npos when the sink did not keep it
    size_t kept = std::numeric_limits<size_t>::max();
    bool decided = false;
    bool selected = false;
};

PathSelection::PathSelection(PathTracker &tracker, SupportedPathSink &sink,
                             const predicate::Bound &bound)
      : tracker_(tracker), sink_(sink), bound_(bound) {
    if (bound.constant())
        throw std::logic_error("pattern: a constant predicate selects no supported path itself");
    const std::vector<LabelRef> &labels = bound.labels();
    label_of_.reserve(labels.size());
    for (LabelId id = 0; id < labels.size(); ++id) {
        label_of_.emplace(labels[id].column, id);
    }
}

PathSelection::~PathSelection() {
    end_pattern();
}

Mode PathSelection::sink_mode(const PathSelectionOptions &options) {
    return options.mode;
}

uint64_t PathSelection::sink_threshold(const PathSelectionOptions &options) {
    return options.either ? options.max_predicate_contexts : options.max_paths;
}

void PathSelection::begin_pattern(const PathSelectionOptions &options, size_t length) {
    static_assert(2 * sizeof(Held) <= kHeldBytes,
                  "a held walk's charge covers its descriptor twice");
    end_pattern();
    options_ = options;
    length_ = length;
    hold_ = options.either;
    // the support only shrinks along a walk: a monotone normal form false on a branch's
    // support is false on every walk below it (with "either" a walk's decision needs its
    // mirror's support too, which only the end of the search gives: no pruning there)
    prune_ = !hold_ && bound_.monotone();
    keeping_ = options.mode != Mode::COUNT;
    stop_ = nullptr;
    supported_ = 0;
    tested_ = 0;
    selected_ = 0;
    output_cut_ = false;
    holding_ = true;
    over_ = false;
    hold_cut_ = nullptr;
    hold_memory_ = false;
    answer_ = PathSelectionAnswer();
}

void PathSelection::end_pattern() {
    // the listed paths' text stays with the answer: only the rest is given back
    tracker_.env().account.release(kept_bytes_ - std::min(answer_bytes_, kept_bytes_));
    kept_bytes_ = 0;
    kept_text_ = 0;
    answer_bytes_ = 0;
    std::vector<std::vector<LabelId>>().swap(kept_labels_);
    drop_held();
}

void PathSelection::drop_held() {
    tracker_.env().account.release(held_bytes_);
    held_bytes_ = 0;
    std::vector<Held>().swap(held_);
}

const char* PathSelection::stop_reason() const {
    if (stop_)
        return stop_;
    if (tracker_.stop_reason())
        return tracker_.stop_reason();
    return sink_.stop_reason();
}

bool PathSelection::stopped_at_threshold() const {
    return stop_ != nullptr;
}

// the predicate's labels among |t|'s top frame's support, as ids into Bound::labels()
static void predicate_labels(const PathTracker &t,
                             const std::unordered_map<Column, LabelId> &label_of,
                             std::vector<LabelId> *frame, std::vector<LabelId> *out) {
    t.frame_support(frame);
    out->clear();
    const std::vector<LabelRef> &dict = t.labels();
    for (LabelId id : *frame) {
        auto it = label_of.find(dict[id].column);
        if (it != label_of.end())
            out->push_back(it->second);
    }
    std::sort(out->begin(), out->end());
}

void PathSelection::support_of_frame(std::vector<LabelId> *out) {
    predicate_labels(tracker_, label_of_, &frame_, out);
}

bool PathSelection::may_hold() {
    support_of_frame(&present_);
    uint64_t units = 0;
    const bool holds = bound_.eval(present_.data(), present_.size(), &units);
    answer_.units += units;
    return holds;
}

bool PathSelection::prepare(const std::vector<Context> &anchors) {
    return tracker_.prepare(anchors);
}

SupportTracker::Verdict PathSelection::open(const SearchState &anchor) {
    const Verdict v = tracker_.open(anchor);
    if (v == Verdict::ALIVE && prune_ && !may_hold()) {
        tracker_.pop();
        ++answer_.pruned;
        ++answer_.pruned_by[anchor.orientation];
        return Verdict::DEAD;
    }
    return v;
}

SupportTracker::Verdict PathSelection::push(const SearchState &child) {
    const Verdict v = tracker_.push(child);
    if (v == Verdict::ALIVE && prune_ && !may_hold()) {
        tracker_.pop();
        ++answer_.pruned;
        ++answer_.pruned_by[child.orientation];
        return Verdict::DEAD;
    }
    return v;
}

void PathSelection::pop() {
    tracker_.pop();
}

SupportTracker::Verdict PathSelection::complete(const PathView &walk) {
    return tracker_.complete(walk);
}

bool PathSelection::accept(const PathView &path, const SupportTracker *support) {
    if (support != this)
        throw std::logic_error("pattern: a supported path from another tracker");
    ++supported_;
    return hold_ ? hold(path) : decide_now(path);
}

bool PathSelection::decide_now(const PathView &path) {
    support_of_frame(&present_);
    uint64_t units = 0;
    const bool selected = bound_.eval(present_.data(), present_.size(), &units);
    answer_.units += units;
    ++tested_;
    if (!selected)
        return true;
    ++selected_;
    if (keeping_) {
        retrieval::Account &account = tracker_.env().account;
        // the sink keeps a selected path under the mode's rule (partial: the first max_paths;
        // all_or_count: all while they are at most max_paths): a path beyond it is never
        // listed, so its selection labels are not charged and cannot end the list on memory
        const bool listed_by_mode = options_.mode == Mode::PARTIAL
                ? sink_.paths().size() < options_.max_paths
                : selected_ <= options_.max_paths;
        uint64_t bytes = 0;
        uint64_t text = 0;
        bool fits = true;
        if (options_.selection_labels && listed_by_mode) {
            const uint64_t names = names_length(bound_.labels(), present_);
            bytes = selection_bytes(present_.size(), names);
            text = selection_text_bytes(present_.size(), names);
            fits = account.charge(bytes);
        }
        if (!fits) {
            // the account cannot hold its selection labels: the list ends before it
            output_cut_ = true;
            keeping_ = false;
        } else {
            const size_t before = sink_.paths().size();
            if (!sink_.accept(path, &tracker_)) {
                account.release(bytes);
                return false;
            }
            const size_t after = sink_.paths().size();
            if (after == before + 1) {
                if (options_.selection_labels)
                    kept_labels_.push_back(present_);
                kept_bytes_ += bytes;
                kept_text_ += text;
            } else {
                account.release(bytes);
                if (after < before) {
                    // all_or_count: more selected than the list holds, none is kept
                    account.release(kept_bytes_);
                    kept_bytes_ = 0;
                    kept_text_ = 0;
                    kept_labels_.clear();
                }
            }
        }
    }
    if (options_.stop_at_threshold && selected_ > options_.max_paths) {
        stop_ = "max_paths";
        return false;
    }
    return true;
}

bool PathSelection::hold(const PathView &path) {
    if (!holding_)
        return true;
    if (held_.size() >= options_.max_predicate_contexts) {
        holding_ = false;
        if (options_.mode == Mode::PARTIAL) {
            hold_cut_ = "max_predicate_contexts";
        } else {
            over_ = true;
            drop_held();
            // the sink's paths can never be listed: given back to the search at once
            sink_.retain({});
        }
        if (options_.stop_at_threshold) {
            stop_ = "max_predicate_contexts";
            return false;
        }
        return true;
    }
    support_of_frame(&present_);
    retrieval::Account &account = tracker_.env().account;
    const uint64_t bytes = held_bytes(length_, present_.size());
    if (!account.charge(bytes)) {
        holding_ = false;
        if (options_.mode == Mode::PARTIAL) {
            hold_cut_ = "max_memory";
        } else {
            hold_memory_ = true;
            drop_held();
            sink_.retain({});
        }
        return true;
    }
    Held h;
    h.sequence = std::string(path.sequence);
    h.orientation = path.anchor.orientation;
    h.support = present_;
    h.bytes = bytes;
    if (options_.mode != Mode::COUNT) {
        const size_t before = sink_.paths().size();
        if (!sink_.accept(path, &tracker_)) {
            account.release(bytes);
            return false;
        }
        const size_t after = sink_.paths().size();
        if (after == before + 1) {
            h.kept = before;
        } else if (after < before) {
            for (Held &x : held_) {
                x.kept = std::numeric_limits<size_t>::max();
            }
        }
    }
    held_bytes_ += bytes;
    held_.push_back(std::move(h));
    return true;
}

const char* PathSelection::read_mirror(const std::string &sequence,
                                       std::vector<LabelId> *out) {
    const PathSupportEnv &env = tracker_.env();
    const size_t k = env.oracle.get_k();
    const size_t n = sequence.size() - k + 1;
    auto stop = [&](const char *reason) {
        if (!answer_.stop)
            answer_.stop = std::make_pair(std::string("selection"), std::string(reason));
        answer_.time_limited |= std::string_view(reason) == "time";
        return reason;
    };
    // the work and the clock before each mirror; its k-mers looked up, one unit a base
    if (!env.budget.check_time())
        return stop("time");
    if (env.units >= env.limits.max_annotation_work)
        return stop("max_annotation_work");
    env.units += sequence.size();
    answer_.mirror_units += sequence.size();
    ++answer_.lookups;
    std::vector<node_index> nodes;
    nodes.reserve(n);
    env.oracle.graph().map_to_nodes_sequentially(sequence, [&](node_index v) {
        nodes.push_back(v);
    });
    out->clear();
    if (nodes.size() != n
            || std::find(nodes.begin(), nodes.end(), DeBruijnGraph::npos) != nodes.end()) {
        // a k-mer of the mirror is not in the graph: no label supports it
        return nullptr;
    }
    if (!mirror_) {
        PathSupportOptions o;
        o.level = tracker_.level();
        mirror_ = std::make_unique<PathTracker>(env, o);
    }
    PathTracker &t = *mirror_;
    t.begin_pattern(sequence.size());
    auto finish = [&]() {
        answer_.rows += t.work().rows;
        answer_.mirror_units += t.work().units;
        // the mirror's refused rows are the selection's
        for (const Json::Value &row : t.rows_refused()) {
            Json::Value r = row;
            r["phase"] = "selection";
            answer_.rows_refused.append(std::move(r));
        }
        const char *reason = t.stop_reason();
        t.end_pattern();
        return reason;
    };
    const Context anchor { Orientation::FORWARD, 0, nodes[0], nodes[0], nodes, sequence };
    if (!t.prepare({ anchor })) {
        const char *reason = finish();
        return stop(reason ? reason : "time");
    }
    SearchState s;
    s.node = nodes[0];
    s.base_node = nodes[0];
    s.stored_reverse_complement = false;
    s.orientation = Orientation::FORWARD;
    s.side = Side::RIGHT;
    s.base = '\0';
    s.position = k;
    s.model = 0;
    s.depth = 0;
    s.spelled = std::string_view(sequence).substr(0, k);
    Verdict v = t.open(s);
    for (size_t i = 1; v == Verdict::ALIVE && i < n; ++i) {
        s.node = nodes[i];
        s.base_node = nodes[i];
        s.base = sequence[k - 1 + i];
        s.position = k + i;
        s.depth = i;
        s.spelled = std::string_view(sequence).substr(0, k + i);
        v = t.push(s);
    }
    if (v == Verdict::STOPPED) {
        const char *reason = finish();
        return stop(reason ? reason : "time");
    }
    if (v == Verdict::ALIVE) {
        const PathView view { anchor, nodes, sequence };
        if (t.complete(view) == Verdict::ALIVE)
            predicate_labels(t, label_of_, &frame_, out);
    }
    finish();
    return nullptr;
}

PathSelectionAnswer PathSelection::finish(const Result &result) {
    const auto t0 = std::chrono::steady_clock::now();
    PathSelectionAnswer &a = answer_;
    retrieval::Account &account = tracker_.env().account;
    Budget &budget = tracker_.env().budget;
    const bool all_or_count = options_.mode == Mode::ALL_OR_COUNT;
    auto done = [&]() {
        a.ms += ms_since(t0);
        return std::move(a);
    };
    auto set_stop = [&](const char *phase, const char *reason) {
        if (!a.stop)
            a.stop = std::make_pair(std::string(phase), std::string(reason));
        a.time_limited |= std::string_view(reason) == "time";
    };
    if (result.refusal || !result.anchors)
        return done();
    const AnchorCounts &anchors = *result.anchors;
    switch (anchors.extension) {
        case Extension::NOT_REQUESTED:
        case Extension::NOT_STARTED:
        case Extension::NOT_ADMITTED:
            // nothing was searched: the engine's withheld release or its stop
            a.pass = SelectionPass::NOT_STARTED;
            return done();
        case Extension::NO_ANCHORS:
            a.pass = SelectionPass::COMPLETED;
            a.tested = Count::exact(Unit::PATHS, 0);
            a.selected = Count::exact(Unit::PATHS, 0);
            a.complete = true;
            return done();
        case Extension::STOPPED:
        case Extension::COMPLETED:
            break;
    }
    // the search ran to its end
    const bool searched = anchors.extension == Extension::COMPLETED;
    const Count &supported = anchors.supported;

    // the release of |chosen| (sink indices, ascending): retained, their selection labels in
    // label order; |text| the bytes of their selection labels' strings, which stay with the
    // answer
    auto release = [&](const std::vector<size_t> &chosen,
                       const std::vector<std::vector<LabelId>> &labels,
                       const std::vector<std::vector<uint8_t>> &rows, uint64_t text) {
        sink_.retain(chosen);
        a.listed = chosen.size();
        a.named = chosen.size();
        if (!options_.selection_labels)
            return;
        answer_bytes_ += text;
        // the label order over the listed paths: paths desc, column asc
        const std::vector<LabelRef> &dict = bound_.labels();
        std::vector<uint64_t> paths(dict.size(), 0);
        for (const auto &list : labels) {
            for (LabelId id : list) {
                paths[id]++;
            }
        }
        std::vector<uint64_t> rank(dict.size(), 0);
        {
            std::vector<LabelId> order;
            for (LabelId id = 0; id < dict.size(); ++id) {
                if (paths[id])
                    order.push_back(id);
            }
            std::sort(order.begin(), order.end(), [&](LabelId x, LabelId y) {
                return paths[x] != paths[y] ? paths[x] > paths[y] : dict[x].name < dict[y].name;
            });
            for (size_t r = 0; r < order.size(); ++r) {
                rank[order[r]] = r;
            }
        }
        a.selection_labels.reserve(labels.size());
        a.selection_rows.reserve(labels.size());
        for (size_t i = 0; i < labels.size(); ++i) {
            std::vector<size_t> by(labels[i].size());
            for (size_t j = 0; j < by.size(); ++j) {
                by[j] = j;
            }
            std::sort(by.begin(), by.end(), [&](size_t x, size_t y) {
                return rank[labels[i][x]] < rank[labels[i][y]];
            });
            std::vector<LabelId> l;
            std::vector<uint8_t> r;
            for (size_t j : by) {
                l.push_back(labels[i][j]);
                r.push_back(rows.empty() ? PathSelectionAnswer::kOnWalk : rows[i][j]);
            }
            a.selection_labels.push_back(std::move(l));
            a.selection_rows.push_back(std::move(r));
        }
    };

    if (!hold_) {
        // decided at completion: the sink kept the selected paths under the request's mode
        a.tested = Count::exact(Unit::PATHS, tested_);
        a.pass = searched ? SelectionPass::COMPLETED : SelectionPass::STOPPED;
        a.selected = searched ? Count::exact(Unit::PATHS, selected_)
                              : Count::at_least(Unit::PATHS, selected_);
        if (options_.mode == Mode::COUNT)
            return done();
        const bool memory = output_cut_ || sink_.output_cut();
        std::vector<size_t> all(sink_.paths().size());
        for (size_t i = 0; i < all.size(); ++i) {
            all[i] = i;
        }
        if (all_or_count) {
            if (!searched)
                return done();      // the engine's stop withholds it (the route)
            if (selected_ > options_.max_paths) {
                a.withheld = "selected_above_threshold";
            } else if (memory) {
                a.withheld = "output_budget";
                set_stop("output", "max_memory");
            } else {
                release(all, kept_labels_, {}, kept_text_);
                a.complete = a.listed == selected_;
            }
            if (a.withheld)
                sink_.retain({});
            return done();
        }
        release(all, kept_labels_, {}, kept_text_);
        if (memory) {
            a.cut = "max_memory";
            a.output_cut = true;
            set_stop("output", "max_memory");
        } else if (selected_ > a.listed) {
            a.cut = "max_paths";
        }
        a.complete = searched && a.listed == selected_;
        return done();
    }

    // ---- "either" on a BASIC graph: the held walks decided now
    if (over_) {
        // more supported walks than the selection holds: not admitted, or, when
        // stop_at_threshold ended the search there, not started (its count unknown)
        a.pass = stop_ ? SelectionPass::NOT_STARTED : SelectionPass::NOT_ADMITTED;
        if (all_or_count)
            a.withheld = "predicate_above_threshold";
        sink_.retain({});
        return done();
    }
    auto bound_selected = [&](uint64_t s, uint64_t t) {
        return supported.relation == Relation::EXACT && supported.value >= t
                ? Count::bounds(Unit::PATHS, s, s + (supported.value - t))
                : Count::at_least(Unit::PATHS, s);
    };
    if (hold_memory_) {
        // all_or_count and count: the set the account could not hold is not decided
        a.pass = SelectionPass::STOPPED;
        a.tested = Count::exact(Unit::PATHS, 0);
        a.selected = bound_selected(0, 0);
        set_stop("selection", "max_memory");
        if (all_or_count)
            a.withheld = "predicate_budget";
        sink_.retain({});
        return done();
    }
    // partial held the first walks only, the account refusing the next one
    if (hold_cut_ && std::string_view(hold_cut_) == "max_memory")
        set_stop("selection", "max_memory");
    // the mirrors are read by a tracker of their own: the search's row cache, which no later
    // step reads, gives its allotment back first, so that the mirrors' replaces it
    if (!options_.mirrors_searched)
        tracker_.release_rows();
    // the held walks by sequence: their mirrors are found among them (when the account refuses
    // the index, no mirror is found among them nor read: every walk is decided on its own
    // support alone, as after a mirror read stopped on the account)
    std::vector<uint32_t> by_sequence;
    const uint64_t index_bytes = kSequenceIndexBytes * held_.size();
    const bool indexed = account.charge(index_bytes);
    bool ended = false;
    if (indexed) {
        by_sequence.resize(held_.size());
        for (uint32_t i = 0; i < by_sequence.size(); ++i) {
            by_sequence[i] = i;
        }
        std::sort(by_sequence.begin(), by_sequence.end(), [&](uint32_t x, uint32_t y) {
            return held_[x].sequence < held_[y].sequence;
        });
        // (n log n comparisons: light work, the clock read once they are done)
        if (!budget.check_time()) {
            set_stop("selection", "time");
            ended = true;
        }
    } else {
        set_stop("selection", "max_memory");
    }
    auto held_with = [&](const std::string &s) -> const Held* {
        auto it = std::lower_bound(by_sequence.begin(), by_sequence.end(), s,
                                   [&](uint32_t i, const std::string &x) {
            return held_[i].sequence < x;
        });
        return it != by_sequence.end() && held_[*it].sequence == s ? &held_[*it] : nullptr;
    };
    // a mirror the search did not hand over has no support only when the search reached every
    // mirror and held every supported walk
    const bool complete_set = options_.mirrors_searched && searched && !hold_cut_;
    std::vector<LabelId> mirror, present, all_labels;
    for (LabelId id = 0; id < bound_.labels().size(); ++id) {
        all_labels.push_back(id);
    }

    // the selected paths the release lists, in answer order, each listed as it is decided
    // (partial: the first max_paths selected; all_or_count: all while they are at most
    // max_paths, none when the engine's stop withholds the release): its selection labels --
    // the labels of the set it was decided on, and beside each the walks carrying it -- are
    // admitted before they are copied; |memory|: the account ended the list
    std::vector<size_t> chosen;
    std::vector<std::vector<LabelId>> labels;
    std::vector<std::vector<uint8_t>> rows;
    std::vector<uint8_t> on;
    bool memory = false;
    uint64_t label_bytes = 0;
    uint64_t label_text = 0;
    auto lists = [&]() {
        if (options_.mode == Mode::COUNT || memory)
            return false;
        // (all_or_count after a stop of the decisions: withheld, none listed)
        return all_or_count ? searched && !a.stop && selected_ <= options_.max_paths
                            : chosen.size() < options_.max_paths;
    };
    auto list = [&](size_t i, const std::vector<LabelId> &evidence) {
        if (held_[i].kept == std::numeric_limits<size_t>::max()) {
            memory = true;
            return;
        }
        if (options_.selection_labels) {
            const uint64_t names = names_length(bound_.labels(), evidence);
            const uint64_t bytes = selection_bytes(evidence.size(), names);
            if (!account.charge(bytes)) {
                memory = true;
                return;
            }
            label_bytes += bytes;
            label_text += selection_text_bytes(evidence.size(), names);
            labels.push_back(evidence);
            rows.push_back(on);
        }
        chosen.push_back(held_[i].kept);
    };
    // all_or_count lists every selected walk or none: once it can list none (more selected
    // than max_paths, the account ended the list, or the decisions stopped), what the list
    // holds is given back at once, before the later mirrors are read, and nothing more is
    // charged for it
    auto unlist = [&]() {
        account.release(label_bytes);
        label_bytes = 0;
        label_text = 0;
        std::vector<size_t>().swap(chosen);
        std::vector<std::vector<LabelId>>().swap(labels);
        std::vector<std::vector<uint8_t>>().swap(rows);
    };
    // the mirrors are read (none without the index) until a read stops on a budget (the work,
    // the account): the walks from there on are decided on their own support alone, as walks
    // whose mirror is not known; a time stop ends the decisions
    bool reads = indexed;
    for (size_t i = 0; i < held_.size() && !ended; ++i) {
        if (!(i % Budget::kReleaseClockStride) && !budget.check_time()) {
            set_stop("selection", "time");
            break;
        }
        Held &h = held_[i];
        const std::string rc = reverse_complement_of(h.sequence);
        uint64_t units = 0;
        std::optional<bool> value;
        // the mirror's support is known (|mirror|)
        bool known = true;
        // the walk was evaluated on its own support (|value|), its units charged
        bool evaluated = false;
        if (rc == h.sequence) {
            mirror = h.support;
        } else if (const Held *m = held_with(rc)) {
            mirror = m->support;
        } else if (options_.mirrors_searched) {
            mirror.clear();
            known = complete_set && indexed;
        } else {
            // the mirror is read only when its support can change the decision or would be
            // listed: the walk's own support is evaluated first (its labels sure, every other
            // label of the predicate maybe). A definite false selects nothing whatever the
            // mirror carries; a definite true selects the walk, whose mirror is then read only
            // for the label evidence a projection lists (selection_labels, selection_strands)
            value = bound_.eval3(h.support.data(), h.support.size(), all_labels.data(),
                                 all_labels.size(), &units);
            // (an evaluation's units are its own: charged here, the decision's below)
            a.units += units;
            units = 0;
            evaluated = true;
            known = false;
            const bool own = value && (!*value || !options_.selection_labels);
            if (!own && reads) {
                if (const char *reason = read_mirror(rc, &mirror)) {
                    if (std::string_view(reason) == "time")
                        break;
                    reads = false;
                } else {
                    known = true;
                }
            }
        }
        if (known) {
            present.clear();
            std::set_union(h.support.begin(), h.support.end(), mirror.begin(), mirror.end(),
                           std::back_inserter(present));
            value = bound_.eval(present.data(), present.size(), &units);
        } else if (!evaluated) {
            // the mirror's support is not known (the search stopped before it, partial held
            // the first walks only, the reads stopped, or the account refused the index):
            // decided only when no mirror support can change it (a walk evaluated so before
            // its read keeps that value)
            value = bound_.eval3(h.support.data(), h.support.size(), all_labels.data(),
                                 all_labels.size(), &units);
        }
        a.units += units;
        if (!value)
            continue;
        h.decided = true;
        h.selected = *value;
        ++tested_;
        selected_ += *value;
        if (*value && lists()) {
            const std::vector<LabelId> &evidence = known ? present : h.support;
            on.clear();
            for (LabelId id : evidence) {
                on.push_back(!known ? PathSelectionAnswer::kOnWalk
                             : (std::binary_search(h.support.begin(), h.support.end(), id)
                                    ? PathSelectionAnswer::kOnWalk : 0)
                               | (std::binary_search(mirror.begin(), mirror.end(), id)
                                    ? PathSelectionAnswer::kOnMirror : 0));
            }
            list(i, evidence);
        }
        if (all_or_count && (memory || a.stop || selected_ > options_.max_paths)
                && !chosen.empty()) {
            unlist();
        }
    }
    const bool all_decided = std::all_of(held_.begin(), held_.end(),
                                         [](const Held &h) { return h.decided; });
    account.release(indexed ? index_bytes : 0);
    const bool completed = searched && !hold_cut_ && all_decided && !a.stop;
    a.pass = completed ? SelectionPass::COMPLETED : SelectionPass::STOPPED;
    a.tested = Count::exact(Unit::PATHS, tested_);
    a.selected = completed ? Count::exact(Unit::PATHS, selected_)
                           : bound_selected(selected_, tested_);
    if (options_.mode == Mode::COUNT)
        return done();

    const char *selection_stop = a.stop ? a.stop->second.c_str() : nullptr;
    if (all_or_count) {
        if (!searched) {
            // the engine's stop withholds it (the route)
        } else if (!completed) {
            a.withheld = selection_stop && std::string_view(selection_stop) == "time"
                    ? "deadline" : selection_stop
                        && std::string_view(selection_stop) == "max_memory" && !indexed
                    ? "predicate_budget" : "annotation_budget";
        } else if (selected_ > options_.max_paths) {
            a.withheld = "selected_above_threshold";
        } else if (memory) {
            a.withheld = "output_budget";
            set_stop("output", "max_memory");
        }
        if (!searched || a.withheld) {
            account.release(label_bytes);
            sink_.retain({});
            return done();
        }
        kept_bytes_ += label_bytes;
        release(chosen, labels, rows, label_text);
        a.complete = true;
        return done();
    }
    // partial: the first max_paths selected among the decided
    kept_bytes_ += label_bytes;
    release(chosen, labels, rows, label_text);
    if (memory) {
        a.cut = "max_memory";
        a.output_cut = true;
        set_stop("output", "max_memory");
    } else if (selection_stop) {
        a.cut = selection_stop;
    } else if (hold_cut_) {
        a.cut = hold_cut_;
    } else if (selected_ > a.listed) {
        a.cut = "max_paths";
    }
    a.complete = completed && a.listed == selected_;
    return done();
}

} // namespace cli
} // namespace mtg
