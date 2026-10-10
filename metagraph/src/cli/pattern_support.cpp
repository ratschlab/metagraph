#include "pattern_support.hpp"
#include "pattern_retrieval_impl.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <limits>
#include <set>
#include <stdexcept>
#include <string_view>
#include <tuple>

#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/coord_to_header.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph::pattern;
using graph::AnnotatedDBG;
using graph::DeBruijnGraph;
using graph::traversal::Column;
using graph::traversal::Coord;
using graph::traversal::FetchRefusal;
using graph::traversal::KeyCost;
using graph::traversal::LabelId;
using graph::traversal::LabelKind;
using graph::traversal::LabelOracle;
using graph::traversal::LabelQuery;
using graph::traversal::LabelRecorder;
using graph::traversal::LabelRef;
using graph::traversal::ReadPacing;
using graph::traversal::Support;
using graph::traversal::kKeyUnits;
using graph::traversal::hits_units;
namespace step = graph::traversal::support_step;
using annot::matrix::DecodeBudget;
using namespace retrieval;

namespace {

using node_index = DeBruijnGraph::node_index;

// a row in the row cache: its key and the map's node beside the runs (RowRuns::bytes)
constexpr uint64_t kRowEntryBytes = 64;

// the memory model of a kept path's labels: the list, and per label and per run of
// occurrences its record, twice its size (the convention of the walker's cost model)
constexpr uint64_t kPathLabelsBytes = 32;
constexpr uint64_t kWalkLabelBytes = 2 * sizeof(PathTracker::WalkLabel);
constexpr uint64_t kOccurrenceRunBytes = 2 * sizeof(PathTracker::Occurrences);
// an anchor's key while the anchors' rows are read: its entry in the list and in the index of
// the keys seen
constexpr uint64_t kAnchorKeyBytes = 64;

double ms_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t)
            .count();
}

} // namespace


// ------------------------------------------------------------------ the tracker

struct PathTracker::Impl {
    explicit Impl(const PathSupportEnv &env) : env(env), k(env.oracle.get_k()) {}

    const PathSupportEnv &env;
    const size_t k;

    // ---- the pattern being searched
    size_t length = 0;
    bool active = false;
    bool prepared = false;
    // the anchors' whole rows: the pattern's dictionary (the permitted set), its names held for
    // the pattern
    std::unique_ptr<LabelRecorder> recorder;
    uint64_t names_bytes = 0;
    // the later rows, over the dictionary (with coordinates when the frames carry chains)
    std::unique_ptr<LabelQuery> query;
    // the row cache: key -> the row as the steps read it, in a fixed allotment of the account
    std::unordered_map<node_index, step::RowRuns> cache;
    uint64_t cache_bytes = 0;
    uint64_t allotment = 0;
    bool path_cache = false;
    // a row the allotment cannot keep, held (charged) for the one step that reads it
    step::RowRuns transient;
    uint64_t transient_bytes = 0;
    // the frames of the DFS: frames[0, depth) are open, each holding frame_bytes[i]
    std::vector<step::Frame> frames;
    std::vector<uint64_t> frame_bytes;
    size_t depth = 0;
    uint64_t frames_bytes = 0;
    // the long merges' clock (support_step), counted across the calls of the pattern
    step::Clock clock;
    step::RecordOf record_of;

    uint64_t statement() const { return statement_bytes(k); }

    // what a read may hold: the account's remainder, less the statement reserved for its row
    DecodeBudget decode_budget() const {
        const uint64_t left = env.account.left();
        DecodeBudget decode(left - std::min(statement(), left));
        if (env.deny_decode)
            decode.deny = env.deny_decode;
        return decode;
    }

    ReadPacing pacing() const {
        ReadPacing p;
        const Deadline *d = &env.budget.deadline();
        p.ms_left = [d]() {
            return d->time_budget_ms() - d->finalize_reserve_ms() - d->elapsed_ms();
        };
        p.stop = [d]() { return d->work_expired(); };
        return p;
    }

    /**
     * The annotation key of a k-mer of the walk: the stored k-mer's on BASIC and wrapped
     * PRIMARY graphs (the engine's base node), the canonical k-mer's on a native CANONICAL one.
     * Off BASIC graphs a k-mer and its reverse complement share that row (on a wrapped PRIMARY
     * graph a node may be stored as its reverse complement:
     * SearchState::stored_reverse_complement), so a frame is told the row is the spelled
     * k-mer's: there is one row for both orientations, and no strand to keep apart.
     */
    node_index key_of(node_index node, node_index base_node, std::string_view kmer) const {
        const node_index key = env.mode == GraphMode::CANONICAL
                ? env.oracle.key_of(node, kmer) : base_node;
        if (key == DeBruijnGraph::npos
                || AnnotatedDBG::graph_to_anno_index(key) >= env.oracle.num_rows()) {
            throw std::runtime_error("pattern: the k-mer " + std::string(kmer)
                                     + " has no annotation row: the annotation does not "
                                     "describe this graph");
        }
        return key;
    }

    // only a wrapped PRIMARY graph stores a node as its reverse complement
    void check_stored(const SearchState &s) const {
        if (s.stored_reverse_complement && env.mode != GraphMode::PRIMARY) {
            throw std::logic_error("pattern: a k-mer stored as its reverse complement off a "
                                   "wrapped primary graph");
        }
    }

    // the frame object at |depth|: the stack grows with the depth reached (an open frame's
    // model counts its object twice, which covers the stack's spare capacity), never to the
    // pattern's n frames at once
    void reserve_frame(size_t depth) {
        if (frames.size() <= depth) {
            frames.resize(depth + 1);
            frame_bytes.resize(depth + 1, 0);
        }
    }

    void release_transient() {
        env.account.release(transient_bytes);
        transient_bytes = 0;
    }

    void clear_cache() {
        cache.clear();
        cache_bytes = 0;
    }
};

uint64_t PathTracker::cache_allotment(uint64_t left, uint64_t ceiling) {
    return std::min(ceiling, left / 4);
}

PathTracker::PathTracker(const PathSupportEnv &env, const PathSupportOptions &options)
      : env_(env), options_(options), impl_(std::make_unique<Impl>(env_)) {
    const AnnotationDescription d = describe_annotation(env.oracle, env.mode);
    const std::string placement = d.placement;
    if (options.level == Support::TRACE && placement != "record") {
        throw std::invalid_argument("pattern: the record level needs a BASIC index with "
                                    "coordinates and the record mapping");
    }
    if (options.require_verified && options.level != Support::TRACE) {
        throw std::invalid_argument("pattern: require_support record_verified needs the record "
                                    "level");
    }
    prune_on_chains_ = options.level == Support::TRACE;
    // the frames carry chains to prune on (the record level), or, without the record mapping,
    // for the occurrences of the labels the answer lists (global placement: never pruned on)
    chains_ = prune_on_chains_
            || (options.labels && env.limits.occurrences && placement == "global");
    if (chains_ && env.oracle.coord_to_header()) {
        const LabelOracle &oracle = env.oracle;
        Impl &m = *impl_;
        m.record_of = [&oracle, &m](LabelId label, Coord c) {
            const LabelOracle::SeqRange r
                    = oracle.sequence_range(m.recorder->labels()[label].column, c);
            return step::RecordRange { r.first, r.last };
        };
    }
    impl_->clock.expired = [this]() { return !env_.budget.check_time(); };
}

PathTracker::~PathTracker() {
    if (impl_->active)
        end_pattern();
}

const std::vector<LabelRef>& PathTracker::labels() const {
    static const std::vector<LabelRef> kNone;
    return impl_->recorder ? impl_->recorder->labels() : kNone;
}

void PathTracker::begin_pattern(size_t length) {
    Impl &m = *impl_;
    if (m.active)
        end_pattern();
    if (length <= m.k)
        throw std::logic_error("pattern: a supported-path search of a pattern not longer than k");
    m.length = length;
    m.active = true;
    m.prepared = false;
    m.depth = 0;
    m.frames_bytes = 0;
    m.clock.since = 0;
    m.clock.stopped = false;
    stop_ = nullptr;
    rows_refused_ = Json::Value(Json::arrayValue);
    work_ = PathSupportWork();
}

void PathTracker::end_pattern() {
    Impl &m = *impl_;
    if (!m.active)
        return;
    // the frames are popped by the engine on every way out; anything left is released
    for (size_t i = 0; i < m.depth; ++i) {
        env_.account.release(m.frame_bytes[i]);
    }
    m.depth = 0;
    m.release_transient();
    m.clear_cache();
    if (m.path_cache) {
        env_.oracle.path_cache().clear();
        env_.oracle.path_cache().set_bound(env_.oracle.path_cache_max());
        m.path_cache = false;
    }
    env_.account.release(m.allotment + m.names_bytes);
    m.allotment = 0;
    m.names_bytes = 0;
    m.query.reset();
    m.recorder.reset();
    m.active = false;
}

void PathTracker::release_rows() {
    Impl &m = *impl_;
    if (!m.active)
        return;
    m.release_transient();
    m.clear_cache();
    if (m.path_cache) {
        env_.oracle.path_cache().clear();
        env_.oracle.path_cache().set_bound(env_.oracle.path_cache_max());
        m.path_cache = false;
    }
    env_.account.release(m.allotment);
    m.allotment = 0;
}

// The rows of the anchors, read whole (every label: never truncated) in answer order of their
// keys' first appearance, one row per read, the time, the work and room for a statement checked
// before each; their labels are the pattern's dictionary. At the label level the rows are kept
// in the row cache for the anchors' frames; at the record level their coordinates are read by
// the frames (the query over the dictionary).
bool PathTracker::prepare(const std::vector<Context> &anchors) {
    Impl &m = *impl_;
    if (!m.active || m.prepared)
        throw std::logic_error("pattern: the supported-path search was prepared out of turn");
    m.prepared = true;
    const auto t0 = std::chrono::steady_clock::now();
    auto stopped = [&](const char *reason) {
        stop_ = reason;
        work_.ms += ms_since(t0);
        return false;
    };

    // the row cache's allotment (a quarter of what the account has left, at most the server's
    // ceiling), and the row-diff path cache in what the rows leave of it
    m.allotment = cache_allotment(env_.account.left(), env_.limits.row_cache_bytes);
    const bool held = env_.account.charge(m.allotment);
    assert(held);
    (void)held;
    if (env_.oracle.row_diff() && m.allotment) {
        const uint64_t max = env_.oracle.path_cache_max();
        env_.oracle.path_cache().set_bound(max ? std::min(max, m.allotment) : m.allotment,
                                           [&m]() {
            return m.cache_bytes < m.allotment ? m.allotment - m.cache_bytes : 0;
        });
        m.path_cache = true;
    }

    m.recorder = std::make_unique<LabelRecorder>(
            env_.oracle, LabelKind::COLUMN, std::max<uint64_t>(1, env_.oracle.num_columns()));
    m.recorder->set_max_cache_bytes(0);

    // the distinct keys, in answer order (the anchors' order), held while the rows are read;
    // spelling a k-mer (CANONICAL) is k - 1 BOSS steps: the clock every 64 anchors
    const uint64_t keys_bytes = kAnchorKeyBytes * anchors.size();
    if (!env_.account.charge(keys_bytes))
        return stopped("max_memory");
    struct Held {
        retrieval::Account &account;
        uint64_t bytes;
        ~Held() { account.release(bytes); }
    } keys_held { env_.account, keys_bytes };
    std::vector<std::pair<node_index, size_t>> keys;
    {
        std::unordered_map<node_index, size_t> seen;
        for (size_t i = 0; i < anchors.size(); ++i) {
            if (i && !(i % Budget::kReleaseClockStride) && !env_.budget.check_time())
                return stopped("time");
            const Context &a = anchors[i];
            std::string kmer;
            if (env_.mode == GraphMode::CANONICAL)
                kmer = env_.oracle.graph().get_node_sequence(a.node);
            const node_index key = m.key_of(a.node, a.base_node, kmer);
            if (seen.emplace(key, i).second)
                keys.emplace_back(key, i);
        }
    }
    const DeBruijnGraph &graph = env_.oracle.graph();
    auto name_bytes = [](std::string_view name) { return label_name_bytes(name); };
    auto refuse = [&](node_index key, const Context &a, const FetchRefusal *r) {
        const bool stated = env_.account.charge(m.statement());
        assert(stated);
        (void)stated;
        Json::Value v;
        v["kmer"] = graph.get_node_sequence(a.node);
        v["row"] = uint_json(AnnotatedDBG::graph_to_anno_index(key));
        v["phase"] = "extension";
        v["reason"] = "max_memory";
        uint64_t need = 0;
        if (r) {
            need = r->cause == FetchRefusal::DECODE ? r->need
                 : r->cause == FetchRefusal::NAMES ? r->demand + r->names_bytes : r->demand;
        }
        v["needed_bytes"] = r ? Json::Value(uint_json(need)) : Json::Value();
        v["available_bytes"] = uint_json(r ? r->left : env_.account.left());
        rows_refused_.append(std::move(v));
    };
    auto charge = [&](uint64_t u) {
        env_.units += u;
        work_.units += u;
    };
    for (const auto &[key, first] : keys) {
        if (!env_.budget.check_time())
            return stopped("time");
        if (env_.units >= env_.limits.max_annotation_work)
            return stopped("max_annotation_work");
        if (env_.account.left() < m.statement())
            return stopped("max_memory");
        ReadPacing pace = m.pacing();
        LabelRecorder::NodeLabels labels;
        uint64_t list_bytes = 0;
        if (env_.budgeted) {
            DecodeBudget decode = m.decode_budget();
            std::vector<LabelRecorder::NodeLabels> out;
            std::vector<KeyCost> costs;
            out.reserve(1);
            costs.reserve(1);
            size_t refused_at = 0;
            if (!m.recorder->fetch(&key, 1, decode, &out, &costs, &refused_at, name_bytes,
                                   &pace)) {
                if (pace.interrupted) {
                    charge(std::max<uint64_t>(kKeyUnits, pace.units));
                    env_.budget.check_time();
                    return stopped("time");
                }
                charge(m.recorder->refusal().units ? m.recorder->refusal().units : kKeyUnits);
                refuse(key, anchors[first], &m.recorder->refusal());
                return stopped("max_memory");
            }
            const bool fits = env_.account.charge(decode.held());
            assert(fits);
            (void)fits;
            // the names stay with the pattern's dictionary; the list and the naming charges
            // are freed once the row is in the cache
            m.names_bytes += m.recorder->last_names_bytes();
            list_bytes = decode.held() - m.recorder->last_names_bytes();
            // the decoded row's whole size (its entries) and its row-diff dependencies
            charge(kKeyUnits + static_cast<uint64_t>(costs[0].entries)
                   + costs[0].dependency_units);
            labels = std::move(out[0]);
        } else {
            const size_t named = m.recorder->labels().size();
            std::vector<LabelRecorder::NodeLabels> out = m.recorder->fetch({ key }, &pace);
            if (pace.interrupted) {
                charge(std::max<uint64_t>(kKeyUnits, pace.units));
                env_.budget.check_time();
                return stopped("time");
            }
            // read without a budget: the names are held whether or not the rest fits (the
            // account may pass its maximum by them, after which nothing fits)
            uint64_t names = 0;
            for (size_t id = named; id < m.recorder->labels().size(); ++id) {
                names += label_name_bytes(m.recorder->labels()[id].name);
            }
            env_.account.force(names);
            m.names_bytes += names;
            charge(kKeyUnits + out[0].total);
            list_bytes = LabelRecorder::held_bytes(out[0]);
            if (!env_.account.charge(list_bytes))
                return stopped("max_memory");
            labels = std::move(out[0]);
        }
        ++work_.rows;
        ++work_.anchor_rows;
        if (labels.truncated())
            throw std::logic_error("pattern: an anchor's row was truncated");
        if (!chains_) {
            // the anchor's frame reads this row: kept in the cache when the allotment holds it
            // (else read again by the frame, over the dictionary)
            step::RowRuns runs(false);
            const uint64_t bytes = kRowEntryBytes
                    + step::RowRuns::model_bytes(false, labels.labels.size(), 0);
            if (m.cache_bytes + bytes <= m.allotment) {
                for (LabelId id : labels.labels) {
                    charge(runs.add(id));
                }
                env_.oracle.make_room(bytes);
                m.cache_bytes += bytes;
                m.cache.emplace(key, std::move(runs));
            }
        }
        env_.account.release(list_bytes);
    }
    work_.permitted_labels = m.recorder->labels().size();
    if (work_.permitted_labels) {
        m.query = std::make_unique<LabelQuery>(env_.oracle, m.recorder->labels(), chains_,
                                               LabelOracle::Access::ROWS);
        m.query->set_max_cache_bytes(0);
    }
    work_.ms += ms_since(t0);
    return true;
}

// The row of |key| (the k-mer |kmer|) as the steps read it: from the row cache, else read
// (the gate before: the time, the work, room for its statement), turned into runs and kept in
// the cache when the allotment holds it (the cache emptied first when it does not hold it
// beside the rows there), else held for this step only (Impl::transient). Null: stopped
// (stop_reason()).
const step::RowRuns* PathTracker::row(uint64_t key, std::string_view kmer) {
    Impl &m = *impl_;
    auto it = m.cache.find(key);
    if (it != m.cache.end()) {
        ++work_.cache_hits;
        return &it->second;
    }
    if (!m.query) {
        // no anchor carries a label: no row holds one of the dictionary's
        m.transient.clear(chains_);
        return &m.transient;
    }
    auto charge = [&](uint64_t u) {
        env_.units += u;
        work_.units += u;
    };
    if (!env_.budget.check_time()) {
        stop_ = "time";
        return nullptr;
    }
    if (env_.units >= env_.limits.max_annotation_work) {
        stop_ = "max_annotation_work";
        return nullptr;
    }
    if (env_.account.left() < m.statement()) {
        stop_ = "max_memory";
        return nullptr;
    }
    ReadPacing pace = m.pacing();
    LabelQuery::NodeHits hits;
    uint64_t hits_bytes = 0;
    const node_index k = key;
    if (env_.budgeted) {
        DecodeBudget decode = m.decode_budget();
        std::vector<LabelQuery::NodeHits> out;
        std::vector<KeyCost> costs;
        out.reserve(1);
        costs.reserve(1);
        size_t refused_at = 0;
        if (!m.query->fetch(&k, 1, decode, &out, &costs, &refused_at, &pace)) {
            if (pace.interrupted) {
                charge(std::max<uint64_t>(kKeyUnits, pace.units));
                env_.budget.check_time();
                stop_ = "time";
                return nullptr;
            }
            const FetchRefusal &r = m.query->refusal();
            charge(r.units ? r.units : kKeyUnits);
            const bool stated = env_.account.charge(m.statement());
            assert(stated);
            (void)stated;
            Json::Value v;
            v["kmer"] = std::string(kmer);
            v["row"] = uint_json(AnnotatedDBG::graph_to_anno_index(key));
            v["phase"] = "extension";
            v["reason"] = "max_memory";
            const uint64_t need = r.cause == FetchRefusal::DECODE ? r.need : r.demand;
            v["needed_bytes"] = uint_json(need);
            v["available_bytes"] = uint_json(r.left);
            rows_refused_.append(std::move(v));
            stop_ = "max_memory";
            return nullptr;
        }
        hits_bytes = decode.held();
        const bool fits = env_.account.charge(hits_bytes);
        assert(fits);
        (void)fits;
        charge(kKeyUnits + static_cast<uint64_t>(costs[0].entries) + costs[0].dependency_units);
        hits = std::move(out[0]);
    } else {
        std::vector<LabelQuery::NodeHits> out = m.query->fetch({ k }, &pace);
        if (pace.interrupted) {
            charge(std::max<uint64_t>(kKeyUnits, pace.units));
            env_.budget.check_time();
            stop_ = "time";
            return nullptr;
        }
        charge(hits_units(out[0]));
        hits_bytes = LabelQuery::held_bytes(out[0]);
        if (!env_.account.charge(hits_bytes)) {
            stop_ = "max_memory";
            return nullptr;
        }
        hits = std::move(out[0]);
    }
    ++work_.rows;
    // the runs: at most one per coordinate, admitted before they are built
    uint64_t coords = 0;
    for (const LabelQuery::Hit &h : hits) {
        coords += h.coords.size();
    }
    const uint64_t bound = kRowEntryBytes
            + step::RowRuns::model_bytes(chains_, hits.size(), chains_ ? coords : 0);
    if (!env_.account.charge(bound)) {
        env_.account.release(hits_bytes);
        stop_ = "max_memory";
        return nullptr;
    }
    step::RowRuns runs(chains_);
    charge(runs.assign(hits, chains_));
    env_.account.release(hits_bytes);
    LabelQuery::NodeHits().swap(hits);
    const uint64_t bytes = kRowEntryBytes + runs.bytes();
    assert(bytes <= bound);
    if (bytes <= m.allotment) {
        // in the allotment from now on: what the row holds is the cache's
        env_.account.release(bound);
        if (m.cache_bytes + bytes > m.allotment) {
            ++work_.evictions;
            m.clear_cache();
        }
        env_.oracle.make_room(bytes);
        m.cache_bytes += bytes;
        return &m.cache.emplace(key, std::move(runs)).first->second;
    }
    // larger than the allotment: held for this step
    ++work_.uncached_rows;
    m.release_transient();
    env_.account.release(bound - bytes);
    m.transient = std::move(runs);
    m.transient_bytes = bytes;
    return &m.transient;
}

SupportTracker::Verdict PathTracker::open(const SearchState &anchor) {
    Impl &m = *impl_;
    if (!m.prepared || m.depth)
        throw std::logic_error("pattern: a supported-path frame opened out of turn");
    if (anchor.side != Side::RIGHT || anchor.spelled.size() != m.k)
        throw std::logic_error("pattern: the supported-path search extends to the right only");
    const auto t0 = std::chrono::steady_clock::now();
    const std::string_view kmer = anchor.spelled;
    m.check_stored(anchor);
    const node_index key = m.key_of(anchor.node, anchor.base_node, kmer);
    m.reserve_frame(0);
    const Verdict v = [&]() {
        if (!m.query)
            return Verdict::DEAD;
        const step::RowRuns *r = row(key, kmer);
        if (!r)
            return Verdict::STOPPED;
        if (env_.units >= env_.limits.max_annotation_work) {
            stop_ = "max_annotation_work";
            return Verdict::STOPPED;
        }
        // the frame's labels and its chains, at least one per coordinate run (a run that
        // crosses into the next record of its column is split there, charged as found)
        const Support level = chains_ ? Support::TRACE : Support::KMER;
        const uint64_t price = step::Frame::model_bytes(level, r->num_labels(), r->total_runs(),
                                                        m.k);
        if (!env_.account.charge(price)) {
            stop_ = "max_memory";
            return Verdict::STOPPED;
        }
        step::Frame &f = m.frames[0];
        const uint64_t units = f.open(step::Orientation::SPELLED, graph::traversal::Arm::RIGHT,
                                      level, kmer, kmer, *r, m.record_of, &m.clock);
        env_.units += units;
        work_.units += units;
        if (m.clock.stopped) {
            env_.account.release(price);
            stop_ = "time";
            return Verdict::STOPPED;
        }
        const uint64_t bytes = f.bytes();
        if (bytes > price && !env_.account.charge(bytes - price)) {
            env_.account.release(price);
            stop_ = "max_memory";
            return Verdict::STOPPED;
        }
        if (bytes < price)
            env_.account.release(price - bytes);
        const bool supported = prune_on_chains_ ? f.total_chains() > 0 : !f.labels().empty();
        if (!supported) {
            env_.account.release(bytes);
            return Verdict::DEAD;
        }
        m.frame_bytes[0] = bytes;
        m.depth = 1;
        m.frames_bytes += bytes;
        work_.frames_peak_bytes = std::max(work_.frames_peak_bytes, m.frames_bytes);
        return Verdict::ALIVE;
    }();
    m.release_transient();
    work_.ms += ms_since(t0);
    return v;
}

SupportTracker::Verdict PathTracker::push(const SearchState &child) {
    Impl &m = *impl_;
    if (!m.depth || child.depth != m.depth || child.side != Side::RIGHT)
        throw std::logic_error("pattern: a supported-path frame pushed out of turn");
    const auto t0 = std::chrono::steady_clock::now();
    const std::string_view kmer = child.spelled.substr(child.spelled.size() - m.k);
    m.check_stored(child);
    const node_index key = m.key_of(child.node, child.base_node, kmer);
    const Verdict v = [&]() {
        const step::RowRuns *r = row(key, kmer);
        if (!r)
            return Verdict::STOPPED;
        if (env_.units >= env_.limits.max_annotation_work) {
            stop_ = "max_annotation_work";
            return Verdict::STOPPED;
        }
        m.reserve_frame(m.depth);
        const step::Frame &top = m.frames[m.depth - 1];
        step::Frame &next = m.frames[m.depth];
        const step::Price price = top.price(*r);
        if (!env_.account.charge(price.bytes)) {
            stop_ = "max_memory";
            return Verdict::STOPPED;
        }
        const uint64_t units = top.step(child.base, kmer, *r, &next, &m.clock);
        env_.units += units;
        work_.units += units;
        if (m.clock.stopped) {
            env_.account.release(price.bytes);
            stop_ = "time";
            return Verdict::STOPPED;
        }
        const uint64_t bytes = next.bytes();
        assert(bytes <= price.bytes);
        env_.account.release(price.bytes - std::min(bytes, price.bytes));
        const bool supported = prune_on_chains_ ? next.total_chains() > 0
                                                : !next.labels().empty();
        if (!supported) {
            env_.account.release(bytes);
            return Verdict::DEAD;
        }
        m.frame_bytes[m.depth] = bytes;
        ++m.depth;
        m.frames_bytes += bytes;
        work_.frames_peak_bytes = std::max(work_.frames_peak_bytes, m.frames_bytes);
        return Verdict::ALIVE;
    }();
    m.release_transient();
    work_.ms += ms_since(t0);
    return v;
}

void PathTracker::pop() {
    Impl &m = *impl_;
    assert(m.depth);
    if (!m.depth)
        return;
    --m.depth;
    env_.account.release(m.frame_bytes[m.depth]);
    m.frames_bytes -= m.frame_bytes[m.depth];
    m.frame_bytes[m.depth] = 0;
}

SupportTracker::Verdict PathTracker::complete(const PathView &walk) {
    const Impl &m = *impl_;
    if (m.depth != m.length - m.k + 1 || walk.path.size() != m.depth)
        throw std::logic_error("pattern: a supported path completed out of turn");
    // every frame pushed was supported: the walk's support is its last frame's
    return Verdict::ALIVE;
}

size_t PathTracker::walk_num_labels() const {
    const Impl &m = *impl_;
    assert(m.depth);
    return m.frames[m.depth - 1].labels().size();
}

size_t PathTracker::walk_num_runs(size_t i) const {
    const Impl &m = *impl_;
    assert(m.depth);
    return chains_ ? m.frames[m.depth - 1].num_chains(i) : 0;
}

void PathTracker::frame_support(std::vector<LabelId> *out) const {
    const Impl &m = *impl_;
    out->clear();
    if (!m.depth)
        return;
    const step::Frame &f = m.frames[m.depth - 1];
    for (size_t i = 0; i < f.labels().size(); ++i) {
        if (!prune_on_chains_ || f.num_chains(i) > 0)
            out->push_back(f.labels()[i]);
    }
}

bool PathTracker::walk_labels(bool occurrences, std::vector<WalkLabel> *out) {
    Impl &m = *impl_;
    if (m.depth != m.length - m.k + 1)
        throw std::logic_error("pattern: the labels of a walk that is not complete");
    const auto t0 = std::chrono::steady_clock::now();
    const step::Frame &f = m.frames[m.depth - 1];
    const size_t n = m.length - m.k + 1;
    const bool records = env_.oracle.coord_to_header() != nullptr;
    std::vector<step::CoordRun> starts;
    out->reserve(out->size() + f.labels().size());
    for (size_t i = 0; i < f.labels().size(); ++i) {
        WalkLabel w;
        w.label = f.labels()[i];
        w.verified = prune_on_chains_ && f.num_chains(i) > 0;
        if (occurrences && chains_ && f.num_chains(i)) {
            starts.clear();
            f.starts(i, &starts);
            const Column column = m.recorder->labels()[w.label].column;
            w.runs.reserve(starts.size());
            for (const step::CoordRun &s : starts) {
                // one mapping per run: a chain run lies in one record (the record level), and
                // its starts are consecutive locals there
                if (!m.clock.tick()) {
                    stop_ = "time";
                    work_.ms += ms_since(t0);
                    return false;
                }
                const uint64_t count = s.last - s.first + 1;
                if (records) {
                    const auto [seq_id, local] = env_.oracle.map_coord(column, s.first);
                    // the last start's walk ends at local + count - 1 + n - 1, inside the record
                    const uint64_t kmers = env_.oracle.num_kmers_in_sequence(column, seq_id);
                    if (local + count - 1 + n - 1 >= kmers) {
                        throw std::runtime_error("pattern: a chain of a supported path passes "
                                                 "the end of its record: the record mapping "
                                                 "does not describe this annotation");
                    }
                    w.runs.push_back(Occurrences { seq_id, local + 1, count });
                } else {
                    w.runs.push_back(Occurrences { s.first, 0, count });
                }
                w.occurrences += count;
            }
        }
        out->push_back(std::move(w));
    }
    work_.ms += ms_since(t0);
    return true;
}


// ------------------------------------------------------------------ the sink

struct SupportedPathSink::Kept {
    // the labels the path lists (all of them, or the verified ones with require_verified),
    // ascending ids
    std::vector<PathTracker::WalkLabel> labels;
    // require_verified: the labels carrying the path that are not verified on it
    std::vector<LabelId> unverified;
    // the labels carrying the path (labels_total)
    uint64_t total = 0;
    // what the account holds for the path: its descriptor (the result object the route builds,
    // the answer's once the path is returned) and its labels' records here
    uint64_t descriptor = 0;
    uint64_t label_bytes = 0;
};

SupportedPathSink::SupportedPathSink(PathTracker &tracker) : tracker_(tracker) {}

SupportedPathSink::~SupportedPathSink() {
    drop_all();
}

void SupportedPathSink::begin_pattern(Mode mode, uint64_t max_paths) {
    drop_all();
    mode_ = mode;
    max_paths_ = max_paths;
    accepted_ = 0;
    output_cut_ = false;
    dropped_ = false;
    keeping_ = mode != Mode::COUNT;
    stop_ = nullptr;
    // half of what is left for the kept paths, the other half for the search's frames and rows
    allowance_ = tracker_.env().account.left() / 2;
}

void SupportedPathSink::end_pattern(size_t returned) {
    // the descriptors of the returned paths stay with the answer (their result objects)
    uint64_t answer = 0;
    for (size_t i = 0; i < returned && i < kept_.size(); ++i) {
        answer += kept_[i].descriptor;
    }
    held_ -= std::min(answer, held_);
    drop_all();
}

void SupportedPathSink::drop_all() {
    tracker_.env().account.release(held_);
    held_ = 0;
    std::vector<Context>().swap(paths_);
    std::vector<Kept>().swap(kept_);
}

void SupportedPathSink::retain(const std::vector<size_t> &indices) {
    std::vector<Context> paths;
    std::vector<Kept> kept;
    paths.reserve(indices.size());
    kept.reserve(indices.size());
    uint64_t bytes = 0;
    for (size_t j = 0; j < indices.size(); ++j) {
        const size_t i = indices[j];
        if (i >= kept_.size() || (j && i <= indices[j - 1]))
            throw std::logic_error("pattern: the supported paths retained are not kept ones");
        bytes += kept_[i].descriptor + kept_[i].label_bytes;
        paths.push_back(std::move(paths_[i]));
        kept.push_back(std::move(kept_[i]));
    }
    tracker_.env().account.release(held_ - bytes);
    held_ = bytes;
    paths_.swap(paths);
    kept_.swap(kept);
}

bool SupportedPathSink::accept(const PathView &path, const SupportTracker *support) {
    if (support != &tracker_)
        throw std::logic_error("pattern: a supported path from another tracker");
    ++accepted_;
    if (!keeping_)
        return true;
    if (mode_ == Mode::ALL_OR_COUNT && accepted_ > max_paths_) {
        // all or none: none, once they are too many
        drop_all();
        dropped_ = true;
        keeping_ = false;
        return true;
    }
    if (mode_ == Mode::PARTIAL && paths_.size() >= max_paths_) {
        keeping_ = false;
        return true;
    }
    const PathSupportEnv &env = tracker_.env();
    const bool labels = tracker_.options().labels;
    const bool occurrences = labels && env.limits.occurrences && tracker_.chains();
    // what the path will hold, from above: its descriptor, every label carrying it, and the
    // runs of the labels with occurrences
    uint64_t bytes = path_descriptor_bytes(env.oracle.get_k(), path.sequence.size());
    if (labels) {
        const size_t n = tracker_.walk_num_labels();
        bytes += kPathLabelsBytes + kWalkLabelBytes * n;
        for (size_t i = 0; occurrences && i < n; ++i) {
            bytes += kOccurrenceRunBytes * tracker_.walk_num_runs(i);
        }
    }
    if (held_ + bytes > allowance_ || !env.account.charge(bytes)) {
        // the list ends here: partial lists the paths before it, all_or_count none
        output_cut_ = true;
        keeping_ = false;
        if (mode_ != Mode::PARTIAL)
            drop_all();
        return true;
    }
    Kept kept;
    kept.descriptor = path_descriptor_bytes(env.oracle.get_k(), path.sequence.size());
    kept.label_bytes = bytes - kept.descriptor;
    if (labels) {
        std::vector<PathTracker::WalkLabel> all;
        if (!tracker_.walk_labels(occurrences, &all)) {
            env.account.release(bytes);
            stop_ = "time";
            return false;
        }
        const bool verified_only = tracker_.options().require_verified;
        const auto &listed = tracker_.options().listed;
        const std::vector<LabelRef> &dict = tracker_.labels();
        for (PathTracker::WalkLabel &w : all) {
            if (listed && !listed(dict[w.label].column))
                continue;
            ++kept.total;
            if (verified_only && !w.verified) {
                kept.unverified.push_back(w.label);
                continue;
            }
            kept.labels.push_back(std::move(w));
        }
    }
    held_ += bytes;
    paths_.push_back(path.copy());
    kept_.push_back(std::move(kept));
    return true;
}


// ------------------------------------------------------------------ the labels of the
// released paths

namespace {

// one occurrence of a label in a path: (seq_id, 1-based start) with the record mapping,
// (kmer_coord, 0) without it; the path's strand
struct Occurrence {
    uint64_t a = 0;
    uint64_t b = 0;
    uint8_t strand = 0;
    bool operator<(const Occurrence &o) const {
        return std::tie(a, b, strand) < std::tie(o.a, o.b, o.strand);
    }
};

// a label listed for a path: its occurrences (partial: the first max_occurrences_per_label of
// the path's), their number, whether it is verified
struct PathLabel {
    LabelId label = 0;
    bool verified = false;
    bool placed = false;
    uint64_t total = 0;
    std::vector<Occurrence> occurrences;
};

// the label order (items desc, column asc) over the labels with |counts|, and the labels the
// results list (partial: the first max_labels, the cut stated); rank[id] the place of a listed
// label, uint64_t's maximum for one not listed
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

Json::Value label_json(const PathLabel &pl, const LabelRef &label, const char *strand,
                       size_t length, const LabelOracle &oracle, bool place, bool records) {
    Json::Value l;
    l["column"] = label.name;
    l["support"] = pl.verified ? "record_verified" : "label_intersection";
    if (!place)
        return l;
    if (!pl.placed) {
        if (records)
            l["occurrences"] = count_json(Relation::UNKNOWN, 0, Unit::PLACED_OCCURRENCES);
        l["occurrence_list"] = Json::Value();
        return l;
    }
    if (records)
        l["occurrences"] = count_json(Relation::EXACT, pl.total, Unit::PLACED_OCCURRENCES);
    const size_t k = oracle.get_k();
    Json::Value occ(Json::arrayValue);
    for (const Occurrence &o : pl.occurrences) {
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

// an occurrence's bytes in the account and its text in the answer
uint64_t occurrence_model(const LabelOracle &oracle, Column column, const Occurrence &o,
                          bool records, uint64_t *text) {
    if (!records) {
        *text += kGlobalOccurrenceText;
        return kGlobalOccurrenceBytes;
    }
    const std::string_view record = oracle.header_name(column, o.a);
    *text += kOccurrenceText + string_text_bytes(record);
    return occurrence_bytes(record);
}

} // namespace

LabelsAnswer SupportedPathSink::labels_answer(const Extraction &x, uint64_t released,
                                              size_t length, const Json::Value &graph_name) {
    if (!tracker_.options().labels)
        throw std::logic_error("pattern: the labels of supported paths that kept none");
    const PathSupportEnv &env = tracker_.env();
    const LabelOracle &oracle = env.oracle;
    const AnnotationDescription d = describe_annotation(oracle, env.mode);
    const std::string index_placement = d.placement;
    const bool partial = mode_ == Mode::PARTIAL;
    // the labels' occurrences are listed where the frames carried chains and the request asks
    // for occurrences: record placement (record_verified) or global (record bounds unknown)
    const bool place = env.limits.occurrences && tracker_.chains();
    const bool records = place && index_placement == "record";
    const bool verified_only = tracker_.options().require_verified;
    const std::vector<LabelRef> &dict = tracker_.labels();
    // what this answer places: nothing at the label level of an index that could place (its
    // search read no coordinates), else the index's placement, or not_requested
    const char *placement = !env.limits.occurrences ? "not_requested"
                          : !place && index_placement == "record" ? "none"
                                                                  : d.placement;

    LabelsAnswer a;
    a.fields["placement"] = placement;
    a.fields["annotation"] = env.budgeted ? "budgeted" : "unbudgeted";
    a.fields["rows_refused"] = tracker_.rows_refused();
    a.fields["anchors_truncated"] = Json::Value(Json::arrayValue);
    a.fields["labels_cut"] = Json::Value();
    a.fields["occurrences_cut"] = Json::Value();
    a.fields["by_label"] = Json::Value();
    if (verified_only)
        a.fields["labels_excluded_unverified"] = count_json(Relation::UNKNOWN, 0, Unit::LABELS);
    a.labels_count = count_json(Relation::UNKNOWN, 0, Unit::LABELS);
    auto support_split = [&](Relation verified_relation, uint64_t verified,
                             Relation intersection_relation, uint64_t intersection) {
        Json::Value v;
        v["record_verified"] = count_json(tracker_.prunes_on_chains() ? verified_relation
                                                                      : Relation::UNKNOWN,
                                          verified, Unit::LABELS);
        v["label_intersection"] = count_json(intersection_relation, intersection, Unit::LABELS);
        return v;
    };
    a.labels_count["by_support"] = support_split(Relation::UNKNOWN, 0, Relation::UNKNOWN, 0);
    a.occurrences_count = count_json(Relation::UNKNOWN, 0, Unit::PLACED_OCCURRENCES);
    if (!env.budgeted)
        a.notes.push_back("annotation_unbudgeted");
    if (place && !records)
        a.notes.push_back("record_bounds_unknown");
    // (at the record level a label is verified by the chains the search carried, placed or
    // not: output.occurrences false leaves its occurrences out, not its support)
    if (!place && !tracker_.prunes_on_chains())
        a.notes.push_back("label_intersection_only");

    std::optional<std::pair<std::string, std::string>> stop;
    bool time_stop = false;
    auto set_stop = [&](const char *reason) {
        if (!stop)
            stop = std::make_pair(std::string("output"), std::string(reason));
        time_stop |= std::string_view(reason) == "time";
    };
    auto finish = [&]() {
        const PathSupportWork &w = tracker_.work();
        a.work["annotation_rows"] = uint_json(w.rows);
        a.work["annotation_units"] = uint_json(w.units);
        a.work["memory_bytes"] = uint_json(env.account.peak());
        a.timing["support_ms"] = w.ms;
        if (env.volume)
            env.volume->add(compact_json_bytes(a.fields["rows_refused"]));
        a.stop = stop;
        a.time_limited = time_stop;
    };

    if (x.withheld) {
        // nothing is released: no label is built
        finish();
        return a;
    }
    const size_t keep = std::min<size_t>(x.returned, paths_.size());
    if (keep < released) {
        // the allowance could not hold every path the release names
        set_stop("max_memory");
        if (!partial) {
            a.withheld = "output_budget";
            finish();
            return a;
        }
        a.cut = "max_memory";
    }

    // the label order of the returned paths: paths desc, column asc
    std::vector<uint64_t> label_paths(dict.size(), 0), label_verified(dict.size(), 0);
    std::vector<bool> label_carries(dict.size(), false);
    for (size_t i = 0; i < keep; ++i) {
        for (const PathTracker::WalkLabel &w : kept_[i].labels) {
            label_carries[w.label] = true;
            label_paths[w.label]++;
            if (w.verified)
                label_verified[w.label]++;
        }
        for (LabelId id : kept_[i].unverified) {
            label_carries[id] = true;
        }
    }
    const LabelOrder labels = order_labels(label_paths, dict, partial, env.limits.max_labels,
                                           &a.fields["labels_cut"]);
    uint64_t excluded_labels = 0;
    for (LabelId id = 0; id < dict.size(); ++id) {
        if (verified_only && label_carries[id] && !label_verified[id])
            excluded_labels++;
    }

    // the text the labels will write, from above and before they are built
    PendingText pending(env.volume);
    std::vector<uint64_t> name_text(dict.size(), 0);
    auto name_bytes = [&](LabelId id) {
        if (!name_text[id])
            name_text[id] = string_text_bytes(dict[id].name);
        return name_text[id];
    };
    const uint64_t graph_text = compact_json_bytes(graph_name);
    uint64_t by_label_text = 0, summary = 0;
    for (size_t r = 0; r < labels.listed; ++r) {
        by_label_text += kPathByLabelText + name_bytes(labels.order[r]) + graph_text;
        summary += by_label_bytes(dict[labels.order[r]].name);
    }
    pending.add(keep * kPathResultLabelsText + by_label_text);
    const bool summary_held = env.account.charge(summary);
    if (!summary_held)
        pending.drop(by_label_text);

    // each path's labels with their occurrences, each label's deduplicated union, the clock
    // read every Budget::kClockStride occurrences
    const uint64_t cap = partial ? env.limits.max_occurrences_per_label
                                 : std::numeric_limits<uint64_t>::max();
    std::vector<std::vector<PathLabel>> lists(keep);
    std::vector<std::set<Occurrence>> unions(dict.size());
    std::vector<bool> output_cut(keep, !summary_held);
    uint64_t output = 0, dedup = 0, unclocked = 0;
    bool output_stopped = !summary_held;
    if (output_stopped)
        set_stop("max_memory");
    auto clocked = [&](uint64_t u) {
        if (unclocked + u > Budget::kClockStride) {
            unclocked = 0;
            if (!env.budget.check_time())
                return false;
        }
        unclocked += u;
        return true;
    };
    for (size_t i = 0; i < keep && !output_stopped; ++i) {
        const Kept &kept = kept_[i];
        const uint8_t strand = strand_rank(paths_[i].orientation);
        std::vector<PathLabel> list;
        uint64_t bytes = 0, new_dedup = 0, text = 0;
        std::vector<std::pair<LabelId, std::set<Occurrence>::iterator>> inserted;
        bool late = false, refused = false;
        for (const PathTracker::WalkLabel &w : kept.labels) {
            const bool listed = labels.rank[w.label] != std::numeric_limits<uint64_t>::max();
            PathLabel pl;
            pl.label = w.label;
            pl.verified = w.verified;
            if (listed) {
                bytes += label_entry_bytes(dict[w.label].name);
                text += kPathLabelText + name_bytes(w.label);
            }
            if (place) {
                pl.placed = true;
                const Column column = dict[w.label].column;
                const uint64_t first = listed ? cap : 0;
                auto hint = unions[w.label].begin();
                for (const PathTracker::Occurrences &run : w.runs) {
                    for (uint64_t t = 0; t < run.count && !late && !refused; ++t) {
                        if (!clocked(1)) {
                            late = true;
                            break;
                        }
                        Occurrence o { records ? run.a : run.a + t, records ? run.b + t : 0,
                                       strand };
                        pl.total++;
                        if (pl.occurrences.size() < first) {
                            pl.occurrences.push_back(o);
                            bytes += occurrence_model(oracle, column, o, records, &text);
                        }
                        const size_t before = unions[w.label].size();
                        auto it = unions[w.label].insert(hint, o);
                        hint = std::next(it);
                        if (unions[w.label].size() > before) {
                            inserted.emplace_back(w.label, it);
                            new_dedup += kDedupBytes;
                        }
                        // (the charge below, the path's records released first, would
                        // be refused: stopped before the union grows past it)
                        if (bytes + new_dedup > env.account.left() + kept.label_bytes)
                            refused = true;
                    }
                    if (late || refused)
                        break;
                }
            }
            if (late || refused)
                break;
            if (listed)
                list.push_back(std::move(pl));
        }
        // the path's label records are no longer needed: its labels are built from here on
        env.account.release(kept.label_bytes);
        held_ -= kept.label_bytes;
        kept_[i].label_bytes = 0;
        bool committed = !late && !refused && env.account.charge(bytes + new_dedup);
        if (committed) {
            pending.add(text);
            if (!env.budget.check_time()) {
                pending.drop(text);
                env.account.release(bytes + new_dedup);
                late = true;
                committed = false;
            }
        }
        if (!committed) {
            set_stop(late ? "time" : "max_memory");
            for (const auto &[id, it] : inserted) {
                unions[id].erase(it);
            }
            for (size_t j = i; j < keep; ++j) {
                output_cut[j] = true;
            }
            output_stopped = true;
            break;
        }
        output += bytes;
        dedup += new_dedup;
        std::sort(list.begin(), list.end(), [&](const PathLabel &x, const PathLabel &y) {
            return labels.rank[x.label] < labels.rank[y.label];
        });
        lists[i] = std::move(list);
    }

    // partial: each label lists the first max_occurrences_per_label occurrences of its union
    uint64_t occurrences_total = 0;
    for (LabelId id : labels.order) {
        occurrences_total += unions[id].size();
    }
    if (partial && place) {
        std::vector<const Occurrence*> bound(dict.size(), nullptr);
        uint64_t cut_labels = 0;
        for (size_t r = 0; r < labels.listed; ++r) {
            const std::set<Occurrence> &u = unions[labels.order[r]];
            if (u.size() <= cap)
                continue;
            cut_labels++;
            if (cap)
                bound[labels.order[r]] = &*std::next(u.begin(), cap - 1);
        }
        if (cut_labels) {
            uint64_t bytes = 0, text = 0;
            for (auto &list : lists) {
                for (PathLabel &pl : list) {
                    if (!bound[pl.label])
                        continue;
                    auto end = std::upper_bound(pl.occurrences.begin(), pl.occurrences.end(),
                                                *bound[pl.label]);
                    for (auto it = end; it != pl.occurrences.end(); ++it) {
                        bytes += occurrence_model(oracle, dict[pl.label].column, *it, records,
                                                  &text);
                    }
                    pl.occurrences.erase(end, pl.occurrences.end());
                }
            }
            env.account.release(bytes);
            output -= std::min(output, bytes);
            pending.drop(text);
            Json::Value cut = reason_json("max_occurrences_per_label");
            cut["labels"] = uint_json(cut_labels);
            a.fields["occurrences_cut"] = std::move(cut);
        }
    }
    // the deduplication state is freed; the labels built for the answer and the summary stay
    env.account.release(dedup);

    const bool all_returned = x.complete && keep == released;
    const bool empty_complete = released == 0 && x.complete;
    const bool exact = all_returned;
    const bool occurrences_exact = all_returned && !output_stopped;
    const Relation fallback = keep ? Relation::AT_LEAST : Relation::UNKNOWN;

    a.complete = all_returned && !output_stopped && !stop && a.fields["labels_cut"].isNull()
                    && a.fields["occurrences_cut"].isNull();
    if (!partial && !a.complete) {
        // all_or_count: all or nothing
        a.withheld = time_stop ? "deadline" : "output_budget";
        env.account.release(output + (summary_held ? summary : 0));
        pending.drop_all();
        finish();
        return a;
    }

    uint64_t verified_labels = 0;
    for (LabelId id : labels.order) {
        verified_labels += label_verified[id] > 0;
    }
    const uint64_t num_labels = labels.order.size();
    const Relation labels_relation = empty_complete || exact ? Relation::EXACT : fallback;
    a.labels_count = count_json(labels_relation, num_labels, Unit::LABELS);
    if (tracker_.prunes_on_chains()) {
        a.labels_count["by_support"] = support_split(
                labels_relation, verified_labels,
                empty_complete || exact ? Relation::EXACT : Relation::UNKNOWN,
                num_labels - verified_labels);
    } else {
        a.labels_count["by_support"] = support_split(Relation::UNKNOWN, 0, labels_relation,
                                                     num_labels);
    }
    if (verified_only) {
        a.fields["labels_excluded_unverified"] = count_json(
                empty_complete || exact ? Relation::EXACT : Relation::UNKNOWN,
                excluded_labels, Unit::LABELS);
    }
    a.occurrences_count = count_json(!records ? Relation::UNKNOWN
                                     : empty_complete || occurrences_exact ? Relation::EXACT
                                     : keep ? Relation::AT_LEAST : Relation::UNKNOWN,
                                     occurrences_total, Unit::PLACED_OCCURRENCES);

    // by_label over the returned paths, in label order
    if (summary_held) {
        Json::Value by_label(Json::arrayValue);
        for (size_t r = 0; r < labels.listed; ++r) {
            const LabelId id = labels.order[r];
            Json::Value b;
            b["graph"] = graph_name;
            b["column"] = dict[id].name;
            b["paths"] = count_json(exact ? Relation::EXACT : Relation::AT_LEAST,
                                    label_paths[id], Unit::PATHS);
            b["paths_record_verified"] = count_json(
                    !tracker_.prunes_on_chains() ? Relation::UNKNOWN
                    : exact ? Relation::EXACT : Relation::AT_LEAST,
                    label_verified[id], Unit::PATHS);
            b["occurrences"] = count_json(!records ? Relation::UNKNOWN
                                          : occurrences_exact ? Relation::EXACT
                                          : keep ? Relation::AT_LEAST : Relation::UNKNOWN,
                                          unions[id].size(), Unit::PLACED_OCCURRENCES);
            by_label.append(std::move(b));
        }
        a.fields["by_label"] = std::move(by_label);
    }

    // the results' labels
    a.result_fields.reserve(keep);
    for (size_t i = 0; i < keep; ++i) {
        Json::Value f;
        f["labels_status"] = output_cut[i] ? "output_budget" : "complete";
        f["labels_total"] = uint_json(kept_[i].total);
        if (verified_only)
            f["labels_excluded_unverified"] = uint_json(kept_[i].unverified.size());
        if (output_cut[i]) {
            f["support"] = Json::Value();
            f["labels"] = Json::Value();
            a.result_fields.push_back(std::move(f));
            continue;
        }
        Json::Value list(Json::arrayValue);
        bool any_verified = false, any_unverified = false;
        const char *strand = strand_of(paths_[i].orientation);
        for (const PathLabel &pl : lists[i]) {
            (pl.verified ? any_verified : any_unverified) = true;
            list.append(label_json(pl, dict[pl.label], strand, length, oracle, place, records));
        }
        f["support"] = !any_verified && !any_unverified ? Json::Value()
                     : any_verified && any_unverified ? Json::Value("mixed")
                     : Json::Value(any_verified ? "record_verified" : "label_intersection");
        f["labels"] = std::move(list);
        a.result_fields.push_back(std::move(f));
    }
    pending.settle();
    finish();
    return a;
}


// ------------------------------------------------------------------ the release

std::optional<Extraction> supported_release(const Result &result, const Request &request,
                                            const SupportedPathSink &sink,
                                            uint64_t *returned) {
    *returned = 0;
    if (request.mode == Mode::COUNT || result.refusal)
        return std::nullopt;
    if (!result.anchors || !request.extend_paths || !request.support || !request.sink)
        throw std::logic_error("pattern: the release of supported paths without their search");
    const AnchorCounts &a = *result.anchors;
    const bool partial = request.mode == Mode::PARTIAL;
    auto withheld_for = [](StopReason reason) {
        switch (reason) {
            case StopReason::MAX_CONTEXTS:
            case StopReason::MAX_ANCHORS:
            case StopReason::MAX_PATHS:
                return Withheld::THRESHOLD_CROSSED;
            case StopReason::TIME:
                return Withheld::DEADLINE;
            case StopReason::MAX_STEPS:
                return Withheld::DISCOVERY_BUDGET;
            case StopReason::EXTERNAL:
                return Withheld::EXTERNAL;
        }
        return Withheld::DISCOVERY_BUDGET;
    };
    Extraction x;
    switch (a.extension) {
        case Extension::NOT_ADMITTED:
            x.withheld = Withheld::ANCHORS_ABOVE_THRESHOLD;
            return x;
        case Extension::NO_ANCHORS:
            x.complete = true;
            return x;
        case Extension::NOT_REQUESTED:
        case Extension::NOT_STARTED:
            // a stop before the extension: an anchor set not known completely is not extended
            if (!result.stop)
                throw std::logic_error("pattern: an extension that did not start without a stop");
            if (partial) {
                x.cut = result.stop->reason;
            } else {
                x.withheld = withheld_for(result.stop->reason);
            }
            return x;
        case Extension::STOPPED:
        case Extension::COMPLETED:
            break;
    }
    const uint64_t supported = a.supported.value;
    if (!partial) {
        if (result.stop) {
            x.withheld = withheld_for(result.stop->reason);
        } else if (supported > request.max_paths) {
            x.withheld = Withheld::COUNT_ABOVE_THRESHOLD;
        } else {
            *returned = supported;
            x.returned = sink.paths().size();
            x.complete = a.supported.relation == Relation::EXACT && x.returned == supported;
        }
        return x;
    }
    // partial: the supported paths completed before a stop, the first max_paths
    *returned = sink.output_cut() ? std::min(sink.accepted(), request.max_paths)
                                  : sink.paths().size();
    x.returned = sink.paths().size();
    if (result.stop) {
        x.cut = result.stop->reason;
    } else if (supported > request.max_paths) {
        x.cut = StopReason::MAX_PATHS;
    }
    x.complete = !result.stop && a.supported.relation == Relation::EXACT
                    && x.returned == supported;
    return x;
}

} // namespace cli
} // namespace mtg
