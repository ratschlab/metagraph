#include "traverse.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <future>
#include <limits>
#include <set>

#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"
#include "common/threads/threading.hpp"
#include "config/config.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using namespace mtg::graph::traversal;

TraverseLimits traverse_limits(const Config &config) {
    TraverseLimits limits;
    limits.max_time_ms = config.traverse_max_time_ms;
    limits.max_seeds = config.traverse_max_seeds;
    return limits;
}

namespace {

/**
 * A JSON object read field by field: every field read is known, and done() refuses any
 * other, so that a misspelled or unsupported field is an error rather than silently ignored.
 */
class Strict {
  public:
    Strict(const Json::Value &value, std::string path) : value_(value), path_(std::move(path)) {
        if (!value_.isObject())
            throw InvalidRequest(path_ + ": expected an object");
    }

    std::string child(const char *key) const { return path_ + "." + key; }

    bool has(const char *key) {
        known_.insert(key);
        return value_.isMember(key);
    }

    const Json::Value& raw(const char *key) {
        known_.insert(key);
        return value_[key];
    }

    std::string str(const char *key, const std::string &def) {
        if (!has(key))
            return def;
        if (!value_[key].isString())
            throw InvalidRequest(child(key) + ": expected a string");
        return value_[key].asString();
    }

    uint64_t uint(const char *key, uint64_t def, uint64_t min = 0,
                  uint64_t max = std::numeric_limits<uint64_t>::max()) {
        if (!has(key))
            return def;
        const Json::Value &v = value_[key];
        if (!v.isUInt64() || (v.isDouble() && v.asDouble() != std::floor(v.asDouble())))
            throw InvalidRequest(child(key) + ": expected a non-negative integer");
        const uint64_t x = v.asUInt64();
        if (x < min || x > max) {
            throw InvalidRequest(child(key) + ": out of range [" + std::to_string(min) + ", "
                                 + std::to_string(max) + "]");
        }
        return x;
    }

    double number(const char *key, double def, double min, double max) {
        if (!has(key))
            return def;
        const Json::Value &v = value_[key];
        if (!v.isNumeric())
            throw InvalidRequest(child(key) + ": expected a number");
        const double x = v.asDouble();
        if (!(x >= min && x <= max))
            throw InvalidRequest(child(key) + ": out of range");
        return x;
    }

    std::vector<std::string> strings(const char *key) {
        std::vector<std::string> out;
        if (!has(key))
            return out;
        const Json::Value &v = value_[key];
        if (!v.isArray())
            throw InvalidRequest(child(key) + ": expected an array of strings");
        for (const Json::Value &x : v) {
            if (!x.isString())
                throw InvalidRequest(child(key) + ": expected an array of strings");
            out.push_back(x.asString());
        }
        return out;
    }

    template <class T>
    T enumeration(const char *key, T def,
                  std::initializer_list<std::pair<const char*, T>> values) {
        if (!has(key))
            return def;
        std::string names;
        if (value_[key].isString()) {
            for (const auto &[name, x] : values) {
                if (value_[key].asString() == name)
                    return x;
            }
        }
        for (const auto &[name, x] : values) {
            names += std::string(names.empty() ? "" : ", ") + "\"" + name + "\"";
        }
        throw InvalidRequest(child(key) + ": expected one of " + names);
    }

    // refuse the fields that were not read
    void done() const {
        for (const std::string &name : value_.getMemberNames()) {
            if (!known_.count(name)) {
                throw InvalidRequest(child(name.c_str()) + ": unknown field");
            }
        }
    }

  private:
    const Json::Value &value_;
    std::string path_;
    std::set<std::string> known_;
};

Json::Value uint_json(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

Json::Value labels_json(const std::vector<LabelId> &labels) {
    Json::Value out(Json::arrayValue);
    for (LabelId l : labels) {
        out.append(l);
    }
    return out;
}

Json::Value interval_json(uint64_t begin, uint64_t end) {
    Json::Value iv(Json::arrayValue);
    iv.append(uint_json(begin));
    iv.append(uint_json(end));
    return iv;
}

Json::Value runs_json(const std::vector<KmerInterval> &runs) {
    Json::Value out(Json::arrayValue);
    for (const KmerInterval &r : runs) {
        out.append(interval_json(r.begin, r.end));
    }
    return out;
}

const char* status_name(ArmResult::Status status) {
    return status == ArmResult::COMPLETE ? "complete" : "truncated";
}

// a label end's reason as stated: the text that qualifies it, if any, else the reason
const char* reason_name(EndReason reason, const char *text) {
    return text && *text ? text : to_string(reason);
}


/********************************* /traverse *********************************/

struct TraverseRequest {
    std::vector<Seed> seeds;
    Strategy strategy;
};

TraverseRequest parse_traverse_request(const Json::Value &json, const TraverseLimits &limits) {
    TraverseRequest req;
    Strict s(json, "request");

    const Json::Value &seeds = s.raw("seeds");
    if (!seeds.isArray() || seeds.empty())
        throw InvalidRequest("request.seeds: expected a non-empty array");
    if (seeds.size() > limits.max_seeds) {
        throw InvalidRequest("request.seeds: " + std::to_string(seeds.size())
                             + " seeds, over the server's limit of "
                             + std::to_string(limits.max_seeds));
    }
    for (Json::ArrayIndex i = 0; i < seeds.size(); ++i) {
        Strict ss(seeds[i], "request.seeds[" + std::to_string(i) + "]");
        Seed seed;
        seed.sequence = ss.str("sequence", "");
        if (seed.sequence.empty())
            throw InvalidRequest(ss.child("sequence") + " is required");
        if (seed.sequence.size() > kMaxSeedBp) {
            throw InvalidRequest(ss.child("sequence") + ": " + std::to_string(seed.sequence.size())
                                 + " bp, over the limit of " + std::to_string(kMaxSeedBp));
        }
        // omitted: derive the permitted set from the seed; an empty array reads like "no
        // label permitted" and is refused
        if (ss.has("labels")) {
            seed.labels = ss.strings("labels");
            if (seed.labels.empty()) {
                throw InvalidRequest(ss.child("labels") + ": expected a non-empty array of "
                                     "labels; omit the field to derive the labels from the seed");
            }
            if (seed.labels.size() > kMaxSeedLabels) {
                throw InvalidRequest(ss.child("labels") + ": " + std::to_string(seed.labels.size())
                                     + " labels, over the limit of "
                                     + std::to_string(kMaxSeedLabels));
            }
        }
        ss.done();
        req.seeds.push_back(std::move(seed));
    }

    Strategy &st = req.strategy;
    st.time_budget_ms = limits.max_time_ms;
    if (s.has("strategy")) {
        Strict t(s.raw("strategy"), "request.strategy");
        st.direction = t.enumeration<Strategy::Direction>("direction", Strategy::BOTH,
            { { "both", Strategy::BOTH }, { "left", Strategy::LEFT }, { "right", Strategy::RIGHT } });
        if (t.has("labels")) {
            Strict l(t.raw("labels"), "request.strategy.labels");
            st.max_seed_labels = l.uint("max_seed_labels", st.max_seed_labels, 1, kMaxSeedLabels);
            using OptionalKind = std::optional<LabelKind>;
            st.seed_label_kind = l.enumeration<OptionalKind>("seed_label_kind", OptionalKind(),
                { { "auto", OptionalKind() },
                  { "column", OptionalKind(LabelKind::COLUMN) },
                  { "header", OptionalKind(LabelKind::HEADER) } });
            l.done();
        }
        if (t.has("branching")) {
            Strict b(t.raw("branching"), "request.strategy.branching");
            // a number in range, or "unlimited"
            auto limit = [&](const char *key, size_t &value) {
                if (!b.has(key))
                    return;
                const Json::Value &v = b.raw(key);
                if (v.isString() && v.asString() == "unlimited") {
                    value = Strategy::kUnlimited;
                } else {
                    // anything else, null included, must be a number in range
                    value = b.uint(key, 0, 0, 1'000'000);
                }
            };
            limit("max_label_branches", st.max_label_branches);
            limit("max_splits_per_path", st.max_splits_per_path);
            b.done();
        }
        if (t.has("bounds")) {
            Strict b(t.raw("bounds"), "request.strategy.bounds");
            st.max_extension_bp = b.uint("max_extension_bp", st.max_extension_bp, 1);
            // the server's bound clamps the request's
            st.time_budget_ms = std::min(b.number("time_budget_ms", limits.max_time_ms, 1,
                                                  std::numeric_limits<double>::max()),
                                         limits.max_time_ms);
            b.done();
        }
        t.done();
    }
    s.done();
    return req;
}

Json::Value strategy_to_json(const Strategy &st) {
    Json::Value j;
    j["direction"] = st.direction == Strategy::BOTH ? "both"
                   : st.direction == Strategy::LEFT ? "left" : "right";
    j["labels"]["max_seed_labels"] = uint_json(st.max_seed_labels);
    j["labels"]["seed_label_kind"] = st.seed_label_kind ? to_string(*st.seed_label_kind) : "auto";
    auto limit_json = [](size_t value) {
        return value == Strategy::kUnlimited ? Json::Value("unlimited") : uint_json(value);
    };
    j["branching"]["max_label_branches"] = limit_json(st.max_label_branches);
    j["branching"]["max_splits_per_path"] = limit_json(st.max_splits_per_path);
    j["bounds"]["max_extension_bp"] = uint_json(st.max_extension_bp);
    j["bounds"]["time_budget_ms"] = st.time_budget_ms;
    return j;
}

Json::Value arm_to_json(const ArmResult &arm) {
    Json::Value j;
    j["status"] = status_name(arm.status);
    j["complete_to_bp"] = uint_json(arm.complete_to_bp);
    j["frontier_remaining"]["live_paths"] = uint_json(arm.frontier_live_paths);
    j["frontier_remaining"]["live_labels"] = uint_json(arm.frontier_live_labels);
    if (arm.cap_trigger) {
        j["cap_trigger"]["reason"] = to_string(arm.cap_trigger->reason);
        j["cap_trigger"]["at_bp"] = uint_json(arm.cap_trigger->at_bp);
        j["cap_trigger"]["demand"] = arm.cap_trigger->demand;
        j["cap_trigger"]["segment"] = uint_json(arm.cap_trigger->segment);
    }

    Json::Value segments(Json::arrayValue);
    for (const Segment &s : arm.segments) {
        Json::Value sj;
        sj["id"] = uint_json(s.id);
        Json::Value parents(Json::arrayValue);
        for (size_t p : s.parents) {
            parents.append(uint_json(p));
        }
        sj["parents"] = std::move(parents);
        sj["from_bp"] = uint_json(s.from_bp);
        sj["length_bp"] = uint_json(s.length_bp);
        sj["sequence"] = s.sequence;
        sj["labels"] = labels_json(s.labels_start);
        sj["labels_at_end"] = labels_json(s.labels_end);
        Json::Value events(Json::arrayValue);
        for (const LabelEnd &e : s.label_ends) {
            Json::Value ej;
            ej["at_bp"] = uint_json(e.at_bp);
            ej["type"] = "label_end";
            ej["label"] = e.label;
            ej["reason"] = reason_name(e.reason, e.text);
            if (e.successors)
                ej["successors"] = uint_json(*e.successors);
            events.append(std::move(ej));
        }
        sj["events"] = std::move(events);
        segments.append(std::move(sj));
    }
    j["segments"] = std::move(segments);

    Json::Value splits(Json::arrayValue);
    for (const Split &s : arm.splits) {
        Json::Value sj;
        sj["at_bp"] = uint_json(s.at_bp);
        sj["segment"] = uint_json(s.segment);
        Json::Value children(Json::arrayValue);
        for (size_t c : s.children) {
            children.append(uint_json(c));
        }
        sj["children"] = std::move(children);
        sj["kind"] = s.ambiguous ? "ambiguous" : "divergence";
        splits.append(std::move(sj));
    }
    j["splits"] = std::move(splits);

    Json::Value paths(Json::arrayValue);
    for (const PathResult &p : arm.paths) {
        Json::Value pj;
        pj["id"] = uint_json(p.id);
        // the path's segments from the root to its leaf
        std::vector<size_t> chain;
        for (size_t s = p.leaf; ; s = arm.segments[s].parents[0]) {
            chain.push_back(s);
            if (arm.segments[s].parents.empty())
                break;
        }
        Json::Value segs(Json::arrayValue);
        for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
            segs.append(uint_json(*it));
        }
        pj["segments"] = std::move(segs);
        pj["length_bp"] = uint_json(p.length_bp);
        if (p.path_reason)
            pj["path_reason"] = to_string(*p.path_reason);
        Json::Value ends(Json::arrayValue);
        for (const PathEnd &e : p.end_labels) {
            Json::Value ej;
            ej["label"] = e.label;
            ej["reason"] = to_string(e.reason);
            ej["branches"] = e.branches;
            ends.append(std::move(ej));
        }
        pj["end_labels"] = std::move(ends);
        paths.append(std::move(pj));
    }
    j["paths"] = std::move(paths);
    j["steps"] = uint_json(arm.steps);
    return j;
}

Json::Value seed_result_to_json(const SeedResult &r) {
    Json::Value j;
    Json::Value &seed = j["seed"];
    seed["length_bp"] = uint_json(r.length_bp);
    seed["num_kmers"] = uint_json(r.num_kmers);
    Json::Value labels(Json::arrayValue);
    for (const LabelRef &l : r.label_dict) {
        labels.append(l.name);
    }
    seed["labels"] = std::move(labels);
    seed["labels_from_seed"] = r.labels_from_seed;
    seed["labels_supporting_total"] = uint_json(r.labels_supporting_total);
    seed["labels_dropped"] = uint_json(r.labels_dropped);
    Json::Value dropped(Json::arrayValue);
    for (const DroppedLabel &d : r.dropped_labels) {
        Json::Value dj;
        dj["label"] = d.name;
        dj["runs"] = runs_json(d.runs);
        dropped.append(std::move(dj));
    }
    seed["dropped_labels"] = std::move(dropped);

    Json::Value dict(Json::arrayValue);
    for (const LabelRef &l : r.label_dict) {
        Json::Value lj;
        lj["name"] = l.name;
        lj["kind"] = to_string(l.kind);
        lj["column"] = uint_json(l.column);
        if (l.kind == LabelKind::HEADER)
            lj["seq_id"] = uint_json(l.seq_id);
        dict.append(std::move(lj));
    }
    j["label_dict"] = std::move(dict);

    Json::Value arms(Json::objectValue);
    for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
        const ArmResult &a = r.arms[static_cast<size_t>(arm)];
        if (a.requested)
            arms[to_string(arm)] = arm_to_json(a);
    }
    j["arms"] = std::move(arms);
    j["elapsed_ms"] = r.elapsed_seconds * 1000;
    return j;
}

Json::Value failed_seed_to_json(const Seed &seed, const SeedDerivationError &e) {
    Json::Value j;
    j["seed"]["length_bp"] = uint_json(seed.sequence.size());
    j["seed"]["labels_from_seed"] = seed.labels.empty();
    j["error"]["cause"] = to_string(e.cause());
    j["error"]["message"] = e.what();
    j["error"]["limit"] = e.limit();
    j["error"]["observed"] = e.observed();
    if (!e.subject().empty())
        j["error"]["subject"] = e.subject();
    return j;
}


} // namespace


Json::Value process_traverse_request(const Json::Value &request,
                                     const AnnotatedDBG &anno_graph,
                                     const TraverseLimits &limits) {
    const TraverseRequest req = parse_traverse_request(request, limits);
    const LabelOracle oracle(anno_graph);
    // refuse a malformed seed before walking any seed
    for (const Seed &seed : req.seeds) {
        check_seed(oracle, seed);
    }
    Json::Value out;
    out["strategy"] = strategy_to_json(req.strategy);
    Json::Value results(Json::arrayValue);
    for (const Seed &seed : req.seeds) {
        try {
            results.append(seed_result_to_json(traverse_seed(oracle, seed, req.strategy)));
        } catch (const SeedDerivationError &e) {
            // a property of the seed and the data: the other seeds are still walked
            results.append(failed_seed_to_json(seed, e));
        }
    }
    out["results"] = std::move(results);
    return out;
}

namespace {

// What one index of a multi-graph server gives for a request
struct IndexWalks {
    // per seed: its walk or its error, or null when the index gives nothing for the seed
    std::vector<Json::Value> walks;
    // per seed and named label: whether the index has the label
    std::vector<std::vector<bool>> known;
    // per seed: whether the seed is shorter than the index's k
    std::vector<bool> too_short;
};

IndexWalks walk_index(const TraverseRequest &req, const TraverseIndex &index) {
    const LabelOracle oracle(*index.anno_graph);
    IndexWalks out;
    for (const Seed &named : req.seeds) {
        out.known.emplace_back();
        Json::Value &walk = out.walks.emplace_back();
        // the index follows the labels it has
        Seed seed;
        seed.sequence = named.sequence;
        for (const std::string &name : named.labels) {
            out.known.back().push_back(oracle.find_column(name) || oracle.find_header(name));
            if (out.known.back().back())
                seed.labels.push_back(name);
        }
        out.too_short.push_back(seed.sequence.size() < oracle.get_k());
        // an index without the seed or its labels gives nothing, like a search without hits
        if ((named.labels.size() && seed.labels.empty()) || out.too_short.back())
            continue;
        try {
            // a malformed seed or an ambiguous label refuses the request
            check_seed(oracle, seed);
            walk = seed_result_to_json(traverse_seed(oracle, seed, req.strategy));
        } catch (const SeedNotInGraph &) {
            continue;
        } catch (const SeedDerivationError &e) {
            if (e.cause() == SeedDerivationError::NO_CARRIER)
                continue;
            walk = failed_seed_to_json(seed, e);
        }
        walk["graph"] = index.graph;
    }
    return out;
}

} // namespace

Json::Value process_traverse_request(const Json::Value &request,
                                     const std::vector<TraverseIndex> &indexes,
                                     const TraverseLimits &limits,
                                     ThreadPool &pool) {
    // the server has read the selection
    Json::Value fields = request;
    if (fields.isObject())
        fields.removeMember("graphs");
    const TraverseRequest req = parse_traverse_request(fields, limits);

    std::vector<IndexWalks> walked(indexes.size());
    // the tasks refer to this frame: they are waited for before it is left, also on a throw
    struct Pending {
        std::vector<std::shared_future<std::exception_ptr>> futures;
        ~Pending() {
            for (auto &future : futures)
                future.wait();
        }
    } pending;
    auto &futures = pending.futures;
    // once a task runs, recording its future must not allocate
    futures.reserve(indexes.size());
    for (size_t i = 0; i < indexes.size(); ++i) {
        futures.push_back(pool.enqueue([&, i]() {
            try {
                walked[i] = walk_index(req, indexes[i]);
            } catch (...) {
                return std::current_exception();
            }
            return std::exception_ptr();
        }));
    }
    std::exception_ptr failed;
    for (auto &future : futures) {
        if (auto e = future.get(); e && !failed)
            failed = e;
    }
    if (failed)
        std::rethrow_exception(failed);

    Json::Value out;
    out["strategy"] = strategy_to_json(req.strategy);
    Json::Value results(Json::arrayValue);
    for (size_t s = 0; s < req.seeds.size(); ++s) {
        if (std::all_of(walked.begin(), walked.end(),
                        [&](const IndexWalks &w) { return w.too_short[s]; }))
            throw std::invalid_argument("Seed shorter than k");
        const std::vector<std::string> &labels = req.seeds[s].labels;
        for (size_t l = 0; l < labels.size(); ++l) {
            if (std::none_of(walked.begin(), walked.end(),
                             [&](const IndexWalks &w) { return w.known[s][l]; })) {
                throw std::invalid_argument("Label not found in the selected graphs: '"
                                            + labels[l] + "'");
            }
        }
        Json::Value walks(Json::arrayValue);
        for (IndexWalks &w : walked) {
            if (!w.walks[s].isNull())
                walks.append(std::move(w.walks[s]));
        }
        results.append(Json::objectValue)["walks"] = std::move(walks);
    }
    out["results"] = std::move(results);
    return out;
}

} // namespace cli
} // namespace mtg
