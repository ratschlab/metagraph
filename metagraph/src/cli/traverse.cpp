#include "traverse.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"
#include "config/config.hpp"
#include "load/load_annotated_graph.hpp"
#include "server_utils.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;
using namespace mtg::graph::traversal;

TraverseLimits traverse_limits(const Config &config) {
    TraverseLimits limits;
    limits.max_time_ms = config.traverse_max_time_ms;
    limits.max_seeds = config.traverse_max_seeds;
    limits.max_seed_bp = config.traverse_max_seed_bp;
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

    bool boolean(const char *key, bool def) {
        if (!has(key))
            return def;
        if (!value_[key].isBool())
            throw InvalidRequest(child(key) + ": expected true or false");
        return value_[key].asBool();
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
                throw InvalidRequest(child(name.c_str()) + ": unknown field (traverse.md lists "
                                     "the fields of the v1 schema)");
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
        seed.seed_id = ss.str("seed_id", "");
        seed.sequence = ss.str("sequence", "");
        if (seed.sequence.empty())
            throw InvalidRequest(ss.child("sequence") + " is required");
        if (seed.sequence.size() > limits.max_seed_bp) {
            throw InvalidRequest(ss.child("sequence") + ": " + std::to_string(seed.sequence.size())
                                 + " bp, over the server's limit of "
                                 + std::to_string(limits.max_seed_bp));
        }
        // omitted: derive the permitted set from the seed; an empty array reads like "no
        // label permitted" and is refused
        if (ss.has("labels")) {
            seed.labels = ss.strings("labels");
            if (seed.labels.empty()) {
                throw InvalidRequest(ss.child("labels") + ": expected a non-empty array of "
                                     "labels; omit the field to derive the labels from the seed");
            }
        }
        ss.done();
        req.seeds.push_back(std::move(seed));
    }

    Strategy &st = req.strategy;
    st.time_budget_ms = limits.max_time_ms;
    if (s.has("strategy")) {
        Strict t(s.raw("strategy"), "request.strategy");
        if (t.uint("schema_version", 1) != 1)
            throw InvalidRequest("request.strategy.schema_version: only version 1 is supported");
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
            const Json::Value &v = b.raw("max_label_branches");
            if (v.isString() && v.asString() == "unlimited") {
                st.max_label_branches = Strategy::kUnlimited;
            } else if (!v.isNull()) {
                st.max_label_branches = b.uint("max_label_branches", 0, 0, 1'000'000);
            }
            b.done();
        }
        if (t.has("bounds")) {
            Strict b(t.raw("bounds"), "request.strategy.bounds");
            st.max_extension_bp = b.uint("max_extension_bp", st.max_extension_bp, 1);
            // the server's bound clamps the request's
            st.time_budget_ms = std::min(b.number("time_budget_ms", limits.max_time_ms, 0,
                                                  std::numeric_limits<double>::max()),
                                         limits.max_time_ms);
            b.done();
        }
        if (t.has("output")) {
            Strict o(t.raw("output"), "request.strategy.output");
            st.sequences = o.boolean("sequences", true);
            o.done();
        }
        t.done();
    }
    s.done();
    return req;
}

Json::Value strategy_to_json(const Strategy &st) {
    Json::Value j;
    j["schema_version"] = 1;
    j["direction"] = st.direction == Strategy::BOTH ? "both"
                   : st.direction == Strategy::LEFT ? "left" : "right";
    j["labels"]["max_seed_labels"] = uint_json(st.max_seed_labels);
    j["labels"]["seed_label_kind"] = st.seed_label_kind ? to_string(*st.seed_label_kind) : "auto";
    j["branching"]["max_label_branches"] = st.max_label_branches == Strategy::kUnlimited
        ? Json::Value("unlimited") : uint_json(st.max_label_branches);
    j["bounds"]["max_extension_bp"] = uint_json(st.max_extension_bp);
    j["bounds"]["time_budget_ms"] = st.time_budget_ms;
    j["output"]["sequences"] = st.sequences;
    return j;
}

Json::Value arm_to_json(const ArmResult &arm, const Strategy &st) {
    Json::Value j;
    j["status"] = status_name(arm.status);
    j["complete_to_bp"] = uint_json(arm.complete_to_bp);
    j["frontier_remaining"]["live_paths"] = uint_json(arm.frontier_live_paths);
    if (arm.cap_trigger) {
        j["cap_trigger"]["reason"] = to_string(arm.cap_trigger->reason);
        j["cap_trigger"]["at_bp"] = uint_json(arm.cap_trigger->at_bp);
        j["cap_trigger"]["demand"] = arm.cap_trigger->demand;
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
        if (st.sequences)
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

Json::Value seed_result_to_json(const SeedResult &r, const Strategy &st) {
    Json::Value j;
    Json::Value &seed = j["seed"];
    seed["seed_id"] = r.seed_id;
    seed["validated_seed_id"] = r.validated_seed_id;
    seed["seed_id_mismatch"] = r.seed_id_mismatch;
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
    seed["labels_dropped_digest"] = r.labels_dropped_digest;
    Json::Value dropped(Json::arrayValue);
    for (const DroppedLabel &d : r.dropped_labels) {
        Json::Value dj;
        dj["label"] = d.name;
        Json::Value runs(Json::arrayValue);
        for (const auto &[a, b] : d.runs) {
            runs.append(interval_json(a, b));
        }
        dj["runs"] = std::move(runs);
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
            arms[to_string(arm)] = arm_to_json(a, st);
    }
    j["arms"] = std::move(arms);
    j["elapsed_ms"] = r.elapsed_seconds * 1000;
    return j;
}

Json::Value failed_seed_to_json(const Seed &seed, const SeedDerivationError &e) {
    Json::Value j;
    j["seed"]["seed_id"] = seed.seed_id;
    j["seed"]["length_bp"] = uint_json(seed.sequence.size());
    j["seed"]["labels_from_seed"] = true;
    j["error"]["cause"] = to_string(e.cause());
    j["error"]["message"] = e.what();
    j["error"]["limit"] = e.limit();
    j["error"]["observed"] = e.observed();
    if (!e.subject().empty())
        j["error"]["subject"] = e.subject();
    return j;
}


/********************************* /resolve **********************************/

struct ResolveRequest {
    std::string sequence;
    ResolveOptions options;
    bool select = false;
    SelectionPolicy policy;
};

ResolveRequest parse_resolve_request(const Json::Value &json) {
    ResolveRequest req;
    Strict s(json, "request");
    req.sequence = s.str("sequence", "");
    if (req.sequence.empty())
        throw InvalidRequest("request.sequence is required");
    if (req.sequence.size() > kMaxResolveQueryBp) {
        throw InvalidRequest("request.sequence: " + std::to_string(req.sequence.size())
                             + " bp, over the limit of " + std::to_string(kMaxResolveQueryBp));
    }
    req.options.labels = s.strings("labels");
    if (req.options.labels.size() > kMaxSeedLabels) {
        throw InvalidRequest("request.labels: " + std::to_string(req.options.labels.size())
                             + " labels, over the limit of " + std::to_string(kMaxSeedLabels));
    }
    if (s.has("discover")) {
        Strict d(s.raw("discover"), "request.discover");
        req.options.discover = true;
        req.options.discover_max_labels = d.uint("max_labels", 1000, 1, kMaxSeedLabels);
        req.options.discover_kind = d.enumeration<LabelKind>("kind", LabelKind::COLUMN,
            { { "column", LabelKind::COLUMN }, { "header", LabelKind::HEADER } });
        d.done();
    }
    if (req.options.labels.empty() == !req.options.discover)
        throw InvalidRequest("request: specify exactly one of 'labels' or 'discover'");
    req.options.min_block_kmers = s.uint("min_block_kmers", 1, 1);
    if (s.has("select")) {
        Strict p(s.raw("select"), "request.select");
        req.select = true;
        req.policy.policy = p.enumeration<SelectionPolicy::Policy>("policy", SelectionPolicy::MAX_SUPPORT,
            { { "max_support", SelectionPolicy::MAX_SUPPORT },
              { "longest_first", SelectionPolicy::LONGEST_FIRST } });
        req.policy.max_seeds = p.uint("max_seeds", 10, 1, 1000);
        req.policy.min_block_bp = p.uint("min_block_bp", 0);
        req.policy.max_labels_per_seed = p.uint("max_labels_per_seed", 1000, 1, kMaxSeedLabels);
        req.policy.merge_overlapping = p.boolean("merge_overlapping", true);
        p.done();
    }
    s.done();
    return req;
}

Json::Value profile_to_json(const SupportProfile &p, const SeedSelection *sel) {
    Json::Value j;
    j["k"] = uint_json(p.k);
    j["regime"] = to_string(p.regime);
    j["num_kmers"] = uint_json(p.num_kmers);
    j["graph_runs"] = runs_json(p.graph_runs);
    Json::Value labels(Json::arrayValue);
    for (const LabelProfile &l : p.labels) {
        Json::Value lj;
        lj["label"] = l.label.name;
        lj["kind"] = to_string(l.label.kind);
        lj["kmers_supported"] = uint_json(l.kmers_supported);
        lj["runs"] = runs_json(l.runs);
        labels.append(std::move(lj));
    }
    j["labels"] = std::move(labels);
    if (p.labels_truncated) {
        Json::Value &t = j["labels_truncated"];
        t["kept"] = uint_json(p.labels_truncated->kept);
        t["total"] = uint_json(p.labels_truncated->total);
        t["min_kept_kmers"] = uint_json(p.labels_truncated->min_kept_kmers);
        t["max_dropped_kmers"] = uint_json(p.labels_truncated->max_dropped_kmers);
        t["dropped_full_length"] = uint_json(p.labels_truncated->dropped_full_length);
    } else {
        j["labels_truncated"] = Json::Value::null;
    }
    Json::Value cands(Json::arrayValue);
    for (const SeedCandidate &c : p.candidates) {
        Json::Value cj;
        cj["kmer_interval"] = interval_json(c.kmers.begin, c.kmers.end);
        const auto [a, b] = p.bp_interval(c.kmers);
        cj["bp_interval"] = interval_json(a, b);
        Json::Value ls(Json::arrayValue);
        for (LabelId l : c.labels) {
            ls.append(p.labels[l].label.name);
        }
        cj["labels"] = std::move(ls);
        cands.append(std::move(cj));
    }
    j["candidates"] = std::move(cands);
    if (sel) {
        Json::Value &s = j["selection"];
        s["policy"]["policy"] = sel->policy.policy == SelectionPolicy::LONGEST_FIRST
                                ? "longest_first" : "max_support";
        s["policy"]["max_seeds"] = uint_json(sel->policy.max_seeds);
        s["policy"]["min_block_bp"] = uint_json(sel->policy.min_block_bp);
        s["policy"]["max_labels_per_seed"] = uint_json(sel->policy.max_labels_per_seed);
        s["policy"]["merge_overlapping"] = sel->policy.merge_overlapping;
        s["counts"]["candidates"] = uint_json(sel->num_candidates);
        s["counts"]["eligible"] = uint_json(sel->num_eligible);
        s["counts"]["selected"] = uint_json(sel->seeds.size());
        Json::Value seeds(Json::arrayValue);
        for (const FrozenSeed &seed : sel->seeds) {
            Json::Value sj;
            sj["seed_id"] = seed.seed_id;
            sj["sequence"] = seed.sequence;
            sj["kmer_interval"] = interval_json(seed.kmers.begin, seed.kmers.end);
            Json::Value ls(Json::arrayValue);
            for (const std::string &l : seed.labels) {
                ls.append(l);
            }
            sj["labels"] = std::move(ls);
            Json::Value &pop = sj["label_population"];
            pop["supporting_total"] = uint_json(seed.population.supporting_total);
            pop["included"] = uint_json(seed.population.included);
            pop["dropped_count"] = uint_json(seed.population.dropped_count);
            if (seed.population.dropped_count) {
                pop["dropped_digest"] = seed.population.dropped_digest;
                Json::Value dl(Json::arrayValue);
                for (const std::string &d : seed.population.dropped) {
                    dl.append(d);
                }
                pop["dropped"] = std::move(dl);
            }
            Json::Value ov(Json::arrayValue);
            for (size_t o : seed.overlaps_with) {
                ov.append(uint_json(o));
            }
            sj["overlaps_with"] = std::move(ov);
            seeds.append(std::move(sj));
        }
        s["seeds"] = std::move(seeds);
    }
    return j;
}

} // namespace


Json::Value process_traverse_request(const Json::Value &request,
                                     const AnnotatedDBG &anno_graph,
                                     const TraverseLimits &limits) {
    const TraverseRequest req = parse_traverse_request(request, limits);
    const LabelOracle oracle(anno_graph);
    Json::Value out;
    out["schema_version"] = 1;
    out["strategy"] = strategy_to_json(req.strategy);
    Json::Value results(Json::arrayValue);
    for (const Seed &seed : req.seeds) {
        try {
            results.append(seed_result_to_json(traverse_seed(oracle, seed, req.strategy),
                                               req.strategy));
        } catch (const SeedDerivationError &e) {
            // a property of the seed and the data: the other seeds are still walked
            results.append(failed_seed_to_json(seed, e));
        }
    }
    out["results"] = std::move(results);
    return out;
}

Json::Value process_resolve_request(const Json::Value &request,
                                    const AnnotatedDBG &anno_graph,
                                    const TraverseLimits &) {
    const ResolveRequest req = parse_resolve_request(request);
    const LabelOracle oracle(anno_graph);
    const SupportProfile profile = resolve_support(oracle, req.sequence, req.options);
    if (!req.select)
        return profile_to_json(profile, nullptr);
    const SeedSelection selection = select_seeds(profile, req.sequence, req.policy,
                                                 oracle.regime() != Regime::BASIC);
    return profile_to_json(profile, &selection);
}

Json::Value traverse_capabilities(const AnnotatedDBG &anno_graph, const TraverseLimits &limits) {
    const LabelOracle oracle(anno_graph);
    Json::Value c;
    c["schema_version"] = 1;
    c["k"] = uint_json(oracle.get_k());
    c["regime"] = to_string(oracle.regime());
    c["alphabet"] = oracle.graph().alphabet();
    c["num_columns"] = uint_json(oracle.num_columns());
    c["has_coordinates"] = oracle.has_coordinates();
    c["has_coord_to_header"] = oracle.coord_to_header() != nullptr;
    c["limits"]["max_time_ms"] = limits.max_time_ms;
    c["limits"]["max_seeds"] = uint_json(limits.max_seeds);
    c["limits"]["max_seed_bp"] = uint_json(limits.max_seed_bp);
    c["limits"]["max_seed_labels"] = uint_json(kMaxSeedLabels);
    c["limits"]["max_resolve_query_bp"] = uint_json(kMaxResolveQueryBp);
    return c;
}

int traverse_graph(Config *config) {
    assert(config);
    assert(config->fnames.size() == 1);

    std::ifstream in(config->fnames[0]);
    if (!in.is_open()) {
        logger->error("Cannot open request file {}", config->fnames[0]);
        return 1;
    }
    std::stringstream content;
    content << in.rdbuf();

    Timer timer;
    std::unique_ptr<AnnotatedDBG> anno_graph = initialize_annotated_dbg(*config);
    logger->trace("Annotated graph loaded in {} sec", timer.elapsed());

    try {
        const Json::Value request = parse_json_string(content.str());
        const TraverseLimits limits = traverse_limits(*config);
        const Json::Value response = config->traverse_resolve
            ? process_resolve_request(request, *anno_graph, limits)
            : process_traverse_request(request, *anno_graph, limits);
        Json::StreamWriterBuilder builder;
        builder["indentation"] = "";
        std::cout << Json::writeString(builder, response) << std::endl;
    } catch (const std::invalid_argument &e) {
        logger->error("Invalid request: {}", e.what());
        return 1;
    }
    return 0;
}

} // namespace cli
} // namespace mtg
