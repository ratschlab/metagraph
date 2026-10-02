#include "traverse.hpp"

#include <algorithm>
#include <fstream>
#include <optional>
#include <set>

#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph::traversal;

namespace {

// ---------------------------------------------------------------- strict JSON access

class Strict {
  public:
    Strict(const Json::Value &v, std::string path) : v_(v), path_(std::move(path)) {
        if (!v_.isObject())
            throw InvalidRequest(path_ + ": expected an object");
    }
    ~Strict() noexcept(false) {
        if (std::uncaught_exceptions())
            return;
        for (const auto &name : v_.getMemberNames()) {
            if (!seen_.count(name))
                throw InvalidRequest(path_ + ": unknown field '" + name + "'");
        }
    }
    bool has(const std::string &k) {
        seen_.insert(k);
        return v_.isMember(k);
    }
    const Json::Value& raw(const std::string &k) { seen_.insert(k); return v_[k]; }
    std::string child_path(const std::string &k) const { return path_ + "." + k; }

    std::string str(const std::string &k, const std::string &def) {
        if (!has(k)) return def;
        if (!v_[k].isString()) throw InvalidRequest(child_path(k) + ": expected a string");
        return v_[k].asString();
    }
    bool boolean(const std::string &k, bool def) {
        if (!has(k)) return def;
        if (!v_[k].isBool()) throw InvalidRequest(child_path(k) + ": expected a boolean");
        return v_[k].asBool();
    }
    uint64_t uint(const std::string &k, uint64_t def, uint64_t min = 0,
                  uint64_t max = std::numeric_limits<uint64_t>::max()) {
        if (!has(k)) return def;
        if (!v_[k].isIntegral() || (v_[k].isInt64() && v_[k].asInt64() < 0))
            throw InvalidRequest(child_path(k) + ": expected a non-negative integer");
        uint64_t x = v_[k].asUInt64();
        if (x < min || x > max)
            throw InvalidRequest(child_path(k) + ": out of range [" + std::to_string(min) + ", "
                                 + std::to_string(max) + "]");
        return x;
    }
    double number(const std::string &k, double def, double min = 0,
                  double max = std::numeric_limits<double>::infinity()) {
        if (!has(k)) return def;
        if (!v_[k].isNumeric()) throw InvalidRequest(child_path(k) + ": expected a number");
        double x = v_[k].asDouble();
        if (!(x >= min && x <= max))
            throw InvalidRequest(child_path(k) + ": out of range");
        return x;
    }
    std::vector<std::string> strings(const std::string &k) {
        std::vector<std::string> out;
        if (!has(k)) return out;
        if (!v_[k].isArray()) throw InvalidRequest(child_path(k) + ": expected an array of strings");
        for (const auto &x : v_[k]) {
            if (!x.isString()) throw InvalidRequest(child_path(k) + ": expected an array of strings");
            out.push_back(x.asString());
        }
        return out;
    }
    template <class E>
    E enumeration(const std::string &k, E def, const std::vector<std::pair<std::string, E>> &values) {
        if (!has(k)) return def;
        std::string s = str(k, "");
        std::string allowed;
        for (const auto &[name, val] : values) {
            if (name == s) return val;
            allowed += (allowed.empty() ? "" : "|") + name;
        }
        throw InvalidRequest(child_path(k) + ": expected one of " + allowed);
    }

  private:
    const Json::Value &v_;
    std::string path_;
    std::set<std::string> seen_;
};

Json::Value uint_json(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

Json::Value labels_json(const std::vector<LabelId> &labels) {
    Json::Value arr(Json::arrayValue);
    for (LabelId l : labels) arr.append(l);
    return arr;
}

// A limit that may also be "unlimited" (Strategy::kUnlimited). |def| is the value for
// an omitted field; a mode in which the knob has no meaning passes |only_unlimited|,
// and a number is then rejected rather than accepted and ignored.
size_t limit_or_unlimited(Strict &s, const std::string &k, size_t def, size_t max,
                          const char *only_unlimited) {
    if (!s.has(k))
        return def;
    const Json::Value &v = s.raw(k);
    if (v.isString()) {
        if (v.asString() == "unlimited")
            return Strategy::kUnlimited;
        throw InvalidRequest(s.child_path(k) + ": expected an integer or \"unlimited\"");
    }
    if (only_unlimited) {
        throw InvalidRequest(s.child_path(k) + ": " + only_unlimited
                             + "; set it to \"unlimited\" or omit it");
    }
    return s.uint(k, def, 0, max);
}

Json::Value limit_json(size_t x) {
    return x == Strategy::kUnlimited ? Json::Value("unlimited") : uint_json(x);
}

// the cost model as the walker sees it, for validating a strategy before any seed
LabelChangeCost representative_cost(const CostSpec &spec) {
    switch (spec.model) {
        case CostSpec::FORBID: return LabelChangeCost::forbid();
        case CostSpec::CONSTANT: return LabelChangeCost::constant(spec.value);
        case CostSpec::TABLE: return LabelChangeCost::table({}, spec.default_cost);
    }
    return LabelChangeCost::forbid();
}

} // namespace


// ---------------------------------------------------------------- parsing

static CostSpec parse_cost(const Json::Value &json, const std::string &path) {
    CostSpec cost;
    Strict s(json, path);
    cost.model = s.enumeration<CostSpec::Model>("model", CostSpec::FORBID,
        { { "forbid", CostSpec::FORBID }, { "constant", CostSpec::CONSTANT }, { "table", CostSpec::TABLE } });
    switch (cost.model) {
        case CostSpec::FORBID:
            break;
        case CostSpec::CONSTANT:
            if (!s.has("value")) throw InvalidRequest(path + ".value is required for model constant");
            cost.value = s.number("value", 0);
            break;
        case CostSpec::TABLE: {
            if (s.has("default")) {
                const auto &d = s.raw("default");
                if (d.isString() && d.asString() == "forbid") {
                    cost.default_cost = kInfiniteLoss;
                } else if (d.isNumeric() && d.asDouble() >= 0) {
                    cost.default_cost = d.asDouble();
                } else {
                    throw InvalidRequest(path + ".default: expected \"forbid\" or a non-negative number");
                }
            }
            const auto &entries = s.raw("entries");
            if (!entries.isArray()) throw InvalidRequest(path + ".entries: expected an array of [from, to, cost]");
            for (const auto &e : entries) {
                if (!e.isArray() || e.size() != 3 || !e[0].isString() || !e[1].isString() || !e[2].isNumeric()
                        || e[2].asDouble() < 0)
                    throw InvalidRequest(path + ".entries: each entry is [from, to, cost >= 0]");
                cost.entries.emplace_back(e[0].asString(), e[1].asString(), e[2].asDouble());
            }
            break;
        }
    }
    return cost;
}

TraverseRequest parse_traverse_request(const Json::Value &json) {
    TraverseRequest req;
    Strict s(json, "request");
    req.release = s.str("release", "");
    // Routing fields: consumed by the server before parsing (resolve_traverse_index).
    // They must still be declared here, or the unknown-field check would reject every
    // multi-graph request.
    req.graph = s.str("graph", "");
    req.graph_path = s.str("graph_path", "");

    const auto &seeds = s.raw("seeds");
    if (!seeds.isArray() || seeds.empty())
        throw InvalidRequest("request.seeds: expected a non-empty array");
    for (Json::ArrayIndex i = 0; i < seeds.size(); ++i) {
        Strict ss(seeds[i], "request.seeds[" + std::to_string(i) + "]");
        Seed seed;
        seed.seed_id = ss.str("seed_id", "");
        seed.sequence = ss.str("sequence", "");
        if (seed.sequence.empty()) throw InvalidRequest(ss.child_path("sequence") + " is required");
        // `labels` is optional: omitting it derives the permitted set from the seed
        // itself (spec §6.1, the design note's `permit: per_hit` default). An
        // explicitly empty array is an error — it reads like "no labels permitted".
        if (ss.has("labels")) {
            seed.labels = ss.strings("labels");
            if (seed.labels.empty()) {
                throw InvalidRequest(ss.child_path("labels") + ": expected a non-empty array of "
                                     "labels; omit the field to derive the permitted set from "
                                     "the seed");
            }
        }
        req.seeds.push_back(std::move(seed));
    }

    if (!s.has("strategy")) throw InvalidRequest("request.strategy is required");
    Strategy &st = req.strategy;
    {
        Strict t(s.raw("strategy"), "strategy");
        uint64_t version = t.uint("schema_version", 1);
        if (version != 1) throw InvalidRequest("strategy.schema_version: only version 1 is supported");
        st.direction = t.enumeration<Strategy::Direction>("direction", Strategy::BOTH,
            { { "both", Strategy::BOTH }, { "left", Strategy::LEFT }, { "right", Strategy::RIGHT } });
        st.support = t.enumeration<Support>("support", Support::KMER,
            { { "kmer", Support::KMER }, { "trace", Support::TRACE } });
        // The preset (spec §6.9) decides the DEFAULTS of the knobs it governs; a knob
        // given explicitly must agree with it, which validate_strategy() checks below.
        st.exhaustive = t.boolean("exhaustive", false);

        if (t.has("labels")) {
            Strict l(t.raw("labels"), "strategy.labels");
            st.label_mode = l.enumeration<LabelMode>("mode", LabelMode::CONSTRAIN,
                { { "constrain", LabelMode::CONSTRAIN }, { "annotate", LabelMode::ANNOTATE } });
            // bounded like its neighbours: a seed of exactly k bases has no intersection
            // to narrow the candidates, so this knob alone decides how large a label
            // dictionary (and with it scratch, stamps, summary and hit maps) the server
            // materialises. The server lowers it further via TraverseLimits.
            st.max_seed_labels = l.uint("max_seed_labels", 1000, 1, 100'000);
            // In annotate mode the per-node cap must admit what the constrain side
            // admits as a derived permitted set, or the structural oracle of a locus
            // with more than 64 labels cuts every node by default and proves nothing.
            const size_t default_cap = st.label_mode == LabelMode::ANNOTATE
                ? std::max<size_t>(64, st.max_seed_labels) : 64;
            st.max_labels_per_node = l.uint("max_labels_per_node", default_cap, 1, 100'000);
            st.extra = l.strings("extra");
            if (l.has("change_cost"))
                req.cost = parse_cost(l.raw("change_cost"), "strategy.labels.change_cost");
            st.loss_budget = l.number("loss_budget", 0);
            st.switch_on_loss_only = l.enumeration<bool>("switch_on", true,
                { { "loss", true }, { "any", false } });
            // unlimited under the preset (a cut source list prunes, see
            // validate_strategy); a number given with it is left for the validation
            // to refuse under a table cost
            st.max_switch_sources = limit_or_unlimited(
                    l, "max_switch_sources", st.exhaustive ? Strategy::kUnlimited : 64,
                    std::numeric_limits<size_t>::max() - 1, nullptr);
            if (st.max_switch_sources == 0)
                throw InvalidRequest("strategy.labels.max_switch_sources: out of range [1, ...]");
            // "auto" (no kind pinned) is what the normalized strategy echoes, so it has
            // to parse back: an echoed strategy must be resubmittable verbatim
            using OptionalKind = std::optional<LabelKind>;
            st.seed_label_kind = l.enumeration<OptionalKind>("seed_label_kind", OptionalKind(),
                { { "auto", OptionalKind() },
                  { "column", OptionalKind(LabelKind::COLUMN) },
                  { "header", OptionalKind(LabelKind::HEADER) } });
        }
        // the preset's default applies to an omitted labels block too
        if (!t.has("labels") && st.exhaustive)
            st.max_switch_sources = Strategy::kUnlimited;
        // Annotate mode has no lineages, so a branch limit and a split limit have no
        // meaning in it; they are "unlimited" there and a number is refused rather than
        // ignored. The exhaustive preset makes "unlimited" the default and keeps a
        // conflicting number for validate_strategy() to refuse.
        const bool annotate = st.label_mode == LabelMode::ANNOTATE;
        const char *no_lineage = annotate
            ? "has no meaning in labels.mode \"annotate\" (no label lineage is tracked, every "
              "structural successor is followed)"
            : nullptr;
        st.max_label_branches = (st.exhaustive || annotate) ? Strategy::kUnlimited : 0;
        st.max_splits_per_path = (st.exhaustive || annotate) ? Strategy::kUnlimited : 64;
        // annotate promises the whole trie, i.e. the per-path walks; merging unites edge
        // histories and drops walks (walk_rule says so when it is on), so there it is
        // opt-in rather than the default
        st.merge_reconverge = !st.exhaustive && !annotate;
        if (t.has("branching")) {
            Strict b(t.raw("branching"), "strategy.branching");
            st.max_label_branches = limit_or_unlimited(b, "max_label_branches",
                                                       st.max_label_branches, 64, no_lineage);
            st.min_successor_labels = b.uint("min_successor_labels", 1, 1);
            st.min_successor_fraction = b.number("min_successor_fraction", 0, 0, 1);
            st.tip_window_bp = b.uint("tip_window_bp", 0);
            st.bubble_window_bp = b.uint("bubble_window_bp", 0);
            st.merge_reconverge = b.enumeration<bool>("on_reconverge", st.merge_reconverge,
                { { "merge", true }, { "keep", false } });
            st.skip_hairpins = b.enumeration<bool>("hairpins", true,
                { { "skip", true }, { "follow", false } });
            st.max_splits_per_path = limit_or_unlimited(b, "max_splits_per_path",
                                                        st.max_splits_per_path,
                                                        std::numeric_limits<size_t>::max() - 1,
                                                        no_lineage);
        }
        if (t.has("bounds")) {
            Strict b(t.raw("bounds"), "strategy.bounds");
            st.max_extension_bp = b.uint("max_extension_bp", 5000, 1);
            st.min_live_labels = b.uint("min_live_labels", 1, 1);
            st.max_steps = b.uint("max_steps", 1'000'000, 1);
            st.max_live_paths = b.uint("max_live_paths", 1000, 1);
            st.max_paths = b.uint("max_paths", 10'000, 1);
            st.max_output_bp = b.uint("max_output_bp", 2'000'000, 1);
            st.time_budget_ms = b.number("time_budget_ms", 30'000);
        }
        if (t.has("frontier")) {
            Strict f(t.raw("frontier"), "strategy.frontier");
            st.order = f.enumeration<Strategy::Order>("order", Strategy::BREADTH_FIRST,
                { { "breadth_first", Strategy::BREADTH_FIRST },
                  { "lowest_loss_first", Strategy::LOWEST_LOSS_FIRST },
                  { "most_supported_first", Strategy::MOST_SUPPORTED_FIRST } });
            st.on_overflow = f.enumeration<Strategy::Overflow>("on_overflow", Strategy::STOP,
                { { "stop", Strategy::STOP }, { "beam", Strategy::BEAM } });
        }
        if (t.has("output")) {
            Strict o(t.raw("output"), "strategy.output");
            req.detail = o.enumeration<std::string>("detail", "full",
                { { "summary", "summary" }, { "tree", "tree" }, { "full", "full" } });
            st.sequences = o.boolean("sequences", true);
            st.profile_bin_bp = o.uint("profile_bin_bp", 100, 1);
            st.max_branch_events = o.uint("max_branch_events", 100);
            st.continuation_bp = o.uint("continuation_bp", 1000, 0);
            req.timing = o.boolean("timing", true);
        }
        if (t.has("annotation")) {
            Strict a(t.raw("annotation"), "strategy.annotation");
            st.batch_kmers = a.uint("batch_kmers", 64, 1, 100'000);
            // The access path is chosen from the annotation representation and the
            // label kinds (header labels require the coordinate path). Overriding it
            // is a C++-API-only facility, so reject it here rather than accept and
            // ignore it.
            const std::string access = a.str("access", "auto");
            if (access != "auto") {
                throw InvalidRequest("strategy.annotation.access: only \"auto\" is supported "
                                     "through this API (the access path follows the annotation "
                                     "representation and the label kinds)");
            }
        }
    }
    // A table is authored by label NAME over the request order (seed labels, then
    // extra). A seed whose labels are derived has no such list at request time, so the
    // entries could not be resolved — rejected rather than silently mapped. (In
    // annotate mode a table conflicts with the mode itself; that message comes first.)
    if (req.cost.model == CostSpec::TABLE && st.label_mode != LabelMode::ANNOTATE) {
        for (size_t i = 0; i < req.seeds.size(); ++i) {
            if (req.seeds[i].labels.empty()) {
                throw InvalidRequest("strategy.labels.change_cost: model \"table\" needs "
                                     "request.seeds[" + std::to_string(i) + "].labels to be given "
                                     "explicitly (table entries are resolved by name against the "
                                     "seed's label list, which is derived from the seed itself "
                                     "when 'labels' is omitted)");
            }
        }
    }
    if (req.cost.model == CostSpec::FORBID && !st.extra.empty())
        throw InvalidRequest("strategy.labels.extra: extra labels are unreachable under change_cost forbid");
    if (req.cost.model == CostSpec::CONSTANT && !st.extra.empty() && req.cost.value > st.loss_budget)
        throw InvalidRequest("strategy.labels.extra: extra labels are unreachable (cost exceeds loss_budget)");
    // the preset / mode conflicts, refused before any seed is looked at
    try {
        validate_strategy(st, representative_cost(req.cost));
    } catch (const std::invalid_argument &e) {
        throw InvalidRequest(std::string("strategy.") + e.what());
    }
    if (st.label_mode == LabelMode::ANNOTATE) {
        for (size_t i = 0; i < req.seeds.size(); ++i) {
            if (!req.seeds[i].labels.empty()) {
                throw InvalidRequest("request.seeds[" + std::to_string(i) + "].labels: labels.mode "
                                     "\"annotate\" records the labels present and filters by "
                                     "none; omit the field (there is no permitted set)");
            }
        }
    }
    return req;
}

ResolveRequest parse_resolve_request(const Json::Value &json) {
    ResolveRequest req;
    Strict s(json, "request");
    // see parse_traverse_request: routing fields are handled by the server
    req.graph = s.str("graph", "");
    req.graph_path = s.str("graph_path", "");
    req.sequence = s.str("sequence", "");
    if (req.sequence.empty()) throw InvalidRequest("request.sequence is required");
    req.options.labels = s.strings("labels");
    if (s.has("discover")) {
        Strict d(s.raw("discover"), "request.discover");
        req.options.discover = true;
        req.options.discover_max_labels = d.uint("max_labels", 1000, 1);
        req.options.discover_kind = d.enumeration<LabelKind>("kind", LabelKind::COLUMN,
            { { "column", LabelKind::COLUMN }, { "header", LabelKind::HEADER } });
    }
    if (req.options.labels.empty() == !req.options.discover)
        throw InvalidRequest("request: specify exactly one of 'labels' or 'discover'");
    req.options.support = s.enumeration<Support>("support", Support::KMER,
        { { "kmer", Support::KMER }, { "trace", Support::TRACE } });
    req.options.min_block_kmers = s.uint("min_block_kmers", 1, 1);
    req.run_format = s.enumeration<std::string>("run_format", "intervals",
        { { "intervals", "intervals" }, { "rle", "rle" } });
    if (s.has("bounds")) {
        Strict b(s.raw("bounds"), "request.bounds");
        req.max_query_bp = b.uint("max_query_bp", 1'000'000, 1);
        // resolve has no cooperative deadline yet, so accepting one would promise
        // something nothing enforces
        if (b.has("time_budget_ms")) {
            throw InvalidRequest("request.bounds.time_budget_ms: not supported for resolve "
                                 "(no deadline is enforced); bound the work with "
                                 "max_query_bp and discover.max_labels instead");
        }
    }
    if (s.has("select")) {
        Strict p(s.raw("select"), "request.select");
        req.select = true;
        req.policy.policy = p.enumeration<SelectionPolicy::Policy>("policy", SelectionPolicy::MAX_SUPPORT,
            { { "longest_first", SelectionPolicy::LONGEST_FIRST },
              { "max_support", SelectionPolicy::MAX_SUPPORT },
              { "explicit", SelectionPolicy::EXPLICIT } });
        req.policy.max_seeds = p.uint("max_seeds", 10, 1);
        req.policy.min_block_bp = p.uint("min_block_bp", 0);
        req.policy.max_labels_per_seed = p.uint("max_labels_per_seed", 1000, 1);
        req.policy.label_order = p.enumeration<SelectionPolicy::LabelOrder>("label_order", SelectionPolicy::HASH,
            { { "hash", SelectionPolicy::HASH }, { "column_id", SelectionPolicy::COLUMN_ID },
              { "kmers_supported", SelectionPolicy::KMERS_SUPPORTED } });
        req.policy.sample_seed = p.uint("sample_seed", 0);
        req.policy.merge_overlapping = p.boolean("merge_overlapping", true);
        if (p.has("seeds")) {
            const auto &ex = p.raw("seeds");
            if (!ex.isArray()) throw InvalidRequest("request.select.seeds: expected an array");
            for (Json::ArrayIndex i = 0; i < ex.size(); ++i) {
                Strict e(ex[i], "request.select.seeds[" + std::to_string(i) + "]");
                ExplicitSeed seed;
                const auto &iv = e.raw("kmer_interval");
                // isIntegral() accepts negatives, and asUInt64() on one aborts, so the
                // sign has to be checked here rather than left to the conversion
                if (!iv.isArray() || iv.size() != 2 || !iv[0].isIntegral() || !iv[1].isIntegral()
                        || iv[0].asInt64() < 0 || iv[1].asInt64() < 0) {
                    throw InvalidRequest(e.child_path("kmer_interval")
                                         + ": expected [begin, end] as non-negative integers");
                }
                seed.kmers = { iv[0].asUInt64(), iv[1].asUInt64() };
                seed.labels = e.strings("labels");
                req.policy.explicit_seeds.push_back(std::move(seed));
            }
        }
        if (req.policy.policy == SelectionPolicy::EXPLICIT && req.policy.explicit_seeds.empty())
            throw InvalidRequest("request.select.seeds is required for policy explicit");
    }
    return req;
}


// ---------------------------------------------------------------- serialization

static Json::Value cost_to_json(const CostSpec &cost) {
    Json::Value j;
    switch (cost.model) {
        case CostSpec::FORBID: j["model"] = "forbid"; break;
        case CostSpec::CONSTANT: j["model"] = "constant"; j["value"] = cost.value; break;
        case CostSpec::TABLE: {
            j["model"] = "table";
            if (cost.default_cost == kInfiniteLoss) j["default"] = "forbid"; else j["default"] = cost.default_cost;
            Json::Value entries(Json::arrayValue);
            for (const auto &[from, to, c] : cost.entries) {
                Json::Value e(Json::arrayValue);
                e.append(from); e.append(to); e.append(c);
                entries.append(e);
            }
            j["entries"] = entries;
        }
    }
    return j;
}

Json::Value strategy_to_json(const Strategy &st, const CostSpec &cost) {
    Json::Value j;
    j["schema_version"] = 1;
    j["direction"] = st.direction == Strategy::BOTH ? "both" : st.direction == Strategy::LEFT ? "left" : "right";
    j["support"] = to_string(st.support);
    j["exhaustive"] = st.exhaustive;
    Json::Value labels;
    labels["mode"] = to_string(st.label_mode);
    labels["max_labels_per_node"] = uint_json(st.max_labels_per_node);
    Json::Value extra(Json::arrayValue);
    for (const auto &e : st.extra) extra.append(e);
    labels["extra"] = extra;
    labels["change_cost"] = cost_to_json(cost);
    labels["loss_budget"] = st.loss_budget;
    labels["switch_on"] = st.switch_on_loss_only ? "loss" : "any";
    labels["max_switch_sources"] = limit_json(st.max_switch_sources);
    labels["max_seed_labels"] = uint_json(st.max_seed_labels);
    // "auto": header when the index has a CoordToHeader, column otherwise
    labels["seed_label_kind"] = st.seed_label_kind ? to_string(*st.seed_label_kind) : "auto";
    j["labels"] = labels;
    Json::Value b;
    b["max_label_branches"] = limit_json(st.max_label_branches);
    b["min_successor_labels"] = uint_json(st.min_successor_labels);
    b["min_successor_fraction"] = st.min_successor_fraction;
    b["tip_window_bp"] = uint_json(st.tip_window_bp);
    b["bubble_window_bp"] = uint_json(st.bubble_window_bp);
    b["on_reconverge"] = st.merge_reconverge ? "merge" : "keep";
    b["hairpins"] = st.skip_hairpins ? "skip" : "follow";
    b["max_splits_per_path"] = limit_json(st.max_splits_per_path);
    j["branching"] = b;
    Json::Value bo;
    bo["max_extension_bp"] = uint_json(st.max_extension_bp);
    bo["min_live_labels"] = uint_json(st.min_live_labels);
    bo["max_steps"] = uint_json(st.max_steps);
    bo["max_live_paths"] = uint_json(st.max_live_paths);
    bo["max_paths"] = uint_json(st.max_paths);
    bo["max_output_bp"] = uint_json(st.max_output_bp);
    bo["time_budget_ms"] = st.time_budget_ms;
    j["bounds"] = bo;
    Json::Value f;
    f["order"] = st.order == Strategy::BREADTH_FIRST ? "breadth_first"
               : st.order == Strategy::LOWEST_LOSS_FIRST ? "lowest_loss_first" : "most_supported_first";
    f["on_overflow"] = st.on_overflow == Strategy::STOP ? "stop" : "beam";
    j["frontier"] = f;
    Json::Value o;
    o["sequences"] = st.sequences;
    o["profile_bin_bp"] = uint_json(st.profile_bin_bp);
    o["max_branch_events"] = uint_json(st.max_branch_events);
    o["continuation_bp"] = uint_json(st.continuation_bp);
    j["output"] = o;
    Json::Value a;
    a["batch_kmers"] = uint_json(st.batch_kmers);
    j["annotation"] = a;
    return j;
}

static Json::Value event_to_json(const Event &ev) {
    Json::Value j;
    j["at_bp"] = uint_json(ev.at_bp);
    j["type"] = to_string(ev.type);
    switch (ev.type) {
        case EventType::LABEL_END:
            j["label"] = ev.label;
            j["reason"] = ev.text.empty() ? to_string(ev.reason) : ev.text.c_str();
            if (ev.reason == EndReason::LOSS_BUDGET) j["needed_budget"] = ev.needed_budget;
            j["structural_successors"] = uint_json(ev.structural_successors);
            break;
        case EventType::SWITCH:
            j["from"] = ev.label; j["to"] = ev.to; j["cost"] = ev.cost;
            break;
        case EventType::BLOCKED:
            j["char"] = std::string(1, ev.ch); j["reason"] = to_string(ev.reason); j["labels"] = labels_json(ev.labels);
            // the labels at the successor not taken appear nowhere else, so a cut list
            // here is reported like any recorded node's
            j["labels_total"] = uint_json(ev.labels_total); j["truncated"] = ev.truncated();
            break;
        case EventType::HAIRPIN:
            j["char"] = std::string(1, ev.ch); j["labels"] = labels_json(ev.labels);
            j["labels_total"] = uint_json(ev.labels_total); j["truncated"] = ev.truncated();
            break;
        case EventType::TIP:
            j["char"] = std::string(1, ev.ch); j["length_bp"] = uint_json(ev.length_bp);
            break;
        case EventType::BUBBLE:
            j["length_bp"] = uint_json(ev.length_bp); j["alleles"] = ev.text;
            break;
        case EventType::RECONVERGE:
        case EventType::REVISIT:
            j["segments"] = labels_json(ev.labels);
            break;
    }
    return j;
}

static Json::Value arm_to_json(const ArmResult &arm, const Strategy &st, const std::string &detail) {
    Json::Value j;
    j["status"] = arm.status == ArmResult::COMPLETE ? "complete"
                : arm.status == ArmResult::TRUNCATED ? "truncated" : "pruned";
    Json::Value fr;
    fr["live_paths"] = uint_json(arm.frontier_live_paths);
    fr["live_labels"] = uint_json(arm.frontier_live_labels);
    // false when a remaining head carried a list cut by labels.max_labels_per_node
    // (annotate mode): live_labels is then a lower bound
    fr["exact"] = arm.frontier_live_labels_exact;
    j["frontier_remaining"] = fr;
    // the completeness guarantee (§6.10): every admissible walk of at most this many
    // bases is present; equals bounds.max_extension_bp exactly when status is complete
    j["complete_to_bp"] = uint_json(arm.complete_to_bp);
    // per-node label lists cut by labels.max_labels_per_node: non-zero means the
    // recorded sets are incomplete and must not be read as "these labels and no other"
    Json::Value lpn;
    lpn["cap"] = uint_json(st.max_labels_per_node);
    lpn["max_seen"] = uint_json(arm.max_labels_at_node);
    lpn["nodes_truncated"] = uint_json(arm.nodes_labels_truncated);
    j["labels_per_node"] = lpn;
    if (arm.cap_trigger) {
        Json::Value c;
        c["reason"] = to_string(arm.cap_trigger->reason);
        c["at_bp"] = uint_json(arm.cap_trigger->at_bp);
        c["segment"] = uint_json(arm.cap_trigger->segment);
        c["live_paths"] = uint_json(arm.cap_trigger->live_paths);
        c["live_labels"] = uint_json(arm.cap_trigger->live_labels);
        c["exact"] = arm.cap_trigger->live_labels_exact;
        j["cap_trigger"] = c;
    }
    if (detail != "summary") {
        Json::Value segs(Json::arrayValue);
        for (const auto &s : arm.segments) {
            Json::Value sj;
            sj["id"] = uint_json(s.id);
            Json::Value parents(Json::arrayValue);
            for (size_t p : s.parents) parents.append(uint_json(p));
            sj["parents"] = parents;
            sj["from_bp"] = uint_json(s.from_bp);
            sj["length_bp"] = uint_json(s.length_bp);
            if (detail == "full" && st.sequences) sj["sequence"] = s.sequence;
            sj["labels"] = labels_json(s.labels_start);
            // the true count at the entry node: for the root (the seed boundary, which
            // is in no run) the only place a cut of its list is reported
            sj["labels_total"] = uint_json(s.labels_start_total);
            sj["labels_truncated"] = s.labels_start_total > s.labels_start.size();
            sj["labels_at_end"] = labels_json(s.labels_end);
            Json::Value evs(Json::arrayValue);
            for (const auto &e : s.events) evs.append(event_to_json(e));
            sj["events"] = evs;
            if (st.label_mode == LabelMode::ANNOTATE) {
                // the labels present along the segment, as runs over its bases
                Json::Value sets(Json::arrayValue);
                for (const auto &r : s.label_sets) {
                    Json::Value rj;
                    rj["from_bp"] = uint_json(r.from_bp);
                    rj["to_bp"] = uint_json(r.to_bp);
                    rj["labels"] = labels_json(r.labels);
                    rj["labels_total"] = uint_json(r.labels_total);
                    rj["truncated"] = r.truncated();
                    sets.append(rj);
                }
                sj["label_sets"] = sets;
            }
            segs.append(sj);
        }
        j["segments"] = segs;
        // the trie view (design note §5.1.1): at_bp is the shared prefix length from the
        // seed boundary; per branch the labels at its first node. A branch's count is
        // not a share of labels_before — one label may follow several branches
        Json::Value splits(Json::arrayValue);
        for (const auto &s : arm.splits) {
            Json::Value sj;
            sj["at_bp"] = uint_json(s.at_bp);
            sj["prefix_bp"] = uint_json(s.at_bp);
            sj["segment"] = uint_json(s.segment);
            Json::Value ch(Json::arrayValue);
            for (size_t c : s.children) ch.append(uint_json(c));
            sj["children"] = ch;
            sj["kind"] = s.ambiguous ? "ambiguous" : "divergence";
            sj["labels_before"] = uint_json(s.labels_before);
            Json::Value branches(Json::arrayValue);
            for (const auto &b : s.branches) {
                Json::Value bj;
                bj["segment"] = uint_json(b.segment);
                bj["char"] = std::string(1, b.ch);
                bj["labels_distinct"] = uint_json(b.labels_distinct);
                bj["labels"] = labels_json(b.labels);
                bj["labels_truncated"] = b.labels_distinct > b.labels.size();
                branches.append(bj);
            }
            sj["branches"] = branches;
            splits.append(sj);
        }
        j["splits"] = splits;
    }
    Json::Value paths(Json::arrayValue);
    for (const auto &p : arm.paths) {
        Json::Value pj;
        pj["id"] = uint_json(p.id);
        if (detail != "summary") {
            Json::Value segs(Json::arrayValue);
            for (size_t s : p.segments) segs.append(uint_json(s));
            pj["segments"] = segs;
        }
        pj["length_bp"] = uint_json(p.length_bp);
        // an object even when empty (annotate mode: no label ends, only path_reason)
        Json::Value reasons(Json::objectValue);
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            if (p.end_reasons[r]) reasons[to_string(static_cast<EndReason>(r))] = p.end_reasons[r];
        }
        pj["end_reasons"] = reasons;
        if (p.path_reason) pj["path_reason"] = to_string(*p.path_reason);
        pj["n_labels"] = uint_json(p.end_labels.size());
        if (detail != "summary") {
            Json::Value els(Json::arrayValue);
            for (const auto &e : p.end_labels) {
                Json::Value ej;
                ej["label"] = e.label; ej["loss"] = e.loss; ej["branches"] = e.branches; ej["run"] = e.run;
                // non-zero: merged in from another route at a reconvergence, so this
                // label does NOT support the spelled bases after route_bp
                ej["route_bp"] = uint_json(e.route_bp);
                els.append(ej);
            }
            pj["end_labels"] = els;
        }
        if (p.continuation) {
            Json::Value c;
            c["sequence"] = p.continuation->sequence;
            c["labels"] = labels_json(p.continuation->labels);
            c["loss_used"] = p.continuation->loss_used;
            c["branches_used"] = p.continuation->branches_used;
            pj["continuation"] = c;
        }
        paths.append(pj);
    }
    j["paths"] = paths;
    if (detail != "summary") {
        Json::Value runs(Json::arrayValue);
        for (const auto &r : arm.runs) {
            Json::Value rj;
            rj["label"] = r.label;
            rj["from_bp"] = uint_json(r.from_bp);
            rj["to_bp"] = uint_json(r.to_bp);
            rj["entered_by"] = r.entered_by_switch ? "switch" : "seed";
            rj["route_bp"] = uint_json(r.route_bp);   // see LabelEnd::route_bp
            if (r.entered_by_switch) { rj["from"] = r.from_label; rj["cost"] = r.switch_cost; }
            if (r.ended) rj["end_reason"] = to_string(r.end_reason);
            runs.append(rj);
        }
        j["runs"] = runs;
    }
    Json::Value growth(Json::arrayValue);
    for (const auto &g : arm.growth) {
        Json::Value gj;
        gj["from_bp"] = uint_json(g.from_bp);
        gj["max_live_paths"] = uint_json(g.max_live_paths);
        gj["distinct_live_labels"] = uint_json(g.max_live_labels);
        gj["live_pairs"] = uint_json(g.max_live_pairs);
        gj["exact"] = g.live_labels_exact;
        gj["steps"] = uint_json(g.steps);
        gj["divergences"] = uint_json(g.divergences);
        gj["ambiguous_branches"] = uint_json(g.ambiguous_branches);
        gj["splits"] = uint_json(g.splits);
        gj["reconvergences"] = uint_json(g.reconvergences);
        gj["bubbles"] = uint_json(g.bubbles);
        gj["tips"] = uint_json(g.tips);
        gj["blocked_repeat"] = uint_json(g.blocked_repeat);
        Json::Value ends;
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            if (g.label_ends[r]) ends[to_string(static_cast<EndReason>(r))] = g.label_ends[r];
        }
        gj["label_ends"] = ends;
        growth.append(gj);
    }
    j["growth"] = growth;
    Json::Value bes(Json::arrayValue);
    for (const auto &b : arm.branch_events) {
        Json::Value bj;
        bj["at_bp"] = uint_json(b.at_bp);
        bj["segment"] = uint_json(b.segment);
        std::string chars(b.chars.begin(), b.chars.end());
        bj["successors"] = chars;
        Json::Value lps(Json::arrayValue);
        for (size_t n : b.labels_per_successor) lps.append(uint_json(n));
        bj["labels_per_successor"] = lps;
        bj["ambiguous"] = labels_json(b.ambiguous);
        bj["dropped"] = labels_json(b.dropped);
        bes.append(bj);
    }
    j["branch_events"] = bes;
    j["branch_events_truncated"] = uint_json(arm.branch_events_total > arm.branch_events.size()
                                             ? arm.branch_events_total - arm.branch_events.size() : 0);
    Json::Value nb(Json::arrayValue);
    for (double x : arm.needed_budgets) nb.append(x);
    j["needed_budgets"] = nb;
    Json::Value counters;
    counters["steps"] = uint_json(arm.steps);
    counters["successor_enumerations"] = uint_json(arm.successor_enumerations);
    counters["output_bp"] = uint_json(arm.output_bp);
    counters["pair_evaluations"] = uint_json(arm.pair_evaluations);
    // the work of the per-path edge-reuse rule and of re-minimising ambiguous nodes
    // (spec §6.4): a pathological locus shows up here, not as silence
    counters["edge_reuse_probes"] = uint_json(arm.edge_reuse_probes);
    counters["reminimisation_rounds"] = uint_json(arm.reminimisation_rounds);
    counters["max_reminimisation_rounds"] = uint_json(arm.max_reminimisation_rounds);
    j["counters"] = counters;
    return j;
}

Json::Value seed_result_to_json(const SeedResult &r, const Strategy &st, const std::string &detail, bool timing) {
    Json::Value j;
    Json::Value seed;
    seed["seed_id"] = r.seed_id;
    seed["validated_seed_id"] = r.validated_seed_id;
    seed["seed_id_mismatch"] = r.seed_id_mismatch;
    seed["length_bp"] = uint_json(r.length_bp);
    seed["num_kmers"] = uint_json(r.num_kmers);
    Json::Value labels(Json::arrayValue);
    for (size_t i = 0; i < r.num_seed_labels; ++i) labels.append(r.label_dict[i].name);
    seed["labels"] = labels;
    // the permitted set was derived from the seed: how many labels support it in full,
    // and what the cap cut (never silently)
    seed["labels_from_seed"] = r.labels_from_seed;
    seed["labels_supporting_total"] = uint_json(r.labels_supporting_total);
    seed["labels_dropped"] = uint_json(r.labels_dropped);
    seed["labels_dropped_digest"] = r.labels_dropped_digest;
    Json::Value dropped(Json::arrayValue);
    for (const auto &d : r.dropped_labels) {
        Json::Value dj;
        dj["label"] = d.name;
        dj["reason"] = d.reason;
        Json::Value runs(Json::arrayValue);
        for (const auto &[a, b] : d.runs) {
            Json::Value iv(Json::arrayValue);
            iv.append(uint_json(a)); iv.append(uint_json(b));
            runs.append(iv);
        }
        dj["runs"] = runs;
        dropped.append(dj);
    }
    seed["dropped_labels"] = dropped;
    j["seed"] = seed;
    j["label_mode"] = to_string(st.label_mode);
    Json::Value dict(Json::arrayValue);
    for (const auto &l : r.label_dict) {
        Json::Value lj;
        lj["name"] = l.name;
        lj["kind"] = to_string(l.kind);
        dict.append(lj);
    }
    j["label_dict"] = dict;
    Json::Value arms;
    for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
        const auto &a = r.arms[static_cast<size_t>(arm)];
        if (!a.requested) continue;
        arms[to_string(arm)] = arm_to_json(a, st, detail);
    }
    j["arms"] = arms;
    Json::Value summary(Json::arrayValue);
    for (size_t l = 0; l < r.label_summary.size(); ++l) {
        Json::Value lj;
        lj["label"] = static_cast<Json::UInt>(l);
        for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
            const auto &s = r.label_summary[l][static_cast<size_t>(arm)];
            if (!r.arms[static_cast<size_t>(arm)].requested) continue;
            Json::Value sj;
            sj["direct_bp"] = uint_json(s.direct_bp);
            sj["reach_bp"] = uint_json(s.reach_bp);
            sj["reentries"] = uint_json(s.reentries);
            Json::Value runs(Json::arrayValue);
            for (uint32_t run : s.runs) runs.append(run);
            sj["runs"] = runs;
            lj[to_string(arm)] = sj;
        }
        summary.append(lj);
    }
    j["label_summary"] = summary;
    Json::Value ann;
    ann["access_path"] = r.access_path;
    ann["keys_mapped"] = uint_json(r.annotation_counters.keys_mapped);
    ann["rows_requested"] = uint_json(r.annotation_counters.rows_requested);
    ann["direct_reads"] = uint_json(r.annotation_counters.direct_reads);
    j["annotation"] = ann;
    if (timing) {
        Json::Value t;
        t["elapsed_ms"] = r.elapsed_seconds * 1000;
        t["cache_hits"] = uint_json(r.annotation_counters.cache_hits);
        // physical fetch work: prefetching along unbranched runs changes these with
        // annotation.batch_kmers while the walk does not (spec §6.8), so they are
        // timing, not part of the invariant result
        t["rows_fetched"] = uint_json(r.annotation_counters.rows_fetched);
        t["tuple_rows_fetched"] = uint_json(r.annotation_counters.tuple_rows_fetched);
        t["coords_mapped"] = uint_json(r.annotation_counters.coords_mapped);
        t["annotation_fetch_ms"] = r.annotation_counters.fetch_seconds * 1000;
        j["timing"] = t;
    }
    return j;
}

static Json::Value interval_json(const KmerInterval &iv) {
    Json::Value j(Json::arrayValue);
    j.append(uint_json(iv.begin));
    j.append(uint_json(iv.end));
    return j;
}

Json::Value profile_to_json(const SupportProfile &p, const SeedSelection *sel, const std::string &run_format) {
    Json::Value j;
    j["k"] = uint_json(p.k);
    j["regime"] = to_string(p.regime);
    j["num_kmers"] = uint_json(p.num_kmers);
    j["support"] = to_string(p.support);
    auto runs_json = [&](const std::vector<KmerInterval> &runs) -> Json::Value {
        if (run_format == "rle")
            return Json::Value(encode_runs(runs, p.num_kmers));
        Json::Value arr(Json::arrayValue);
        for (const auto &r : runs) arr.append(interval_json(r));
        return arr;
    };
    j["graph_runs"] = runs_json(p.graph_runs);
    Json::Value labels(Json::arrayValue);
    for (const auto &l : p.labels) {
        Json::Value lj;
        lj["label"] = l.label.name;
        lj["kind"] = to_string(l.label.kind);
        lj["kmers_supported"] = uint_json(l.kmers_supported);
        lj["runs"] = runs_json(l.runs);
        if (p.support == Support::TRACE) {
            Json::Value tb(Json::arrayValue);
            for (uint64_t x : l.trace_breaks) tb.append(uint_json(x));
            lj["trace_breaks"] = tb;
        }
        labels.append(lj);
    }
    j["labels"] = labels;
    if (p.labels_truncated) {
        Json::Value t;
        t["kept"] = uint_json(p.labels_truncated->kept);
        t["total"] = uint_json(p.labels_truncated->total);
        t["min_kept_kmers"] = uint_json(p.labels_truncated->min_kept_kmers);
        t["max_dropped_kmers"] = uint_json(p.labels_truncated->max_dropped_kmers);
        t["dropped_full_length"] = uint_json(p.labels_truncated->dropped_full_length);
        j["labels_truncated"] = t;
    } else {
        j["labels_truncated"] = Json::Value::null;
    }
    Json::Value cands(Json::arrayValue);
    for (const auto &c : p.candidates) {
        Json::Value cj;
        cj["kmer_interval"] = interval_json(c.kmers);
        auto [a, b] = p.bp_interval(c.kmers);
        cj["bp_interval"] = interval_json({ a, b });
        Json::Value ls(Json::arrayValue);
        for (LabelId l : c.labels) ls.append(p.labels[l].label.name);
        cj["labels"] = ls;
        cands.append(cj);
    }
    j["candidates"] = cands;
    if (sel) {
        Json::Value s;
        Json::Value pol;
        pol["policy"] = sel->policy.policy == SelectionPolicy::LONGEST_FIRST ? "longest_first"
                      : sel->policy.policy == SelectionPolicy::MAX_SUPPORT ? "max_support" : "explicit";
        pol["max_seeds"] = uint_json(sel->policy.max_seeds);
        pol["min_block_bp"] = uint_json(sel->policy.min_block_bp);
        pol["max_labels_per_seed"] = uint_json(sel->policy.max_labels_per_seed);
        pol["label_order"] = sel->policy.label_order == SelectionPolicy::HASH ? "hash"
                           : sel->policy.label_order == SelectionPolicy::COLUMN_ID ? "column_id" : "kmers_supported";
        pol["sample_seed"] = uint_json(sel->policy.sample_seed);
        pol["merge_overlapping"] = sel->policy.merge_overlapping;
        s["policy"] = pol;
        Json::Value counts;
        counts["candidates"] = uint_json(sel->num_candidates);
        counts["eligible"] = uint_json(sel->num_eligible);
        counts["selected"] = uint_json(sel->seeds.size());
        s["counts"] = counts;
        Json::Value seeds(Json::arrayValue);
        for (const auto &seed : sel->seeds) {
            Json::Value sj;
            sj["seed_id"] = seed.seed_id;
            sj["sequence"] = seed.sequence;
            sj["kmer_interval"] = interval_json(seed.kmers);
            Json::Value ls(Json::arrayValue);
            for (const auto &l : seed.labels) ls.append(l);
            sj["labels"] = ls;
            Json::Value pop;
            pop["supporting_total"] = uint_json(seed.population.supporting_total);
            pop["included"] = uint_json(seed.population.included);
            pop["dropped_count"] = uint_json(seed.population.dropped_count);
            if (seed.population.dropped_count) {
                pop["dropped_digest"] = seed.population.dropped_digest;
                Json::Value dl(Json::arrayValue);
                for (const auto &d : seed.population.dropped) dl.append(d);
                pop["dropped"] = dl;
            }
            sj["label_population"] = pop;
            Json::Value ov(Json::arrayValue);
            for (size_t o : seed.overlaps_with) ov.append(uint_json(o));
            sj["overlaps_with"] = ov;
            if (!seed.labels_not_covering.empty()) {
                Json::Value nc(Json::arrayValue);
                for (const auto &l : seed.labels_not_covering) nc.append(l);
                sj["labels_not_covering"] = nc;
            }
            seeds.append(sj);
        }
        s["seeds"] = seeds;
        j["selection"] = s;
    }
    return j;
}

Json::Value capabilities_to_json(const LabelOracle &oracle, const std::string &release) {
    Json::Value c;
    c["schema_version"] = 1;
    c["k"] = uint_json(oracle.get_k());
    c["regime"] = to_string(oracle.regime());
    c["alphabet"] = oracle.graph().alphabet();
    c["num_labels"] = uint_json(oracle.num_columns());
    c["has_coordinates"] = oracle.has_coordinates();
    c["has_coord_to_header"] = oracle.coord_to_header() != nullptr;
    c["supports_trace"] = oracle.has_coordinates() && oracle.regime() == Regime::BASIC;
    Json::Value models(Json::arrayValue);
    models.append("forbid"); models.append("constant"); models.append("table");
    c["cost_models_available"] = models;
    Json::Value modes(Json::arrayValue);
    modes.append("constrain"); modes.append("annotate");
    c["label_modes"] = modes;
    c["direct_access"] = oracle.supports_direct();
    c["release"] = release;
    return c;
}


// ---------------------------------------------------------------- processing

static LabelChangeCost make_cost(const CostSpec &spec, const std::vector<std::string> &dict) {
    switch (spec.model) {
        case CostSpec::FORBID:
            return LabelChangeCost::forbid();
        case CostSpec::CONSTANT:
            return LabelChangeCost::constant(spec.value);
        case CostSpec::TABLE: {
            std::map<std::pair<LabelId, LabelId>, double> entries;
            auto id_of = [&](const std::string &name) -> std::optional<LabelId> {
                auto it = std::find(dict.begin(), dict.end(), name);
                if (it == dict.end()) return std::nullopt;
                return static_cast<LabelId>(it - dict.begin());
            };
            for (const auto &[from, to, c] : spec.entries) {
                auto a = id_of(from), b = id_of(to);
                if (a && b) entries[{ *a, *b }] = c;   // entries for labels outside P are irrelevant
            }
            return LabelChangeCost::table(std::move(entries), spec.default_cost);
        }
    }
    return LabelChangeCost::forbid();
}

Json::Value process_traverse_request(const Json::Value &json,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const std::string &release,
                                     const TraverseLimits &limits) {
    TraverseRequest req = parse_traverse_request(json);
    if (!req.release.empty() && !release.empty() && req.release != release)
        throw InvalidRequest("request.release '" + req.release + "' does not match the loaded index release '"
                             + release + "'");
    if (limits.max_seeds && req.seeds.size() > limits.max_seeds) {
        throw InvalidRequest("request.seeds: at most " + std::to_string(limits.max_seeds)
                             + " seeds per request");
    }
    if (limits.max_seed_bp) {
        // The derivation reads one full annotation row per seed k-mer, so the seed length
        // is itself a cost the caller declares; /resolve bounds its query the same way.
        for (const auto &seed : req.seeds) {
            if (seed.sequence.size() > limits.max_seed_bp) {
                throw InvalidRequest("request.seeds: a seed is longer than the server limit of "
                                     + std::to_string(limits.max_seed_bp) + " bp");
            }
        }
    }
    Json::Value clamped(Json::arrayValue);
    auto clamp = [&](const char *field, double requested, double effective) {
        Json::Value c;
        c["field"] = field;
        c["requested"] = requested;
        c["effective"] = effective;
        clamped.append(c);
    };
    const bool derives = std::any_of(req.seeds.begin(), req.seeds.end(),
                                     [](const Seed &s) { return s.labels.empty(); });
    if (limits.max_time_ms > 0 && req.strategy.time_budget_ms > limits.max_time_ms) {
        clamp("bounds.time_budget_ms", req.strategy.time_budget_ms, limits.max_time_ms);
        req.strategy.time_budget_ms = limits.max_time_ms;
    } else if (limits.max_time_ms > 0 && derives && req.strategy.time_budget_ms <= 0) {
        // A zero budget means "do not extend", and the walk reads it that way. But this
        // budget is also what bounds the DERIVATION of a permitted set the caller did not
        // name, which runs before any extension: zero would leave that phase unbounded.
        clamp("bounds.time_budget_ms", req.strategy.time_budget_ms, limits.max_time_ms);
        req.strategy.time_budget_ms = limits.max_time_ms;
    }
    if (limits.max_seed_labels && req.strategy.max_seed_labels > limits.max_seed_labels) {
        clamp("labels.max_seed_labels", static_cast<double>(req.strategy.max_seed_labels),
              static_cast<double>(limits.max_seed_labels));
        req.strategy.max_seed_labels = limits.max_seed_labels;
    }

    LabelOracle oracle(anno_graph);
    Timer timer;
    Json::Value out;
    out["release"] = release;
    out["capabilities"] = capabilities_to_json(oracle, release);
    Json::Value strategy = strategy_to_json(req.strategy, req.cost);
    strategy["clamped"] = clamped;
    out["strategy"] = strategy;
    // what every arm's complete_to_bp quantifies over: the per-path edge-reuse rule is
    // what makes "all walks" a finite set in a graph with cycles
    out["walk_rule"] = walk_rule_statement(req.strategy, oracle);
    out["algorithm_version"] = "traverse-0.2";

    Json::Value results(Json::arrayValue);
    std::set<std::string> seen_ids;
    // The oracle's counters accumulate over the whole request, but each seed must report
    // its own work, or the same seed would read differently depending on its position
    // (and the per-seed result would not be invariant under reordering).
    LabelOracle::Counters before = oracle.counters();
    auto per_seed = [&](const LabelOracle::Counters &now) {
        LabelOracle::Counters d;
        d.keys_mapped = now.keys_mapped - before.keys_mapped;
        d.rows_requested = now.rows_requested - before.rows_requested;
        d.cache_hits = now.cache_hits - before.cache_hits;
        d.rows_fetched = now.rows_fetched - before.rows_fetched;
        d.tuple_rows_fetched = now.tuple_rows_fetched - before.tuple_rows_fetched;
        d.direct_reads = now.direct_reads - before.direct_reads;
        d.coords_mapped = now.coords_mapped - before.coords_mapped;
        d.fetch_seconds = now.fetch_seconds - before.fetch_seconds;
        before = now;
        return d;
    };
    for (const auto &seed : req.seeds) {
        std::vector<std::string> dict = seed.labels;
        dict.insert(dict.end(), req.strategy.extra.begin(), req.strategy.extra.end());
        LabelChangeCost cost = make_cost(req.cost, dict);
        try {
            SeedResult r = traverse_seed(oracle, seed, req.strategy, cost, release);
            r.annotation_counters = per_seed(r.annotation_counters);
            Json::Value rj = seed_result_to_json(r, req.strategy, req.detail, req.timing);
            if (!seen_ids.insert(r.validated_seed_id).second)
                rj["duplicate"] = true;
            results.append(rj);
        } catch (const SeedDerivationError &e) {
            // The permitted set could not be derived from THIS seed (nothing carries it
            // in full, the budget ran out, the names are ambiguous). The caller named no
            // labels, so it had no lever on that and nothing to fix in the request:
            // report it against the seed and keep the other seeds' traversals, instead of
            // discarding a 100-seed batch because seed 57 spans a recombination point.
            Json::Value rj;
            Json::Value sj;
            sj["seed_id"] = seed.seed_id;
            sj["length_bp"] = uint_json(seed.sequence.size());
            sj["labels_from_seed"] = true;
            rj["seed"] = sj;
            rj["error"] = e.what();
            results.append(rj);
            per_seed(oracle.counters());   // do not bill this seed's reads to the next
        } catch (const std::invalid_argument &e) {
            throw InvalidRequest(std::string("seed '") + (seed.seed_id.empty() ? seed.sequence.substr(0, 32) : seed.seed_id)
                                 + "': " + e.what());
        }
    }
    out["results"] = results;
    if (req.timing) {
        Json::Value t;
        t["elapsed_ms"] = timer.elapsed() * 1000;
        out["timing"] = t;
    }
    return out;
}

Json::Value process_resolve_request(const Json::Value &json,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const std::string &release,
                                    uint64_t max_query_bp) {
    ResolveRequest req = parse_resolve_request(json);
    if (max_query_bp && req.sequence.size() > max_query_bp)
        throw InvalidRequest("request.sequence: longer than the server limit of "
                             + std::to_string(max_query_bp) + " bp");
    if (req.sequence.size() > req.max_query_bp)
        throw InvalidRequest("request.sequence: longer than bounds.max_query_bp");

    LabelOracle oracle(anno_graph);
    Timer timer;
    SupportProfile profile;
    std::optional<SeedSelection> selection;
    try {
        profile = resolve_support(oracle, req.sequence, req.options);
        if (req.select) {
            req.policy.release_id = release;
            selection = select_seeds(profile, req.sequence, req.policy, oracle.regime() != Regime::BASIC);
        }
    } catch (const std::invalid_argument &e) {
        throw InvalidRequest(e.what());
    }
    Json::Value out = profile_to_json(profile, selection ? &*selection : nullptr, req.run_format);
    out["release"] = release;
    out["capabilities"] = capabilities_to_json(oracle, release);
    Json::Value t;
    t["elapsed_ms"] = timer.elapsed() * 1000;
    out["timing"] = t;
    return out;
}


// ---------------------------------------------------------------- CLI

int traverse_graph(Config *config) {
    assert(config);
    assert(config->infbase_annotators.size() == 1);

    auto loaded = load_graph_with_async_annotation(*config);
    auto graph = loaded.first.get();
    auto anno_graph = loaded.second.get();
    assert(anno_graph);
    graph.reset();

    Json::StreamWriterBuilder builder;
    builder["indentation"] = config->output_json ? "" : "  ";
    int status = 0;
    for (const auto &file : config->fnames) {
        std::ifstream in(file);
        if (!in.good()) {
            logger->error("Cannot open request file {}", file);
            return 1;
        }
        Json::Value json;
        Json::CharReaderBuilder reader;
        std::string errs;
        if (!Json::parseFromStream(reader, in, &json, &errs)) {
            logger->error("Invalid JSON in {}: {}", file, errs);
            return 1;
        }
        try {
            Json::Value out = config->traverse_resolve
                ? process_resolve_request(json, *anno_graph, config->index_release)
                : process_traverse_request(json, *anno_graph, config->index_release);
            std::cout << Json::writeString(builder, out) << std::endl;
        } catch (const InvalidRequest &e) {
            logger->error("Invalid request in {}: {}", file, e.what());
            Json::Value err;
            err["error"] = e.what();
            std::cout << Json::writeString(builder, err) << std::endl;
            status = 1;
        }
    }
    return status;
}

} // namespace cli
} // namespace mtg
