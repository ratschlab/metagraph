#include "traverse.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <tuple>

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

// The arm counters by name, in the walker's order: one list for the JSON `counters` and
// the graphlet's A record, so that a counter added here reaches both (an MGT reader
// keeps a name it does not know).
std::vector<std::pair<const char*, uint64_t>> arm_counters(const ArmResult &arm) {
    return {
        { "steps", arm.steps },
        { "successor_enumerations", arm.successor_enumerations },
        { "output_bp", arm.output_bp },
        { "pair_evaluations", arm.pair_evaluations },
        // the work of the per-path edge-reuse rule, of re-minimising ambiguous nodes
        // (spec §6.4) and of recording the refusals (§7.2): a pathological locus shows
        // up here, not as silence
        { "edge_reuse_probes", arm.edge_reuse_probes },
        { "reminimisation_rounds", arm.reminimisation_rounds },
        { "max_reminimisation_rounds", arm.max_reminimisation_rounds },
        { "refusal_scans", arm.refusal_scans },
        // derivations a max_switch_sources cut may have changed (the switch_sources
        // limitation's observed)
        { "switch_sources_cut", arm.switch_sources_cut },
    };
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
            // Not implemented (spec §6.5): a window would be accepted and then walked as
            // 0, a different question than the one asked, so only 0 (the echo's value,
            // which must stay resubmittable) is accepted.
            for (const char *window : { "tip_window_bp", "bubble_window_bp" }) {
                if (b.uint(window, 0)) {
                    throw InvalidRequest(b.child_path(window) + ": not implemented (tip and bubble "
                                         "windows are not supported yet); set it to 0 or omit it");
                }
            }
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
                { { "summary", "summary" }, { "tree", "tree" }, { "full", "full" },
                  { "graphlet", "graphlet" } });
            st.sequences = o.boolean("sequences", true);
            st.profile_bin_bp = o.uint("profile_bin_bp", 100, 1);
            // refusals are the evidence for a successor not taken, so the cap on them
            // must be removable: "unlimited" keeps every event (§7.2)
            st.max_branch_events = limit_or_unlimited(o, "max_branch_events", 100,
                                                      std::numeric_limits<size_t>::max() - 1,
                                                      nullptr);
            // 0: no continuation sequence (labels and loss still reported); 1 .. k - 1 is
            // refused in process_traverse_request, where k is known: a continuation
            // shorter than k could not be resubmitted as a seed
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
                entries.append(std::move(e));
            }
            j["entries"] = std::move(entries);
        }
    }
    return j;
}

Json::Value strategy_to_json(const Strategy &st, const CostSpec &cost, const std::string &detail,
                             bool timing) {
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
    labels["extra"] = std::move(extra);
    labels["change_cost"] = cost_to_json(cost);
    labels["loss_budget"] = st.loss_budget;
    labels["switch_on"] = st.switch_on_loss_only ? "loss" : "any";
    labels["max_switch_sources"] = limit_json(st.max_switch_sources);
    labels["max_seed_labels"] = uint_json(st.max_seed_labels);
    // "auto": header when the index has a CoordToHeader, column otherwise
    labels["seed_label_kind"] = st.seed_label_kind ? to_string(*st.seed_label_kind) : "auto";
    j["labels"] = std::move(labels);
    Json::Value b;
    b["max_label_branches"] = limit_json(st.max_label_branches);
    b["min_successor_labels"] = uint_json(st.min_successor_labels);
    b["min_successor_fraction"] = st.min_successor_fraction;
    b["tip_window_bp"] = uint_json(st.tip_window_bp);
    b["bubble_window_bp"] = uint_json(st.bubble_window_bp);
    b["on_reconverge"] = st.merge_reconverge ? "merge" : "keep";
    b["hairpins"] = st.skip_hairpins ? "skip" : "follow";
    b["max_splits_per_path"] = limit_json(st.max_splits_per_path);
    j["branching"] = std::move(b);
    Json::Value bo;
    bo["max_extension_bp"] = uint_json(st.max_extension_bp);
    bo["min_live_labels"] = uint_json(st.min_live_labels);
    bo["max_steps"] = uint_json(st.max_steps);
    bo["max_live_paths"] = uint_json(st.max_live_paths);
    bo["max_paths"] = uint_json(st.max_paths);
    bo["max_output_bp"] = uint_json(st.max_output_bp);
    bo["time_budget_ms"] = st.time_budget_ms;
    j["bounds"] = std::move(bo);
    Json::Value f;
    f["order"] = st.order == Strategy::BREADTH_FIRST ? "breadth_first"
               : st.order == Strategy::LOWEST_LOSS_FIRST ? "lowest_loss_first" : "most_supported_first";
    f["on_overflow"] = st.on_overflow == Strategy::STOP ? "stop" : "beam";
    j["frontier"] = std::move(f);
    Json::Value o;
    // request fields, not walker knobs; echoed so that the echo is the whole request
    o["detail"] = detail;
    o["timing"] = timing;
    o["sequences"] = st.sequences;
    o["profile_bin_bp"] = uint_json(st.profile_bin_bp);
    o["max_branch_events"] = limit_json(st.max_branch_events);
    o["continuation_bp"] = uint_json(st.continuation_bp);
    j["output"] = std::move(o);
    Json::Value a;
    a["batch_kmers"] = uint_json(st.batch_kmers);
    j["annotation"] = std::move(a);
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
            // followed (hairpins: follow, the child is present) or skipped (the only
            // record of a successor not taken); a checker needs to know which
            j["followed"] = ev.text == "followed";
            break;
        case EventType::TIP:
            j["char"] = std::string(1, ev.ch); j["length_bp"] = uint_json(ev.length_bp);
            break;
        case EventType::BUBBLE:
            j["length_bp"] = uint_json(ev.length_bp); j["alleles"] = ev.text;
            break;
        case EventType::RECONVERGE:
            j["segments"] = labels_json(ev.labels);
            break;
        case EventType::REVISIT:
            j["segments"] = labels_json(ev.labels);
            // the distance difference to the first arrival; same_distance: a keep-mode
            // join of two walks at one node and depth (no difference)
            j["length_bp"] = uint_json(ev.length_bp);
            j["same_distance"] = ev.text == "same_distance";
            break;
    }
    return j;
}

// ---------------------------------------------------------------- stated limitations

// One entry of `limitations` (spec §7.0): what limited the result, the request field
// that controls it (relative to `strategy`, as in strategy.clamped), its value, what the
// run met against it, and what is missing because of it — so that an agent turns the
// knob instead of parsing prose, and nothing is cut without saying so.
static Json::Value limitation(const char *kind, const std::string &knob, Json::Value limit,
                              Json::Value observed, const std::string &effect) {
    Json::Value j;
    j["kind"] = kind;
    j["knob"] = knob;
    j["limit"] = std::move(limit);
    j["observed"] = std::move(observed);
    j["effect"] = effect;
    return j;
}

// The request field a resource stop answers to (§6.7) and its value; a beam's width is
// bounds.max_live_paths.
static std::pair<std::string, Json::Value> cap_knob(EndReason reason, const Strategy &st) {
    switch (reason) {
        case EndReason::MAX_STEPS: return { "bounds.max_steps", uint_json(st.max_steps) };
        case EndReason::MAX_LIVE_PATHS:
        case EndReason::BEAM_PRUNED: return { "bounds.max_live_paths", uint_json(st.max_live_paths) };
        case EndReason::MAX_PATHS: return { "bounds.max_paths", uint_json(st.max_paths) };
        case EndReason::MAX_OUTPUT: return { "bounds.max_output_bp", uint_json(st.max_output_bp) };
        case EndReason::TIME_BUDGET: return { "bounds.time_budget_ms", Json::Value(st.time_budget_ms) };
        default: return { "", Json::Value() };
    }
}

static bool ended_by(const ArmResult &arm, EndReason reason) {
    if (arm.cap_trigger && arm.cap_trigger->reason == reason)
        return true;
    return std::any_of(arm.paths.begin(), arm.paths.end(),
                       [&](const PathResult &p) { return p.path_reason == reason; });
}

// Every cap that limited this arm, each emitted only when it did (spec §7.0).
// With a finite change cost and a finite branch limit, the sources left after a branch-limit
// exclusion are re-minimised greedily (spec §6.4), so a reported loss can exceed the optimum.
// pair_evaluations > 0 means switch costs were actually priced (forbid prices none, and its
// losses are all zero); max_reminimisation_rounds > 0 means an exclusion re-ran derive().
static bool greedy_losses(const ArmResult &arm, const Strategy &st) {
    return arm.max_reminimisation_rounds > 0 && arm.pair_evaluations > 0
        && st.max_label_branches != Strategy::kUnlimited;
}

static Json::Value arm_limitations(const ArmResult &arm, const Strategy &st) {
    Json::Value out(Json::arrayValue);
    // The walk domain: the cap that set complete_to_bp, then any other cap that ended
    // walks after it — a beam prunes a level and a step cap trips later, and raising
    // only the first knob would run into the second.
    if (arm.status != ArmResult::COMPLETE && arm.cap_trigger) {
        std::vector<EndReason> caps { arm.cap_trigger->reason };
        std::map<EndReason, size_t> walks_cut;
        for (const PathResult &p : arm.paths) {
            if (!p.path_reason || !is_resource_stop(*p.path_reason))
                continue;
            walks_cut[*p.path_reason]++;
            if (std::find(caps.begin(), caps.end(), *p.path_reason) == caps.end())
                caps.push_back(*p.path_reason);
        }
        for (EndReason r : caps) {
            const bool trigger = r == arm.cap_trigger->reason;
            std::string effect = "every walk of at most complete_to_bp bp is present, longer walks may be "
                                 "missing: ";
            effect += r == EndReason::BEAM_PRUNED
                ? "the beam kept the best-supported max_live_paths heads of a level and pruned the rest"
                : std::string("the exploration stopped at this cap (") + to_string(r) + ")";
            if (!trigger)
                effect += " after the cap that set complete_to_bp, so raising only that knob stops here";
            effect += "; raise the knob";
            auto [knob, limit] = cap_knob(r, st);
            // at the trigger what the cap compared, which exceeded the limit; for a later
            // cap (no trigger of its own) the walks it ended
            const double demand = arm.cap_trigger->demand;
            Json::Value observed = !trigger ? uint_json(walks_cut[r])
                                 : r == EndReason::TIME_BUDGET ? Json::Value(demand)
                                 : uint_json(static_cast<uint64_t>(demand));
            Json::Value j = limitation("walk_domain", knob, limit, observed, effect);
            j["complete_to_bp"] = uint_json(arm.complete_to_bp);
            out.append(std::move(j));
        }
    }
    if (arm.branch_events_complete_to_bp != std::numeric_limits<uint64_t>::max()) {
        Json::Value j = limitation("branch_events", "output.max_branch_events",
                                   limit_json(st.max_branch_events), uint_json(arm.branch_events_total),
                                   "branch decisions and refusals at or beyond complete_to_bp are not "
                                   "reported; more were produced: raise the knob or set it to "
                                   "\"unlimited\"");
        j["complete_to_bp"] = uint_json(arm.branch_events_complete_to_bp);
        out.append(std::move(j));
    }
    if (arm.nodes_labels_truncated) {
        out.append(limitation("label_lists", "labels.max_labels_per_node",
                              uint_json(st.max_labels_per_node), uint_json(arm.max_labels_at_node),
                              std::to_string(arm.nodes_labels_truncated) + " recorded label list(s) are "
                              "cut at the cap; every list carries its true count, and what is derived "
                              "from the lists (annotate mode's label_summary, an oracle filtered from "
                              "the recorded sets) is a lower bound; raise the knob"));
    }
    // A cut source list (§6.3) is stated by the derivations it may have changed, not only
    // by the label ends it caused: a cut source that goes on along another successor
    // does not end, yet a target it was the cheapest way into was entered at a higher
    // loss or not at all.
    size_t cut_ends = 0;
    for (const Segment &s : arm.segments) {
        for (const Event &ev : s.events) {
            cut_ends += ev.type == EventType::LABEL_END && ev.text == "switch_sources";
        }
    }
    if (arm.switch_sources_cut || cut_ends) {
        std::string effect = "the switch sources of " + std::to_string(arm.switch_sources_cut)
                           + " successor derivation(s) were cut to the cheapest by (loss, branches, "
                             "label) while a cut source could switch into a target there: losses "
                             "may be overestimated and switch entries missed because cheaper "
                             "sources were cut";
        if (cut_ends) {
            effect += "; " + std::to_string(cut_ends) + " label(s) ended label_lost "
                      "(switch_sources) because they were not priced";
        }
        effect += "; raise the knob or set it to \"unlimited\"";
        Json::Value j = limitation("switch_sources", "labels.max_switch_sources",
                                   limit_json(st.max_switch_sources),
                                   uint_json(arm.switch_sources_cut), effect);
        // the label ends the cut caused (label_lost, text switch_sources); a subset of
        // the derivations' effects, since a cut source that goes on does not end
        j["label_ends"] = uint_json(cut_ends);
        out.append(std::move(j));
    }
    size_t inexact = !arm.frontier_live_labels_exact
                   + (arm.cap_trigger && !arm.cap_trigger->live_labels_exact);
    for (const GrowthBin &g : arm.growth) {
        inexact += !g.live_labels_exact;
    }
    if (inexact) {
        out.append(limitation("inexact_counts", "labels.max_labels_per_node",
                              uint_json(st.max_labels_per_node), uint_json(inexact),
                              "the live-label counts flagged \"exact\": false (frontier_remaining, "
                              "cap_trigger, growth bins) counted heads carrying a cut list and are lower "
                              "bounds; raise the knob"));
    }
    if (std::string(arm.completeness_scope) == "united_history") {
        // observed: the merges that actually united histories (0: none was reduced)
        size_t merges = 0;
        for (const GrowthBin &g : arm.growth) {
            merges += g.reconvergences;
        }
        out.append(limitation("scope", "branching.on_reconverge", "merge", uint_json(merges),
                              "completeness holds for the united-history rule, not per path; use "
                              "\"keep\" for the per-path guarantee"));
    }
    if (greedy_losses(arm, st)) {
        out.append(limitation("greedy_losses", "branching.max_label_branches",
                              uint_json(st.max_label_branches), uint_json(arm.reminimisation_rounds),
                              "after a branch-limit exclusion the switch losses were re-minimised "
                              "greedily over the remaining sources (observed: re-minimisation rounds): a "
                              "reported loss, and an entry decided by the loss budget, may differ from the "
                              "optimum; set the knob to \"unlimited\" (or use the forbid cost) for exact "
                              "losses"));
    }
    return out;
}

// The per-seed outcome (spec §7.0), one axis per guarantee, read off the stated limitations
// of the result itself (seed level and every arm): conservative by construction
// (DESIGN-traverse-graphlet.md §14, v5.2) — an axis is complete only when no limitation of
// its class applies, so a reader never finds a limitation whose axis still reads complete.
//   walks               partial: walk_domain (an arm stopped at a cap), seed_labels (a carrier
//                       of the seed was not taken, so the walks only it carries are missing),
//                       or scope with at least one merge (a history was united, so the walks
//                       are complete for the united-history rule, not per path); failed: no
//                       traversal (the derivation failed)
//   branch_diagnostics  cut: branch_events (decisions at and beyond evidence.complete_to_bp
//                       are not reported)
//   label_evidence      lower_bound: evidence may be missing or understated — label_lists,
//                       inexact_counts, seed_labels, switch_sources, greedy_losses;
//                       qualified: something reported may be overstated —
//                       trace_record_boundaries; qualified wins when both apply
//   delivery            inline (the whole result is in this response); spooled / paged
//                       are reserved for the graphlet delivery path
static Json::Value outcome_of(const Json::Value &result, bool failed) {
    bool partial = false, cut = false, lower = false, qualified = false;
    auto classify = [&](const Json::Value &lims) {
        for (const Json::Value &l : lims) {
            const std::string kind = l["kind"].asString();
            partial |= kind == "walk_domain" || kind == "seed_labels"
                    || (kind == "scope" && l["observed"].asUInt64() > 0);
            cut |= kind == "branch_events";
            lower |= kind == "label_lists" || kind == "inexact_counts" || kind == "seed_labels"
                  || kind == "switch_sources" || kind == "greedy_losses";
            qualified |= kind == "trace_record_boundaries";
        }
    };
    classify(result["limitations"]);
    if (result.isMember("arms")) {
        for (const std::string &side : result["arms"].getMemberNames()) {
            const Json::Value &arm = result["arms"][side];
            classify(arm["limitations"]);
            // every non-complete arm states a walk_domain limitation; the status is read
            // too so that an unstated stop could never read as complete
            partial |= arm["status"].asString() != "complete";
        }
    }
    Json::Value o;
    o["walks"] = failed ? "failed" : partial ? "partial" : "complete";
    o["branch_diagnostics"] = cut ? "cut" : "complete";
    o["label_evidence"] = qualified ? "qualified" : lower ? "lower_bound" : "complete";
    o["delivery"] = "inline";
    return o;
}

// The per-arm fields every detail level carries: the certificate (status, complete_to_bp
// and what it quantifies over), the frontier, the label-list cuts, the evidence boundary,
// the stated limitations and the cap trigger.
static void arm_certificate_json(Json::Value *j, const ArmResult &arm, const Strategy &st) {
    (*j)["status"] = arm.status == ArmResult::COMPLETE ? "complete"
                   : arm.status == ArmResult::TRUNCATED ? "truncated" : "pruned";
    Json::Value fr;
    fr["live_paths"] = uint_json(arm.frontier_live_paths);
    fr["live_labels"] = uint_json(arm.frontier_live_labels);
    // false when a remaining head carried a list cut by labels.max_labels_per_node
    // (annotate mode): live_labels is then a lower bound
    fr["exact"] = arm.frontier_live_labels_exact;
    (*j)["frontier_remaining"] = std::move(fr);
    // the completeness guarantee (§6.10): every admissible walk of at most this many
    // bases is present; equals bounds.max_extension_bp exactly when status is complete
    (*j)["complete_to_bp"] = uint_json(arm.complete_to_bp);
    // what it quantifies over: the per-path edge history (on_reconverge keep) or the
    // united history of merged routes (merge), see §6.10
    (*j)["completeness_scope"] = arm.completeness_scope;
    // per-node label lists cut by labels.max_labels_per_node: non-zero means the
    // recorded sets are incomplete and must not be read as "these labels and no other"
    Json::Value lpn;
    lpn["cap"] = uint_json(st.max_labels_per_node);
    lpn["max_seen"] = uint_json(arm.max_labels_at_node);
    lpn["nodes_truncated"] = uint_json(arm.nodes_labels_truncated);
    (*j)["labels_per_node"] = std::move(lpn);
    // the evidence boundary (§7.0): below complete_to_bp every branch decision and
    // refusal is reported; null when nothing was cut
    const bool evidence_complete
        = arm.branch_events_complete_to_bp == std::numeric_limits<uint64_t>::max();
    Json::Value evidence;
    evidence["complete"] = evidence_complete;
    evidence["complete_to_bp"] = evidence_complete ? Json::Value()
                                                   : uint_json(arm.branch_events_complete_to_bp);
    (*j)["evidence"] = std::move(evidence);
    // every cap that limited this arm, with the knob to turn; empty when none did
    (*j)["limitations"] = arm_limitations(arm, st);
    if (arm.cap_trigger) {
        Json::Value c;
        c["reason"] = to_string(arm.cap_trigger->reason);
        c["at_bp"] = uint_json(arm.cap_trigger->at_bp);
        c["segment"] = uint_json(arm.cap_trigger->segment);
        c["live_paths"] = uint_json(arm.cap_trigger->live_paths);
        c["live_labels"] = uint_json(arm.cap_trigger->live_labels);
        c["exact"] = arm.cap_trigger->live_labels_exact;
        (*j)["cap_trigger"] = std::move(c);
    }
}

static Json::Value counters_json(const ArmResult &arm) {
    Json::Value counters(Json::objectValue);
    for (const auto &[name, value] : arm_counters(arm)) {
        counters[name] = uint_json(value);
    }
    return counters;
}

static Json::Value arm_to_json(const ArmResult &arm, const Strategy &st, const std::string &detail) {
    Json::Value j;
    arm_certificate_json(&j, arm, st);
    if (detail != "summary") {
        Json::Value segs(Json::arrayValue);
        for (const auto &s : arm.segments) {
            Json::Value sj;
            sj["id"] = uint_json(s.id);
            Json::Value parents(Json::arrayValue);
            for (size_t p : s.parents) parents.append(uint_json(p));
            sj["parents"] = std::move(parents);
            if (s.parents.size() > 1) {
                // per parent, the labels whose kept entry came through it: the per-label
                // routes through the merge (empty lists in annotate mode, no lineages)
                Json::Value via(Json::arrayValue);
                for (const auto &l : s.labels_via_parent) via.append(labels_json(l));
                sj["labels_via_parent"] = std::move(via);
            }
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
            sj["events"] = std::move(evs);
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
                    sets.append(std::move(rj));
                }
                sj["label_sets"] = std::move(sets);
            }
            segs.append(std::move(sj));
        }
        j["segments"] = std::move(segs);
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
            sj["children"] = std::move(ch);
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
                branches.append(std::move(bj));
            }
            sj["branches"] = std::move(branches);
            splits.append(std::move(sj));
        }
        j["splits"] = std::move(splits);
    }
    Json::Value paths(Json::arrayValue);
    for (const auto &p : arm.paths) {
        Json::Value pj;
        pj["id"] = uint_json(p.id);
        if (detail != "summary") {
            Json::Value segs(Json::arrayValue);
            for (size_t s : p.segments) segs.append(uint_json(s));
            pj["segments"] = std::move(segs);
        }
        pj["length_bp"] = uint_json(p.length_bp);
        // an object even when empty (annotate mode: no label ends, only path_reason)
        Json::Value reasons(Json::objectValue);
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            if (p.end_reasons[r]) reasons[to_string(static_cast<EndReason>(r))] = p.end_reasons[r];
        }
        pj["end_reasons"] = std::move(reasons);
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
                els.append(std::move(ej));
            }
            pj["end_labels"] = std::move(els);
        }
        if (p.continuation) {
            Json::Value c;
            c["sequence"] = p.continuation->sequence;
            c["labels"] = labels_json(p.continuation->labels);
            c["loss_used"] = p.continuation->loss_used;
            c["branches_used"] = p.continuation->branches_used;
            pj["continuation"] = std::move(c);
        }
        paths.append(std::move(pj));
    }
    j["paths"] = std::move(paths);
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
            // where the run ended or was closed (a label end is there at to_bp; a
            // closed run's segment is a parent of the merge), and the lineage's branch
            // count and loss at that point — what a claim inside a walk reports
            rj["segment"] = r.segment == UINT32_MAX ? Json::Value() : Json::Value(r.segment);
            rj["branches"] = r.branches;
            rj["loss"] = r.loss;
            runs.append(std::move(rj));
        }
        j["runs"] = std::move(runs);
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
        // an object even when no label ended in the bin (was null)
        Json::Value ends(Json::objectValue);
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            if (g.label_ends[r]) ends[to_string(static_cast<EndReason>(r))] = g.label_ends[r];
        }
        gj["label_ends"] = std::move(ends);
        growth.append(std::move(gj));
    }
    j["growth"] = std::move(growth);
    Json::Value bes(Json::arrayValue);
    for (const auto &b : arm.branch_events) {
        Json::Value bj;
        bj["at_bp"] = uint_json(b.at_bp);
        bj["segment"] = uint_json(b.segment);
        std::string chars(b.chars.begin(), b.chars.end());
        bj["successors"] = chars;
        Json::Value lps(Json::arrayValue);
        for (size_t n : b.labels_per_successor) lps.append(uint_json(n));
        bj["labels_per_successor"] = std::move(lps);
        bj["ambiguous"] = labels_json(b.ambiguous);
        bj["dropped"] = labels_json(b.dropped);
        // per successor not followed for labels present on it: why, and for which
        Json::Value refused(Json::arrayValue);
        for (const auto &r : b.refused) {
            Json::Value rj;
            rj["char"] = std::string(1, r.ch);
            rj["cause"] = r.cause;
            rj["labels"] = labels_json(r.labels);
            refused.append(std::move(rj));
        }
        bj["refused"] = std::move(refused);
        bes.append(std::move(bj));
    }
    j["branch_events"] = std::move(bes);
    j["branch_events_truncated"] = uint_json(arm.branch_events_total > arm.branch_events.size()
                                             ? arm.branch_events_total - arm.branch_events.size() : 0);
    Json::Value nb(Json::arrayValue);
    for (double x : arm.needed_budgets) nb.append(x);
    j["needed_budgets"] = std::move(nb);
    j["counters"] = counters_json(arm);
    return j;
}

// The arm of the `detail: graphlet` summary (DESIGN-traverse-graphlet.md §3): the
// certificate, the counters and the counts that validate the body, nothing per segment,
// leaf or label (those are in the MGT text).
static Json::Value arm_summary_json(const ArmResult &arm, const Strategy &st) {
    Json::Value j;
    arm_certificate_json(&j, arm, st);
    j["counters"] = counters_json(arm);
    j["branch_events_total"] = uint_json(arm.branch_events_total);
    Json::Value counts;
    size_t merges = 0;
    uint64_t bases = 0, max_bp = 0;
    for (const Segment &s : arm.segments) {
        merges += s.parents.size() > 1;
        bases += s.length_bp;
    }
    std::array<uint32_t, kNumEndReasons> ends {};
    for (const GrowthBin &g : arm.growth) {
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            ends[r] += g.label_ends[r];
        }
    }
    // leaves by their path reason; "semantic" when every label ended for its own reason
    std::map<std::string, uint64_t> by_reason;
    for (const PathResult &p : arm.paths) {
        max_bp = std::max(max_bp, p.length_bp);
        by_reason[p.path_reason ? to_string(*p.path_reason) : "semantic"]++;
    }
    counts["segments"] = uint_json(arm.segments.size());
    counts["leaves"] = uint_json(arm.paths.size());
    counts["splits"] = uint_json(arm.splits.size());
    counts["merges"] = uint_json(merges);
    counts["runs"] = uint_json(arm.runs.size());
    counts["bases"] = uint_json(bases);
    counts["max_bp"] = uint_json(max_bp);
    // label (run) ends by reason, silent switch-source ends included: the growth bins'
    Json::Value label_ends(Json::objectValue);
    for (size_t r = 0; r < kNumEndReasons; ++r) {
        if (ends[r]) label_ends[to_string(static_cast<EndReason>(r))] = ends[r];
    }
    counts["label_ends"] = std::move(label_ends);
    Json::Value leaves(Json::objectValue);
    for (const auto &[reason, n] : by_reason) {
        leaves[reason] = uint_json(n);
    }
    counts["leaves_by_reason"] = std::move(leaves);
    j["counts"] = std::move(counts);
    return j;
}

Json::Value seed_result_to_json(const SeedResult &r, const Strategy &st, const std::string &detail, bool timing) {
    const bool graphlet = detail == "graphlet";
    Json::Value j;
    Json::Value seed;
    seed["seed_id"] = r.seed_id;
    seed["validated_seed_id"] = r.validated_seed_id;
    seed["seed_id_mismatch"] = r.seed_id_mismatch;
    seed["length_bp"] = uint_json(r.length_bp);
    seed["num_kmers"] = uint_json(r.num_kmers);
    if (!graphlet) {
        // the graphlet carries these in its S, X and L records
        Json::Value labels(Json::arrayValue);
        for (size_t i = 0; i < r.num_seed_labels; ++i) labels.append(r.label_dict[i].name);
        seed["labels"] = std::move(labels);
    }
    // the permitted set was derived from the seed: how many labels support it in full,
    // and what the cap cut (never silently)
    seed["labels_from_seed"] = r.labels_from_seed;
    seed["labels_supporting_total"] = uint_json(r.labels_supporting_total);
    seed["labels_dropped"] = uint_json(r.labels_dropped);
    seed["labels_dropped_digest"] = r.labels_dropped_digest;
    if (graphlet) {
        seed["num_labels"] = uint_json(r.label_dict.size());
        seed["num_seed_labels"] = uint_json(r.num_seed_labels);
    } else {
        Json::Value dropped(Json::arrayValue);
        for (const auto &d : r.dropped_labels) {
            Json::Value dj;
            dj["label"] = d.name;
            dj["reason"] = d.reason;
            Json::Value runs(Json::arrayValue);
            for (const auto &[a, b] : d.runs) {
                Json::Value iv(Json::arrayValue);
                iv.append(uint_json(a)); iv.append(uint_json(b));
                runs.append(std::move(iv));
            }
            dj["runs"] = std::move(runs);
            // the runs are k-mer PRESENCE runs on the seed, also under support: trace, where
            // the label was dropped for lacking a coordinate-consecutive occurrence
            dj["runs_kind"] = "presence";
            dropped.append(std::move(dj));
        }
        seed["dropped_labels"] = std::move(dropped);
    }
    j["seed"] = std::move(seed);
    // the caps that limited the seed itself (§7.0); process_traverse_request adds the
    // server clamps this seed ran into
    Json::Value lims(Json::arrayValue);
    if (r.labels_dropped) {
        lims.append(limitation("seed_labels", "labels.max_seed_labels", uint_json(st.max_seed_labels),
                               uint_json(r.labels_supporting_total),
                               std::to_string(r.labels_dropped) + " label(s) carrying the whole seed "
                               "were not taken (labels_dropped_digest identifies them): walks only they "
                               "carry are missing; raise the knob, or traverse the complement with an "
                               "explicit list"));
    }
    if (st.support == Support::TRACE) {
        size_t columns = 0;
        for (const auto &l : r.label_dict) {
            columns += l.kind == LabelKind::COLUMN;
        }
        if (columns) {
            lims.append(limitation("trace_record_boundaries", "labels.seed_label_kind", "column",
                                   uint_json(columns),
                                   "a column label's trace follows consecutive coordinates of the whole "
                                   "column, so it cannot detect the boundary between two records whose "
                                   "coordinates happen to be adjacent; use header labels (seed_label_kind: "
                                   "\"header\") where the index has a CoordToHeader"));
        }
    }
    j["limitations"] = std::move(lims);
    j["label_mode"] = to_string(st.label_mode);
    if (!graphlet) {
        Json::Value dict(Json::arrayValue);
        for (const auto &l : r.label_dict) {
            Json::Value lj;
            lj["name"] = l.name;
            lj["kind"] = to_string(l.kind);
            // the stable cross-retrieval key (LabelRef): a name alone can be ambiguous
            // (one accession in two columns); seq_id only means something for a header
            lj["column"] = uint_json(l.column);
            if (l.kind == LabelKind::HEADER)
                lj["seq_id"] = uint_json(l.seq_id);
            dict.append(std::move(lj));
        }
        j["label_dict"] = std::move(dict);
    }
    Json::Value arms(Json::objectValue);
    for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
        const auto &a = r.arms[static_cast<size_t>(arm)];
        if (!a.requested) continue;
        arms[to_string(arm)] = graphlet ? arm_summary_json(a, st) : arm_to_json(a, st, detail);
    }
    j["arms"] = std::move(arms);
    // §7.0: the guarantees are independent, so the outcome states each on its own axis
    // instead of folding them into one value that would read "partial" for a complete
    // walk with cut diagnostics, or "complete" for walks whose label lists were cut
    j["outcome"] = outcome_of(j, false);
    if (!graphlet) {
        // derived from the runs (constrain) or the recorded sets (annotate) by the
        // graphlet's reader, so not in its summary
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
                sj["runs"] = std::move(runs);
                lj[to_string(arm)] = std::move(sj);
            }
            summary.append(std::move(lj));
        }
        j["label_summary"] = std::move(summary);
    }
    Json::Value ann;
    ann["access_path"] = r.access_path;
    ann["keys_mapped"] = uint_json(r.annotation_counters.keys_mapped);
    ann["rows_requested"] = uint_json(r.annotation_counters.rows_requested);
    ann["direct_reads"] = uint_json(r.annotation_counters.direct_reads);
    j["annotation"] = std::move(ann);
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
        j["timing"] = std::move(t);
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
            lj["trace_breaks"] = std::move(tb);
        }
        labels.append(std::move(lj));
    }
    j["labels"] = std::move(labels);
    if (p.labels_truncated) {
        Json::Value t;
        t["kept"] = uint_json(p.labels_truncated->kept);
        t["total"] = uint_json(p.labels_truncated->total);
        t["min_kept_kmers"] = uint_json(p.labels_truncated->min_kept_kmers);
        t["max_dropped_kmers"] = uint_json(p.labels_truncated->max_dropped_kmers);
        t["dropped_full_length"] = uint_json(p.labels_truncated->dropped_full_length);
        j["labels_truncated"] = std::move(t);
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
        cj["labels"] = std::move(ls);
        cands.append(std::move(cj));
    }
    j["candidates"] = std::move(cands);
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
        s["policy"] = std::move(pol);
        Json::Value counts;
        counts["candidates"] = uint_json(sel->num_candidates);
        counts["eligible"] = uint_json(sel->num_eligible);
        counts["selected"] = uint_json(sel->seeds.size());
        s["counts"] = std::move(counts);
        Json::Value seeds(Json::arrayValue);
        for (const auto &seed : sel->seeds) {
            Json::Value sj;
            sj["seed_id"] = seed.seed_id;
            sj["sequence"] = seed.sequence;
            sj["kmer_interval"] = interval_json(seed.kmers);
            Json::Value ls(Json::arrayValue);
            for (const auto &l : seed.labels) ls.append(l);
            sj["labels"] = std::move(ls);
            Json::Value pop;
            pop["supporting_total"] = uint_json(seed.population.supporting_total);
            pop["included"] = uint_json(seed.population.included);
            pop["dropped_count"] = uint_json(seed.population.dropped_count);
            if (seed.population.dropped_count) {
                pop["dropped_digest"] = seed.population.dropped_digest;
                Json::Value dl(Json::arrayValue);
                for (const auto &d : seed.population.dropped) dl.append(d);
                pop["dropped"] = std::move(dl);
            }
            sj["label_population"] = std::move(pop);
            Json::Value ov(Json::arrayValue);
            for (size_t o : seed.overlaps_with) ov.append(uint_json(o));
            sj["overlaps_with"] = std::move(ov);
            if (!seed.labels_not_covering.empty()) {
                Json::Value nc(Json::arrayValue);
                for (const auto &l : seed.labels_not_covering) nc.append(l);
                sj["labels_not_covering"] = std::move(nc);
            }
            seeds.append(std::move(sj));
        }
        s["seeds"] = std::move(seeds);
        j["selection"] = std::move(s);
    }
    return j;
}

static IndexIdentity identity_or_default(const IndexIdentity *given, const LabelOracle &oracle);

Json::Value capabilities_to_json(const LabelOracle &oracle, const std::string &release,
                                 const IndexIdentity *identity) {
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
    c["cost_models_available"] = std::move(models);
    Json::Value modes(Json::arrayValue);
    modes.append("constrain"); modes.append("annotate");
    c["label_modes"] = std::move(modes);
    c["direct_access"] = oracle.supports_direct();
    c["release"] = release;
    // the retrieval format (spec §7.5): `detail: graphlet` embeds MGT text of this version
    c["graphlet_format"] = kGraphletFormatVersion;
    Json::Value details(Json::arrayValue);
    for (const char *d : { "summary", "tree", "full", "graphlet" }) details.append(d);
    c["detail_levels"] = std::move(details);
    // which index this is (DESIGN-traverse-graphlet.md §3.1): labels are joined across
    // retrievals only on an equal index_fp; null fp = no manifest, joins unverifiable;
    // index_meta_fp is a negative check only
    const IndexIdentity id = identity_or_default(identity, oracle);
    c["index_ns"] = id.name.empty() ? Json::Value() : Json::Value(id.name);
    c["index_fp"] = id.fp.empty() ? Json::Value() : Json::Value(id.fp);
    c["index_meta_fp"] = id.meta_fp;
    return c;
}


// ---------------------------------------------------------------- MGT v1 codec

namespace mgt {

namespace {

[[noreturn]] void fail(const std::string &what, std::string_view token) {
    throw FormatError(what + ": '" + std::string(token) + "'");
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// a canonical unsigned decimal: no sign, no leading zero, fits 64 bits
bool parse_uint(std::string_view s, uint64_t *x) {
    if (s.empty() || (s.size() > 1 && s[0] == '0'))
        return false;
    uint64_t v = 0;
    for (char c : s) {
        if (!is_digit(c))
            return false;
        const uint64_t d = c - '0';
        if (v > (std::numeric_limits<uint64_t>::max() - d) / 10)
            return false;
        v = v * 10 + d;
    }
    *x = v;
    return true;
}

void append_uint(std::string &out, uint64_t x) {
    char buf[24];
    auto res = std::to_chars(buf, buf + sizeof(buf), x);
    out.append(buf, res.ptr);
}

bool is_continuation(unsigned char c) { return (c & 0xC0) == 0x80; }

// well-formed UTF-8 (no overlong forms, no surrogates, at most U+10FFFF)
bool valid_utf8(std::string_view s) {
    size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = s[i];
        size_t n;
        uint32_t cp;
        if (c < 0x80) { i++; continue; }
        if ((c & 0xE0) == 0xC0) { n = 1; cp = c & 0x1F; }
        else if ((c & 0xF0) == 0xE0) { n = 2; cp = c & 0x0F; }
        else if ((c & 0xF8) == 0xF0) { n = 3; cp = c & 0x07; }
        else return false;
        if (i + n >= s.size())
            return false;
        for (size_t j = 1; j <= n; ++j) {
            const unsigned char d = s[i + j];
            if (!is_continuation(d))
                return false;
            cp = (cp << 6) | (d & 0x3F);
        }
        if ((n == 1 && cp < 0x80) || (n == 2 && cp < 0x800) || (n == 3 && cp < 0x10000)
                || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF))
            return false;
        i += n + 1;
    }
    return true;
}

int hex_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;   // lowercase hex is not canonical
}

const char kHexUpper[] = "0123456789ABCDEF";

} // namespace

std::string encode_float(double x) {
    if (x == 0)
        return "0";     // -0 too
    if (std::isinf(x) && x > 0)
        return "inf";
    if (!(x > 0))
        throw std::logic_error("MGT floats are non-negative (got a negative number or NaN)");
    // the shortest round-trip digits d1..dn and the exponent E of d1.d2..dn x 10^E
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof(buf), x, std::chars_format::scientific);
    std::string_view sci(buf, res.ptr - buf);
    const size_t e = sci.find('e');
    std::string digits;
    for (char c : sci.substr(0, e)) {
        if (c != '.')
            digits.push_back(c);
    }
    while (digits.size() > 1 && digits.back() == '0')
        digits.pop_back();
    long exponent = std::strtol(std::string(sci.substr(e + 1)).c_str(), nullptr, 10);
    // the decimal point sits after |point| digits (point <= 0: before leading zeros)
    const long point = exponent + 1;
    const long n = digits.size();
    std::string out;
    if (point <= 0) {
        out = "0.";
        out.append(-point, '0');
        out += digits;
    } else if (point >= n) {
        out = digits;
        out.append(point - n, '0');
    } else {
        out = digits.substr(0, point) + "." + digits.substr(point);
    }
    return out;
}

double decode_float(std::string_view token) {
    if (token == "inf")
        return std::numeric_limits<double>::infinity();
    // (0|[1-9][0-9]*)(\.[0-9]*[1-9])?
    size_t i = 0;
    if (token.empty() || !is_digit(token[0]))
        fail("not a canonical float", token);
    if (token[0] == '0') {
        i = 1;
    } else {
        while (i < token.size() && is_digit(token[i])) i++;
    }
    if (i < token.size()) {
        if (token[i] != '.' || i + 1 == token.size() || token.back() == '0')
            fail("not a canonical float", token);
        for (size_t j = i + 1; j < token.size(); ++j) {
            if (!is_digit(token[j]))
                fail("not a canonical float", token);
        }
    }
    // the C locale's strtod is correctly rounded; whether the token is THE spelling of
    // the double it denotes is decided by re-encoding it (rejects the exact expansion of
    // 0.1, 1e23 written by to_chars(fixed), an overflow to inf and an underflow to 0)
    const std::string s(token);
    const double x = std::strtod(s.c_str(), nullptr);
    if (encode_float(x) != token)
        fail("not the shortest spelling of its float", token);
    return x;
}

std::string encode_ranges(const std::vector<uint64_t> &ids) {
    std::string out;
    if (ids.empty())
        return ".";
    for (size_t i = 0; i < ids.size(); ) {
        if (i && ids[i] <= ids[i - 1])
            throw std::logic_error("RANGES: ids must be strictly ascending");
        size_t j = i;
        while (j + 1 < ids.size() && ids[j + 1] == ids[j] + 1) j++;
        if (!out.empty())
            out.push_back(',');
        append_uint(out, ids[i]);
        if (j > i) {
            out.push_back('-');
            append_uint(out, ids[j]);
        }
        i = j + 1;
    }
    return out;
}

std::vector<uint64_t> decode_ranges(std::string_view token, std::optional<uint64_t> bound) {
    std::vector<uint64_t> out;
    if (token == ".")
        return out;
    if (token.empty())
        fail("empty RANGES", token);
    size_t pos = 0;
    while (true) {
        size_t comma = token.find(',', pos);
        std::string_view item = token.substr(pos, comma == std::string_view::npos
                                                      ? std::string_view::npos : comma - pos);
        const size_t dash = item.find('-');
        uint64_t a, b;
        if (dash == std::string_view::npos) {
            if (!parse_uint(item, &a))
                fail("malformed RANGES", token);
            b = a;
        } else if (!parse_uint(item.substr(0, dash), &a) || !parse_uint(item.substr(dash + 1), &b)
                       || b <= a) {
            fail("malformed RANGES", token);
        }
        // ascending, and every maximal run collapsed: an item may not continue the last
        if (!out.empty() && a <= out.back() + 1)
            fail("RANGES not ascending or not collapsed", token);
        if (bound && b >= *bound)
            fail("RANGES id beyond the id space", token);
        for (uint64_t x = a; ; ++x) {
            out.push_back(x);
            if (x == b) break;
        }
        if (comma == std::string_view::npos)
            break;
        pos = comma + 1;
    }
    return out;
}

std::string encode_setexpr(const std::vector<uint64_t> &ids, const std::vector<uint64_t> *base) {
    std::string explicit_form = encode_ranges(ids);
    if (!base)
        return explicit_form;
    std::vector<uint64_t> removed, added;
    std::set_difference(base->begin(), base->end(), ids.begin(), ids.end(),
                        std::back_inserter(removed));
    std::set_difference(ids.begin(), ids.end(), base->begin(), base->end(),
                        std::back_inserter(added));
    std::string delta = "!";
    if (!removed.empty())
        delta += encode_ranges(removed);
    if (!added.empty())
        delta += "+" + encode_ranges(added);
    return delta.size() < explicit_form.size() ? delta : explicit_form;
}

std::vector<uint64_t> decode_setexpr(std::string_view token, const std::vector<uint64_t> *base,
                                     std::optional<uint64_t> bound) {
    if (token.empty() || token[0] != '!')
        return decode_ranges(token, bound);
    if (!base)
        fail("a delta SETEXPR needs a base", token);
    std::string_view rest = token.substr(1);
    const size_t plus = rest.find('+');
    std::string_view r = rest.substr(0, plus);
    std::string_view a = plus == std::string_view::npos ? std::string_view() : rest.substr(plus + 1);
    // an empty part is omitted, never written '.'
    if (r == "." || (plus != std::string_view::npos && (a.empty() || a == ".")))
        fail("malformed SETEXPR", token);
    std::vector<uint64_t> removed = r.empty() ? std::vector<uint64_t>() : decode_ranges(r, bound);
    std::vector<uint64_t> added = a.empty() ? std::vector<uint64_t>() : decode_ranges(a, bound);
    // a mismatch with the base means the reader derived another base than the writer
    if (!std::includes(base->begin(), base->end(), removed.begin(), removed.end()))
        fail("SETEXPR removes ids its base does not hold", token);
    std::vector<uint64_t> kept, out;
    std::set_difference(base->begin(), base->end(), removed.begin(), removed.end(),
                        std::back_inserter(kept));
    for (uint64_t x : added) {
        if (std::binary_search(base->begin(), base->end(), x))
            fail("SETEXPR adds ids its base holds", token);
    }
    std::set_union(kept.begin(), kept.end(), added.begin(), added.end(), std::back_inserter(out));
    return out;
}

std::string pct_escape(std::string_view raw) {
    std::string out;
    out.reserve(raw.size());
    for (char c : raw) {
        switch (c) {
            case '%': out += "%25"; break;
            case '\n': out += "%0A"; break;
            case '\r': out += "%0D"; break;
            default: out.push_back(c);
        }
    }
    return out;
}

std::string pct_unescape(std::string_view token) {
    std::string out;
    out.reserve(token.size());
    for (size_t i = 0; i < token.size(); ++i) {
        const char c = token[i];
        if (c == '\n' || c == '\r')
            fail("raw line break in a free-text field", token);
        if (c != '%') {
            out.push_back(c);
            continue;
        }
        std::string_view esc = token.substr(i + 1, 2);
        if (esc == "25") out.push_back('%');
        else if (esc == "0A") out.push_back('\n');
        else if (esc == "0D") out.push_back('\r');
        else fail("not one of %25 %0A %0D", token);
        i += 2;
    }
    return out;
}

std::string front_code(std::string_view previous, std::string_view name) {
    size_t p = 0;
    while (p < previous.size() && p < name.size() && previous[p] == name[p]) p++;
    // the suffix must start at a code point of its own (valid UTF-8 alone)
    while (p > 0 && p < name.size() && is_continuation(name[p])) p--;
    std::string out;
    append_uint(out, p);
    out.push_back(' ');
    out += pct_escape(name.substr(p));
    return out;
}

std::string to_valid_utf8(std::string_view s) {
    if (valid_utf8(s))
        return std::string(s);
    std::string out;
    out.reserve(s.size() + 8);
    size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = s[i];
        if (c < 0x80) {
            out.push_back(c);
            ++i;
            continue;
        }
        // the lead byte's sequence length and the range of its FIRST continuation byte,
        // which is what excludes overlong forms, surrogates and code points > U+10FFFF
        // (Unicode Table 3-7); C0, C1, F5..FF and a stray continuation lead nothing
        size_t n = 0;
        unsigned char lo = 0x80, hi = 0xBF;
        if (c >= 0xC2 && c <= 0xDF) { n = 1; }
        else if (c == 0xE0) { n = 2; lo = 0xA0; }
        else if ((c >= 0xE1 && c <= 0xEC) || c == 0xEE || c == 0xEF) { n = 2; }
        else if (c == 0xED) { n = 2; hi = 0x9F; }
        else if (c == 0xF0) { n = 3; lo = 0x90; }
        else if (c >= 0xF1 && c <= 0xF3) { n = 3; }
        else if (c == 0xF4) { n = 3; hi = 0x8F; }
        // the maximal subpart: the lead byte and the continuation bytes valid so far
        size_t j = i + 1;
        for (size_t m = 0; m < n && j < s.size(); ++m) {
            const unsigned char d = s[j];
            if (d < (m == 0 ? lo : 0x80) || d > (m == 0 ? hi : 0xBF))
                break;
            ++j;
        }
        if (n && j - i == n + 1) {
            out.append(s.substr(i, n + 1));
        } else {
            out += "\xEF\xBF\xBD";   // U+FFFD for the whole ill-formed subpart
        }
        i = j;
    }
    return out;
}

std::string front_decode(std::string_view previous, std::string_view token) {
    const size_t space = token.find(' ');
    uint64_t p;
    if (space == std::string_view::npos || !parse_uint(token.substr(0, space), &p))
        fail("malformed front-coded name", token);
    if (p > previous.size() || (p < previous.size() && is_continuation(previous[p])))
        fail("prefix beyond the previous name or inside a character", token);
    return std::string(previous.substr(0, p)) + pct_unescape(token.substr(space + 1));
}

std::string encode_kvalue(const KValue &value) {
    switch (value.type) {
        case KValue::INTEGER: {
            std::string out = "i:";
            append_uint(out, value.integer);
            return out;
        }
        case KValue::FLOAT:
            return "f:" + encode_float(value.real);
        case KValue::STRING: {
            // a fixed field is ASCII without spaces, and ',' separates <extra> items
            std::string out = "s:";
            for (unsigned char c : value.string) {
                if (c < 0x21 || c > 0x7E || c == '%' || c == ',') {
                    out.push_back('%');
                    out.push_back(kHexUpper[c >> 4]);
                    out.push_back(kHexUpper[c & 15]);
                } else {
                    out.push_back(c);
                }
            }
            return out;
        }
        case KValue::UNLIMITED:
            return "u";
    }
    return "u";
}

KValue decode_kvalue(std::string_view token) {
    KValue v;
    if (token == "u") {
        v.type = KValue::UNLIMITED;
        return v;
    }
    if (token.size() < 2 || token[1] != ':')
        fail("malformed K value", token);
    std::string_view body = token.substr(2);
    switch (token[0]) {
        case 'i':
            v.type = KValue::INTEGER;
            if (!parse_uint(body, &v.integer))
                fail("malformed K integer", token);
            return v;
        case 'f':
            v.type = KValue::FLOAT;
            v.real = decode_float(body);
            return v;
        case 's':
            v.type = KValue::STRING;
            for (size_t i = 0; i < body.size(); ++i) {
                const unsigned char c = body[i];
                if (c != '%') {
                    if (c < 0x21 || c > 0x7E || c == ',')
                        fail("a K string byte that must be escaped", token);
                    v.string.push_back(c);
                    continue;
                }
                const bool whole = i + 2 < body.size();
                const int hi = whole ? hex_value(body[i + 1]) : -1;
                const int lo = whole ? hex_value(body[i + 2]) : -1;
                if (hi < 0 || lo < 0)
                    fail("malformed escape in a K string", token);
                const unsigned char b = hi * 16 + lo;
                // an escape is canonical only for a byte that cannot be written raw
                if (b >= 0x21 && b <= 0x7E && b != '%' && b != ',')
                    fail("escape of a byte that is written raw", token);
                v.string.push_back(b);
                i += 2;
            }
            if (!valid_utf8(v.string))
                fail("a K string that is not UTF-8", token);
            return v;
    }
    fail("unknown K value type", token);
}

} // namespace mgt

// ---------------------------------------------------------------- index identity

std::string sha256_hex(std::string_view data) {
    static const uint32_t K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    };
    uint32_t h[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
    auto rotr = [](uint32_t x, int n) { return (x >> n) | (x << (32 - n)); };
    std::string msg(data);
    const uint64_t bits = static_cast<uint64_t>(data.size()) * 8;
    msg.push_back(static_cast<char>(0x80));
    while (msg.size() % 64 != 56)
        msg.push_back('\0');
    for (int i = 7; i >= 0; --i)
        msg.push_back(static_cast<char>((bits >> (8 * i)) & 0xff));
    for (size_t chunk = 0; chunk < msg.size(); chunk += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            const auto *p = reinterpret_cast<const unsigned char*>(msg.data() + chunk + 4 * i);
            w[i] = (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3];
        }
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t t1 = hh + S1 + ch + K[i] + w[i];
            const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2 = S0 + maj;
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    static const char hex[] = "0123456789abcdef";
    std::string out;
    for (uint32_t x : h) {
        for (int i = 28; i >= 0; i -= 4)
            out.push_back(hex[(x >> i) & 15]);
    }
    return out;
}

bool valid_index_name(const std::string &name) {
    if (name.empty())
        return false;
    for (unsigned char c : name) {
        if (!std::isalnum(c) && c != '.' && c != '_' && c != '-')
            return false;
    }
    return true;
}

// FNV-1a-64 over what the loaded index says about itself. Two indexes that differ here
// are different; two that agree may still differ (swapped column memberships keep every
// count and name), which is why this is never treated as identity.
std::string index_meta_fingerprint(const LabelOracle &oracle) {
    auto field = [](const std::string &s, uint64_t h) {
        // length-prefixed, so that no column name can imitate a field boundary
        return fnv1a64(std::to_string(s.size()) + ":" + s + "\n", h);
    };
    uint64_t h = fnv1a64("mgt-index-meta 1\n");
    h = field("k=" + std::to_string(oracle.get_k()), h);
    h = field(std::string("regime=") + to_string(oracle.regime()), h);
    h = field("alphabet=" + oracle.graph().alphabet(), h);
    h = field("rows=" + std::to_string(oracle.num_rows()), h);
    h = field("columns=" + std::to_string(oracle.num_columns()), h);
    for (uint64_t c = 0; c < oracle.num_columns(); ++c) {
        h = field(oracle.column_name(c), h);
    }
    h = field(std::string("coordinates=") + (oracle.has_coordinates() ? "1" : "0"), h);
    h = field(std::string("coord_to_header=") + (oracle.coord_to_header() ? "1" : "0"), h);
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return buf;
}

std::string index_manifest_fingerprint(const std::string &manifest_path,
                                       const std::vector<std::string> &loaded) {
    std::ifstream in(manifest_path);
    if (!in.good())
        throw std::runtime_error("index manifest " + manifest_path + ": cannot be read");
    Json::Value json;
    Json::CharReaderBuilder reader;
    std::string errs;
    if (!Json::parseFromStream(reader, in, &json, &errs))
        throw std::runtime_error("index manifest " + manifest_path + ": invalid JSON: " + errs);
    auto bad = [&](const std::string &what) {
        return std::runtime_error("index manifest " + manifest_path + ": " + what);
    };
    if (!json.isObject() || !json["files"].isArray() || json["files"].empty())
        throw bad("expected an object with a non-empty array 'files' of {path, size, sha256}");
    std::map<std::string, std::pair<uint64_t, std::string>> files;
    for (const Json::Value &f : json["files"]) {
        if (!f.isObject() || !f["path"].isString() || !f["size"].isUInt64() || !f["sha256"].isString())
            throw bad("every entry of 'files' is {path: string, size: integer, sha256: string}");
        const std::string path = f["path"].asString();
        const std::string digest = f["sha256"].asString();
        if (path.empty() || path.find_first_of("\t\n\r") != std::string::npos)
            throw bad("a file path is empty or holds a tab or line break");
        if (digest.size() != 64 || digest.find_first_not_of("0123456789abcdef") != std::string::npos)
            throw bad("the sha256 of " + path + " is not 64 lowercase hex digits");
        if (!files.emplace(path, std::make_pair(f["size"].asUInt64(), digest)).second)
            throw bad("the file " + path + " is listed twice");
    }
    // A manifest of another bundle must not lend its identity: every loaded file that
    // exists as named is listed (by base name) with its size. Hashing it is what the
    // manifest saves the server from (hundreds of GB at start-up).
    auto base_name = [](const std::string &p) {
        const size_t slash = p.find_last_of('/');
        return slash == std::string::npos ? p : p.substr(slash + 1);
    };
    for (const std::string &path : loaded) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f.good())
            continue;   // named without its extension: the loader resolves it, not us
        const uint64_t size = static_cast<uint64_t>(f.tellg());
        bool listed = false;
        for (const auto &[p, entry] : files) {
            listed |= base_name(p) == base_name(path) && entry.first == size;
        }
        if (!listed) {
            throw bad("the loaded file " + path + " (" + std::to_string(size)
                      + " bytes) is not listed with that size: the manifest describes another index");
        }
    }
    // the canonical file list, in byte order of path (std::map's order)
    std::string canonical;
    for (const auto &[path, entry] : files) {
        canonical += path + "\t" + std::to_string(entry.first) + "\t" + entry.second + "\n";
    }
    const std::string fp = sha256_hex(canonical);
    // a manifest may state its own fingerprint for humans; it must be the one computed
    if (json.isMember("index_fp") && json["index_fp"].asString() != fp)
        throw bad("its index_fp " + json["index_fp"].asString() + " is not the digest of its files ("
                  + fp + "): the file list was edited");
    return fp;
}

// the identity a response states: |given|, or none with the meta fingerprint computed now
static IndexIdentity identity_or_default(const IndexIdentity *given, const LabelOracle &oracle) {
    if (given && !given->meta_fp.empty())
        return *given;
    IndexIdentity id = given ? *given : IndexIdentity();
    id.meta_fp = index_meta_fingerprint(oracle);
    return id;
}

IndexIdentity index_identity(const Config &config, const graph::AnnotatedDBG &anno_graph) {
    IndexIdentity id;
    if (!config.index_name.empty() && !valid_index_name(config.index_name)) {
        throw std::runtime_error("--index-name '" + config.index_name + "': expected "
                                 "[A-Za-z0-9._-]+");
    }
    id.name = config.index_name;
    if (!config.index_manifest.empty()) {
        std::vector<std::string> loaded { config.infbase };
        loaded.insert(loaded.end(), config.infbase_annotators.begin(),
                      config.infbase_annotators.end());
        id.fp = index_manifest_fingerprint(config.index_manifest, loaded);
    }
    id.meta_fp = index_meta_fingerprint(LabelOracle(anno_graph));
    return id;
}

// ---------------------------------------------------------------- the graphlet (MGT v1)

namespace {

// End reasons as MGT codes (DESIGN-traverse-graphlet.md §2.1), in EndReason order. 'Y'
// (resource_limit) is reserved for the budgets of §14, which no walk produces yet.
const char kReasonCodes[] = "DLBRUVJTXSPNOMW";
static_assert(sizeof(kReasonCodes) == kNumEndReasons + 1, "one MGT code per EndReason");

char reason_code(EndReason reason) { return kReasonCodes[static_cast<size_t>(reason)]; }

// The walker's text qualifier of a label end, kept as the code's second letter (today's
// JSON replaces the enum by the text; the graphlet keeps both). Each text belongs to one
// reason, which is checked rather than trusted: a new text must get a letter here.
char qualifier_code(const std::string &text, EndReason reason) {
    struct Qualifier { const char *text; EndReason reason; char code; };
    static const Qualifier table[] = {
        { "minority", EndReason::BRANCH, 'm' },
        { "below_min_labels", EndReason::BRANCH, 'b' },
        { "split_limit", EndReason::BRANCH, 's' },
        { "hairpin", EndReason::DEAD_END, 'h' },
        { "superseded", EndReason::LABEL_LOST, 's' },
        { "switch_sources", EndReason::LABEL_LOST, 'x' },
    };
    if (text.empty())
        return 0;
    for (const Qualifier &q : table) {
        if (text == q.text && reason == q.reason)
            return q.code;
    }
    throw std::logic_error("graphlet: the label end text '" + text + "' of reason "
                           + to_string(reason) + " has no MGT qualifier");
}

[[noreturn]] void unrepresentable(const std::string &what) {
    throw std::logic_error("graphlet: " + what);
}

void put_uint(std::string &out, uint64_t x) {
    char buf[24];
    auto res = std::to_chars(buf, buf + sizeof(buf), x);
    out.append(buf, res.ptr);
}

// RANGES straight into |out| (the codec's rule, without a temporary vector<uint64_t>)
template <class T>
void put_ranges(std::string &out, const std::vector<T> &ids) {
    if (ids.empty()) {
        out.push_back('.');
        return;
    }
    for (size_t i = 0; i < ids.size(); ) {
        if (i && ids[i] <= ids[i - 1])
            unrepresentable("a label or segment list is not strictly ascending");
        size_t j = i;
        while (j + 1 < ids.size() && ids[j + 1] == ids[j] + 1) j++;
        if (i)
            out.push_back(',');
        put_uint(out, ids[i]);
        if (j > i) {
            out.push_back('-');
            put_uint(out, ids[j]);
        }
        i = j + 1;
    }
}

// SETEXPR against |base| (null: none): the shorter of explicit and delta, a tie explicit
template <class T>
void put_setexpr(std::string &out, const std::vector<T> &ids, const std::vector<T> *base) {
    if (!base) {
        put_ranges(out, ids);
        return;
    }
    // the differences below need a sorted base; a bad one must not yield a wrong delta
    if (std::adjacent_find(base->begin(), base->end(), std::greater_equal<T>()) != base->end())
        unrepresentable("a SETEXPR base is not strictly ascending");
    std::string explicit_form, delta = "!";
    put_ranges(explicit_form, ids);
    std::vector<T> removed, added;
    std::set_difference(base->begin(), base->end(), ids.begin(), ids.end(),
                        std::back_inserter(removed));
    std::set_difference(ids.begin(), ids.end(), base->begin(), base->end(),
                        std::back_inserter(added));
    if (!removed.empty())
        put_ranges(delta, removed);
    if (!added.empty()) {
        delta.push_back('+');
        put_ranges(delta, added);
    }
    out += delta.size() < explicit_form.size() ? delta : explicit_form;
}

template <class T>
std::vector<T> set_union_of(const std::vector<std::vector<T>> &sets) {
    std::vector<T> out, next;
    for (const auto &s : sets) {
        next.clear();
        std::set_union(out.begin(), out.end(), s.begin(), s.end(), std::back_inserter(next));
        out.swap(next);
    }
    return out;
}

class GraphletWriter {
  public:
    GraphletWriter(const SeedResult &r, const Seed &seed, const Strategy &st,
                   const GraphletContext &ctx, const Json::Value &result_json)
          : r_(r), seed_(seed), st_(st), ctx_(ctx), json_(result_json),
            annotate_(st.label_mode == LabelMode::ANNOTATE) {}

    std::string write(size_t *lines) {
        header();
        seed_record();
        dropped();
        dictionary();
        outcome();
        limitations();
        for (Arm a : { Arm::LEFT, Arm::RIGHT }) {
            const ArmResult &arm = r_.arms[static_cast<size_t>(a)];
            if (arm.requested)
                write_arm(arm);
        }
        out_ += "Z ";
        put_uint(out_, lines_ + 1);
        end();
        if (lines)
            *lines = lines_;
        return std::move(out_);
    }

  private:
    // what the records of one arm are derived from, computed once per arm
    struct ArmIndex {
        std::vector<int64_t> path_of;          // leaf segment -> path, -1 elsewhere
        std::vector<int> split_of;             // segment -> Split::ambiguous, -1: no split
        std::vector<char> first_char;          // split child -> its branch's base
        // the runs anchored on each segment (ids into ArmResult::runs), so that a segment's
        // end set and a leaf's end labels cost its own runs, not all of them
        std::vector<std::vector<uint32_t>> anchored;
        // LABEL_END events by (segment, label, at): each is claimed by exactly one run
        std::map<std::tuple<size_t, LabelId, uint64_t>, std::pair<const Event*, bool>> label_ends;
    };

    // ---- output primitives
    void sp() { out_.push_back(' '); }
    void end() { out_.push_back('\n'); lines_++; }
    void num(uint64_t x) { put_uint(out_, x); }
    void real(double x) { out_ += mgt::encode_float(x); }
    void flag(bool x) { out_.push_back(x ? '1' : '0'); }
    // a fixed field: non-empty printable ASCII without spaces
    void token(std::string_view t, const char *what) {
        if (t.empty())
            unrepresentable(std::string("empty ") + what);
        for (unsigned char c : t) {
            if (c < 0x21 || c > 0x7E)
                unrepresentable(std::string(what) + " '" + std::string(t) + "' is not a token");
        }
        out_ += t;
    }
    void base(char c) {
        if (static_cast<unsigned char>(c) < 0x21 || static_cast<unsigned char>(c) > 0x7E)
            unrepresentable("a base outside printable ASCII");
        out_.push_back(c);
    }

    void header() {
        out_ += "H mgt ";
        put_uint(out_, kGraphletFormatVersion);
        sp(); num(ctx_.k);
        sp(); token(ctx_.regime, "regime");
        sp(); token(ctx_.alphabet, "alphabet");
        sp(); out_.push_back(annotate_ ? 'a' : 'c');
        sp(); out_.push_back(st_.support == Support::TRACE ? 't' : 'k');
        sp(); out_.push_back(st_.merge_reconverge ? 'm' : 'k');
        sp(); num(st_.max_labels_per_node);
        sp(); num(st_.continuation_bp);
        sp(); num(ctx_.seed_index);
        // the orientation rule: every G stores its bases in walking order (§2.4)
        out_ += " walk ";
        if (ctx_.identity.name.empty()) out_.push_back('*'); else token(ctx_.identity.name, "index_ns");
        sp();
        if (ctx_.identity.fp.empty()) out_.push_back('*'); else token(ctx_.identity.fp, "index_fp");
        sp(); token(ctx_.identity.meta_fp, "index_meta_fp");
        end();
    }

    void seed_record() {
        // the seed as the walker validated it: case-mapped, request orientation
        std::string sequence = seed_.sequence;
#if ! _DNA_CASE_SENSITIVE_GRAPH
        for (char &c : sequence) {
            c = std::toupper(static_cast<unsigned char>(c));
        }
#endif
        if (sequence.size() != r_.length_bp)
            unrepresentable("the seed sequence does not match the result's length");
        out_ += "S ";
        token(r_.validated_seed_id, "validated_seed_id");
        sp(); num(r_.length_bp);
        sp(); num(r_.num_kmers);
        sp(); num(r_.num_seed_labels);
        sp(); token(sequence, "seed sequence");
        end();
    }

    void dropped() {
        for (const DroppedLabel &d : r_.dropped_labels) {
            out_ += "X ";
            token(d.reason, "dropped label reason");
            sp();
            if (d.runs.empty())
                out_.push_back('.');
            for (size_t i = 0; i < d.runs.size(); ++i) {
                if (i)
                    out_.push_back(',');
                num(d.runs[i].first);
                out_.push_back('-');
                num(d.runs[i].second);
            }
            sp();
            // process_traverse_request made every name valid UTF-8; a name that is not
            // would be replaced inside the JSON string, so graphlet_bytes and the front
            // coding (byte counts) would no longer describe the transported body
            if (!mgt::valid_utf8(d.name))
                unrepresentable("a dropped label name that is not UTF-8");
            out_ += mgt::pct_escape(d.name);
            end();
        }
    }

    void dictionary() {
        const std::string *previous = nullptr;
        static const std::string kNone;
        for (const LabelRef &l : r_.label_dict) {
            out_ += "L ";
            out_.push_back(l.kind == LabelKind::HEADER ? 'h' : 'c');
            sp(); num(l.column);
            sp();
            if (l.kind == LabelKind::HEADER) {
                num(l.seq_id);
            } else if (l.seq_id) {
                unrepresentable("a column label with a sequence id");
            } else {
                out_.push_back('*');
            }
            sp();
            if (!mgt::valid_utf8(l.name))
                unrepresentable("a label name that is not UTF-8");
            out_ += mgt::front_code(previous ? *previous : kNone, l.name);
            end();
            previous = &l.name;
        }
    }

    void outcome() {
        const Json::Value &o = json_["outcome"];
        auto code = [&](const char *axis, std::initializer_list<std::pair<const char*, char>> values) {
            const std::string v = o[axis].asString();
            for (const auto &[name, c] : values) {
                if (v == name) {
                    out_.push_back(c);
                    return;
                }
            }
            unrepresentable(std::string("outcome.") + axis + " '" + v + "'");
        };
        out_ += "O ";
        code("walks", { { "complete", 'c' }, { "partial", 'p' }, { "failed", 'f' } });
        sp();
        code("branch_diagnostics", { { "complete", 'c' }, { "cut", 'x' } });
        sp();
        code("label_evidence", { { "complete", 'c' }, { "lower_bound", 'l' }, { "qualified", 'q' } });
        sp();
        code("delivery", { { "inline", 'i' }, { "spooled", 's' }, { "paged", 'p' } });
        end();
    }

    // a JSON limitation value as a typed K value (§2.2 v5.1): the knob's own type
    mgt::KValue kvalue(const Json::Value &v, const std::string &what) {
        mgt::KValue k;
        if (v.isString()) {
            // "unlimited" is how every limit knob spells "no limit" (limit_json)
            k.type = v.asString() == "unlimited" ? mgt::KValue::UNLIMITED : mgt::KValue::STRING;
            k.string = v.asString();
        } else if (v.type() == Json::uintValue || (v.type() == Json::intValue && v.asInt64() >= 0)) {
            k.type = mgt::KValue::INTEGER;
            k.integer = v.asUInt64();
        } else if (v.type() == Json::realValue) {
            k.type = mgt::KValue::FLOAT;
            k.real = v.asDouble();
        } else {
            unrepresentable("limitation field " + what + " is not a K value");
        }
        return k;
    }

    void limitation(char arm, const Json::Value &l) {
        static const std::set<std::string> kFixed {
            "kind", "knob", "limit", "observed", "complete_to_bp", "effect"
        };
        out_ += "K ";
        out_.push_back(arm);
        sp(); token(l["kind"].asString(), "limitation kind");
        sp(); token(l["knob"].asString(), "limitation knob");
        sp(); out_ += mgt::encode_kvalue(kvalue(l["limit"], "limit"));
        sp(); out_ += mgt::encode_kvalue(kvalue(l["observed"], "observed"));
        sp();
        if (!l.isMember("complete_to_bp")) {
            out_.push_back('*');
        } else if (l["complete_to_bp"].isUInt64()) {
            num(l["complete_to_bp"].asUInt64());
        } else {
            unrepresentable("limitation complete_to_bp is not an integer");
        }
        sp();
        // further fields (server_limit, label_ends, cause, ...), in byte order of name
        bool any = false;
        for (const std::string &name : l.getMemberNames()) {
            if (kFixed.count(name))
                continue;
            for (char c : name) {
                if (!(std::islower(static_cast<unsigned char>(c)) || c == '_'))
                    unrepresentable("limitation field name '" + name + "'");
            }
            out_ += any ? "," : "";
            out_ += name;
            out_.push_back('=');
            out_ += mgt::encode_kvalue(kvalue(l[name], name));
            any = true;
        }
        if (!any)
            out_.push_back('.');
        sp();
        if (!l["effect"].isString())
            unrepresentable("a limitation without an effect");
        out_ += mgt::pct_escape(l["effect"].asString());
        end();
    }

    // seed-level ones first, then each requested arm's, each list in its JSON order
    void limitations() {
        for (const Json::Value &l : json_["limitations"]) {
            limitation('*', l);
        }
        for (Arm a : { Arm::LEFT, Arm::RIGHT }) {
            if (!r_.arms[static_cast<size_t>(a)].requested)
                continue;
            for (const Json::Value &l : json_["arms"][to_string(a)]["limitations"]) {
                limitation(a == Arm::LEFT ? 'l' : 'r', l);
            }
        }
    }

    ArmIndex index_arm(const ArmResult &arm) {
        const auto &segs = arm.segments;
        ArmIndex ix;
        ix.path_of.assign(segs.size(), -1);
        ix.split_of.assign(segs.size(), -1);
        ix.first_char.assign(segs.size(), 0);
        ix.anchored.resize(segs.size());
        for (size_t s = 0; s < segs.size(); ++s) {
            if (segs[s].id != s)
                unrepresentable("segment ids are not ordinals");
            for (size_t p : segs[s].parents) {
                // parents first (new_segment appends), so ids order the DAG
                if (p >= s)
                    unrepresentable("a parent created after its child");
            }
        }
        // paths: leaf ordinal in segment-id order, chain through parents[0] (finalize)
        size_t last_leaf = 0;
        for (size_t i = 0; i < arm.paths.size(); ++i) {
            const PathResult &p = arm.paths[i];
            if (p.id != i || p.segments.empty())
                unrepresentable("path ids are not ordinals");
            const size_t leaf = p.segments.back();
            if (leaf >= segs.size() || ix.path_of[leaf] >= 0 || (i && leaf <= last_leaf))
                unrepresentable("paths are not one per leaf in segment order");
            last_leaf = leaf;
            ix.path_of[leaf] = i;
            size_t cur = leaf;
            for (size_t j = p.segments.size(); j-- > 0; ) {
                if (p.segments[j] != cur)
                    unrepresentable("a path is not its leaf's first-parent chain");
                if (j)
                    cur = segs[cur].parents.empty() ? SIZE_MAX : segs[cur].parents[0];
            }
            if (!segs[cur].parents.empty()
                    || p.length_bp != segs[leaf].from_bp + segs[leaf].length_bp)
                unrepresentable("a path does not start at the root or has another length");
        }
        // splits: the children with one parent, grouped by parent, in (at_bp, first child)
        std::pair<uint64_t, size_t> last_split { 0, 0 };
        for (size_t i = 0; i < arm.splits.size(); ++i) {
            const Split &sp = arm.splits[i];
            if (sp.segment >= segs.size() || ix.split_of[sp.segment] >= 0 || sp.children.size() < 2)
                unrepresentable("a split that its children do not reproduce");
            const Segment &parent = segs[sp.segment];
            if (sp.at_bp != parent.from_bp + parent.length_bp)
                unrepresentable("a split not at its segment's end");
            std::pair<uint64_t, size_t> key { sp.at_bp, sp.children[0] };
            if (i && !(last_split < key))
                unrepresentable("splits not ordered by (at_bp, first child)");
            last_split = key;
            ix.split_of[sp.segment] = sp.ambiguous;
            // labels_before: the count at the split node
            size_t before = annotate_ ? (parent.label_sets.empty() ? parent.labels_start_total
                                                                  : parent.label_sets.back().labels_total)
                                      : parent.labels_end.size();
            if (sp.labels_before != before || sp.branches.size() != sp.children.size())
                unrepresentable("a split's labels_before or branches are not derivable");
            for (size_t j = 0; j < sp.children.size(); ++j) {
                const size_t c = sp.children[j];
                const SplitBranch &br = sp.branches[j];
                if (c >= segs.size() || br.segment != c || segs[c].parents.size() != 1
                        || segs[c].parents[0] != sp.segment || (j && c <= sp.children[j - 1]))
                    unrepresentable("a split's children are not its single-parent children");
                const Segment &child = segs[c];
                const size_t cut = std::min(child.labels_start.size(), st_.max_labels_per_node);
                if (br.labels_distinct != child.labels_start_total
                        || br.labels != std::vector<LabelId>(child.labels_start.begin(),
                                                             child.labels_start.begin() + cut))
                    unrepresentable("a branch's labels are not its child's entry labels");
                ix.first_char[c] = br.ch;
            }
        }
        size_t children = 0;
        for (const Segment &s : segs) {
            children += s.parents.size() == 1;
        }
        size_t split_children = 0;
        for (const Split &sp : arm.splits) {
            split_children += sp.children.size();
        }
        if (children != split_children)
            unrepresentable("a single-parent segment that no split lists");
        // label ends and reconvergences, derived from R and G
        for (const Segment &s : segs) {
            size_t joins = 0;
            for (const Event &ev : s.events) {
                if (ev.type == EventType::LABEL_END) {
                    auto [it, fresh] = ix.label_ends.emplace(
                            std::make_tuple(s.id, ev.label, ev.at_bp), std::make_pair(&ev, false));
                    if (!fresh)
                        unrepresentable("two label ends of one label at one position");
                } else if (ev.type == EventType::RECONVERGE) {
                    joins++;
                    if (ev.at_bp != s.from_bp || ev.labels.size() != s.parents.size()
                            || !std::equal(ev.labels.begin(), ev.labels.end(), s.parents.begin()))
                        unrepresentable("a reconverge event that its segment does not reproduce");
                }
            }
            if (joins != (s.parents.size() > 1))
                unrepresentable("a merge without exactly one reconverge event");
        }
        for (size_t i = 0; i < arm.runs.size(); ++i) {
            if (arm.runs[i].segment >= segs.size())
                unrepresentable("a run without an anchor");
            ix.anchored[arm.runs[i].segment].push_back(i);
        }
        return ix;
    }

    void write_arm(const ArmResult &arm) {
        ArmIndex ix = index_arm(arm);
        const char side = arm.arm == Arm::LEFT ? 'l' : 'r';
        // ---- A
        size_t merges = 0;
        uint64_t bases = 0;
        for (const Segment &s : arm.segments) {
            merges += s.parents.size() > 1;
            bases += s.length_bp;
        }
        out_ += "A ";
        out_.push_back(side);
        sp(); out_.push_back(arm.status == ArmResult::COMPLETE ? 'c'
                             : arm.status == ArmResult::TRUNCATED ? 't' : 'p');
        sp(); num(arm.complete_to_bp);
        sp(); out_.push_back(std::string(arm.completeness_scope) == "united_history" ? 'u' : 'p');
        sp(); num(arm.frontier_live_paths);
        sp(); num(arm.frontier_live_labels);
        sp(); flag(arm.frontier_live_labels_exact);
        sp(); num(arm.max_labels_at_node);
        sp(); num(arm.nodes_labels_truncated);
        sp();
        bool first = true;
        for (const auto &[name, value] : arm_counters(arm)) {
            out_ += first ? "" : ",";
            out_ += name;
            out_.push_back('=');
            num(value);
            first = false;
        }
        sp();
        if (arm.cap_trigger) {
            const CapTrigger &c = *arm.cap_trigger;
            out_.push_back(reason_code(c.reason));
            out_.push_back(','); num(c.at_bp);
            out_.push_back(','); num(c.segment);
            out_.push_back(','); num(c.live_paths);
            out_.push_back(','); num(c.live_labels);
            out_.push_back(','); flag(c.live_labels_exact);
            out_.push_back(','); real(c.demand);
        } else {
            out_.push_back('*');
        }
        sp(); num(arm.branch_events_total);
        sp();
        if (arm.branch_events_complete_to_bp == std::numeric_limits<uint64_t>::max()) {
            out_.push_back('*');
        } else {
            num(arm.branch_events_complete_to_bp);
        }
        sp(); num(arm.segments.size());
        sp(); num(arm.runs.size());
        sp(); num(arm.paths.size());
        sp(); num(arm.splits.size());
        sp(); num(merges);
        sp(); num(bases);
        end();
        // ---- B: the growth bins
        for (const GrowthBin &g : arm.growth) {
            out_ += "B ";
            num(g.from_bp);
            sp(); num(g.max_live_paths);
            sp(); num(g.max_live_labels);
            sp(); num(g.max_live_pairs);
            sp(); flag(g.live_labels_exact);
            sp(); num(g.steps);
            sp(); num(g.divergences);
            sp(); num(g.ambiguous_branches);
            sp(); num(g.splits);
            sp(); num(g.reconvergences);
            sp(); num(g.bubbles);
            sp(); num(g.tips);
            sp(); num(g.blocked_repeat);
            sp();
            bool any = false;
            for (size_t r = 0; r < kNumEndReasons; ++r) {
                if (!g.label_ends[r])
                    continue;
                out_ += any ? "," : "";
                out_.push_back(kReasonCodes[r]);
                out_.push_back(':');
                num(g.label_ends[r]);
                any = true;
            }
            if (!any)
                out_.push_back('.');
            end();
        }
        // ---- V: the branch events as stored (the only record of a successor not taken)
        for (const BranchEvent &b : arm.branch_events) {
            out_ += "V ";
            num(b.at_bp);
            sp(); num(b.segment);
            sp();
            if (b.chars.empty())
                out_.push_back('.');
            for (char c : b.chars) {
                base(c);
            }
            sp();
            if (b.labels_per_successor.empty())
                out_.push_back('.');
            for (size_t i = 0; i < b.labels_per_successor.size(); ++i) {
                out_ += i ? "," : "";
                num(b.labels_per_successor[i]);
            }
            sp(); put_ranges(out_, b.ambiguous);
            sp(); put_ranges(out_, b.dropped);
            sp();
            if (b.refused.empty())
                out_.push_back('.');
            for (size_t i = 0; i < b.refused.size(); ++i) {
                const BranchEvent::Refusal &rf = b.refused[i];
                out_ += i ? ";" : "";
                base(rf.ch);
                out_.push_back(':');
                for (const char *c = rf.cause; *c; ++c) {
                    if (!(std::islower(static_cast<unsigned char>(*c)) || *c == '_'))
                        unrepresentable(std::string("refusal cause '") + rf.cause + "'");
                }
                out_ += rf.cause;
                out_.push_back(':');
                put_ranges(out_, rf.labels);
            }
            end();
        }
        // ---- per segment: G P* E* T? C?
        for (const Segment &s : arm.segments) {
            write_segment(arm, s, ix);
        }
        // ---- R: the runs, in ArmResult::runs order (run ids keep their meaning)
        std::vector<double> needed;
        for (size_t i = 0; i < arm.runs.size(); ++i) {
            write_run(arm.runs[i], ix, &needed);
        }
        for (const auto &[key, ev] : ix.label_ends) {
            if (!ev.second)
                unrepresentable("a label end that no run claims");
        }
        // needed_budgets is derived from the B runs (a histogram: order is not information)
        std::vector<double> stored = arm.needed_budgets;
        std::sort(needed.begin(), needed.end());
        std::sort(stored.begin(), stored.end());
        if (needed != stored)
            unrepresentable("needed_budgets that the loss_budget ends do not reproduce");
    }

    // the walking-order bases of |s| (the left arm's Segment::sequence is natural)
    std::string walk_bases(const ArmResult &arm, const Segment &s) const {
        std::string w = s.sequence;
        if (arm.arm == Arm::LEFT)
            std::reverse(w.begin(), w.end());
        return w;
    }

    std::vector<LabelId> chronological_end(const ArmResult &arm, const Segment &s,
                                           const ArmIndex &ix) const {
        // ends before switch-ins at one position; ends anchored at the segment's last node
        // stay (labels_end is taken there)
        const uint64_t last = s.from_bp + s.length_bp;
        std::vector<std::tuple<uint64_t, int, LabelId>> changes;
        for (uint32_t i : ix.anchored[s.id]) {
            if (arm.runs[i].to_bp < last)
                changes.emplace_back(arm.runs[i].to_bp, 0, arm.runs[i].label);
        }
        for (const Event &ev : s.events) {
            if (ev.type == EventType::SWITCH)
                changes.emplace_back(ev.at_bp, 1, ev.to);
        }
        std::sort(changes.begin(), changes.end());
        std::vector<LabelId> cur = s.labels_start;
        for (const auto &[at, kind, label] : changes) {
            auto it = std::lower_bound(cur.begin(), cur.end(), label);
            const bool present = it != cur.end() && *it == label;
            if (kind == 0 && present) {
                cur.erase(it);
            } else if (kind == 1 && !present) {
                cur.insert(it, label);
            }
        }
        return cur;
    }

    void write_segment(const ArmResult &arm, const Segment &s, const ArmIndex &ix) {
        const auto &segs = arm.segments;
        out_ += "G ";
        if (s.parents.empty())
            out_.push_back('*');
        for (size_t j = 0; j < s.parents.size(); ++j) {
            out_ += j ? "," : "";
            num(s.parents[j]);
        }
        sp(); num(s.from_bp);
        sp(); num(s.length_bp);
        sp();
        // ---- entry, against parent[0]'s end set (a merge: the union of its parents')
        std::vector<LabelId> entry_base;
        if (s.parents.size() == 1) {
            entry_base = segs[s.parents[0]].labels_end;
        } else if (s.parents.size() > 1) {
            std::vector<std::vector<LabelId>> ends;
            for (size_t p : s.parents) {
                ends.push_back(segs[p].labels_end);
            }
            entry_base = set_union_of(ends);
        }
        bool entry_rule = false;
        if (!annotate_ && s.parents.empty()) {
            // the root holds every seed label (walker init_arm)
            std::vector<LabelId> seeds(r_.num_seed_labels);
            std::iota(seeds.begin(), seeds.end(), 0);
            entry_rule = s.labels_start == seeds;
        } else if (!annotate_ && s.parents.size() > 1) {
            entry_rule = s.labels_start == set_union_of(s.labels_via_parent);
        }
        if (entry_rule) {
            out_.push_back('*');
        } else {
            put_setexpr(out_, s.labels_start, s.parents.empty() ? nullptr : &entry_base);
        }
        sp();
        if (s.labels_start_total == s.labels_start.size()) {
            out_.push_back('*');
        } else {
            num(s.labels_start_total);
        }
        sp();
        // ---- end, against the entry set
        const std::vector<LabelId> end_rule = !annotate_ ? chronological_end(arm, s, ix)
            : s.label_sets.empty() ? s.labels_start : s.label_sets.back().labels;
        if (end_rule == s.labels_end) {
            out_.push_back('*');
        } else {
            put_setexpr(out_, s.labels_end, &s.labels_start);
        }
        sp();
        // ---- partition: none without a merge; one empty list per parent in annotate mode
        std::vector<std::vector<LabelId>> partition_rule;
        if (s.parents.size() > 1 && annotate_)
            partition_rule.resize(s.parents.size());
        if (s.labels_via_parent == partition_rule) {
            out_.push_back('*');
        } else {
            if (s.parents.size() < 2 || s.labels_via_parent.size() != s.parents.size())
                unrepresentable("labels_via_parent without one list per parent of a merge");
            for (size_t j = 0; j < s.labels_via_parent.size(); ++j) {
                out_ += j ? "|" : "";
                put_ranges(out_, s.labels_via_parent[j]);
            }
        }
        sp();
        if (ix.split_of[s.id] < 0) out_.push_back('*'); else flag(ix.split_of[s.id]);
        sp();
        // ---- first_base (only without bases) and the bases in walking order
        if (st_.sequences) {
            if (s.sequence.size() != s.length_bp)
                unrepresentable("a segment's sequence is not length_bp bases");
            const std::string w = walk_bases(arm, s);
            if (ix.first_char[s.id] && w[0] != ix.first_char[s.id])
                unrepresentable("a branch's base is not its child's first base");
            out_ += "* ";
            if (w.empty())
                out_.push_back('.');
            for (char c : w) {
                base(c);
            }
        } else {
            if (ix.first_char[s.id]) base(ix.first_char[s.id]); else out_.push_back('*');
            out_ += " *";
        }
        end();
        // ---- P: annotate presence runs, each against the previous one (the first: entry)
        const std::vector<LabelId> *previous = &s.labels_start;
        for (const LabelSetRun &run : s.label_sets) {
            out_ += "P ";
            num(run.from_bp);
            sp(); num(run.to_bp);
            sp(); num(run.labels_total);
            sp(); put_setexpr(out_, run.labels, previous);
            end();
            previous = &run.labels;
        }
        // ---- E: the stored events; label ends come from R, reconvergences from G
        for (const Event &ev : s.events) {
            write_event(ev);
        }
        if (ix.path_of[s.id] >= 0)
            write_leaf(arm, s, arm.paths[ix.path_of[s.id]], ix);
    }

    void write_event(const Event &ev) {
        auto labels = [&]() {
            num(ev.labels_total);
            sp(); put_ranges(out_, ev.labels);
        };
        switch (ev.type) {
            case EventType::LABEL_END:
            case EventType::RECONVERGE:
                return;
            case EventType::SWITCH:
                out_ += "E "; num(ev.at_bp);
                out_ += " s "; num(ev.label);
                sp(); num(ev.to);
                sp(); real(ev.cost);
                break;
            case EventType::BLOCKED:
                if (!ev.text.empty())
                    unrepresentable("a blocked event with a text");
                out_ += "E "; num(ev.at_bp);
                out_ += " b "; base(ev.ch);
                sp(); out_.push_back(reason_code(ev.reason));
                sp(); labels();
                break;
            case EventType::HAIRPIN:
                if (!ev.text.empty() && ev.text != "followed")
                    unrepresentable("a hairpin event text '" + ev.text + "'");
                out_ += "E "; num(ev.at_bp);
                out_ += " h "; base(ev.ch);
                sp(); labels();
                sp(); out_.push_back(ev.text.empty() ? 's' : 'f');
                break;
            case EventType::REVISIT:
                if (ev.labels.size() != 1 || (ev.text == "same_distance") != (ev.length_bp == 0)
                        || (!ev.text.empty() && ev.text != "same_distance"))
                    unrepresentable("a revisit event the record cannot hold");
                out_ += "E "; num(ev.at_bp);
                out_ += " v "; num(ev.labels[0]);
                sp();
                if (ev.length_bp) num(ev.length_bp); else out_.push_back('=');
                break;
            case EventType::TIP:
                out_ += "E "; num(ev.at_bp);
                out_ += " t "; base(ev.ch);
                sp(); num(ev.length_bp);
                break;
            case EventType::BUBBLE:
                out_ += "E "; num(ev.at_bp);
                out_ += " u "; num(ev.length_bp);
                sp(); token(ev.text, "bubble alleles");
                break;
        }
        end();
    }

    void write_leaf(const ArmResult &arm, const Segment &s, const PathResult &p,
                    const ArmIndex &ix) {
        // end_labels: the runs anchored here ending at the leaf's last node, ascending
        // label; end_reasons their reasons
        std::vector<std::pair<LabelId, uint32_t>> derived;
        std::array<uint32_t, kNumEndReasons> reasons {};
        for (uint32_t i : ix.anchored[s.id]) {
            const LabelRun &run = arm.runs[i];
            if (run.to_bp == p.length_bp) {
                derived.emplace_back(run.label, i);
                if (!run.ended)
                    unrepresentable("a merge closure at a leaf");
                reasons[static_cast<size_t>(run.end_reason)]++;
            }
        }
        std::sort(derived.begin(), derived.end());
        if (derived.size() != p.end_labels.size() || reasons != p.end_reasons)
            unrepresentable("a leaf whose end labels its runs do not reproduce");
        for (size_t i = 0; i < derived.size(); ++i) {
            if (derived[i].first != p.end_labels[i].label || derived[i].second != p.end_labels[i].run)
                unrepresentable("a leaf whose end labels its runs do not reproduce");
        }
        out_ += "T ";
        if (p.path_reason) out_.push_back(reason_code(*p.path_reason)); else out_.push_back('*');
        sp();
        bool any = false;
        for (const LabelEnd &e : p.end_labels) {
            if (!e.loss && !e.branches && !e.route_bp)
                continue;
            out_ += any ? "," : "";
            num(e.label);
            out_.push_back(':'); real(e.loss);
            out_.push_back(':'); num(e.branches);
            out_.push_back(':'); num(e.route_bp);
            any = true;
        }
        if (!any)
            out_.push_back('.');
        end();
        if (p.continuation)
            write_continuation(arm, s, p, *p.continuation);
    }

    // C: the walker's labels are primary, the spelling is derived from the bases
    //   right arm: the last n bases of seed + natural(right flank)
    //   left arm:  the first n bases of natural(left flank) + seed
    void write_continuation(const ArmResult &arm, const Segment &leaf, const PathResult &p,
                            const Continuation &c) {
        const uint64_t n = c.sequence.size();
        std::string derived;
        bool derivable = n == 0;
        if (st_.sequences && n <= r_.length_bp + p.length_bp) {
            // the last min(n, flank) walked bases, leaf -> root, as make_continuation does:
            // O(n) per leaf, never the whole chain (a comb-shaped trie has long ones)
            std::string rev;
            const uint64_t from_walk = std::min<uint64_t>(n, p.length_bp);
            for (size_t i = p.segments.size(); i-- > 0 && rev.size() < from_walk; ) {
                // walking order is the stored order on the right, reversed on the left
                const std::string &w = arm.segments[p.segments[i]].sequence;
                auto take = [&](auto begin, auto end) {
                    for (auto it = begin; it != end && rev.size() < from_walk; ++it) {
                        rev.push_back(*it);
                    }
                };
                if (arm.arm == Arm::RIGHT) {
                    take(w.rbegin(), w.rend());
                } else {
                    take(w.begin(), w.end());
                }
            }
            const uint64_t from_seed = n - from_walk;
            std::string seq = seed_.sequence;
#if ! _DNA_CASE_SENSITIVE_GRAPH
            for (char &ch : seq) {
                ch = std::toupper(static_cast<unsigned char>(ch));
            }
#endif
            if (arm.arm == Arm::RIGHT) {
                derived = seq.substr(seq.size() - from_seed);
                derived.append(rev.rbegin(), rev.rend());
            } else {
                derived = rev + seq.substr(0, from_seed);
            }
            derivable = true;
        }
        out_ += "C ";
        num(n);
        sp(); real(c.loss_used);
        sp(); num(c.branches_used);
        sp(); put_setexpr(out_, c.labels, &leaf.labels_end);
        if (!derivable || derived != c.sequence) {
            sp();
            token(c.sequence, "continuation sequence");
        }
        end();
    }

    void write_run(const LabelRun &run, ArmIndex &ix, std::vector<double> *needed) {
        out_ += "R ";
        num(run.segment);
        sp(); num(run.label);
        sp(); num(run.from_bp);
        sp(); num(run.to_bp);
        sp();
        const Event *ev = nullptr;
        if (run.ended) {
            auto it = ix.label_ends.find({ run.segment, run.label, run.to_bp });
            if (it != ix.label_ends.end()) {
                ev = it->second.first;
                it->second.second = true;
                if (ev->reason != run.end_reason)
                    unrepresentable("a run whose label end has another reason");
                out_.push_back(reason_code(run.end_reason));
                if (char q = qualifier_code(ev->text, ev->reason))
                    out_.push_back(q);
            } else if (run.end_reason == EndReason::LABEL_LOST) {
                // a switch source that went on only under other names: no event
                out_ += "Lw";
            } else {
                unrepresentable("a run ended without a label end event");
            }
        } else {
            out_.push_back('m');   // closed by the merge at to_bp
        }
        sp(); num(run.route_bp);
        sp();
        if (run.entered_by_switch) {
            num(run.from_label);
            out_.push_back(':');
            real(run.switch_cost);
        } else if (run.from_label || run.switch_cost) {
            unrepresentable("a seed-entered run with a switch source");
        } else {
            out_.push_back('*');
        }
        sp();
        if (run.prev_run == UINT32_MAX) out_.push_back('*'); else num(run.prev_run);
        sp();
        if (ev) num(ev->structural_successors); else out_.push_back('*');
        sp(); num(run.branches);
        sp(); real(run.loss);
        if (ev && run.end_reason == EndReason::LOSS_BUDGET) {
            sp(); real(ev->needed_budget);
            needed->push_back(ev->needed_budget);
        }
        end();
    }

    const SeedResult &r_;
    const Seed &seed_;
    const Strategy &st_;
    const GraphletContext &ctx_;
    const Json::Value &json_;
    const bool annotate_;
    std::string out_;
    size_t lines_ = 0;
};

} // namespace

std::string graphlet_text(const SeedResult &result, const Seed &seed, const Strategy &strategy,
                          const GraphletContext &context, const Json::Value &result_json,
                          size_t *lines) {
    return GraphletWriter(result, seed, strategy, context, result_json).write(lines);
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

// The request-level clamps (strategy.clamped) that bound THIS seed's result, added to
// its `limitations` as kind server_clamp: a request cannot raise past them, so the agent
// must know which ones it ran into — a lowered derived-set cap that cut the derived set,
// a lowered time budget that tripped, and a budget raised from zero (the walk ran
// although the request asked for none). The seed_labels entry also gets the server's
// maximum, beyond which raising its knob does nothing.
static void state_server_clamps(Json::Value *rj, const SeedResult &r, const Json::Value &clamped) {
    Json::Value &lims = (*rj)["limitations"];
    for (const Json::Value &c : clamped) {
        const std::string field = c["field"].asString();
        const bool lowered = c["effective"].asDouble() < c["requested"].asDouble();
        bool affected = false;
        if (field == "labels.max_seed_labels") {
            affected = r.labels_dropped > 0;
            for (Json::Value &l : lims) {
                if (affected && l["kind"].asString() == "seed_labels")
                    l["server_limit"] = uint_json(c["effective"].asUInt64());
            }
        } else if (field == "bounds.time_budget_ms") {
            affected = !lowered;
            for (const ArmResult &a : r.arms) {
                affected |= a.requested && ended_by(a, EndReason::TIME_BUDGET);
            }
        }
        if (!affected)
            continue;
        // strategy.clamped already carries each value in its knob's type
        lims.append(limitation("server_clamp", field, c["effective"], c["requested"],
                               lowered ? "the server lowered the requested value to its maximum and "
                                         "this seed ran into it; a request cannot raise it further"
                                       : "the server raised the requested value (it also bounds the "
                                         "derivation of a permitted set); the walk ran under it"));
    }
    // the outcome is read off the final limitations (server_clamp itself is in no class;
    // what a clamp caused is stated by the walk_domain or seed_labels it led to)
    (*rj)["outcome"] = outcome_of(*rj, false);
}

// The result of a seed whose permitted set could not be DERIVED (§6.1 step 4): no arms,
// `outcome.walks: failed`, the message as `error`, and the cause as a `derivation` limitation
// naming the request field that would get past it (§7.0) — so that an agent acts on the
// knob instead of parsing the message. |clamped| supplies `server_limit` when the server
// clamped that knob: raising it beyond the server's value does nothing.
static Json::Value failed_seed_to_json(const Seed &seed, const SeedDerivationError &e,
                                       const Strategy &st, const Json::Value &clamped) {
    auto server_limit = [&](const std::string &knob, Json::Value *l) {
        for (const Json::Value &c : clamped) {
            if (c["field"].asString() == knob)
                (*l)["server_limit"] = c["effective"];
        }
    };
    std::string knob;
    Json::Value limit, observed;
    std::string effect;
    switch (e.cause()) {
        case SeedDerivationError::NO_CARRIER:
            knob = "seeds[].sequence";
            limit = uint_json(static_cast<uint64_t>(e.limit()));
            observed = uint_json(static_cast<uint64_t>(e.observed()));
            effect = "no label carries every k-mer of the seed (limit: its k-mers; observed: the "
                     "k-mers read when no candidate was left), so no permitted set exists to "
                     "traverse under: shorten the seed to a stretch one label carries (resolve "
                     "reports each label's runs) or name the labels explicitly";
            break;
        case SeedDerivationError::NO_TRACE_CARRIER:
            knob = "support";
            limit = "trace";
            observed = uint_json(static_cast<uint64_t>(e.observed()));
            effect = "the observed number of labels carry every k-mer of the seed, none of them "
                     "as one coordinate-consecutive occurrence: set support to \"kmer\" (presence), "
                     "shorten the seed, or name the labels explicitly";
            break;
        case SeedDerivationError::TOO_WIDE:
            knob = "seeds[].sequence";
            limit = uint_json(static_cast<uint64_t>(e.limit()));
            observed = uint_json(static_cast<uint64_t>(e.observed()));
            effect = "the narrowest of the seed's first 64 k-mers has more annotation entries "
                     "(columns plus coordinates) than the derivation materialises (64 x "
                     "labels.max_seed_labels, at least 65536): start the seed in a less repetitive "
                     "k-mer, raise labels.max_seed_labels, or name the labels explicitly";
            break;
        case SeedDerivationError::TIME_BUDGET:
            knob = "bounds.time_budget_ms";
            limit = e.limit();
            observed = e.observed();
            effect = "the time budget ran out while deriving the permitted set, before any "
                     "extension: raise the budget, shorten the seed, or name the labels explicitly "
                     "(an explicit list reads only its own columns)";
            break;
        case SeedDerivationError::AMBIGUOUS_HEADER:
            knob = "labels.seed_label_kind";
            limit = "header";
            observed = e.subject();
            effect = "the derived header label (observed) occurs in more than one annotation "
                     "column or names a column, so the derived list could not be resubmitted as "
                     "explicit labels: set the knob to \"column\" or name the labels explicitly";
            break;
        case SeedDerivationError::OVER_SEED_LABEL_CAP:
            knob = "labels.max_seed_labels";
            limit = uint_json(static_cast<uint64_t>(e.limit()));
            observed = uint_json(static_cast<uint64_t>(e.observed()));
            effect = "under `exhaustive` the derived set is refused rather than cut (every walk "
                     "of a dropped carrier would be missing): raise the knob to at least observed "
                     "or name the labels explicitly";
            break;
    }
    Json::Value rj;
    Json::Value sj;
    sj["seed_id"] = seed.seed_id;
    sj["length_bp"] = uint_json(seed.sequence.size());
    sj["labels_from_seed"] = true;
    rj["seed"] = std::move(sj);
    rj["error"] = e.what();
    Json::Value lims(Json::arrayValue);
    Json::Value d = limitation("derivation", knob, limit, observed, effect);
    d["cause"] = to_string(e.cause());
    server_limit(knob, &d);
    lims.append(std::move(d));
    if (e.labels_cut()) {
        // the trace check runs after the cap: a carrier the cap cut was never checked
        Json::Value s = limitation(
                "seed_labels", "labels.max_seed_labels", uint_json(st.max_seed_labels),
                uint_json(static_cast<uint64_t>(e.observed())),
                std::to_string(e.labels_cut()) + " label(s) carrying every k-mer of the seed were "
                "cut before the trace check and never checked; one of them may carry the seed "
                "as one occurrence: raise the knob");
        server_limit("labels.max_seed_labels", &s);
        lims.append(std::move(s));
    }
    rj["limitations"] = std::move(lims);
    // no walk was made, so nothing was cut on the other axes — except the carriers a cap
    // cut before the trace check (seed_labels), whose evidence is then missing
    rj["outcome"] = outcome_of(rj, true);
    return rj;
}

Json::Value process_traverse_request(const Json::Value &json,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const std::string &release,
                                     const TraverseLimits &limits,
                                     const IndexIdentity *identity) {
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
    // the values keep the type of their knob: an integer knob echoes as an integer (as a
    // double it printed 10000.0, which a client reads as a float where the field is an
    // integer), a time budget as the number it was parsed as
    auto clamp = [&](const char *field, Json::Value requested, Json::Value effective) {
        Json::Value c;
        c["field"] = field;
        c["requested"] = std::move(requested);
        c["effective"] = std::move(effective);
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
        clamp("labels.max_seed_labels", uint_json(req.strategy.max_seed_labels),
              uint_json(limits.max_seed_labels));
        req.strategy.max_seed_labels = limits.max_seed_labels;
    }

    LabelOracle oracle(anno_graph);
    // checked here, where k is known: a continuation shorter than k is not a valid seed,
    // so the promise that continuations are resubmittable (§7.1) would not hold
    const uint64_t k = oracle.get_k();
    if (req.strategy.continuation_bp > 0 && req.strategy.continuation_bp < k) {
        throw InvalidRequest("strategy.output.continuation_bp: "
                             + std::to_string(req.strategy.continuation_bp) + " is between 1 and "
                             "k - 1 (k = " + std::to_string(k) + "), which would give continuations "
                             "shorter than k that are not valid traverse input; use 0 (no "
                             "continuation sequence) or at least " + std::to_string(k));
    }
    Timer timer;
    // computed once per request: the H record of every graphlet repeats it
    GraphletContext graphlet;
    graphlet.identity = identity_or_default(identity, oracle);
    graphlet.k = k;
    graphlet.regime = to_string(oracle.regime());
    graphlet.alphabet = oracle.graph().alphabet();
    Json::Value out;
    out["release"] = release;
    out["capabilities"] = capabilities_to_json(oracle, release, &graphlet.identity);
    Json::Value strategy = strategy_to_json(req.strategy, req.cost, req.detail, req.timing);
    strategy["clamped"] = clamped;
    out["strategy"] = std::move(strategy);
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
    // building the response (the JSON tree and the graphlet text), the part a request
    // can make slow without walking more: reported apart from the walk
    double serialize_seconds = 0;
    for (size_t i = 0; i < req.seeds.size(); ++i) {
        const Seed &seed = req.seeds[i];
        std::vector<std::string> dict = seed.labels;
        dict.insert(dict.end(), req.strategy.extra.begin(), req.strategy.extra.end());
        LabelChangeCost cost = make_cost(req.cost, dict);
        try {
            SeedResult r = traverse_seed(oracle, seed, req.strategy, cost, release);
            r.annotation_counters = per_seed(r.annotation_counters);
            // Names come from FASTA headers and file names, which need not be UTF-8. The
            // JSON writer would replace an ill-formed sequence its own way (and differently
            // for different bytes), so detail: full, the summary and the MGT body, whose
            // graphlet_bytes and L prefix lengths count bytes, would disagree. Made valid
            // once, here, every output carries the same bytes.
            for (LabelRef &l : r.label_dict)
                l.name = mgt::to_valid_utf8(l.name);
            for (DroppedLabel &d : r.dropped_labels)
                d.name = mgt::to_valid_utf8(d.name);
            Timer serialize;
            Json::Value rj = seed_result_to_json(r, req.strategy, req.detail, req.timing);
            state_server_clamps(&rj, r, clamped);
            if (!seen_ids.insert(r.validated_seed_id).second)
                rj["duplicate"] = true;
            if (req.detail == "graphlet") {
                // after the clamps: the K records are the final limitations
                graphlet.seed_index = i;
                size_t lines = 0;
                std::string text = graphlet_text(r, seed, req.strategy, graphlet, rj, &lines);
                rj["graphlet_bytes"] = uint_json(text.size());
                rj["graphlet_lines"] = uint_json(lines);
                rj["graphlet"] = std::move(text);
            }
            const double seconds = serialize.elapsed();
            serialize_seconds += seconds;
            if (req.timing)
                rj["timing"]["serialize_ms"] = seconds * 1000;
            results.append(std::move(rj));
        } catch (const SeedDerivationError &e) {
            // The permitted set could not be derived from THIS seed (nothing carries it
            // in full, the budget ran out, the names are ambiguous). The caller named no
            // labels, so it had no lever on that and nothing to fix in the request:
            // report it against the seed and keep the other seeds' traversals, instead of
            // discarding a 100-seed batch because seed 57 spans a recombination point.
            results.append(failed_seed_to_json(seed, e, req.strategy, clamped));
            per_seed(oracle.counters());   // do not bill this seed's reads to the next
        } catch (const std::invalid_argument &e) {
            throw InvalidRequest(std::string("seed '") + (seed.seed_id.empty() ? seed.sequence.substr(0, 32) : seed.seed_id)
                                 + "': " + e.what());
        }
    }
    out["results"] = std::move(results);
    if (req.timing) {
        Json::Value t;
        t["elapsed_ms"] = timer.elapsed() * 1000;
        t["serialize_ms"] = serialize_seconds * 1000;
        out["timing"] = std::move(t);
    }
    return out;
}

Json::Value process_resolve_request(const Json::Value &json,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const std::string &release,
                                    uint64_t max_query_bp,
                                    const IndexIdentity *identity) {
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
    out["capabilities"] = capabilities_to_json(oracle, release, identity);
    Json::Value t;
    t["elapsed_ms"] = timer.elapsed() * 1000;
    out["timing"] = std::move(t);
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

    // the index identity every response states (§3.1), computed once for all requests
    IndexIdentity identity;
    try {
        identity = index_identity(*config, *anno_graph);
    } catch (const std::exception &e) {
        logger->error("{}", e.what());
        return 1;
    }

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
                ? process_resolve_request(json, *anno_graph, config->index_release, 0, &identity)
                : process_traverse_request(json, *anno_graph, config->index_release, {}, &identity);
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
