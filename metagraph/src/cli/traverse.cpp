#include "traverse.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>

#include "cli/config/config.hpp"
#include "cli/json_helpers.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse_attempts.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph::traversal;

namespace {

// The rules of the capabilities' coordinates and resolve blocks, stated by reference to the
// SPEC section that holds each (the documents a service returns in one piece have a ceiling of
// 32 KiB; the numbers a client computes with are fields beside them). ASCII only
constexpr char kCoordinatesRule[]
        = "SPEC-labeled-traversal-core.md section 7.1, record coordinates";
constexpr char kCoordinatesOutputBound[]
        = "SPEC-labeled-traversal-core.md section 7.1, record coordinates: what bounds the block";
constexpr char kResolveTimeBudgetRule[] = "SPEC-labeled-traversal-core.md section 4.5";

// ---------------------------------------------------------------- the attempt's delivery check

// The check of the server's attempt while the results are built (traverse_attempts.hpp): its
// client (gone: AttemptAborted) and its bound (reached: AttemptAtBound), every 4096 objects of
// the JSON tree and of the MGT text, so that an attempt stops at its bound wherever it is
// rather than outliving the lease that bound is (DESIGN-traverse-graphlet.md §14). Set by
// process_traverse_request for its own thread (the builders take no attempt); null otherwise,
// where a tick is one test of a thread-local pointer.
thread_local Attempt *t_delivery = nullptr;
thread_local uint32_t t_delivery_ticks = 0;

inline void delivery_tick() {
    if (t_delivery && !(++t_delivery_ticks & 4095))
        t_delivery->check_delivery();
}

class DeliveryScope {
  public:
    explicit DeliveryScope(Attempt *attempt) : previous_(t_delivery) { t_delivery = attempt; }
    ~DeliveryScope() { t_delivery = previous_; }
    DeliveryScope(const DeliveryScope&) = delete;
    DeliveryScope& operator=(const DeliveryScope&) = delete;

  private:
    Attempt *previous_;
};

// ---------------------------------------------------------------- strict JSON access

struct RefuseRequest {
    InvalidRequest operator()(const std::string &message) const { return InvalidRequest(message); }
};

// A request object read strictly (StrictObject): an unknown field is refused when it goes out
// of scope, unless an exception is already on its way
class Strict : public StrictObject<RefuseRequest> {
  public:
    using StrictObject::StrictObject;
    ~Strict() noexcept(false) {
        if (!std::uncaught_exceptions())
            finish();
    }
    std::string child_path(const std::string &k) const { return path(k); }
};

Json::Value labels_json(const std::vector<LabelId> &labels) {
    Json::Value arr(Json::arrayValue);
    for (LabelId l : labels) arr.append(l);
    return arr;
}

// A limit that may also be "unlimited" (Strategy::kUnlimited). |def| is the value for
// an omitted field; a mode in which the knob has no meaning passes |only_unlimited|,
// and a number is then rejected rather than accepted and ignored.
// |min| > 0: a knob 0 means nothing for. Every number out of [min, max] is then refused with
// that whole range, "unlimited" included, so that a client following the refusal is not
// refused again (a refusal naming [0, max] would misstate it, and 0 would need a refusal of
// its own). With |min| 0, every number out of [0, max] is refused.
size_t limit_or_unlimited(Strict &s, const std::string &k, size_t def, size_t max,
                          const char *only_unlimited, size_t min = 0) {
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
    // a non-negative integer, as Strict::uint takes it (anything else is its refusal)
    if (min && v.isIntegral() && !(v.isInt64() && v.asInt64() < 0)
            && (v.asUInt64() < min || v.asUInt64() > max)) {
        throw InvalidRequest(s.child_path(k) + ": out of range [" + std::to_string(min) + ", "
                             + std::to_string(max) + "] (or \"unlimited\")");
    }
    return s.uint(k, def, min, max);
}

Json::Value limit_json(size_t x) {
    return x == Strategy::kUnlimited ? Json::Value("unlimited") : uint_json(x);
}

// The arm counters by name, in the walker's order: one list for the JSON `counters` and
// the graphlet's A record, so that a counter added here reaches both (an MGT reader
// keeps a name it does not know). work_units (the work budget's meter) is listed when
// the request set a budget, so that a response without one is unchanged.
std::vector<std::pair<const char*, uint64_t>> arm_counters(const ArmResult &arm,
                                                           const Strategy &st) {
    std::vector<std::pair<const char*, uint64_t>> counters {
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
    if (st.max_memory_bytes || st.max_work_units)
        counters.emplace_back("work_units", arm.work_units);
    return counters;
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
    // Routing fields: consumed by the server before parsing (resolve_traverse_index). They
    // are still checked and declared here, or the unknown-field check would reject every
    // multi-graph request. `graphs` ([name], /search's form of `graph`) and `in_ram` (load the
    // index into RAM for the request, as /search does) likewise
    s.str("graph", "");
    s.str("graph_path", "");
    s.strings("graphs");
    s.boolean("in_ram", false);
    // The attempt's fields (the frozen wire contract, DESIGN-traverse-graphlet.md §14.1): the
    // server reads them before the request is parsed, to register the attempt (attempt_ids,
    // the same rule); declared here so that strict parsing accepts them
    const AttemptIds ids = attempt_ids(json);
    for (const char *field : kAttemptFields) {
        s.has(field);
    }
    req.attempt_id = ids.attempt_id;
    req.budget_id = ids.budget_id;
    req.locus_id = ids.locus_id;
    req.not_after_ms = ids.not_after_ms;

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
            // the request's budgets (DESIGN-traverse-graphlet.md §14); omitted = none
            if (b.has("max_memory_mb"))
                st.max_memory_bytes = b.uint("max_memory_mb", 0, 1, 1'048'576) << 20;
            if (b.has("max_work_units"))
                st.max_work_units = b.uint("max_work_units", 0, 1);
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
            // Record coordinates (DESIGN-traverse-graphlet.md §18; opt-in). The cap is refused
            // without coordinates: true, where it would change nothing (§7.0's rule for knobs),
            // and accepted with it whatever the support — under support kmer it is inert, and
            // the result says why no coordinates are reported, so that a request switching its
            // support stays valid
            st.coordinates = o.boolean("coordinates", false);
            if (o.has("max_coordinate_occurrences")) {
                if (!st.coordinates) {
                    throw InvalidRequest("strategy.output.max_coordinate_occurrences: only with "
                                         "output.coordinates true, where it bounds the occurrence "
                                         "lists; without it, it would change nothing: omit it "
                                         "or set output.coordinates to true");
                }
                // 0 would record no occurrence while asking for them: refused with the range
                st.max_coordinate_occurrences = limit_or_unlimited(
                        o, "max_coordinate_occurrences", 16,
                        std::numeric_limits<size_t>::max() - 1, nullptr, 1);
            }
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
    // see parse_traverse_request: routing fields and in_ram are handled by the server
    s.str("graph", "");
    s.str("graph_path", "");
    s.strings("graphs");
    s.boolean("in_ram", false);
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
        // The deadline (the search service runs /resolve as a job): opt-in, so that a
        // request without it is answered without a time limit. Any finite non-negative
        // number here; the route refuses a budget not above its finalisation reserve and
        // lowers one above the server's cap, which it knows
        if (b.has("time_budget_ms")) {
            req.time_budget_ms = b.number("time_budget_ms", 0, 0,
                                          std::numeric_limits<double>::max());
            req.time_budget_given = b.raw("time_budget_ms");
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
    // echoed only when given: an omitted budget is no budget, and the echo of a request
    // without one stays what it always was (still resubmittable either way)
    if (st.max_memory_bytes)
        bo["max_memory_mb"] = uint_json(st.max_memory_bytes >> 20);
    if (st.max_work_units)
        bo["max_work_units"] = uint_json(st.max_work_units);
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
    // echoed only when asked for, with the cap, default included, so that the echo is the
    // whole request; a request without them echoes neither
    if (st.coordinates) {
        o["coordinates"] = true;
        o["max_coordinate_occurrences"] = limit_json(st.max_coordinate_occurrences);
    }
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

// What a walk_domain of a budget at the server's maximum (SPEC §10.3; state_budget_clamp)
// says in place of "raise the knob"
static const char *const kAtServerMaximum
    = "; the knob is at the server's maximum (server_limit), which a request cannot raise";

// A stop requested from outside the walk (AttemptControl): the attempt was cancelled or
// reached the duration bound the server enforces for it. No budget of the request ran out,
// so its statements name the attempt, never a budget knob to raise
static bool is_external_stop(const ResourceStop &q) {
    return q.resource == ResourceStop::CANCELLED || q.resource == ResourceStop::ATTEMPT_DEADLINE;
}

// what an external stop's statements say happened, in their words
static const char* external_cause(const ResourceStop &q) {
    return q.resource == ResourceStop::CANCELLED
        ? "the attempt was cancelled (POST /traverse/cancel)"
        // the bound as enforced: capped at attempts.hard_cap_ms, the content timeout less one
        // second. An effect is priced at 640 bytes (DeliveryCosts), so the cap is named by its
        // capabilities field rather than spelled out
        : "the attempt reached the time at which the server stops walking it (the duration "
          "bound it enforces, the seeds' time budgets plus its allowance (at most hard_cap_ms), "
          "less the larger of half the allowance and the delivery reserve, kept for the "
          "delivery; see usage.bound.walk_until_ms)";
}

// The request field a resource stop answers to (§6.7) and its value; a beam's width is
// bounds.max_live_paths; a head a budget did not admit answers to the budget that refused
// it (|stop|; a memory refusal without a budget is a test hook's, stated as unlimited).
static std::pair<std::string, Json::Value> cap_knob(EndReason reason, const Strategy &st,
                                                    const ResourceStop *stop) {
    switch (reason) {
        case EndReason::RESOURCE_LIMIT:
            if (stop && stop->resource == ResourceStop::WORK)
                return { "bounds.max_work_units", uint_json(st.max_work_units) };
            // a stop from outside the walk answers to no budget: the field that names the
            // attempt, with the attempt's bound (ms) as its limit
            if (stop && is_external_stop(*stop))
                return { "attempt_id", uint_json(static_cast<uint64_t>(stop->limit)) };
            return { "bounds.max_memory_mb", st.max_memory_bytes
                                                 ? uint_json(st.max_memory_bytes >> 20)
                                                 : Json::Value("unlimited") };
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

// A stop by a budget-aware annotation read that did not fit
// (DESIGN-traverse-graphlet.md §14.1): its statements name the decoding, its levers the seed
static bool is_decode_stop(const ResourceStop &q) {
    return std::string(q.phase) == "annotation_decode";
}

// A seed-level memory stop's walk_domain observed, in the knob's unit: the smallest budget
// (whole MiB) that holds its demand beside the caches' allotments of THAT budget — the demand
// includes this budget's allotments, which grow with the knob, so the demand itself, raised
// to, would fail again
static uint64_t knob_mib(const ResourceStop &q) {
    return memory_budget_holding(static_cast<uint64_t>(std::ceil(q.demand)), q.allotted) >> 20;
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

static Json::Value arm_limitations(const ArmResult &arm, const Strategy &st,
                                   const ResourceStop *stop) {
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
            const bool decode = r == EndReason::RESOURCE_LIMIT && stop && is_decode_stop(*stop);
            const ResourceStop::Cause cause = r == EndReason::RESOURCE_LIMIT && stop
                ? stop->cause : ResourceStop::HEAD;
            const bool external = r == EndReason::RESOURCE_LIMIT && stop && is_external_stop(*stop);
            effect += r == EndReason::BEAM_PRUNED
                ? "the beam kept the best-supported max_live_paths heads of a level and pruned the rest"
                : external
                ? std::string(external_cause(*stop)) + " and the exploration stopped at the next "
                  "checkpoint (resource_limit; see resource_stop; "
                  + (trigger ? "limit: the attempt's bound, observed: its elapsed time, ms"
                             : "after the cap that set complete_to_bp; limit: the attempt's "
                               "bound, ms, observed: the walks it ended")
                  + "); no budget of the request ran out: a new attempt can continue from the "
                    "leaves"
                : decode && stop->injected
                ? "an injected refusal of an annotation read (a test hook, not the budget) stopped "
                  "the exploration there (resource_limit; see resource_stop)"
                : decode
                ? "the request's memory budget did not admit decoding the next level's annotation "
                  "rows with their row-diff dependency rows, and the exploration stopped there "
                  "(resource_limit; see resource_stop)"
                : cause == ResourceStop::LABEL_NAMES
                ? "the request's memory budget did not admit the dictionary labels the next "
                  "level's annotation rows name first (with their delivery in the requested "
                  "detail), and the exploration stopped there (resource_limit; see resource_stop)"
                : cause == ResourceStop::LEVEL_LISTS
                ? "the next level's key and successor lists, held beside the request's account, "
                  "left none of the memory budget for its annotation rows, and the exploration "
                  "stopped there (resource_limit; see resource_stop)"
                : r == EndReason::RESOURCE_LIMIT && stop && stop->injected
                ? "an injected allocation refusal (a test hook, not the budget) did not admit the "
                  "next head, and the exploration stopped there (resource_limit; see resource_stop)"
                : r == EndReason::RESOURCE_LIMIT
                ? "the request's budget did not admit the next head, and the exploration stopped "
                  "there (resource_limit; see resource_stop)"
                : std::string("the exploration stopped at this cap (") + to_string(r) + ")";
            if (!trigger && !external)
                effect += " after the cap that set complete_to_bp, so raising only that knob stops here";
            // A refused level read and a level whose lists left no room state the least the
            // stop needed (the exact demand would depend on annotation.batch_kmers, §6.8), so
            // the value is a lower bound, said as such, and raising the knob to it is no
            // promise
            if (trigger && r == EndReason::RESOURCE_LIMIT && stop && !stop->injected
                    && stop->lower_bound) {
                effect += " (observed: at least what admitting it needed, MiB rounded up to the "
                          "smallest budget that admits it with the caches' allotments, which grow "
                          "with the budget; raising the knob to it may still fail: unread rows or "
                          "later state may need more)";
            } else if (trigger && r == EndReason::RESOURCE_LIMIT && stop && !stop->injected
                    && stop->resource == ResourceStop::MEMORY && stop->allotted) {
                // the need at this budget includes the caches' allotments of this budget, which
                // grow with it: the budget that admits the head is stated
                // (memory_budget_holding)
                effect += " (observed: what admitting it needed, MiB rounded up to the smallest "
                          "budget that admits it with the caches' allotments, which grow with the "
                          "budget)";
            }
            effect += external ? ""
                : r == EndReason::RESOURCE_LIMIT && stop && stop->injected
                ? "; no request knob caused it" : "; raise the knob";
            auto [knob, limit] = cap_knob(r, st, stop);
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
// (DESIGN-traverse-graphlet.md §14) — an axis is complete only when no limitation of
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
//                       trace_record_boundaries, and a walked result's derivation (its
//                       permitted set derived from part of the seed); qualified wins when
//                       both apply. The coordinates limitation (a cut occurrence list) is in
//                       no class: its block states complete: false
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
            // a walked result's derivation limitation: its permitted set was derived from part
            // of the seed, a superset of the whole seed's carriers. A failed seed's states
            // why there is no result at all and qualifies nothing
            qualified |= kind == "trace_record_boundaries" || (kind == "derivation" && !failed);
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
static void arm_certificate_json(Json::Value *j, const ArmResult &arm, const Strategy &st,
                                 const ResourceStop *stop) {
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
    (*j)["limitations"] = arm_limitations(arm, st, stop);
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

static Json::Value counters_json(const ArmResult &arm, const Strategy &st) {
    Json::Value counters(Json::objectValue);
    for (const auto &[name, value] : arm_counters(arm, st)) {
        counters[name] = uint_json(value);
    }
    return counters;
}

static Json::Value arm_to_json(const ArmResult &arm, const Strategy &st, const std::string &detail,
                               const ResourceStop *stop) {
    Json::Value j;
    arm_certificate_json(&j, arm, st, stop);
    if (detail != "summary") {
        Json::Value segs(Json::arrayValue);
        for (const auto &s : arm.segments) {
            delivery_tick();
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
    // a path's chain is produced here from parent pointers (the walker keeps only the
    // leaf), through one buffer reused across paths
    std::vector<size_t> chain;
    for (const auto &p : arm.paths) {
        delivery_tick();
        Json::Value pj;
        pj["id"] = uint_json(p.id);
        if (detail != "summary") {
            chain.clear();
            walk_path_leaf_first(arm, p, [&](size_t s) { chain.push_back(s); return true; });
            Json::Value segs(Json::arrayValue);
            for (auto it = chain.rbegin(); it != chain.rend(); ++it) segs.append(uint_json(*it));
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
                // non-zero: merged in from another route at a reconvergence (the latest
                // such merge), so this label does NOT support the spelled bases BEFORE
                // route_bp; its displayed support at the leaf starts there
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
            delivery_tick();
            Json::Value rj;
            rj["label"] = r.label;
            rj["from_bp"] = uint_json(r.from_bp);
            rj["to_bp"] = uint_json(r.to_bp);
            rj["entered_by"] = r.entered_by_switch ? "switch" : "seed";
            // the run's EARLIEST non-first-parent merge stamp (LabelRun::route_bp), not where
            // its displayed support begins: that is derived (evidence_from, spec §7.1)
            rj["route_bp"] = uint_json(r.route_bp);
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
    j["counters"] = counters_json(arm, st);
    return j;
}

// The arm of the `detail: graphlet` summary (DESIGN-traverse-graphlet.md §3): the
// certificate, the counters and the counts that validate the body, nothing per segment,
// leaf or label (those are in the MGT text).
static Json::Value arm_summary_json(const ArmResult &arm, const Strategy &st,
                                    const ResourceStop *stop) {
    Json::Value j;
    arm_certificate_json(&j, arm, st, stop);
    j["counters"] = counters_json(arm, st);
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

// A resource stop is stated (`resource_stop`) when a budget of §14 ended the walk, and for
// a time stop only when the request set such a budget: a time stop has its walk_domain
// limitation either way, and a response without a budget keeps its form.
static const ResourceStop* stated_stop(const SeedResult &r, const Strategy &st) {
    if (!r.resource_stop)
        return nullptr;
    if (r.resource_stop->resource == ResourceStop::TIME && !st.max_memory_bytes
            && !st.max_work_units)
        return nullptr;
    return &*r.resource_stop;
}

// `resource_stop` (DESIGN-traverse-graphlet.md §14, the MGT Q record): the amounts in the
// unit of the knob that controls them — whole MiB for memory (used rounded up, remaining
// down), work units, milliseconds — the next actions, and what the stop means. A property
// of the walk, so the same in every detail: the graphlet's Q record must rebuild the
// detail: full rendering of the same walk. |failed|: the budget did not hold the seed
// itself (no walk, no leaves to continue from), whose actions are the levers on the seed.
static Json::Value resource_stop_json(const ResourceStop &q, const Strategy &st,
                                      const SeedBudgetError *failed = nullptr,
                                      const ResourceAccount *account = nullptr) {
    constexpr double kMiB = 1 << 20;
    Json::Value j;
    // a stop from outside the walk stops the attempt, not a budget of the locus
    j["scope"] = is_external_stop(q) ? "attempt" : "locus";
    Json::Value actions(Json::arrayValue);
    std::string knob, unit;
    switch (q.resource) {
        case ResourceStop::MEMORY: {
            j["resource"] = "memory";
            knob = "bounds.max_memory_mb";
            const Json::Value limit = q.limit > 0
                ? uint_json(static_cast<uint64_t>(q.limit) >> 20) : Json::Value("unlimited");
            j["requested"] = limit;
            j["effective"] = limit;
            j["used"] = uint_json(static_cast<uint64_t>(std::ceil(q.used / kMiB)));
            j["remaining"] = q.limit > 0
                ? uint_json(q.limit > q.used ? static_cast<uint64_t>((q.limit - q.used) / kMiB) : 0)
                : Json::Value("unlimited");
            // the output is charged in the requested detail: a graphlet (or no bases)
            // delivers more walk within the same budget. No lever helps against a refusal
            // a test hook injected (C++ callers only): only continuing from the leaves
            const bool annotate = st.label_mode == LabelMode::ANNOTATE;
            // A seed-phase read (the derivation's window, the validation) holds annotation
            // rows, not output: the detail does not change what it needs. An annotate root's
            // read competes with the account, which holds the other arm's root and its labels
            // once that was read: the detail and the label cap shrink those (they turn such a
            // failure into a walk)
            const bool seed_read = failed && q.cause == ResourceStop::READ_ROW
                && q.where != ResourceStop::ROOT;
            const bool root_read = failed && q.cause == ResourceStop::READ_ROW
                && q.where == ResourceStop::ROOT;
            // the other root's state competes with a root's read when the row (its demand, or
            // the least its read alone was seen to need) would fit without it
            const bool competes = root_read && q.label_bytes
                && q.demand - q.used <= static_cast<double>(q.left + q.label_bytes);
            if (!q.injected) {
                actions.append("raise_memory_budget");
                if (!seed_read && !(root_read && !competes)) {
                    actions.append("use_graphlet");
                    if (st.sequences)
                        actions.append("drop_sequences");
                    // Wherever coordinates were asked for, the output's coordinate part is in the
                    // account: the recorded block (a run's coordinates can cost 4-20 times the
                    // run itself) or the null form with its reason, which no request can avoid
                    // but by dropping them (an index without coordinates, support kmer: about
                    // 1.6 KB a seed, enough to move a stop). drop_coordinates is offered
                    // for both
                    if (st.coordinates)
                        actions.append("drop_coordinates");
                }
            }
            if (q.injected) {
                // no lever helps against a test hook
            } else if (q.cause == ResourceStop::READ_ROW) {
                // A refused annotation read: on a row-diff annotation every read
                // decodes a whole row with its dependency rows, so what reads fewer rows is
                // the lever — a more selective seed, or in annotate mode (which reads every
                // node's row) a label-constrained query; naming fewer labels or a smaller
                // radius does not shrink a row
                if (competes && q.labels)
                    actions.append("lower_max_labels_per_node");
                actions.append("more_selective_seed");
                if (annotate)
                    actions.append("label_constrained_query");
            } else if (q.cause == ResourceStop::LABEL_NAMES) {
                // the labels a level's rows name first (annotate mode): fewer recorded per node,
                // or a label-constrained query, which names none
                actions.append("lower_max_labels_per_node");
                actions.append("label_constrained_query");
            } else if (failed) {
                // a depth-0 state costs per label at the roots: fewer labels fit. The lever
                // is the knob that sets them — the derived set's cap, the named list, or
                // (annotate mode, which permits no set) the cap on a node's recorded list
                actions.append(st.label_mode == LabelMode::ANNOTATE ? "lower_max_labels_per_node"
                               : failed->labels_from_seed() ? "lower_max_seed_labels"
                                                            : "name_fewer_labels");
            }
            break;
        }
        case ResourceStop::WORK:
            j["resource"] = "work";
            knob = "bounds.max_work_units";
            j["requested"] = uint_json(static_cast<uint64_t>(q.limit));
            j["effective"] = uint_json(static_cast<uint64_t>(q.limit));
            j["used"] = uint_json(static_cast<uint64_t>(q.used));
            j["remaining"] = uint_json(q.limit > q.used ? static_cast<uint64_t>(q.limit - q.used) : 0);
            actions.append("raise_work_budget");
            // the seed phase reads one annotation row per seed k-mer
            if (failed)
                actions.append("shorten_seed");
            break;
        case ResourceStop::TIME:
            j["resource"] = "time";
            knob = "bounds.time_budget_ms";
            j["requested"] = q.limit;
            j["effective"] = q.limit;
            j["used"] = q.used;
            j["remaining"] = q.limit > q.used ? q.limit - q.used : 0.0;
            actions.append("raise_time_budget");
            break;
        case ResourceStop::CANCELLED:
        case ResourceStop::ATTEMPT_DEADLINE:
            // the attempt's bound and its elapsed time, whole milliseconds; no budget to raise:
            // what was walked is delivered, and a new attempt does the rest
            j["resource"] = q.resource == ResourceStop::CANCELLED ? "cancelled" : "attempt_deadline";
            knob = "attempt_id";
            j["requested"] = uint_json(static_cast<uint64_t>(q.limit));
            j["effective"] = uint_json(static_cast<uint64_t>(q.limit));
            j["used"] = uint_json(static_cast<uint64_t>(q.used));
            j["remaining"] = uint_json(q.limit > q.used ? static_cast<uint64_t>(q.limit - q.used) : 0);
            if (failed)
                actions.append("retry_attempt");
            break;
    }
    if (!failed)
        actions.append("continue_from_leaves");
    j["phase"] = q.phase;
    j["actions"] = std::move(actions);
    if (failed && is_external_stop(q)) {
        j["message"] = std::string(failed->what()) + " (the seed is failed: its seed phase was "
                       "stopped before any traversal, so no result exists, not even one complete "
                       "to 0 bp; no budget of the request ran out: a new attempt can retry it)";
        return j;
    }
    if (failed) {
        std::string message = std::string(failed->what()) + " (the seed is failed: no valid "
                              "traversal exists within the budget, not even one complete to 0 bp)";
        if (q.resource == ResourceStop::WORK) {
            // how far the failed seed phase ran past the budget, with the number that bounds
            // it, as a walk's work stop states it
            message += "; the seed phase is compared with the budget once every "
                     + std::to_string(kWorkCheckInterval) + " units and fails at a comparison "
                       "finding it at least " + std::to_string(kWorkCheckInterval)
                     + " units over budget (the most this seed charged between two comparisons: "
                     + std::to_string(account ? account->largest_charge
                                              : failed->account().largest_charge)
                     + " units)";
        }
        j["message"] = message;
        return j;
    }
    if (q.resource == ResourceStop::MEMORY && q.cause != ResourceStop::HEAD) {
        // A level's budget-aware annotation read did not fit: the level is censored from its
        // first head, as at a refused head. The message names what did not fit — the row, the
        // labels it would name first, or the level's own lists — with its bytes
        const std::string where = std::to_string(q.at_bp) + " bp on the " + to_string(q.arm)
                                + " arm";
        const std::string present = "every walk up to each arm's complete_to_bp is present and "
                                    "the heads not expanded end with resource_limit";
        std::string message;
        if (q.injected) {
            message = "an injected refusal of an annotation read (a test hook, not " + knob
                    + ": the budget, if any, would have admitted the read) stopped the walk at "
                    + where + " while decoding annotation: " + present;
        } else if (q.cause == ResourceStop::LEVEL_LISTS) {
            message = "the memory budget (" + knob + ") stopped the walk at " + where
                    + " before it read the next level's annotation: the level's key and "
                      "successor lists, held until its heads are processed beside the walk's "
                      "account (" + std::to_string(q.held) + " bytes together), left none of the "
                    + std::to_string(static_cast<uint64_t>(q.limit)) + " bytes for its rows; "
                    + present + "; the lists grow with the level's width, the account with the "
                      "output's detail";
        } else if (q.cause == ResourceStop::LABEL_NAMES && q.names_after_read) {
            // a format whose reads are not budget-aware: no row was refused, so there is no
            // row demand or bytes left to state, and the need is a lower bound
            message = "the memory budget (" + knob + ") stopped the walk at " + where
                    + " after reading the next level's annotation: its rows named "
                    + std::to_string(q.labels) + " new dictionary label(s), whose entries and "
                      "delivery in the requested detail ("
                    + std::to_string(q.label_bytes) + " bytes) put the walk's account over the "
                      "budget (" + std::to_string(q.held) + " bytes held with the level's "
                      "lists); " + present + "; on this annotation's format a level's rows are "
                      "read whole and the labels they name are charged after the read, so the "
                      "level needs at least the account with them, its heads more; a lower "
                      "labels.max_labels_per_node, detail graphlet or a label-constrained query "
                      "names fewer or cheaper labels";
        } else if (q.cause == ResourceStop::LABEL_NAMES) {
            message = "the memory budget (" + knob + ") stopped the walk at " + where
                    + " while reading the next level's annotation: a row that fits ("
                    + std::to_string(q.row_demand) + " bytes, read alone) names "
                    + std::to_string(q.labels) + " new dictionary label(s) whose entries, "
                      "delivery in the requested detail and provisional naming come to "
                    + std::to_string(q.label_bytes) + " bytes, more than the "
                    + std::to_string(q.left) + " bytes it had left beside the walk's account, "
                      "the level's lists and the rows read before it (" + std::to_string(q.held)
                    + " bytes); " + present + "; labels are named and charged with the rows that "
                      "first carry them, a refused read names none; a lower "
                      "labels.max_labels_per_node, detail graphlet or a label-constrained query "
                      "names fewer or cheaper labels";
        } else {
            message = "the memory budget (" + knob + ") stopped the walk at " + where
                    + " while decoding annotation: a row of the next level read alone with its "
                      "row-diff dependency rows needs more than the "
                    + std::to_string(q.left) + " bytes the walk had left beside "
                      "its account, the level's lists and the rows read before it ("
                    + std::to_string(q.held) + " bytes); " + present + "; each dependency row and "
                      "coordinate tuple is charged before it is held, every row a fetch returns is "
                      "admitted against its standalone demand, and a read that does not fit is "
                      "refused whole, so the stop does not depend on annotation.batch_kmers; "
                    + (st.label_mode == LabelMode::ANNOTATE
                        ? "a more selective seed or a label-constrained query reads fewer rows"
                        : "a more selective seed reads fewer rows")
                    + ", a smaller radius does not";
        }
        j["message"] = message + "; a continuation from a leaf is a new traversal";
        return j;
    }
    if (is_external_stop(q)) {
        // when the walk stopped, on the attempt's clock: what a ledger needs to release the
        // attempt's capacity (DESIGN §14: a cancellation must be acknowledged)
        j["message"] = std::string(external_cause(q)) + " and the walk stopped at the next "
            "checkpoint, " + std::to_string(static_cast<uint64_t>(q.used)) + " ms after the "
            "request was received (the bound: " + std::to_string(static_cast<uint64_t>(q.limit))
            + " ms), at " + std::to_string(q.at_bp) + " bp on the " + to_string(q.arm)
            + " arm: every walk up to each arm's complete_to_bp is present and the heads not "
              "expanded end with resource_limit; no budget of the request ran out; a "
              "continuation from a leaf is a new traversal";
        return j;
    }
    // A refusal injected by a test hook (WalkerHooks::deny, C++ callers only) is reported as
    // a memory stop so that it stays representable in Q and K, but no budget refused the
    // head (without one the limit reads "unlimited"), so the message must not say a budget did
    const std::string cause = q.injected
        ? "an injected allocation refusal (a test hook, not " + knob + ": the budget, if any, "
          "would have admitted the head)"
        : "the " + j["resource"].asString() + " budget (" + knob + ")";
    std::string message = cause + " stopped the walk at " + std::to_string(q.at_bp) + " bp on the "
        + to_string(q.arm) + " arm: every walk up to each arm's complete_to_bp is present and the "
          "heads not expanded end with "
        + (q.resource == ResourceStop::TIME ? "time_budget" : "resource_limit");
    if (q.resource == ResourceStop::MEMORY && !q.injected && account && account->decode_charged) {
        message += "; memory is modelled (the walker's state and the output in the "
                   "requested detail, charged per head; detail graphlet costs the least), "
                   "annotation reads are charged inside the decoder (what a level's fetch holds "
                   "until its heads are processed is not: memory_bound_soft)";
    } else if (q.resource == ResourceStop::MEMORY && !q.injected) {
        message += "; memory is modelled (the walker's state and the output in the "
                   "requested detail, charged per head; detail graphlet costs the least), "
                   "annotation decoding is not charged yet (memory_bound_soft)";
    } else if (q.resource == ResourceStop::WORK) {
        // How far used can exceed the budget, stated with the number that bounds it rather
        // than as a fixed maximum or a kind of charge that another kind could break (two
        // roots' rows, a row's coordinates and a label-state scan can each exceed "the
        // widest row"): every comparison follows one that passed, so the overrun is at
        // most what was charged since, and the walker records the most it charged between
        // two comparisons. Work is the walk's, the same in every detail (charging delivery
        // would make the stop depend on the detail), so it also says what bounds the
        // output.
        message += "; work is compared with the budget after every charge, so used exceeds it "
                   "by at most what was charged since the previous comparison";
        if (account && account->largest_charge) {
            message += " (the most this seed charged between two comparisons: "
                     + std::to_string(account->largest_charge) + " units)";
        }
        // with the budget-aware reads a row is charged with its row-diff dependency rows; a
        // row-diff annotation without them says that they are not counted
        const char *weights = account && account->decode_charged
            ? "(8 per key and per dependency row, 1 per entry and coordinate; near the budget a "
              "call reads one key)"
            : account && account->row_diff_uncounted
            ? "(8 per key, 1 per entry and coordinate, none per dependency row; near the budget a "
              "call reads one key)"
            : "(8 per key, 1 per entry and coordinate; near the budget a call reads one key)";
        // W is the threshold at which a seed phase is failed, not a ceiling on how far it ran
        // past the budget
        message += std::string(", one indivisible charge: a fetch call's rows, decoded whole with "
                   "their coordinates ") + weights + ", a label-state scan, or the roots' rows "
                   "with the end of the seed phase (one charge, so that the result complete to 0 "
                   "bp is delivered, as wide as the index makes them; the seed phase fails at a "
                   "comparison finding it at least " + std::to_string(kWorkCheckInterval)
                 + " units over budget); work bounds the walk, not its delivery: the output's "
                   "size (and the time to write it) is bounded by bounds.max_memory_mb";
    }
    message += "; a continuation from a leaf is a new traversal";
    j["message"] = message;
    return j;
}

// memory_bound_soft (§7.0): stated by every response under a memory budget — a failed seed's
// too — because something is always held beyond the admitted account. With the budget-aware
// reads (account.decode_charged) the annotation reads are charged before they are held, and
// the statement names only what is still uncharged; otherwise the decoded rows are uncharged
// too, and the statement says so. Both name a failed result's echo of the request's seed_id,
// which only the request bounds (failed_soft prices it).
static Json::Value memory_bound_soft(const Strategy &st, const ResourceAccount &account) {
    return limitation("memory_bound_soft", "bounds.max_memory_mb",
                      uint_json(st.max_memory_bytes >> 20),
                      uint_json((account.soft_overshoot + (1 << 20) - 1) >> 20),
                      account.decode_charged
                      ? "the memory budget is enforced on the walker's modelled state, this "
                        "response's output and the annotation reads (dependency rows and "
                        "coordinate tuples are charged before they are held; a read that does not "
                        "fit is refused whole); held beyond the admitted account, so able to "
                        "exceed the budget: a level's keys and fetched rows until its heads are "
                        "processed, the seed phase's intersection and hits, a failed seed's "
                        "depth-0 dictionary and a failed result's echo of seed_id (observed: the "
                        "largest excess seen, MiB, rounded up), and, not observed, an index-wide "
                        "header lookup and the label dictionary's first table and growth copy "
                        "(about 1.5 KB)"
                      : "the memory budget is enforced on the walker's modelled state and on "
                        "this response's output, admitted per head; the annotation rows the "
                        "seed phase and each level decode (with an annotate dictionary's "
                        "growth and a cache beyond its allotment), and a failed result's echo of "
                        "seed_id, are held before they can be charged, so the peak can exceed the "
                        "budget by them (observed: the excess seen, MiB, rounded up)");
}

// One occurrence, the half-open base interval [start, end) in record (header labels) or column
// (column labels) coordinates
static Json::Value occurrence_json(uint64_t start, uint64_t end) {
    Json::Value j(Json::arrayValue);
    j.append(uint_json(start));
    j.append(uint_json(end));
    return j;
}

/**
 * The `coordinates` block of a seed whose coordinates were recorded (DESIGN-traverse-graphlet.md
 * §18.2), the same in every detail: the kind (record: every label a header, positions within its
 * record; column: every label a column, positions global in its column; mixed: each label's kind
 * says which), k, the cap, whether every list is whole; per seed label its occurrences of the
 * seed; per requested arm one entry per run, in R order, with the occurrences of the run's own
 * bases [from_bp, to_bp) — from a chain's k-mer coordinate c at the run's last node and the
 * run's length L: [c + k - L, c + k) on the right arm, [c, c + L) on the left — the true count
 * where the list was cut, the chains that ended before the last node, and lower_bound where the
 * run's chains are a lower bound (a switch into a label whose own lineage was live). A cut list
 * is stated by the `coordinates` limitation appended to |lims| (no outcome class); lower-bound
 * runs by the block alone. A column's coordinates number its k-mers (record i's k-mer j is
 * offset_i + j), so a column interval numbers base p of record i as offset_i + p: record i's
 * last k - 1 bases share their numbers with record i + 1's first k - 1, and an interval touching
 * them is attributed to one record only with the record lengths, which this block does not use.
 * The contract and the docs state it, beside trace_record_boundaries (a column's trace can run
 * across two records whose coordinates are adjacent).
 * delivery_tick() per entry and occurrence: a large block is built under the attempt's check.
 */
static Json::Value coordinates_json(const SeedResult &r, const Strategy &st, Json::Value *lims) {
    const uint64_t k = r.k;
    size_t headers = 0, columns = 0;
    for (const LabelRef &l : r.label_dict) {
        (l.kind == LabelKind::HEADER ? headers : columns)++;
    }
    Json::Value block;
    block["kind"] = !columns ? "record" : !headers ? "column" : "mixed";
    block["k"] = uint_json(k);
    block["max_occurrences"] = limit_json(st.max_coordinate_occurrences);
    size_t lists_cut = 0, runs_lower_bound = 0;
    uint64_t largest_cut = 0;
    // the true count beside a cut list, never a silent cut
    auto state_total = [&](Json::Value *entry, uint64_t total, size_t listed) {
        if (total > listed) {
            (*entry)["occurrences_total"] = uint_json(total);
            lists_cut++;
            largest_cut = std::max<uint64_t>(largest_cut, total);
        }
    };
    Json::Value seed(Json::arrayValue);
    for (size_t l = 0; l < r.seed_coordinates.size(); ++l) {
        const SeedCoordinates &sc = r.seed_coordinates[l];
        Json::Value e;
        e["label"] = static_cast<Json::UInt>(l);
        Json::Value occurrences(Json::arrayValue);
        for (Coord c : sc.starts) {
            occurrences.append(occurrence_json(c, c + r.length_bp));
            delivery_tick();
        }
        e["occurrences"] = std::move(occurrences);
        state_total(&e, sc.total, sc.starts.size());
        seed.append(std::move(e));
        delivery_tick();
    }
    block["seed"] = std::move(seed);
    Json::Value arms(Json::objectValue);
    for (Arm side : { Arm::LEFT, Arm::RIGHT }) {
        const ArmResult &a = r.arms[static_cast<size_t>(side)];
        if (!a.requested)
            continue;
        assert(a.run_coordinates.size() == a.runs.size());
        Json::Value list(Json::arrayValue);
        for (size_t i = 0; i < a.runs.size() && i < a.run_coordinates.size(); ++i) {
            const LabelRun &run = a.runs[i];
            const RunCoordinates &rc = a.run_coordinates[i];
            const uint64_t length = run.to_bp - run.from_bp;
            Json::Value e;
            e["run"] = static_cast<Json::UInt64>(i);
            e["label"] = run.label;
            e["from_bp"] = uint_json(run.from_bp);
            e["to_bp"] = uint_json(run.to_bp);
            Json::Value occurrences(Json::arrayValue);
            for (Coord c : rc.ends) {
                // a chain carries the run's whole length: on the right arm its k-mer at the last
                // node ends at c + k, the run's own bases end there
                assert(side == Arm::LEFT || c + k >= length);
                occurrences.append(side == Arm::RIGHT ? occurrence_json(c + k - length, c + k)
                                                      : occurrence_json(c, c + length));
                delivery_tick();
            }
            e["occurrences"] = std::move(occurrences);
            state_total(&e, rc.total, rc.ends.size());
            if (rc.chains_ended)
                e["chains_ended"] = uint_json(rc.chains_ended);
            if (rc.lower_bound) {
                e["lower_bound"] = true;
                runs_lower_bound++;
            }
            list.append(std::move(e));
            delivery_tick();
        }
        arms[to_string(side)] = std::move(list);
    }
    block["arms"] = std::move(arms);
    // false iff some list was cut or some run's occurrences are a lower bound
    block["complete"] = !lists_cut && !runs_lower_bound;
    if (runs_lower_bound)
        block["runs_lower_bound"] = uint_json(runs_lower_bound);
    if (lists_cut) {
        Json::Value l = limitation("coordinates", "output.max_coordinate_occurrences",
                                   limit_json(st.max_coordinate_occurrences),
                                   uint_json(largest_cut),
                                   std::to_string(lists_cut) + " occurrence list(s) (seed labels "
                                   "and runs) are cut at the cap to their first occurrences by "
                                   "start; each states its true count (occurrences_total; "
                                   "observed: the largest): raise the knob or set it to "
                                   "\"unlimited\"");
        l["lists_cut"] = uint_json(lists_cut);
        lims->append(std::move(l));
    }
    return block;
}

// compact_json_size with |check| (the attempt's delivery check) called every 4096 values, so
// that counting a large block is under the attempt's bound like writing it; |values| counts them
static uint64_t compact_size(const Json::Value &v, const std::function<void()> &check,
                             uint64_t *values) {
    if (check && !(++*values & 4095))
        check();
    switch (v.type()) {
        case Json::nullValue:
            return 4;
        case Json::booleanValue:
            return v.asBool() ? 4 : 5;
        case Json::uintValue:
            return decimal_digits(v.asUInt64());
        case Json::stringValue: {
            const char *begin = nullptr, *end = nullptr;
            v.getString(&begin, &end);
            return 2 + json_escaped_size(std::string_view(begin, end - begin));
        }
        case Json::intValue:
        case Json::realValue:
            // the writer's own text (a sign and digits; 17 significant digits or its special
            // values): rare in what this counts
            return json_text(v, true).size();
        case Json::arrayValue: {
            uint64_t n = 2 + (v.size() ? v.size() - 1 : 0);
            for (const Json::Value &e : v) {
                n += compact_size(e, check, values);
            }
            return n;
        }
        case Json::objectValue: {
            uint64_t n = 2 + (v.size() ? v.size() - 1 : 0);
            for (auto it = v.begin(); it != v.end(); ++it) {
                const char *end = nullptr;
                const char *key = it.memberName(&end);
                n += 2 + json_escaped_size(std::string_view(key, end - key)) + 1
                   + compact_size(*it, check, values);
            }
            return n;
        }
    }
    return json_text(v, true).size();
}

uint64_t compact_json_size(const Json::Value &v) {
    uint64_t values = 0;
    return compact_size(v, {}, &values);
}

uint64_t coordinates_text_bytes(const Json::Value &result, const std::function<void()> &check) {
    if (!result.isObject())
        return 0;
    uint64_t values = 0;
    // a member of an object that holds others (a seed's result always does): its quoted key,
    // the colon, its value and one comma
    auto member = [&](const char *key) -> uint64_t {
        if (!result.isMember(key))
            return 0;
        return 2 + std::strlen(key) + 1 + compact_size(result[key], check, &values) + 1;
    };
    uint64_t bytes = member("coordinates") + member("coordinates_reason");
    // a cut list's limitation (seed level, at most one), an element of the list: with a comma
    // when the list holds another
    bool cut = false;
    const Json::Value &lims = result["limitations"];
    if (lims.isArray()) {
        for (const Json::Value &l : lims) {
            if (l.isObject() && l["kind"] == "coordinates") {
                bytes += compact_json_size(l) + (lims.size() > 1 ? 1 : 0);
                cut = true;
            }
        }
    }
    // drop_coordinates: offered right after use_graphlet and drop_sequences, so never first
    // (resource_stop_json): the quoted token and its comma
    bool drop = false;
    if (result.isMember("resource_stop")) {
        for (const Json::Value &a : result["resource_stop"]["actions"]) {
            if (a == "drop_coordinates") {
                bytes += 2 + std::strlen("drop_coordinates") + 1;
                drop = true;
            }
        }
    }
    if (result.isMember("graphlet") && (cut || drop)) {
        // in the MGT body (a JSON string, escaped once more): the K record of the cut list and
        // the Q record's token (",drop_coordinates": printable, nothing to escape), and the Z
        // record's count, one record more with the K
        const char *begin = nullptr, *end = nullptr;
        result["graphlet"].getString(&begin, &end);
        const std::string_view body(begin, end - begin);
        uint64_t mgt = 0, escaped = 0;
        if (cut) {
            // K records follow the H, S, X, L, O and Q records, so never the first line
            const size_t at = body.find("\nK * coordinates ");
            if (at != std::string_view::npos) {
                const size_t eol = body.find('\n', at + 1);
                const std::string_view line = body.substr(at + 1, eol == std::string_view::npos
                                                                      ? std::string_view::npos
                                                                      : eol - at);
                mgt += line.size();
                escaped += json_escaped_size(line);
            }
        }
        if (drop) {
            mgt += 1 + std::strlen("drop_coordinates");
            escaped += 1 + std::strlen("drop_coordinates");
        }
        // Z states the line count, as graphlet_lines does: the K record is one line
        const uint64_t lines = result["graphlet_lines"].asUInt64();
        const uint64_t z = cut && lines ? decimal_digits(lines) - decimal_digits(lines - 1) : 0;
        mgt += z;
        escaped += z;
        bytes += escaped + z;
        // graphlet_bytes states the body's length
        const uint64_t body_bytes = result["graphlet_bytes"].asUInt64();
        if (body_bytes >= mgt)
            bytes += decimal_digits(body_bytes) - decimal_digits(body_bytes - mgt);
    }
    return bytes;
}

Json::Value seed_result_to_json(const SeedResult &r, const Strategy &st, const std::string &detail, bool timing) {
    const bool graphlet = detail == "graphlet";
    const ResourceStop *stop = stated_stop(r, st);
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
    if (r.derivation_partial) {
        // The time budget ran out while the permitted set was derived, after part of the seed (a
        // partial derivation, SPEC §7.0): the set of the k-mers read is a superset of the whole
        // seed's carriers, so a label here may not carry the whole seed — label evidence
        // qualified (outcome_of) — and the walk stopped at the seed. Observed: j, the k-mers
        // read (an integer), of num_kmers — NOT the elapsed ms that the same (kind, cause, knob)
        // observes on a seed whose derivation failed (derivation_out_of_time; a number). One
        // triple, two units, told apart by the result's shape (walked: arms, no `error`);
        // SPEC §7.0
        const uint64_t read = r.derivation_partial->kmers_read;
        Json::Value d = limitation(
                "derivation", "bounds.time_budget_ms", Json::Value(st.time_budget_ms),
                uint_json(read),
                "the time budget ran out while deriving the permitted set from the seed, after "
                + std::to_string(read) + " of " + std::to_string(r.num_kmers) + " k-mers: the "
                "labels carrying the k-mers read were taken, a superset of those carrying the "
                "whole seed (a label the whole seed would exclude may be among them"
                + (st.support == Support::TRACE ? ", and none was checked for one coordinate-"
                                                  "consecutive occurrence of the seed"
                                                : "")
                + "), and the walk stopped at the seed (complete_to_bp 0): raise the budget, "
                  "shorten the seed, or name the labels explicitly");
        // the fields a failed derivation's limitation has (n is the seed's num_kmers, in the
        // seed object and the S record): no field the contract does not already know
        d["cause"] = to_string(SeedDerivationError::TIME_BUDGET);
        lims.append(std::move(d));
    }
    if (r.labels_dropped) {
        // A set derived from part of the seed is a superset of the whole seed's carriers: the
        // cap cut it in index order, so what it dropped need not carry the whole seed and what
        // it kept may not either — the true carriers may be among either. No walk is missing
        // for them (the walk stopped at the seed); the levers are the cap, the time budget and
        // an explicit list
        lims.append(limitation("seed_labels", "labels.max_seed_labels", uint_json(st.max_seed_labels),
                               uint_json(r.labels_supporting_total),
                               r.derivation_partial
                               ? std::to_string(r.labels_dropped) + " of the "
                                 + std::to_string(r.labels_supporting_total) + " label(s) carrying "
                                 "the k-mers read were not taken (labels_dropped_digest identifies "
                                 "them): the set derived from part of the seed is a superset of the "
                                 "labels carrying the whole seed, so the ones dropped may include "
                                 "labels carrying the whole seed and the ones kept labels it would "
                                 "exclude; raise the knob or the time budget, or name the labels "
                                 "explicitly"
                               : std::to_string(r.labels_dropped) + " label(s) carrying the whole seed "
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
    if (st.max_memory_bytes) {
        // the budget is enforced on the modelled state and output (§14.1); a read that is
        // not budget-aware decodes a level's annotation rows before they can be charged,
        // so the bound is soft, and every response under it says so
        lims.append(memory_bound_soft(st, r.account));
    }
    if (st.coordinates) {
        // asked for: the block where they were recorded, otherwise null with the reason
        // (§18.1); without the request nothing
        if (r.coordinates_recorded) {
            j["coordinates"] = coordinates_json(r, st, &lims);
        } else {
            assert(r.coordinates_reason);
            j["coordinates"] = Json::Value();
            j["coordinates_reason"] = r.coordinates_reason ? r.coordinates_reason
                                                           : kCoordinatesNoTraversal;
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
        arms[to_string(arm)] = graphlet ? arm_summary_json(a, st, stop)
                                        : arm_to_json(a, st, detail, stop);
    }
    j["arms"] = std::move(arms);
    // §7.0: the guarantees are independent, so the outcome states each on its own axis
    // instead of folding them into one value that would read "partial" for a complete
    // walk with cut diagnostics, or "complete" for walks whose label lists were cut
    j["outcome"] = outcome_of(j, false);
    if (stop)
        j["resource_stop"] = resource_stop_json(*stop, st, nullptr, &r.account);
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
        // the seed phase (validation or derivation): its time, the part of it spent
        // resolving label names (a first header name builds the header index) and the part
        // spent reading (in annotation_fetch_ms too)
        const LabelOracle::Counters &c = r.annotation_counters;
        t["seed_phase_ms"] = c.seed_phase_seconds * 1000;
        t["seed_fetch_ms"] = c.seed_fetch_seconds * 1000;
        t["label_resolve_ms"] = c.label_resolve_seconds * 1000;
        // the deadline record (SPEC §6.8): the seed's longest uninterruptible piece, and what
        // stopped its walk how long after its deadline — which piece made a stop late
        Json::Value dl;
        const DeadlineRecord &d = r.deadline;
        Json::Value piece;
        piece["ms"] = d.longest.ms;
        piece["kind"] = d.longest.kind;
        piece["rows"] = uint_json(d.longest.rows);
        piece["coordinates"] = uint_json(d.longest.coordinates);
        dl["longest_piece"] = std::move(piece);
        if (!d.stopped_by.empty())
            dl["stopped_by"] = d.stopped_by;
        if (d.after_deadline_ms)
            dl["stop_after_deadline_ms"] = *d.after_deadline_ms;
        t["deadline"] = std::move(dl);
        if (c.path_cache_hits || c.path_cache_stored_rows || c.path_cache_peak_bytes) {
            // the row-diff path cache's physical work (only here: what the cache keeps
            // never changes a result, a charge or an admission)
            Json::Value pc;
            pc["hits"] = uint_json(c.path_cache_hits);
            pc["stored_rows_read"] = uint_json(c.path_cache_stored_rows);
            pc["rows_kept"] = uint_json(c.path_cache_rows_kept);
            pc["bytes_kept"] = uint_json(c.path_cache_bytes_kept);
            pc["peak_bytes"] = uint_json(c.path_cache_peak_bytes);
            t["path_cache"] = std::move(pc);
        }
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

Json::Value profile_to_json(const SupportProfile &p, const SeedSelection *sel, const std::string &run_format,
                            const std::function<void()> &check) {
    // the answer's deadline, read every kResolveCheckLabels objects (labels, candidates, seeds)
    size_t objects = 0;
    auto built = [&]() {
        if (check && ++objects % kResolveCheckLabels == 0)
            check();
    };
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
        built();
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
        built();
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
            built();
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

Json::Value coordinates_capabilities_json(const LabelOracle &oracle) {
    Json::Value c;
    // recorded under support trace only (§18.1), so exactly where trace is supported
    const bool supported = oracle.has_coordinates() && oracle.regime() == Regime::BASIC;
    c["supported"] = supported;
    c["knob"] = "output.coordinates";
    c["cap_knob"] = "output.max_coordinate_occurrences";
    c["max_occurrences_default"] = uint_json(Strategy().max_coordinate_occurrences);
    // what this index can report: record and mixed need header labels (a CoordToHeader)
    Json::Value kinds(Json::arrayValue);
    if (supported) {
        if (oracle.coord_to_header())
            kinds.append("record");
        kinds.append("column");
        if (oracle.coord_to_header())
            kinds.append("mixed");
    }
    c["kinds"] = std::move(kinds);
    c["limitation"] = "coordinates";
    c["action"] = "drop_coordinates";
    // what bounds the block's size (no server limit does on its own: --traverse-max-memory-mb
    // is 0 by default) and the coordinates' rule, both the SPEC's
    c["output_bound"] = kCoordinatesOutputBound;
    c["rule"] = kCoordinatesRule;
    return c;
}

Json::Value capabilities_to_json(const LabelOracle &oracle, const std::string &release,
                                 const IndexIdentity *identity) {
    Json::Value c;
    // the REQUEST schema this server accepts (strategy.schema_version): not a feature level
    c["schema_version"] = 1;
    // What the server offers beyond the base contract, monotonic: a client states a feature as
    // "feature_level >= n" (fields are only ever added; SPEC §10.3 states what each level adds;
    // kTraverseFeatureLevel). In the probe and in every response, so that a client reading only
    // responses states it too.
    c["feature_level"] = kTraverseFeatureLevel;
    c["k"] = uint_json(oracle.get_k());
    c["regime"] = to_string(oracle.regime());
    c["alphabet"] = oracle.graph().alphabet();
    c["num_labels"] = uint_json(oracle.num_columns());
    c["has_coordinates"] = oracle.has_coordinates();
    c["has_coord_to_header"] = oracle.coord_to_header() != nullptr;
    c["supports_trace"] = oracle.has_coordinates() && oracle.regime() == Regime::BASIC;
    c["cost_models_available"] = strings_json({ "forbid", "constant", "table" });
    c["label_modes"] = strings_json({ "constrain", "annotate" });
    c["direct_access"] = oracle.supports_direct();
    c["release"] = release;
    // the retrieval format (spec §7.5): `detail: graphlet` embeds MGT text of this version
    c["graphlet_format"] = kGraphletFormatVersion;
    c["detail_levels"] = strings_json({ "summary", "tree", "full", "graphlet" });
    // which index this is (DESIGN-traverse-graphlet.md §3.1): labels are joined across
    // retrievals only on an equal index_fp; null fp = no manifest, joins unverifiable;
    // index_meta_fp is a negative check only
    const IndexIdentity id = identity_or_default(identity, oracle);
    c["index_ns"] = string_or_null(id.name);
    c["index_fp"] = string_or_null(id.fp);
    c["index_meta_fp"] = id.meta_fp;
    return c;
}


// ---------------------------------------------------------------- MGT v1 codec

namespace mgt {

namespace {

void append_uint(std::string &out, uint64_t x) {
    char buf[24];
    auto res = std::to_chars(buf, buf + sizeof(buf), x);
    out.append(buf, res.ptr);
}

bool is_continuation(unsigned char c) { return (c & 0xC0) == 0x80; }


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

static bool ends_with(const std::string &s, std::string_view suffix) {
    return s.size() >= suffix.size()
        && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

const std::vector<IndexAnnotationKind>& index_annotation_kinds() {
    using namespace annot;
    // parse_annotation_type's extensions, in its order (the first suffix that matches decides),
    // with what the loader does beyond the file itself for each: RowDiff<ColumnMajor> reads
    // the row-diff anchors and fork successors beside the GRAPH (build_annotated_dbg), and a
    // MultiIntMatrix (coordinate) annotation the sequence headers (load_coord_to_header). A
    // unit test checks both against the types initialize_annotation constructs
    static const std::vector<IndexAnnotationKind> kinds = {
        { ColumnCompressed<>::kExtension, false, false },
        { ColumnCoordAnnotator::kExtension, false, true },
        { MultiBRWTCoordAnnotator::kExtension, false, true },
        { RowDiffCoordAnnotator::kExtension, false, true },
        { RowDiffBRWTCoordAnnotator::kExtension, false, true },
        { RowDiffColumnAnnotator::kExtension, true, false },
        { RowCompressed<>::kExtension, false, false },
        { MultiBRWTAnnotator::kExtension, false, false },
        { RowDiffBRWTAnnotator::kExtension, false, false },
        { BinRelWTAnnotator::kExtension, false, false },
        { RowFlatAnnotator::kExtension, false, false },
        { RowSparseAnnotator::kExtension, false, false },
        { RowDiffRowFlatAnnotator::kExtension, false, false },
        { RowDiffRowSparseAnnotator::kExtension, false, false },
        { RowDiffDiskAnnotator::kExtension, false, false },
        { IntRowDiffDiskAnnotator::kExtension, false, false },
        { RowDiffDiskCoordAnnotator::kExtension, false, true },
        { RainbowfishAnnotator::kExtension, false, false },
        { RbBRWTAnnotator::kExtension, false, false },
        { IntMultiBRWTAnnotator::kExtension, false, false },
        { IntRowDiffBRWTAnnotator::kExtension, false, false },
    };
    return kinds;
}

std::vector<IndexFile> index_load_inventory(const std::string &graph,
                                            const std::string &annotation,
                                            bool coord_mapping) {
    std::vector<IndexFile> files;
    // The paths are derived from the LISTED spelling, as the loaders derive them: a sidecar is
    // looked up next to a symlinked main file, not next to its target (otherwise two symlinks
    // to one graph and annotation could hide different .seqs files)
    files.push_back({ graph, "graph", true });
    // the dummy-edge mask and the Bloom filter DBGSuccinct::load reads beside the graph are
    // its derived data, not part of the identity: index_derived_files
    files.push_back({ annotation, "annotation", true });
    for (const IndexAnnotationKind &kind : index_annotation_kinds()) {
        if (!utils::ends_with(annotation, kind.extension))
            continue;
        if (kind.row_diff_anchors) {
            // config.infbase + kRowDiffAnchorExt: the graph's spelling, required (the loader
            // exits without them)
            files.push_back({ graph + annot::matrix::kRowDiffAnchorExt, "row_diff_anchors",
                              true });
            files.push_back({ graph + annot::matrix::kRowDiffForkSuccExt,
                              "row_diff_fork_succ", true });
        }
        if (kind.coordinates && coord_mapping) {
            const std::string seqs = utils::remove_suffix(annotation, kind.extension)
                                   + annot::CoordToHeader::kExtension;
            if (std::filesystem::exists(seqs))
                files.push_back({ seqs, "coord_to_header", false });
        }
        break;
    }
    return files;
}

std::vector<std::string> index_bundle_files(const std::string &graph,
                                            const std::string &annotation,
                                            bool coord_mapping) {
    std::vector<std::string> files;
    for (const IndexFile &f : index_load_inventory(graph, annotation, coord_mapping)) {
        files.push_back(f.path);
    }
    return files;
}

const char *const kIndexDerivedDataRule
        = "the dummy-edge mask (.edgemask) and the Bloom filter (.bloom) are derived data of "
          "the graph, not part of index_fp: they decide which counts are exact and how fast "
          "k-mers are looked up, never what an exact answer is, so adding, removing or "
          "rebuilding one leaves index_fp unchanged, and a manifest lists neither";

const char* index_derived_role(const std::string &path) {
    using graph::DBGSuccinct;
    const size_t slash = path.find_last_of('/');
    const std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    if (ends_with(name, DBGSuccinct::kDummyMaskExtension))
        return "graph_mask";
    if (ends_with(name, DBGSuccinct::kBloomFilterExtension))
        return "graph_bloom";
    return nullptr;
}

std::vector<IndexDerivedFile> index_derived_files(const std::string &graph) {
    using graph::DBGSuccinct;
    std::vector<IndexDerivedFile> files;
    if (!utils::ends_with(graph, DBGSuccinct::kExtension))
        return files;
    // DBGSuccinct::load: the dummy-edge mask when it opens, and only then the Bloom filter
    // when it exists (a mask that opens but does not fit the graph fails the load)
    const std::string prefix = utils::remove_suffix(graph, DBGSuccinct::kExtension);
    const std::string mask = prefix + DBGSuccinct::kDummyMaskExtension;
    const std::string bloom = prefix + DBGSuccinct::kBloomFilterExtension;
    const bool mask_loaded = std::ifstream(mask).good();
    const bool bloom_exists = std::filesystem::exists(bloom);
    files.push_back({ mask, "graph_mask", std::filesystem::exists(mask), mask_loaded });
    files.push_back({ bloom, "graph_bloom", bloom_exists, mask_loaded && bloom_exists });
    return files;
}

std::vector<std::string> index_unloaded_optional_files(const std::string &graph,
                                                       const std::string &annotation,
                                                       bool coord_mapping) {
    // every optional identity file the inventory could hold for the pair, as it derives the
    // paths (the graph's mask and Bloom filter are derived data: a manifest lists them never)
    std::vector<std::string> candidates;
    for (const IndexAnnotationKind &kind : index_annotation_kinds()) {
        if (!utils::ends_with(annotation, kind.extension))
            continue;
        if (kind.coordinates) {
            candidates.push_back(utils::remove_suffix(annotation, kind.extension)
                                 + annot::CoordToHeader::kExtension);
        }
        break;
    }
    const std::vector<std::string> loaded = index_bundle_files(graph, annotation, coord_mapping);
    std::vector<std::string> unloaded;
    for (const std::string &c : candidates) {
        if (std::find(loaded.begin(), loaded.end(), c) == loaded.end())
            unloaded.push_back(c);
    }
    return unloaded;
}

Json::Value index_inventory_json(const std::string &graph, const std::string &annotation,
                                 bool coord_mapping) {
    Json::Value j;
    j["graph"] = graph;
    j["annotation"] = annotation;
    j["coord_mapping"] = coord_mapping;
    Json::Value files(Json::arrayValue);
    for (const IndexFile &f : index_load_inventory(graph, annotation, coord_mapping)) {
        Json::Value e;
        e["path"] = f.path;
        e["role"] = f.role;
        e["required"] = f.required;
        e["exists"] = std::filesystem::exists(f.path);
        files.append(std::move(e));
    }
    j["files"] = std::move(files);
    // what the loader reads beside the graph that is not part of the identity: shown so that
    // an operator sees what is loaded, never listed in a manifest
    Json::Value derived(Json::arrayValue);
    for (const IndexDerivedFile &f : index_derived_files(graph)) {
        Json::Value e;
        e["path"] = f.path;
        e["role"] = f.role;
        e["exists"] = f.exists;
        e["loaded"] = f.loaded;
        derived.append(std::move(e));
    }
    j["derived"] = std::move(derived);
    j["derived_rule"] = kIndexDerivedDataRule;
    // the table the paths are derived by (scripts/traversal/index_manifest.py mirrors it, and
    // an integration test compares the two)
    Json::Value kinds(Json::arrayValue);
    for (const IndexAnnotationKind &kind : index_annotation_kinds()) {
        Json::Value k;
        k["extension"] = kind.extension;
        k["row_diff_anchors"] = kind.row_diff_anchors;
        k["coordinates"] = kind.coordinates;
        kinds.append(std::move(k));
    }
    j["annotation_kinds"] = std::move(kinds);
    return j;
}

std::string index_manifest_fingerprint(const std::string &manifest_path,
                                       const std::vector<std::string> &loaded,
                                       std::string *stated_name,
                                       const std::vector<std::string> &not_loaded) {
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
    // The graph's derived data is not part of the identity: a manifest that lists a mask or a
    // Bloom filter — any entry named *.edgemask
    // or *.bloom, whichever graph it belongs to — is refused rather than its entry skipped,
    // since skipping it would change the index_fp the manifest states without a word
    for (const auto &[p, entry] : files) {
        if (index_derived_role(p)) {
            const size_t slash = p.find_last_of('/');
            throw bad("it lists " + p + ": " + kIndexDerivedDataRule + " (write it again "
                      "without " + (slash == std::string::npos ? p : p.substr(slash + 1))
                      + ": scripts/traversal/index_manifest.py leaves them out; the new "
                        "manifest's index_fp differs from this one's)");
        }
    }
    // A manifest of another bundle must not lend its identity: every loaded file that
    // exists as named is listed (by base name) with its size. Hashing it is what the
    // manifest saves the server from (hundreds of GB at start-up).
    auto base_name = [](const std::string &p) {
        const size_t slash = p.find_last_of('/');
        return slash == std::string::npos ? p : p.substr(slash + 1);
    };
    // Base names unique: loaded files are matched to entries by base name, so a manifest of a
    // directory listing several bundles that share names (A/graph.dbg, B/graph.dbg,
    // A/annotation.seqs, B/annotation.seqs) would pass for every one of them, which would then
    // state one index_fp for different indexes (index_manifest.py writes bare base names, and
    // --verify refuses such a manifest too)
    std::map<std::string, std::string> by_name;
    for (const auto &[p, entry] : files) {
        auto [it, inserted] = by_name.emplace(base_name(p), p);
        if (!inserted) {
            throw bad("its entries " + it->second + " and " + p + " share the base name "
                      + it->first + ": loaded files are matched by base name, so a manifest "
                        "describes one bundle with distinct base names (index_manifest.py "
                        "writes one per pair)");
        }
    }
    // The optional files of this pair's inventory that it does not load (missing beside the
    // listed spelling, --no-coord-mapping) must not be listed either: index_fp would describe
    // a file set that is not the loaded one, while an index that does load the file states the
    // same index_fp. Files outside the inventory (--extra: a
    // column annotation's .coords, the .weights, anchors beside another annotation type) stay
    // allowed
    for (const std::string &path : not_loaded) {
        auto it = by_name.find(base_name(path));
        if (it != by_name.end()) {
            throw bad("it lists " + it->second + ", which the server does not load for this "
                      "index (" + path + ": missing beside the listed file, or "
                      "--no-coord-mapping): a manifest lists exactly the optional files of the "
                      "loader inventory that are loaded (traverse --index-inventory), so that "
                      "index_fp identifies the loaded files");
        }
    }
    for (const std::string &path : loaded) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f.good())
            continue;   // a required file that is missing: the loader names it
        const uint64_t size = static_cast<uint64_t>(f.tellg());
        bool named = false, listed = false;
        for (const auto &[p, entry] : files) {
            named |= base_name(p) == base_name(path);
            listed |= base_name(p) == base_name(path) && entry.first == size;
        }
        if (!named) {
            // an identity sidecar the loader reads (a .seqs, the row-diff anchors) that the
            // manifest does not cover could change answers under an unchanged fingerprint
            throw bad("it does not cover the file " + path + " (" + std::to_string(size)
                      + " bytes), which the server loads for this index: the manifest must "
                        "list every file of the loader inventory (traverse --index-inventory; "
                        "scripts/traversal/index_manifest.py lists them)");
        }
        if (!listed) {
            throw bad("the loaded file " + path + " (" + std::to_string(size)
                      + " bytes) is not listed with that size: the manifest describes another index");
        }
    }
    // One graph with one annotation: a manifest that also lists another graph or annotation
    // (one written for a directory, or with --extra) would lend one fingerprint to every index
    // of the directory, two annotations of one graph included (two annotations with
    // swapped memberships would state one index_fp and one index_meta_fp, and compare as
    // the same index)
    std::set<std::string> loaded_names;
    for (const std::string &path : loaded) {
        loaded_names.insert(base_name(path));
    }
    for (const auto &[p, entry] : files) {
        const std::string name = base_name(p);
        const bool annotation = ends_with(name, ".annodbg");
        if ((annotation || ends_with(name, "dbg")) && !loaded_names.count(name)) {
            throw bad("it lists the " + std::string(annotation ? "annotation " : "graph ") + p
                      + ", which this index does not load: a manifest describes one graph with "
                        "one annotation (scripts/traversal/index_manifest.py writes one per pair)");
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
    if (stated_name)
        *stated_name = json["index_ns"].isString() ? json["index_ns"].asString() : "";
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
    // --index-name is checked by the flag parser
    IndexIdentity id;
    id.name = config.index_name;
    if (!config.index_manifest.empty()) {
        if (config.infbase_annotators.size() != 1)
            throw std::runtime_error("--index-manifest: one annotation (-a) is loaded with it");
        id.fp = index_manifest_fingerprint(
                config.index_manifest,
                index_bundle_files(config.infbase, config.infbase_annotators[0],
                                   !config.no_coord_mapping),
                nullptr,
                index_unloaded_optional_files(config.infbase, config.infbase_annotators[0],
                                              !config.no_coord_mapping));
    }
    id.meta_fp = index_meta_fingerprint(LabelOracle(anno_graph));
    return id;
}

// ---------------------------------------------------------------- the graphlet (MGT v1)

namespace {

// End reasons as MGT codes (DESIGN-traverse-graphlet.md §2.1), in EndReason order; 'Y'
// is resource_limit, a head the memory or work budget of §14 did not admit.
const char kReasonCodes[] = "DLBRUVJTXSPNOMWY";
static_assert(sizeof(kReasonCodes) == kNumEndReasons + 1, "one MGT code per EndReason");

char reason_code(EndReason reason) { return kReasonCodes[static_cast<size_t>(reason)]; }

// The walker's text qualifier of a label end, kept as the code's second letter (the
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
        // adjacency without overflow (T's maximum + 1 wraps to 0)
        while (j + 1 < ids.size() && ids[j] != std::numeric_limits<T>::max()
                   && ids[j + 1] == static_cast<T>(ids[j] + 1)) j++;
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
        resource_stop();
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
            // process_traverse_request refuses a seed with a name that is not UTF-8; one
            // here would be replaced inside the JSON string, so graphlet_bytes and the
            // front coding (byte counts) would no longer describe the transported body
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

    // a JSON limitation value as a typed K value (§2.2): the knob's own type
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

    // Q: the result's resource_stop, the JSON object itself (DESIGN §14), when there is one
    void resource_stop() {
        if (!json_.isMember("resource_stop"))
            return;
        const Json::Value &q = json_["resource_stop"];
        out_ += "Q ";
        token(q["scope"].asString(), "resource_stop scope");
        sp(); token(q["resource"].asString(), "resource_stop resource");
        sp(); token(q["phase"].asString(), "resource_stop phase");
        for (const char *amount : { "requested", "effective", "used", "remaining" }) {
            sp();
            if (q[amount].isNull()) {
                out_.push_back('*');
            } else {
                out_ += mgt::encode_kvalue(kvalue(q[amount], amount));
            }
        }
        sp();
        if (q["actions"].empty())
            out_.push_back('.');
        for (Json::ArrayIndex i = 0; i < q["actions"].size(); ++i) {
            const std::string action = q["actions"][i].asString();
            for (char c : action) {
                if (!(std::islower(static_cast<unsigned char>(c)) || c == '_'))
                    unrepresentable("resource_stop action '" + action + "'");
            }
            out_ += i ? "," : "";
            out_ += action;
        }
        sp();
        out_ += mgt::pct_escape(q["message"].asString());
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
        // paths: leaf ordinal in segment-id order (finalize). The chain through parents[0]
        // holds by construction: a path stores only its leaf
        size_t last_leaf = 0;
        for (size_t i = 0; i < arm.paths.size(); ++i) {
            const PathResult &p = arm.paths[i];
            if (p.id != i)
                unrepresentable("path ids are not ordinals");
            const size_t leaf = p.leaf;
            if (leaf >= segs.size() || ix.path_of[leaf] >= 0 || (i && leaf <= last_leaf))
                unrepresentable("paths are not one per leaf in segment order");
            last_leaf = leaf;
            ix.path_of[leaf] = i;
            if (p.length_bp != segs[leaf].from_bp + segs[leaf].length_bp)
                unrepresentable("a path whose length is not its leaf's end");
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
        for (const auto &[name, value] : arm_counters(arm, st_)) {
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
            delivery_tick();
            write_segment(arm, s, ix);
        }
        // ---- R: the runs, in ArmResult::runs order (run ids keep their meaning)
        std::vector<double> needed;
        for (size_t i = 0; i < arm.runs.size(); ++i) {
            delivery_tick();
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
            walk_path_leaf_first(arm, p, [&](size_t s) {
                // walking order is the stored order on the right, reversed on the left
                const std::string &w = arm.segments[s].sequence;
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
                return rev.size() < from_walk;
            });
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

// The length jsoncpp's writer gives |s| inside a JSON string, quotes excluded (its
// escaping with emitUTF8 off, which every writer here uses): '"' '\\' and \b \f \n \r \t
// take two bytes, any other control character six (\u00XX), a code point above U+007F six
// (twelve as a surrogate pair), and a malformed sequence the six of U+FFFD for the bytes
// the writer consumes with it. GraphletCodec.JsonEscapedSizeIsTheWriters checks it against
// the writer itself.
uint64_t json_escaped_size(std::string_view s) {
    uint64_t n = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == '"' || c == '\\' || c == '\b' || c == '\f' || c == '\n' || c == '\r' || c == '\t') {
            n += 2;
        } else if (c < 0x20) {
            n += 6;
        } else if (c < 0x80) {
            n += 1;
        } else if (c >= 0xF8) {
            n += 6;     // no lead byte: one replacement character
        } else {
            // a lead byte claims up to three more bytes, as jsoncpp's utf8ToCodepoint does
            // (without checking that they are continuation bytes); a sequence cut by the
            // end of the string is one replacement character for its lead byte
            const size_t width = c < 0xE0 ? 2 : c < 0xF0 ? 3 : 4;
            if (s.size() - i < width) {
                n += 6;
                continue;
            }
            n += width == 4 ? 12 : 6;
            i += width - 1;
        }
    }
    return n;
}

// the decimal exponent of x > 0, finite: x = d.ddd x 10^E (exact, from the shortest digits)
static int decimal_exponent(double x) {
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof(buf), x, std::chars_format::scientific);
    const std::string_view sci(buf, res.ptr - buf);
    return static_cast<int>(std::strtol(std::string(sci.substr(sci.find('e') + 1)).c_str(),
                                        nullptr, 10));
}

// The widest mgt::encode_float of any value in [lo, 1) (0 when lo >= 1): positional, with
// lo's leading zeros after the point and at most 17 significant digits, "0." + zeros + 17
static uint64_t width_from(double lo) {
    if (!(lo > 0) || lo >= 1)
        return 0;
    return static_cast<uint64_t>(2 + (-decimal_exponent(lo) - 1) + 17);
}

// The widest encoding of any value in [1, hi] (0 when hi < 1): the integer digits of hi, or
// 17 significant digits and the point; every finite double has at most 309 integer digits
static uint64_t width_upto(double hi) {
    if (!(hi >= 1))
        return 0;
    if (std::isinf(hi))
        return 309;
    return std::max<uint64_t>(18, decimal_exponent(hi) + 1);
}

uint64_t mgt_float_width(const Strategy &st, const LabelChangeCost &cost,
                         double requested_time_ms, double max_time_ms) {
    uint64_t width = kMgtFloatWidth;
    // A result's floats are costs (E, R), losses (R, T, C) and the budget a loss-budget end
    // needed (R), the time a time stop compared (A, Q, K) and the time knobs (Q, K). A loss is
    // a sum of the costs of the switches taken (positive ones at least the smallest positive
    // cost) and at most the loss budget; a needed budget is a loss plus one cost: so every
    // such value is 0, +inf, or in [smallest positive cost, loss budget + largest cost], and
    // the sums are rounded monotonically. Canonical MGT writes them positionally: 1e-300 is
    // 302 characters (switches of 1e-300 priced at 24 characters a float would deliver
    // 22 MB within an account of 16 MiB).
    std::vector<double> costs;
    switch (cost.model()) {
        case LabelChangeCost::FORBID:
            break;
        case LabelChangeCost::CONSTANT:
            costs.push_back(cost.default_cost());
            break;
        case LabelChangeCost::TABLE:
            costs.push_back(cost.default_cost());
            for (const auto &entry : cost.table()) {
                costs.push_back(entry.second);
            }
            break;
    }
    double smallest = kInfiniteLoss, largest = 0;
    for (double c : costs) {
        if (c > 0 && c != kInfiniteLoss) {
            smallest = std::min(smallest, c);
            largest = std::max(largest, c);
        }
    }
    if (smallest != kInfiniteLoss) {
        width = std::max(width, width_from(smallest));
        width = std::max(width, width_upto(st.loss_budget + largest));
    }
    // A time stop compares the elapsed time, at least the budget, with it; an elapsed time
    // is measured in clock ticks (1 ns: at least 0.000001 ms, so at most 24 characters below
    // 1 ms) and below 10^17 ms; a budget of +inf is written "inf"
    for (double t : { st.time_budget_ms, requested_time_ms, max_time_ms }) {
        if (t > 0 && !std::isinf(t))
            width = std::max({ width, width_from(t), width_upto(t) });
    }
    return width;
}

/**
 * What one object costs this response to deliver in |detail| (bytes, upper bounds checked
 * against the serialisers by Graphlet.DeliveryCostsBoundTheOutput, adversarial names
 * included), for the memory budget (DESIGN-traverse-graphlet.md §14: "delivery is accounted
 * per expansion"). Every bound is a worst case, not an average (delivery margins come from
 * demonstrated upper bounds, escaping included).
 *
 * JSON (summary / tree / full, and the summary of a graphlet). jsoncpp's tree costs per
 * value at most: an object member a std::map node (rb-tree links and colour, the CZString key,
 * the Value: 96 B with malloc's 16-byte quantum and 8-byte header) plus its key's copy (48 B:
 * every key here is at most 25 characters); an array element a node (96 B); an object or
 * array value its map (64 B); a string value its length prefix, NUL and rounding (48 B for
 * the short ones, 32 B plus its bytes for a long one). Its text at most: in the CLI's
 * indented form (two spaces a level, at most ten levels deep in a traversal response) a
 * member writes a line break, its indent, the quoted key, " : ", its value (24 B for a number
 * at 17 significant digits or 20 digits, 30 B for the longest fixed string) and a comma: 80 B;
 * an element 48 B; an object's or array's closing line 24 B; the server's compact form is
 * shorter. Json::writeString holds the text three times at the peak: its stream's buffer,
 * which grows by doubling (less than twice the text), and the copy it returns, while the tree
 * is still alive. So a member costs 144 + 3 * 80, an element 96 + 3 * 48, a container
 * 64 + 3 * 24; a base 1 + 3; a name its escaped length three times plus its bytes and 32
 * (DeliveryCosts::name), per place it is written. Compressing the text later holds the text
 * and the deflate output (no larger than the text plus 0.03%) beside zlib's 256 KiB state,
 * less than the writeString peak plus that state, which the fixed part carries.
 *
 * The graphlet: the MGT writer's text grows by doubling (less than twice the body p) and is
 * copied into the JSON value, 3p at once; then the JSON value is held while writeString holds
 * the escaped body E three times, p + 3E. E is p for every record field (printable ASCII, no
 * quote) plus one byte for the line break that ends a record, so a record of b bytes costs
 * 4b + 3, a name p + 3E with p its percent-escaped length; plus the writer's per-arm index
 * (per segment, run and label end). Record fields are bounded at their widest: ids, counts
 * and positions below 10^10 (10 digits; a budget of at most 1 TiB charges more than 110 B
 * per object, so no count reaches 10^10), floats at |float_width| characters (24 unless the
 * request's costs, loss budget or time budget can be written wider: mgt_float_width).
 *
 * The fixed part holds zlib's deflate state (256 KiB), the envelope (the seed object, its
 * outcome, annotation, timing and resource_stop, both arms' certificates, counters and
 * empty lists; a graphlet's arm summaries and its H S O Q A Z records) and every limitation
 * a result can state: at most 5 at the seed level (seed_labels, trace_record_boundaries,
 * memory_bound_soft, two server_clamp) and 8 per arm (a stopping cap and a beam's
 * walk_domain, branch_events, label_lists, switch_sources, inexact_counts, scope,
 * greedy_losses), each with an effect of at most 640 bytes; with recorded coordinates one more
 * (a cut list's `coordinates`), and a permitted set derived from part of the seed states its
 * `derivation`, which the walker charges where it applies (DeliveryCosts::extra_limitation). A path chain costs only where
 * JSON spells it (tree, full).
 */
DeliveryCosts delivery_costs(const std::string &detail, bool sequences, uint64_t float_width,
                             CoordinatesOutput coordinates) {
    DeliveryCosts d;
    // the characters every float can take beyond the 24 the record bounds below assume
    const uint64_t wide = float_width > kMgtFloatWidth ? float_width - kMgtFloatWidth : 0;
    constexpr uint64_t kTextCopies = 3;
    constexpr uint64_t kMember = 96 + 48 + kTextCopies * 80;     // 384
    constexpr uint64_t kElement = 96 + kTextCopies * 48;         // 240
    constexpr uint64_t kMap = 64 + kTextCopies * 24;             // 136
    constexpr uint64_t kShort = 48;                              // a short string's buffer
    constexpr uint64_t kLongString = 32;                         // a long string's, beyond its bytes
    constexpr uint64_t kLimitations = 5 + 2 * 8;
    constexpr uint64_t kEffect = 640;
    // a limitation: its element and object, at most 8 members (kind, knob, limit, observed,
    // complete_to_bp, effect, one extra), 3 short strings and the effect (escaped: its two
    // quotes around "unlimited" add two bytes)
    constexpr uint64_t kLimitationJson = kElement + kMap + 8 * kMember + 3 * kShort
                                       + kEffect + kLongString + kTextCopies * (kEffect + 4);
    // the envelope: the seed object (11 members, 3 lists), the result's own members (11, 7
    // maps), outcome (4), annotation (4), timing (9), resource_stop (9 members, its actions,
    // a message of at most 1 KiB), and per arm its member, its certificate (21 members, 5
    // maps), lists and counters (21 members, 9 maps)
    constexpr uint64_t kMessage = 1024;
    constexpr uint64_t kEnvelopeJson = (11 + 11 + 4 + 4 + 9 + 9 + 2 * 43) * kMember
                                     + (3 + 7 + 3 + 2 * 14) * kMap + 8 * kElement + 24 * kShort
                                     + kMessage + kLongString + kTextCopies * (kMessage + 8);
    constexpr uint64_t kDeflate = 1 << 18;
    // the name of a dictionary label, a dropped label or the seed: how often the detail
    // writes it, and what each written copy costs
    auto json_name = [](std::string_view name, uint64_t copies) {
        return copies * (name.size() + kLongString + kTextCopies * (json_escaped_size(name) + 2));
    };
    // Record coordinates (§18.2; opt-in, so nothing here without them), JSON in every detail
    // (a graphlet's in its summary). Fixed: the block's member, its object and 7 members (kind,
    // k, max_occurrences, complete, runs_lower_bound, seed, arms), the seed list, the arms
    // object with 2 members and their lists, 2 short strings, and the two echo members of
    // strategy.output (charged per seed) — which bounds the null form with its reason (2
    // members, a short string) — and, where a block can be recorded, one limitation more (a cut
    // list's). A run's entry: an element, its object, 8 members and its occurrence list; a seed
    // label's: an element, its object, 3 members and its list; an occurrence: an element, its
    // pair's array and two elements (numbers of at most 20 digits, the 48 B an element's text
    // takes). A seed-level limitation beyond the fixed part's (a partial derivation's): one more
    if (coordinates != CoordinatesOutput::NONE) {
        d.coordinate_run = kElement + kMap + 8 * kMember + kMap;
        d.coordinate_seed = kElement + kMap + 3 * kMember + kMap;
        d.occurrence = 3 * kElement + kMap;
    }
    const uint64_t coordinate_fixed = coordinates == CoordinatesOutput::BLOCK
        ? 12 * kMember + 5 * kMap + 3 * kShort + kLimitationJson
        : coordinates == CoordinatesOutput::REASON ? 4 * kMember + kShort : 0;
    d.extra_limitation = kLimitationJson;
    if (detail == "graphlet") {
        constexpr uint64_t kCopies = 4;    // p + 3E per record byte (E = p, see above)
        constexpr uint64_t kLine = 3;      // the escaped line break of a record, three times
        auto record = [](uint64_t widest) { return kCopies * widest + kLine; };
        // the JSON summary: the envelope with per arm a summary (85 members, 9 maps) in
        // place of its lists, every limitation in JSON, and in the body H S O Q A Z and a K
        // record per limitation (at most 240 B of fields with its effect)
        constexpr uint64_t kEnvelope = (11 + 11 + 4 + 4 + 9 + 9 + 2 * 85) * kMember
                                     + (3 + 7 + 3 + 2 * 9) * kMap + 8 * kElement + 24 * kShort
                                     + kMessage + kLongString + kTextCopies * (kMessage + 8);
        constexpr uint64_t kRecords = 4 * 1024 + kMessage + 2 * 1024
                                    + kLimitations * (240 + kEffect);
        // the floats of the fixed records: Q's four amounts, each A's cap trigger, each K's
        // limit and observed value
        const uint64_t fixed_floats = 4 + 2 + 2 * kLimitations;
        d.fixed = kDeflate + kEnvelope + kLimitations * kLimitationJson
                + kCopies * (kRecords + fixed_floats * wide) + 16 * kLine;
        // a limitation more also writes its K record (240 B of fields with its effect, its
        // limit and observed floats)
        const uint64_t k_record = kCopies * (240 + kEffect + 2 * wide) + kLine;
        d.coordinate_fixed = coordinate_fixed
                           + (coordinates == CoordinatesOutput::BLOCK ? k_record : 0);
        d.fixed += d.coordinate_fixed;
        d.extra_limitation += k_record;
        // L <c|h> <column> <seq_id> <prefix_len> <suffix>: 38 B without the name
        d.label = record(40);
        // G: parents (the first), from_bp, length_bp, the set codes and entry_total, split,
        // first_base, the bases' marker: 64 B; plus the writer's index per segment
        d.segment = record(64) + 40;
        d.segment_label = kCopies * 11;
        d.merge_parent = kCopies * 13;     // ",<id>" and "|" with its list
        d.base = sequences ? kCopies : 0;
        d.continuation_base = kCopies;
        // R: seven ids and positions, the end code, from:cost, branches, loss, needed: 172 B
        // with three floats; its entry in the writer's index of anchored runs
        d.run = record(172 + 3 * wide) + 8;
        // E: at most 64 B of fields (a switch's cost one float); a label end's node in the
        // writer's index instead
        d.event = record(64 + wide) + 96;
        d.event_label = kCopies * 11;
        // T, and C without its labels and bases (its loss_used one float)
        d.leaf = record(32) + record(56 + wide);
        d.leaf_label = kCopies * (70 + wide);   // label:loss:branches:route_bp, a label in C
        d.split_branch = kCopies;          // its first base
        d.branch_event = record(60);
        d.branch_event_entry = kCopies * 12;
        d.refusal = kCopies * 24;
        d.presence_run = record(48);
        d.bin = record(14 * 11 + kNumEndReasons * 13 + 2);
        d.dropped = record(24);            // X <reason> without its runs and name
        d.dropped_run = kCopies * 22;
        d.name = [](std::string_view name, DeliveryCosts::Name use) -> uint64_t {
            if (use == DeliveryCosts::Name::SEED_ID) {
                // the JSON summary's seed.seed_id: the body holds no seed_id
                return name.size() + kLongString + kTextCopies * (json_escaped_size(name) + 2);
            }
            // one L or X record each: the percent-escaped name p, then p + 3E
            const std::string escaped = mgt::pct_escape(name);
            return escaped.size() + kTextCopies * json_escaped_size(escaped);
        };
        return d;
    }
    const bool segments = detail != "summary";
    d.coordinate_fixed = coordinate_fixed;
    d.fixed = kDeflate + kEnvelopeJson + kLimitations * kLimitationJson + coordinate_fixed;
    // label_dict: element, object, 4 members, its kind; label_summary: element, object,
    // "label" and per arm a member with an object of 4 members and the runs list; a seed
    // label's element in seed.labels
    d.label = (kElement + kMap + 4 * kMember + kShort) + (kElement + 5 * kMap + 11 * kMember)
            + kElement;
    // element, object and 12 members (id, parents, labels_via_parent, from_bp, length_bp,
    // sequence, labels, labels_total, labels_truncated, labels_at_end, events, label_sets),
    // 6 lists, its first parent's id, the sequence's buffer
    d.segment = segments ? 2 * kElement + 7 * kMap + 12 * kMember + kLongString : 0;
    d.segment_label = segments ? kElement : 0;
    d.merge_parent = segments ? 2 * kElement + kMap : 0;
    d.base = segments && detail == "full" && sequences ? 1 + kTextCopies : 0;
    d.continuation_base = 1 + kTextCopies;
    // a run: element, object, 11 members, entered_by and end_reason (tree, full), and its id
    // in label_summary's runs (every detail)
    d.run = (segments ? kElement + kMap + 11 * kMember + 2 * kShort : 0) + kElement;
    // an event: element, object, its list, 7 members, 3 short strings (tree, full); a
    // loss-budget end's needed_budgets element (every detail)
    d.event = (segments ? kElement + 2 * kMap + 7 * kMember + 3 * kShort : 0) + kElement;
    d.event_label = segments ? kElement : 0;
    // a path: element, object, 12 members (with its continuation's), 6 lists and maps,
    // path_reason and the continuation's buffer
    d.leaf = kElement + 6 * kMap + 12 * kMember + kShort + kLongString;
    // an end label: element, object, 5 members (tree, full); its reason in the path's
    // end_reasons (every detail) and, charged alike, a continuation label's element
    d.leaf_label = (segments ? kElement + kMap + 5 * kMember : kElement) + kMember;
    d.chain_entry = segments ? kElement : 0;
    // element, object, children, branches, 7 members, kind
    d.split = segments ? kElement + 3 * kMap + 7 * kMember + kShort : 0;
    // element, object, labels, 5 members, char; its id in children
    d.split_branch = segments ? 2 * kElement + 2 * kMap + 5 * kMember + kShort : 0;
    // element, object, 4 lists, 7 members, successors' buffer (every detail)
    d.branch_event = kElement + 5 * kMap + 7 * kMember + kShort;
    // an element, and its successor's character (1 + 3)
    d.branch_event_entry = kElement + 1 + kTextCopies;
    // element, object, labels, 3 members, char and cause
    d.refusal = kElement + 2 * kMap + 3 * kMember + 2 * kShort;
    d.presence_run = segments ? kElement + 2 * kMap + 5 * kMember : 0;
    // element, object, label_ends, 14 members and one per end reason (every detail)
    d.bin = kElement + 2 * kMap + (14 + kNumEndReasons) * kMember;
    // element, object, runs, 4 members, reason and runs_kind; per run its pair
    d.dropped = kElement + 2 * kMap + 4 * kMember + 2 * kShort;
    d.dropped_run = 3 * kElement + kMap;
    d.name = [json_name](std::string_view name, DeliveryCosts::Name use) -> uint64_t {
        // a seed label's name is in label_dict and seed.labels, every other once
        return json_name(name, use == DeliveryCosts::Name::SEED_LABEL ? 2 : 1);
    };
    return d;
}

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

// The strategy.clamped entry of the server's maximum (SPEC §10.3) of the budget a stop ran
// into — bounds.max_memory_mb for a memory stop, bounds.max_work_units for a work stop — or
// null: the budget was the request's own (a smaller one is kept), or the stop is no budget's
static const Json::Value* budget_clamp(const Json::Value &clamped, ResourceStop::Resource resource) {
    const char *field = resource == ResourceStop::MEMORY ? "bounds.max_memory_mb"
                      : resource == ResourceStop::WORK ? "bounds.max_work_units" : nullptr;
    if (!field)
        return nullptr;
    for (const Json::Value &c : clamped) {
        if (c["field"].asString() == field)
            return &c;
    }
    return nullptr;
}

// What a seed stopped by a budget at the server's maximum (SPEC §10.3) states, walked or failed
// in its seed phase alike: the server_clamp limitation; the stop's `requested`, what the request
// asked for ("unlimited": it gave no budget); no action raising the budget, and the knob's
// walk_domain entries carry server_limit (as a derivation's and seed_labels' do) and say that a
// request cannot raise it — an arm's "; raise the knob" is replaced by that statement, a failed
// seed's effect is written with it (budget_failed_seed_to_json)
static void state_budget_clamp(Json::Value *rj, const Json::Value &c) {
    const std::string field = c["field"].asString();
    Json::Value &lims = (*rj)["limitations"];
    if (rj->isMember("resource_stop")) {
        Json::Value &q = (*rj)["resource_stop"];
        q["requested"] = c["requested"];
        const char *raise = field == "bounds.max_memory_mb" ? "raise_memory_budget"
                                                            : "raise_work_budget";
        Json::Value actions(Json::arrayValue);
        for (const Json::Value &a : q["actions"]) {
            if (a.asString() != raise)
                actions.append(a);
        }
        q["actions"] = std::move(actions);
    }
    static const std::string kRaise = "; raise the knob";
    auto at_max = [&](Json::Value &l) {
        if (l["kind"].asString() != "walk_domain" || l["knob"].asString() != field)
            return;
        l["server_limit"] = c["effective"];
        std::string effect = l["effect"].asString();
        if (effect.size() >= kRaise.size()
                && effect.compare(effect.size() - kRaise.size(), kRaise.size(), kRaise) == 0) {
            effect.resize(effect.size() - kRaise.size());
            l["effect"] = effect + kAtServerMaximum;
        }
    };
    for (Json::Value &l : lims) {
        at_max(l);
    }
    // (looked up, not indexed: indexing would add the member to an arm that has none)
    if (rj->isMember("arms")) {
        for (const char *arm : { "left", "right" }) {
            if ((*rj)["arms"].isMember(arm) && (*rj)["arms"][arm].isMember("limitations")) {
                for (Json::Value &l : (*rj)["arms"][arm]["limitations"]) {
                    at_max(l);
                }
            }
        }
    }
    lims.append(limitation("server_clamp", field, c["effective"], c["requested"],
                           c["requested"].isString()
                             ? "the request gave no budget and the server set its maximum, which "
                               "this seed ran into; a request cannot raise it further"
                             : "the server lowered the requested value to its maximum and this "
                               "seed ran into it; a request cannot raise it further"));
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
        // ("unlimited": an omitted budget the server set to its maximum)
        const bool lowered = c["requested"].isNumeric()
                && c["effective"].asDouble() < c["requested"].asDouble();
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
            // a set derived from part of the seed names the clamped knob too, with the
            // server's value, as a failed derivation's entry does (failed_seed_to_json)
            for (Json::Value &l : lims) {
                if (l["kind"].asString() == "derivation" && l["knob"].asString() == field)
                    l["server_limit"] = c["effective"];
            }
            // a time stop states what the request asked for beside what bound it
            if (rj->isMember("resource_stop")
                    && (*rj)["resource_stop"]["resource"].asString() == "time")
                (*rj)["resource_stop"]["requested"] = c["requested"];
        } else if (field == "bounds.max_memory_mb" || field == "bounds.max_work_units") {
            // the server's maximum of a budget bound this seed when the budget stopped
            // its walk (state_budget_clamp)
            const ResourceStop::Resource resource = field == "bounds.max_memory_mb"
                ? ResourceStop::MEMORY : ResourceStop::WORK;
            if (r.resource_stop && r.resource_stop->resource == resource)
                state_budget_clamp(rj, c);
            continue;
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

// What a result failed or refused without an admitted result holds until the response is
// written: the fixed part of a delivered result (an upper bound of this shorter one) and its
// echo of seed_id in every copy the serialisers hold (failed_soft's price; also the usage's
// final_bytes of a failed or never started seed)
static uint64_t failed_result_bytes(const Strategy &st, const Seed &seed) {
    const DeliveryCosts &d = st.delivery;
    return d.fixed + seed.seed_id.size()
        + (d.name ? d.name(seed.seed_id, DeliveryCosts::Name::SEED_ID) : 0);
}

// memory_bound_soft's observed excess (bytes) for a seed failed or refused without an
// admitted result. Such a result is not admitted, and it echoes the request's seed_id,
// which only the request's size bounds: a seed_id that the depth-0 admission refused still
// comes back whole (100,000 emoji: a 1.2 MB result under 1 MiB). It is priced as a delivered
// result prices it — the fixed part, an
// upper bound of this shorter result, and the seed_id in every copy the serialisers hold —
// and what exceeds the budget is stated with what the walk observed. Index-supplied names
// in its messages are cut under a budget (Walker::echoed), so they add no more than a
// bounded prefix.
static uint64_t failed_soft(const Strategy &st, const Seed &seed, uint64_t observed) {
    if (!st.max_memory_bytes)
        return observed;
    const uint64_t held = failed_result_bytes(st, seed);
    return std::max(observed, held > st.max_memory_bytes ? held - st.max_memory_bytes : 0);
}

// The head of a result without a walk (a failed, refused or never started seed): its seed
// (seed_id, length_bp, labels_from_seed) and |error|
static Json::Value failed_seed_head(const Seed &seed, bool labels_from_seed, std::string error) {
    Json::Value rj;
    Json::Value &sj = rj["seed"];
    sj["seed_id"] = seed.seed_id;
    sj["length_bp"] = uint_json(seed.sequence.size());
    sj["labels_from_seed"] = labels_from_seed;
    rj["error"] = std::move(error);
    return rj;
}

// What such a result states under a memory budget (§7.0): memory_bound_soft, with what the seed
// was seen to hold beyond the budget (|observed|) and the result's echo of seed_id
static void append_failed_memory(Json::Value *lims, const Strategy &st, const Seed &seed,
                                 uint64_t observed, bool decode_charged) {
    if (!st.max_memory_bytes)
        return;
    ResourceAccount account;
    account.memory_limit = st.max_memory_bytes;
    account.soft_overshoot = failed_soft(st, seed, observed);
    account.decode_charged = decode_charged;
    lims->append(memory_bound_soft(st, account));
}

// A failed, refused or never started seed of a request that asked for coordinates states that
// it has none, and why: |reason| is computed once per request — the index's or the support's
// when either rules them out whatever the seed, otherwise "no traversal" (§18.1). Null: not
// asked for, nothing is added
static void state_no_coordinates(Json::Value *rj, const char *reason) {
    if (!reason)
        return;
    (*rj)["coordinates"] = Json::Value();
    (*rj)["coordinates_reason"] = reason;
}

// The result of a seed whose permitted set could not be DERIVED (§6.1 step 4): no arms,
// `outcome.walks: failed`, the message as `error`, and the cause as a `derivation` limitation
// naming the request field that would get past it (§7.0) — so that an agent acts on the
// knob instead of parsing the message. |clamped| supplies `server_limit` when the server
// clamped that knob: raising it beyond the server's value does nothing.
static Json::Value failed_seed_to_json(const Seed &seed, const SeedDerivationError &e,
                                       const Strategy &st, const Json::Value &clamped,
                                       bool decode_charged, const char *coordinates_reason) {
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
        case SeedDerivationError::UNREPRESENTABLE_LABEL_NAME:
            // stated by unrepresentable_seed_to_json, which knows where the name came from
            knob = "seeds[].labels";
            limit = Json::Value();
            observed = uint_json(static_cast<uint64_t>(e.observed()));
            effect = e.what();
            break;
    }
    // a derivation is made only for such a seed: labels_from_seed true, by the rule every failed
    // writer states
    Json::Value rj = failed_seed_head(seed, labels_derived_from_seed(seed, st.label_mode),
                                      e.what());
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
    // every response under a memory budget states it (§7.0), a failed derivation's too: it
    // decoded whole annotation rows that no admission charged, and observed is what it was seen
    // to hold beyond the budget
    append_failed_memory(&lims, st, seed, e.soft_overshoot(), decode_charged);
    rj["limitations"] = std::move(lims);
    // no walk was made, so nothing was cut on the other axes — except the carriers a cap
    // cut before the trace check (seed_labels), whose evidence is then missing
    rj["outcome"] = outcome_of(rj, true);
    state_no_coordinates(&rj, coordinates_reason);
    return rj;
}

// A seed whose label dictionary (or dropped labels) would hold a name that is not valid
// UTF-8 (spec §6.1 step 4, §7.0): refused per seed, in both label modes, in the shape of a
// failed derivation -- no arms, `outcome.walks: failed`, an `error`, and a `derivation`
// limitation with cause `unrepresentable_label_name`. Names come from FASTA headers and
// file names, which need not be UTF-8; no output (JSON, MGT) carries such a name verbatim,
// and a REPLACED name (U+FFFD) can be another label's name anywhere in the index: a
// continuation that resubmitted it would go on under that other label. So nothing of
// the seed is delivered, and neither the
// error nor the effect echoes the name's bytes: they name the label by its column (and
// sequence id). The knob is the one that avoids recording the label where one exists.
// Returns a null value when every name is valid.
static Json::Value unrepresentable_seed_to_json(const Seed &seed, const SeedResult &r,
                                                const Strategy &st,
                                                const char *coordinates_reason) {
    size_t bad = 0;
    std::optional<size_t> first_dict;       // index into label_dict
    std::optional<size_t> first_dropped;    // index into dropped_labels
    for (size_t i = 0; i < r.label_dict.size(); ++i) {
        if (!mgt::valid_utf8(r.label_dict[i].name)) {
            bad++;
            if (!first_dict)
                first_dict = i;
        }
    }
    for (size_t i = 0; i < r.dropped_labels.size(); ++i) {
        if (!mgt::valid_utf8(r.dropped_labels[i].name)) {
            bad++;
            if (!first_dropped)
                first_dropped = i;
        }
    }
    if (!bad)
        return Json::Value();
    std::string where;
    std::string knob;
    Json::Value limit;
    std::string lever;
    if (first_dict) {
        const LabelRef &l = r.label_dict[*first_dict];
        where = "column " + std::to_string(l.column);
        if (l.kind == LabelKind::HEADER)
            where += ", sequence " + std::to_string(l.seq_id);
        if (st.label_mode == LabelMode::ANNOTATE) {
            knob = "labels.mode";
            limit = "annotate";
            lever = "traverse in constrain mode with an explicit list of labels whose names "
                    "are UTF-8 (resolve lists the labels on the seed)";
        } else if (*first_dict >= r.num_seed_labels) {
            knob = "labels.extra";
            limit = uint_json(static_cast<uint64_t>(r.label_dict.size() - r.num_seed_labels));
            lever = "leave that label out of labels.extra";
        } else if (r.labels_from_seed && l.kind == LabelKind::HEADER) {
            knob = "labels.seed_label_kind";
            limit = "header";
            lever = "set the knob to \"column\" (labels are then the annotation columns) or "
                    "name the labels explicitly, leaving that one out";
        } else {
            knob = "seeds[].labels";
            limit = uint_json(static_cast<uint64_t>(r.num_seed_labels));
            lever = r.labels_from_seed
                ? "name the labels explicitly, leaving that one out"
                : "leave that label out of seeds[].labels";
        }
    } else {
        where = "dropped seed label " + std::to_string(*first_dropped + 1) + " of "
              + std::to_string(r.dropped_labels.size());
        knob = "seeds[].labels";
        limit = uint_json(static_cast<uint64_t>(r.num_seed_labels + r.dropped_labels.size()));
        lever = "leave that label out of seeds[].labels";
    }
    const std::string count = std::to_string(bad) + " label name" + (bad > 1 ? "s" : "");
    const std::string effect
        = count + " of this seed (observed; the first: " + where + ") " + (bad > 1 ? "are" : "is")
        + " not valid UTF-8: no output can carry such a name verbatim, and a replaced one could "
          "name another label of the index, so the seed is refused and nothing of it is "
          "delivered: " + lever + ", or rename the label in the index";
    Json::Value rj = failed_seed_head(seed, r.labels_from_seed,
        "The seed's labels include " + count + " that " + (bad > 1 ? "are" : "is")
        + " not valid UTF-8 (the first: " + where + "): refused, since no output "
          "carries such a name verbatim and a replaced name could resolve to another "
          "label");
    Json::Value lims(Json::arrayValue);
    Json::Value d = limitation("derivation", knob, std::move(limit),
                               uint_json(static_cast<uint64_t>(bad)), effect);
    d["cause"] = to_string(SeedDerivationError::UNREPRESENTABLE_LABEL_NAME);
    lims.append(std::move(d));
    rj["limitations"] = std::move(lims);
    rj["outcome"] = outcome_of(rj, true);
    // nothing of the walk is delivered, its coordinates neither
    state_no_coordinates(&rj, coordinates_reason);
    return rj;
}

// A seed a request budget does not hold (SeedBudgetError, DESIGN-traverse-graphlet.md §14):
// no valid traversal exists within the budget — not even one complete to 0 bp — so the seed
// is failed, in the shape of a failed derivation (no arms, an error, outcome.walks: failed)
// that every client already handles. The structured reason is in the limitations and the
// resource_stop (DESIGN §14: "failed: no valid traversal exists (structured reason in
// limitations / resource_stop)"): a seed-level walk_domain naming the budget's knob, and
// the resource_stop with what the budget held and needed and the levers on the seed.
// DECISION (the conservative rule of DESIGN §14: a stated limitation is never
// omitted, and no axis is non-complete without the limitation that explains it): the
// walk_domain is stated although it is otherwise an arm's — the walk domain was cut to
// nothing by a cap, which is exactly what a walk_domain states, and a reader that only
// reads limitations must still find the knob (a work-budget failure has no other entry).
// It is the same kind, knob and typed values as an arm's resource_limit walk_domain, so
// no record, field or token is new. The other seeds of the request are traversed;
// nothing per label is delivered.
static Json::Value budget_failed_seed_to_json(const Seed &seed, const SeedBudgetError &e,
                                              const Strategy &st, const Json::Value &clamped,
                                              const char *coordinates_reason) {
    const ResourceStop &q = e.stop();
    // the budget at the server's maximum, which a request cannot raise: its lever is
    // not offered (state_budget_clamp states the rest)
    const Json::Value *clamp = q.injected ? nullptr : budget_clamp(clamped, q.resource);
    const bool at_max = clamp != nullptr;
    Json::Value rj = failed_seed_head(seed, e.labels_from_seed(), e.what());
    Json::Value lims(Json::arrayValue);
    if (q.resource == ResourceStop::MEMORY && is_decode_stop(q)) {
        // A seed-phase or root read did not fit. Observed: what admitting the refused row needed
        // — its standalone demand beside what was held — or, where its read alone was refused,
        // the least it was seen to need (the budget's MiB plus one would tell an agent to raise
        // the budget to a value that fails again)
        const bool root = q.where == ResourceStop::ROOT;
        lims.append(limitation("walk_domain", "bounds.max_memory_mb",
                               st.max_memory_bytes ? uint_json(st.max_memory_bytes >> 20)
                                                   : Json::Value("unlimited"),
                               uint_json(knob_mib(q)),
                               q.injected
                               ? "an injected refusal of an annotation read (a test hook, not the "
                                 "budget) failed the seed before any traversal; no request knob "
                                 "caused it"
                               : std::string("the memory budget did not admit reading ")
                                 + (root ? "an annotate root's annotation row"
                                         : q.where == ResourceStop::DERIVATION
                                         ? "the seed's annotation rows to derive its labels (a "
                                           "window of up to 64 held at once)"
                                         : "the seed's annotation rows to validate its labels")
                                 + " with their row-diff dependency rows, so no traversal was made "
                                   "(observed: what admitting the refused row needed beside what "
                                   "was held, MiB rounded up"
                                 + (q.allotted ? " to the smallest budget that admits it with the "
                                                 "caches' allotments, which grow with the budget"
                                               : "")
                                 + (q.lower_bound ? ", at least: its read alone was refused; "
                                                    "raising the knob to it may still fail: unread "
                                                    "roots, labels or later state may need more"
                                                  : "")
                                 + (at_max ? "); start from a more selective seed"
                                           : "); raise the knob or start from a more selective seed")
                                 + (root && q.label_bytes
                                         && q.demand - q.used
                                                <= static_cast<double>(q.left + q.label_bytes)
                                     ? "; detail graphlet or a lower labels.max_labels_per_node "
                                       "shrinks the other root's state it competes with"
                                     : "")
                                 + (at_max ? kAtServerMaximum : "")));
    } else if (q.resource == ResourceStop::MEMORY) {
        lims.append(limitation("walk_domain", "bounds.max_memory_mb",
                               uint_json(st.max_memory_bytes >> 20),
                               uint_json(knob_mib(q)),
                               std::string("the memory budget does not hold the seed's depth-0 state (its "
                               "label dictionary and both arms' roots, with what ending and "
                               "delivering them costs in the requested detail), so no traversal "
                               "was made (observed: the MiB it needs, rounded up")
                               + (q.allotted ? " to the smallest budget that holds it with the "
                                               "caches' allotments, which grow with the budget"
                                             : "")
                               + (q.lower_bound ? ", at least: a root's row or labels were not "
                                                  "built once the budget was reached; raising the "
                                                  "knob to it may still fail: unread roots, labels "
                                                  "or later state may need more" : "")
                               + (at_max ? "); " : "); raise the knob, ")
                               + "use detail graphlet, or start from fewer labels (fewer permitted; "
                                 "in annotate mode a lower labels.max_labels_per_node)"
                               + (at_max ? kAtServerMaximum : "")));
    } else if (is_external_stop(q)) {
        // stopped from outside in its seed phase: no result exists to deliver, and no budget
        // of the request is to blame (the knob names the attempt)
        lims.append(limitation("walk_domain", "attempt_id",
                               uint_json(static_cast<uint64_t>(q.limit)),
                               uint_json(static_cast<uint64_t>(q.used)),
                               std::string(external_cause(q)) + " while the seed was read "
                               "(validated against its labels, or its permitted set derived), so "
                               "no traversal was made (limit: the attempt's bound, observed: its "
                               "elapsed time, ms); no budget of the request ran out: a new attempt "
                               "can retry the seed"));
    } else {
        lims.append(limitation("walk_domain", "bounds.max_work_units", uint_json(st.max_work_units),
                               uint_json(static_cast<uint64_t>(q.used)),
                               std::string("the work budget ran out while the seed was read "
                               "(validated against its labels, or its permitted set derived: one "
                               "annotation row per seed k-mer), so no traversal was made (observed: "
                               "the work units spent); ")
                               + (at_max ? std::string("shorten the seed") + kAtServerMaximum
                                         : std::string("raise the knob or shorten the seed"))));
    }
    // every response under a memory budget states it (§7.0), with what this result's echo
    // of the request holds beyond the budget (failed_soft)
    if (st.max_memory_bytes) {
        ResourceAccount account = e.account();
        account.soft_overshoot = failed_soft(st, seed, account.soft_overshoot);
        lims.append(memory_bound_soft(st, account));
    }
    rj["limitations"] = std::move(lims);
    rj["resource_stop"] = resource_stop_json(q, st, &e);
    if (clamp)
        state_budget_clamp(&rj, *clamp);
    rj["outcome"] = outcome_of(rj, true);
    state_no_coordinates(&rj, coordinates_reason);
    return rj;
}

// the token a resource takes in resource_stop.resource (and the usage's stopped_by)
static const char* resource_name(ResourceStop::Resource resource) {
    switch (resource) {
        case ResourceStop::MEMORY: return "memory";
        case ResourceStop::WORK: return "work";
        case ResourceStop::TIME: return "time";
        case ResourceStop::CANCELLED: return "cancelled";
        case ResourceStop::ATTEMPT_DEADLINE: return "attempt_deadline";
    }
    return "";
}

// A seed the attempt never started (DESIGN-traverse-graphlet.md §14): the attempt was
// cancelled, or reached the duration bound the server enforces for it, before this seed's walk
// began. Failed per seed in the shape of a failed derivation (no arms, an error,
// outcome.walks: failed) with a seed-level walk_domain naming the attempt and the
// resource_stop, phase not_started, whose action is a new attempt: nothing of the seed was
// read, and no budget of the request ran out. Only free tokens take new values (MGT v1).
static Json::Value not_started_seed_to_json(const Seed &seed, ExternalStop stop, double bound_ms,
                                            double elapsed_ms, const Strategy &st,
                                            bool decode_charged, const char *coordinates_reason) {
    ResourceStop q;
    q.resource = stop == ExternalStop::CANCELLED ? ResourceStop::CANCELLED
                                                 : ResourceStop::ATTEMPT_DEADLINE;
    q.phase = "not_started";
    q.limit = std::ceil(bound_ms);
    q.used = std::ceil(elapsed_ms);
    q.demand = q.used;
    const uint64_t limit = static_cast<uint64_t>(q.limit);
    const uint64_t used = static_cast<uint64_t>(q.used);
    // labels_from_seed as a walked or budget-failed seed of the same request states it: false
    // in annotate mode (seed.labels.empty() would say true for every annotate seed)
    Json::Value rj = failed_seed_head(seed, labels_derived_from_seed(seed, st.label_mode),
        std::string("not started: ") + external_cause(q) + " before this seed began, "
        + std::to_string(used) + " ms after the request was received; no budget of the "
          "request ran out: a new attempt can traverse it");
    Json::Value lims(Json::arrayValue);
    lims.append(limitation("walk_domain", "attempt_id", uint_json(limit), uint_json(used),
                           std::string(external_cause(q)) + " before this seed's walk began, so "
                           "no traversal was made (limit: the attempt's bound, observed: its "
                           "elapsed time, ms); no budget of the request ran out: a new attempt "
                           "can traverse the seed"));
    // every response under a memory budget states it (§7.0): this result echoes seed_id
    append_failed_memory(&lims, st, seed, 0, decode_charged);
    rj["limitations"] = std::move(lims);
    rj["outcome"] = outcome_of(rj, true);
    Json::Value j;
    j["scope"] = "attempt";
    j["resource"] = q.resource == ResourceStop::CANCELLED ? "cancelled" : "attempt_deadline";
    j["phase"] = q.phase;
    j["requested"] = uint_json(limit);
    j["effective"] = uint_json(limit);
    j["used"] = uint_json(used);
    j["remaining"] = uint_json(limit > used ? limit - used : 0);
    j["actions"] = strings_json({ "retry_attempt" });
    j["message"] = std::string(external_cause(q)) + " before this seed's walk began, "
                 + std::to_string(used) + " ms after the request was received (the bound: "
                 + std::to_string(limit) + " ms): the seed is failed, nothing of it was read; no "
                   "budget of the request ran out: a new attempt can traverse it";
    rj["resource_stop"] = std::move(j);
    state_no_coordinates(&rj, coordinates_reason);
    return rj;
}

// A seed with a k-mer the graph does not have, on a server whose graphs need not hold the seeds
// sent to them (TraverseLimits::not_in_graph_per_seed: a multi-graph server's chunks): no walk,
// outcome.walks "not_in_graph", the walker's message as `error` (the graph runs of its k-mers)
// and not_in_graph {kmers, kmers_present} — 0 present: none of it is in this graph; some: only
// part of it, which is not a seed here (§6.1 step 2) — in the shape of a failed seed, without
// a limitation: no knob of the request would get past it. Under a memory budget it states
// memory_bound_soft, as every result does
static Json::Value not_in_graph_seed_to_json(const Seed &seed, const SeedNotInGraph &e,
                                             const Strategy &st, bool decode_charged,
                                             const char *coordinates_reason) {
    Json::Value rj = failed_seed_head(seed, labels_derived_from_seed(seed, st.label_mode),
                                      e.what());
    Json::Value absent;
    absent["kmers"] = uint_json(e.kmers());
    absent["kmers_present"] = uint_json(e.present());
    rj["not_in_graph"] = std::move(absent);
    Json::Value lims(Json::arrayValue);
    append_failed_memory(&lims, st, seed, 0, decode_charged);
    rj["limitations"] = std::move(lims);
    rj["outcome"] = outcome_of(rj, true);
    rj["outcome"]["walks"] = "not_in_graph";
    state_no_coordinates(&rj, coordinates_reason);
    return rj;
}

// A request with attempt_id outside the server (the CLI, a test): its usage is stated as the
// server states it, with the bound computed but not enforced (no lease to protect)
static const std::string& local_instance() {
    static const std::string instance = AttemptRegistry(AttemptSettings()).server_instance();
    return instance;
}

static std::unique_ptr<Attempt> local_attempt(const AttemptIds &ids) {
    AttemptSettings settings;
    settings.hard_cap_ms = 0;
    auto attempt = std::make_unique<Attempt>(0, std::chrono::system_clock::now(), settings,
                                             local_instance(), false);
    attempt->set_ids(ids);
    return attempt;
}

Json::Value process_traverse_request(const Json::Value &json,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const std::string &release,
                                     const TraverseLimits &limits,
                                     const IndexIdentity *identity,
                                     Attempt *attempt,
                                     ResultTexts *texts) {
    TraverseRequest req = parse_traverse_request(json);
    std::unique_ptr<Attempt> local;
    if (!attempt && !req.attempt_id.empty()) {
        local = local_attempt({ req.attempt_id, req.budget_id, req.locus_id, req.not_after_ms, {} });
        attempt = local.get();
    }
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
        note_clamped(&clamped, field, std::move(requested), std::move(effective));
    };
    const bool derives = std::any_of(req.seeds.begin(), req.seeds.end(),
                                     [](const Seed &s) { return s.labels.empty(); });
    // the time budget as requested: a clamp states it beside the effective one (K, Q)
    const double requested_time_ms = req.strategy.time_budget_ms;
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
    // The server's maxima of the budgets (SPEC §10.3): a larger budget is lowered to the
    // maximum and an omitted one (none: unbounded) set to it, echoed with `requested`
    // "unlimited" (how a knob spells no limit, which the graphlet's K and Q records can hold,
    // unlike null), so that no request on such a server runs without the bound the operator
    // chose
    if (limits.max_memory_mb) {
        const uint64_t requested_mb = req.strategy.max_memory_bytes >> 20;
        if (!req.strategy.max_memory_bytes || requested_mb > limits.max_memory_mb) {
            clamp("bounds.max_memory_mb",
                  req.strategy.max_memory_bytes ? uint_json(requested_mb)
                                                : Json::Value("unlimited"),
                  uint_json(limits.max_memory_mb));
            req.strategy.max_memory_bytes = limits.max_memory_mb << 20;
        }
    }
    if (limits.max_work_units) {
        if (!req.strategy.max_work_units || req.strategy.max_work_units > limits.max_work_units) {
            clamp("bounds.max_work_units",
                  req.strategy.max_work_units ? uint_json(req.strategy.max_work_units)
                                              : Json::Value("unlimited"),
                  uint_json(limits.max_work_units));
            req.strategy.max_work_units = limits.max_work_units;
        }
    }
    // the attempt's duration bound: n_seeds x the effective per-seed time budget (after the
    // server's clamp) plus the server's allowance, under the HTTP server's cap
    if (attempt) {
        attempt->set_bound(req.seeds.size(), req.strategy.time_budget_ms,
                           req.strategy.max_memory_bytes, limits.load_ms);
        attempt->set_delivery_detail(req.detail);
    }

    LabelOracle oracle(anno_graph);
    // Record coordinates (opt-in): the reason every seed of the request reports none when the
    // index or the support rules them out, otherwise the one a seed without a traversal states
    // (failed, refused, not started; §18.1); null when not asked for. And what they cost the
    // output: none, the null form, or the block where they can be recorded
    const char *const coordinates_ruled_out
        = coordinates_reason(req.strategy, oracle.has_coordinates());
    const char *const no_coordinates = !req.strategy.coordinates ? nullptr
                                     : coordinates_ruled_out ? coordinates_ruled_out
                                                             : kCoordinatesNoTraversal;
    const CoordinatesOutput coordinates_output
        = !req.strategy.coordinates ? CoordinatesOutput::NONE
        : coordinates_ruled_out ? CoordinatesOutput::REASON : CoordinatesOutput::BLOCK;
    // what the requested output costs per object, for the memory budget (§14); re-priced per
    // seed where its floats can be wider (mgt_float_width)
    req.strategy.delivery = delivery_costs(req.detail, req.strategy.sequences, kMgtFloatWidth,
                                           coordinates_output);
    uint64_t priced_width = kMgtFloatWidth;

    // the chunked deadlines (spec §6.8): every annotation read of this request is decoded in
    // pieces of about this duration under a deadline, the deadline checked between them
    oracle.pacer().target_ms = limits.chunk_target_ms;
    // the row-diff path cache of this request's reads (the walker bounds it per seed under a
    // memory budget): what a read returns does not depend on it, the decoding work does
    oracle.set_path_cache_max(limits.path_cache_bytes);
    if (limits.path_cache_retention) {
        const auto [checkpoint, successors, narrow] = *limits.path_cache_retention;
        oracle.path_cache().set_retention(checkpoint, successors, narrow);
    }
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
    out["algorithm_version"] = kTraverseAlgorithmVersion;

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
        d.path_cache_hits = now.path_cache_hits - before.path_cache_hits;
        d.path_cache_stored_rows = now.path_cache_stored_rows - before.path_cache_stored_rows;
        d.path_cache_rows_kept = now.path_cache_rows_kept - before.path_cache_rows_kept;
        d.path_cache_bytes_kept = now.path_cache_bytes_kept - before.path_cache_bytes_kept;
        d.path_cache_peak_bytes = now.path_cache_peak_bytes;   // the request's so far
        d.seed_phase_seconds = now.seed_phase_seconds - before.seed_phase_seconds;
        d.label_resolve_seconds = now.label_resolve_seconds - before.label_resolve_seconds;
        d.seed_fetch_seconds = now.seed_fetch_seconds - before.seed_fetch_seconds;
        before = now;
        return d;
    };
    // building the response (the JSON tree and the graphlet text), the part a request
    // can make slow without walking more: reported apart from the walk
    double serialize_seconds = 0;
    // The attempt's control of the walk: its stop (a cancel, its bound, a client gone) polled
    // at the walker's checkpoints, its clock and bound for what such a stop states, and the
    // seed's usage (the meter); the results are built under its delivery check
    AttemptControl control;
    if (attempt) {
        control.poll = [attempt]() { return attempt->poll(); };
        control.elapsed_ms = [attempt]() { return attempt->elapsed_ms(); };
        control.bound_ms = attempt->bound_ms();
        // a paced read's chunks: sized by the time left to the walk-until, and preceded by a
        // poll that reads the clock and the client (the stop flag alone is read elsewhere)
        control.ms_left = [attempt]() { return attempt->ms_left(); };
        control.poll_now = [attempt]() { return attempt->poll(/* force */ true); };
        // the walked seed's account, from which its output is estimated (the delivery reserve),
        // with its record coordinates' share apart (opt-in; 0 without them): an occurrence's
        // account is 20-73 times its text (872 bytes for 12-44), the rest of the output's
        // 115-129, so the reserve prices that share at its own bound
        // (kCoordinateAccountPerTextByte)
        control.progress = [attempt](uint64_t account, uint64_t coordinates) {
            attempt->progress(account, coordinates);
        };
    }
    DeliveryScope delivery(attempt);
    // A seed's result as the response holds it: in the tree, or (|texts|, the server) written
    // as compact text now, under the attempt's delivery check, so that its tree is freed at
    // once and the attempt's delivery reserve counts its exact bytes. |built_seconds|: what
    // building its tree took (the build rate measured on large results)
    const std::function<void()> text_check = attempt
        ? std::function<void()>([attempt]() { attempt->check_delivery(); })
        : std::function<void()>();
    // |account|: the walk's final modelled account (0: no walk), whose ratio to the text the
    // server measures for its next attempts' estimates; |coordinates|: its record coordinates'
    // share, which the ratio leaves out with the exact text it wrote (coordinates_text_bytes,
    // counted before the tree is freed)
    auto append = [&](Json::Value &&rj, double built_seconds, uint64_t account = 0,
                      uint64_t coordinates = 0) {
        if (!texts) {
            results.append(std::move(rj));
            return;
        }
        Timer written;
        const uint64_t coordinate_text = attempt && account && req.strategy.coordinates
                                       ? coordinates_text_bytes(rj, text_check) : 0;
        double gap = 0;
        texts->texts.push_back(json_text(rj, true, text_check, attempt ? &gap : nullptr));
        if (attempt)
            attempt->note_delivery_gap_ms(gap);
        rj = Json::Value();
        if (attempt) {
            attempt->note_delivered(texts->texts.back().size(), built_seconds + written.elapsed(),
                                    account, coordinates, coordinate_text);
        }
    };
    // a seed never started reads its annotation as one started would have
    const bool decode_charged = (req.strategy.max_memory_bytes || req.strategy.max_work_units)
                              && oracle.decode_charged();
    for (size_t i = 0; i < req.seeds.size(); ++i) {
        const Seed &seed = req.seeds[i];
        if (attempt) {
            // between seeds: a stop that came after the previous seed's walk (or while it was
            // built) leaves the rest of the seeds unstarted, each failed with the stop
            const ExternalStop stop = attempt->poll(/* force */ true);
            if (stop == ExternalStop::CLIENT_GONE) {
                throw AttemptAborted("the client closed its connection: the request was "
                                     "abandoned before seed " + std::to_string(i) + " of "
                                     + std::to_string(req.seeds.size()));
            }
            if (stop != ExternalStop::NONE) {
                for (size_t j = i; j < req.seeds.size(); ++j) {
                    append(not_started_seed_to_json(req.seeds[j], stop, attempt->bound_ms(),
                                                    attempt->elapsed_ms(), req.strategy,
                                                    decode_charged, no_coordinates), 0);
                    // its usage is what its result states and holds
                    SeedUsage usage;
                    usage.meter.memory_final = failed_result_bytes(req.strategy, req.seeds[j]);
                    usage.meter.soft_excess = failed_soft(req.strategy, req.seeds[j], 0);
                    attempt->seed_not_started(j, usage);
                }
                break;
            }
            attempt->seed_started(i);
        }
        Timer seed_timer;
        AttemptMeter meter;
        control.meter = &meter;
        // What the seed's walk consumed and what stopped it, recorded when its walk is over
        // (before its result is built, which the attempt's bound can still interrupt). A
        // failed or refused seed's result is not the walk's: its usage states what that result
        // holds and what its memory_bound_soft states, the echo of seed_id included.
        // |refused|: the demand a memory budget refused, when one stopped or failed the seed
        auto walked = [&](const std::string &outcome, const std::string &stopped_by,
                          bool failed_result = false,
                          std::optional<uint64_t> refused = std::nullopt) {
            if (!attempt)
                return;
            // the longest read or head piece so far (usage, and the server's deadline_check)
            attempt->note_max_uninterruptible_ms(oracle.pacer().max_uninterruptible_ms());
            SeedUsage usage;
            usage.outcome = outcome;
            usage.stopped_by = stopped_by;
            usage.meter = meter;
            if (failed_result) {
                usage.meter.memory_final = failed_result_bytes(req.strategy, seed);
                usage.meter.soft_excess = failed_soft(req.strategy, seed, meter.soft_excess);
            }
            usage.refused_bytes = refused;
            usage.elapsed_ms = seed_timer.elapsed() * 1000;
            attempt->seed_walked(i, usage);
        };
        auto refused_of = [](const ResourceStop &q) -> std::optional<uint64_t> {
            if (q.resource != ResourceStop::MEMORY || q.injected)
                return std::nullopt;
            return static_cast<uint64_t>(std::ceil(q.demand));
        };
        // the seed's result is appended: its outcome, and its walk and building in its time
        auto delivered = [&](const std::string &outcome) {
            if (attempt)
                attempt->seed_delivered(i, outcome, seed_timer.elapsed() * 1000);
        };
        std::vector<std::string> dict = seed.labels;
        dict.insert(dict.end(), req.strategy.extra.begin(), req.strategy.extra.end());
        LabelChangeCost cost = make_cost(req.cost, dict);
        // the output's floats are priced at the widest the seed's costs and the time budgets
        // can be written (24 characters for every usual request)
        const uint64_t float_width = mgt_float_width(req.strategy, cost, requested_time_ms,
                                                     limits.max_time_ms);
        if (float_width != priced_width) {
            req.strategy.delivery = delivery_costs(req.detail, req.strategy.sequences, float_width,
                                                   coordinates_output);
            priced_width = float_width;
        }
        try {
            SeedResult r = traverse_seed(oracle, seed, req.strategy, cost, release, nullptr,
                                         attempt ? &control : nullptr);
            if (attempt) {
                // provisional until the result is built (outcome_of reads its limitations)
                bool partial = false;
                for (const ArmResult &arm : r.arms) {
                    partial |= arm.requested && arm.status != ArmResult::COMPLETE;
                }
                walked(partial || r.resource_stop ? "partial" : "complete",
                       r.resource_stop ? resource_name(r.resource_stop->resource) : "", false,
                       r.resource_stop ? refused_of(*r.resource_stop) : std::nullopt);
            }
            r.annotation_counters = per_seed(r.annotation_counters);
            // Names come from FASTA headers and file names, which need not be UTF-8. No
            // output carries such a name verbatim, and a replaced one (U+FFFD) can be the
            // name of another label anywhere in the index -- resubmitted, it resolves to
            // that label. The seed is refused instead (in both modes; checked after the
            // walk, since an annotate dictionary is complete only then), like a failed
            // derivation; the writers' own refusals stay as the backstop.
            Json::Value refused = unrepresentable_seed_to_json(seed, r, req.strategy,
                                                               no_coordinates);
            if (!refused.isNull() && r.derivation_partial) {
                // A partial derivation's set is delivered as a walk or not at all (Walker::run):
                // such a name among labels the whole seed may exclude fails the seed with the
                // time budget that cut its derivation (the handler below states it)
                SeedDerivationError e = derivation_out_of_time(
                        r.derivation_partial->kmers_read, r.num_kmers, req.strategy.time_budget_ms,
                        r.derivation_partial->elapsed_ms);
                e.set_soft_overshoot(r.account.soft_overshoot);
                throw e;
            }
            if (!refused.isNull()) {
                if (req.strategy.max_memory_bytes) {
                    // every response under a memory budget states it (§7.0), with what
                    // this result's echo of the request holds beyond it (failed_soft)
                    ResourceAccount account = r.account;
                    account.soft_overshoot = failed_soft(req.strategy, seed, account.soft_overshoot);
                    refused["limitations"].append(memory_bound_soft(req.strategy, account));
                    refused["outcome"] = outcome_of(refused, true);
                }
                // its result is the refusal, not the walk
                walked("failed", r.resource_stop ? resource_name(r.resource_stop->resource) : "",
                       true, r.resource_stop ? refused_of(*r.resource_stop) : std::nullopt);
                const std::string outcome = refused["outcome"]["walks"].asString();
                append(std::move(refused), 0);
                delivered(outcome);
                continue;
            }
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
            const std::string outcome = rj["outcome"]["walks"].asString();
            append(std::move(rj), seconds, meter.memory_final, meter.memory_coordinates);
            delivered(outcome);
        } catch (const SeedDerivationError &e) {
            // The permitted set could not be derived from THIS seed (nothing carries it
            // in full, the budget ran out, the names are ambiguous). The caller named no
            // labels, so it had no lever on that and nothing to fix in the request:
            // report it against the seed and keep the other seeds' traversals, instead of
            // discarding a 100-seed batch because seed 57 spans a recombination point.
            walked("failed", e.cause() == SeedDerivationError::TIME_BUDGET ? "time" : "", true);
            Json::Value failed = failed_seed_to_json(seed, e, req.strategy, clamped,
                                                     oracle.decode_charged(), no_coordinates);
            const std::string outcome = failed["outcome"]["walks"].asString();
            append(std::move(failed), 0);
            delivered(outcome);
            oracle.sync_path_cache_counters();
            per_seed(oracle.counters());   // do not bill this seed's reads to the next
        } catch (const SeedBudgetError &e) {
            // a request budget does not hold this seed (§14): failed per seed, like a
            // derivation, since the budget is per seed and the other seeds may fit (or the
            // attempt stopped it in its seed phase)
            walked("failed", resource_name(e.stop().resource), true, refused_of(e.stop()));
            Json::Value failed = budget_failed_seed_to_json(seed, e, req.strategy, clamped,
                                                            no_coordinates);
            const std::string outcome = failed["outcome"]["walks"].asString();
            append(std::move(failed), 0);
            delivered(outcome);
            oracle.sync_path_cache_counters();
            per_seed(oracle.counters());
        } catch (const AttemptAborted &) {
            // the client is gone: nothing is written, but what the seed consumed until it was
            // abandoned is the attempt's (GET /traverse/attempt states it)
            walked("abandoned", "client_gone");
            throw;
        } catch (const SeedNotInGraph &e) {
            if (!limits.not_in_graph_per_seed) {
                // the only graph of this server (or the CLI's): the whole request fails (400),
                // as every malformed seed fails it
                walked("failed", "");
                throw InvalidRequest(std::string("seed '") + (seed.seed_id.empty() ? seed.sequence.substr(0, 32) : seed.seed_id)
                                     + "': " + e.what());
            }
            // a chunk of a multi-graph server that does not hold the seed: stated per seed, the
            // other seeds traversed
            walked("not_in_graph", "", true);
            Json::Value absent = not_in_graph_seed_to_json(seed, e, req.strategy,
                                                           oracle.decode_charged(),
                                                           no_coordinates);
            append(std::move(absent), 0);
            delivered("not_in_graph");
            oracle.sync_path_cache_counters();
            per_seed(oracle.counters());
        } catch (const std::invalid_argument &e) {
            // the whole request fails (400); what the seed consumed is still the attempt's
            walked("failed", "");
            throw InvalidRequest(std::string("seed '") + (seed.seed_id.empty() ? seed.sequence.substr(0, 32) : seed.seed_id)
                                 + "': " + e.what());
        }
    }
    if (texts) {
        texts->active = true;
    } else {
        out["results"] = std::move(results);
    }
    if (attempt) {
        attempt->walk_ended();
        // every response to a ledger-managed request states what it consumed (DESIGN §14.1)
        if (attempt->managed())
            out["usage"] = attempt->usage_json(attempt->walk_reason());
    }
    if (req.timing) {
        Json::Value t;
        t["elapsed_ms"] = timer.elapsed() * 1000;
        t["serialize_ms"] = serialize_seconds * 1000;
        // a request with in_ram: the time before its work began, waiting for and loading its
        // index into RAM (0: served by the index the server holds)
        if (limits.load_ms)
            t["load_ms"] = *limits.load_ms;
        out["timing"] = std::move(t);
    }
    return out;
}

Json::Value resolve_deadline_body(const ResolveDeadline &e) {
    Json::Value b;
    b["error"] = e.what();
    b["code"] = "deadline";
    return b;
}

Json::Value resolve_capabilities_json(const ResolveTimeLimits &limits) {
    Json::Value t;
    t["accepted"] = true;
    t["knob"] = "bounds.time_budget_ms";
    // no deadline without the field, whatever the cap: such a request has no time limit
    t["default"] = Json::Value();
    // the cap (--traverse-max-time-ms, as /traverse's); 0: none
    t["max_time_ms"] = number_json(limits.max_time_ms);
    t["finalize_reserve_ms"] = number_json(limits.finalize_ms);
    t["finalize_reserve_configurable"] = false;
    t["check_kmers"] = uint_json(kResolveCheckKmers);
    t["check_labels"] = uint_json(kResolveCheckLabels);
    t["stop_phases"] = strings_json({ "rows", "support" });
    // the rule (opt-in, the clamp, where the deadline is read, the stop and the prefix it
    // answers, the reserve and its 503, what is not polled) is the SPEC's
    t["rule"] = kResolveTimeBudgetRule;
    Json::Value r;
    r["time_budget"] = std::move(t);
    return r;
}

// Under a stop, why an explicit selection cannot be made from the prefix (empty: it can, and
// is the prefix's): an interval ending past the k-mers resolved, or, on a discovery, a label
// the prefix did not discover — whether either holds on the whole query is unknown, so the
// selection is not made rather than refused (a 400 would blame the request for the deadline).
// What is known whatever the stop is refused as without one, every seed checked first, in
// select_seeds' order: an interval empty or out of range, an
// interval not fully in the graph (the whole query's k-mers are mapped before the deadline is
// read), and, with explicit labels (every one of them profiled whatever the stop), a seed
// label that is not one of them
static std::string explicit_selection_blocked(const SupportProfile &profile,
                                              const SelectionPolicy &policy, bool discover) {
    for (size_t i = 0; i < policy.explicit_seeds.size(); ++i) {
        const ExplicitSeed &ex = policy.explicit_seeds[i];
        if (ex.kmers.begin >= ex.kmers.end || ex.kmers.end > profile.stop->query_kmers)
            throw InvalidRequest("Explicit seed interval out of range");
        const bool in_one_run = std::any_of(profile.stop->query_graph_runs.begin(),
                                            profile.stop->query_graph_runs.end(),
            [&](const KmerInterval &run) {
                return run.begin <= ex.kmers.begin && run.end >= ex.kmers.end;
            });
        if (!in_one_run)
            throw InvalidRequest("Explicit seed interval is not fully present in the graph");
        if (discover)
            continue;
        for (const std::string &name : ex.labels) {
            const bool profiled = std::any_of(profile.labels.begin(), profile.labels.end(),
                [&](const LabelProfile &lp) { return lp.label.name == name; });
            if (!profiled)
                throw InvalidRequest("Explicit seed label was not profiled: '" + name + "'");
        }
    }
    for (size_t i = 0; i < policy.explicit_seeds.size(); ++i) {
        const ExplicitSeed &ex = policy.explicit_seeds[i];
        const std::string seed = "select.seeds[" + std::to_string(i) + "]";
        if (ex.kmers.end > profile.num_kmers) {
            return seed + ".kmer_interval [" + std::to_string(ex.kmers.begin) + ", "
                    + std::to_string(ex.kmers.end) + ") ends past the k-mers resolved";
        }
        for (const std::string &name : ex.labels) {
            const bool profiled = std::any_of(profile.labels.begin(), profile.labels.end(),
                [&](const LabelProfile &lp) { return lp.label.name == name; });
            if (!profiled)
                return seed + " names a label the resolved prefix did not discover";
        }
    }
    return "";
}

static Json::Value resolve_stop_json(const ResolveStop &stop, size_t k, double budget_ms,
                                      double reserve_ms, const std::string &not_made) {
    Json::Value s;
    s["phase"] = stop.phase == ResolveStop::ROWS ? "rows" : "support";
    s["reason"] = "time";
    s["resolved_kmers"] = uint_json(stop.resolved_kmers);
    s["query_kmers"] = uint_json(stop.query_kmers);
    // the bases the resolved k-mers cover (none for none), and where the rest starts: k-mer i
    // starts at base i, so the sequence from remainder_from_bp has the k-mers not resolved
    s["resolved_bp"] = uint_json(stop.resolved_kmers ? stop.resolved_kmers + k - 1 : 0);
    s["query_bp"] = uint_json(stop.query_kmers + k - 1);
    s["remainder_from_bp"] = uint_json(stop.resolved_kmers);
    std::string message = "bounds.time_budget_ms " + ms_text(budget_ms)
        + " ms less the finalisation reserve of " + ms_text(reserve_ms) + " ms ran out while "
        + (stop.phase == ResolveStop::ROWS ? "the annotation rows of the query's k-mers were read"
                                           : "the explicit labels' hits were read")
        + ": this answer is exactly the resolve of the sequence's first "
        + std::to_string(stop.resolved_kmers) + " of " + std::to_string(stop.query_kmers)
        + " k-mers" + (stop.resolved_kmers ? " (a run ending there may continue past it)" : "")
        + "; resolve the rest from base " + std::to_string(stop.resolved_kmers)
        + ", or raise bounds.time_budget_ms";
    if (!not_made.empty())
        message += "; the explicit selection was not made (selection null): " + not_made;
    s["message"] = std::move(message);
    return s;
}

Json::Value process_resolve_request(
        const Json::Value &json,
        const graph::AnnotatedDBG &anno_graph,
        const std::string &release,
        uint64_t max_query_bp,
        const IndexIdentity *identity,
        const std::function<bool()> &client_gone,
        const ResolveTimeLimits &time,
        ResolveDelivery *delivery,
        const std::function<std::chrono::steady_clock::time_point()> &clock) {
    using graph::pattern::Deadline;
    // the deadline of a request with bounds.time_budget_ms starts here, its body parsed (as
    // /pattern's, DESIGN-pattern-search.md §5.3), before its fields are read; a request without
    // the field reads no clock for it
    const std::function<Deadline::Clock::time_point()> now
            = clock ? clock : std::function<Deadline::Clock::time_point()>(&Deadline::Clock::now);
    std::optional<Deadline::Clock::time_point> start;
    if (json.isObject() && json.isMember("bounds") && json["bounds"].isObject()
            && json["bounds"].isMember("time_budget_ms")) {
        start = now();
    }
    ResolveRequest req = parse_resolve_request(json);
    // a client that is gone is not answered: the request is abandoned at the next phase
    auto abandon = []() {
        throw AttemptAborted("the client closed its connection: the resolve request was abandoned");
    };
    if (client_gone) {
        req.options.stop = client_gone;
        req.options.abandon = abandon;
    }
    if (max_query_bp && req.sequence.size() > max_query_bp)
        throw InvalidRequest("request.sequence: longer than the server limit of "
                             + std::to_string(max_query_bp) + " bp");
    if (req.sequence.size() > req.max_query_bp)
        throw InvalidRequest("request.sequence: longer than bounds.max_query_bp");

    // the deadline, only when asked for: without it no deadline is set or read (the clock of
    // timing.elapsed_ms is read as it always was), and the answer is the one this route
    // always gave
    std::optional<Deadline> deadline;
    Json::Value clamped(Json::arrayValue);
    std::string late;
    if (req.time_budget_ms) {
        double budget = *req.time_budget_ms;
        // the reserve is inside the budget: a budget not above it leaves no time to work
        if (!(budget > time.finalize_ms)) {
            throw InvalidRequest("request.bounds.time_budget_ms: expected more than the "
                                 "finalisation reserve of " + ms_text(time.finalize_ms) + " ms");
        }
        if (time.max_time_ms > 0 && budget > time.max_time_ms) {
            // lowered to the server's cap and stated, as /traverse states its clamps
            note_clamped(&clamped, "bounds.time_budget_ms", req.time_budget_given,
                         number_json(time.max_time_ms));
            budget = time.max_time_ms;
        }
        // (the parse admits the field only where the peek above saw it)
        if (!start)
            start = now();
        deadline.emplace(*start, budget, time.finalize_ms, now);
        late = "resolve: the answer could not be built and written within "
               "bounds.time_budget_ms (" + ms_text(budget) + " ms, the finalisation reserve of "
               + ms_text(time.finalize_ms) + " ms included): nothing partial is sent";
        const Deadline *d = &*deadline;
        req.options.time_up = [d]() { return d->work_expired(); };
        req.options.finish_check = [d, &late]() {
            if (d->respond_expired())
                throw ResolveDeadline(late);
        };
        // the writing and the compression of the answer, after this function returned: the
        // check holds its own copy of the deadline
        if (delivery) {
            delivery->set_check([copy = *deadline, late]() {
                if (copy.respond_expired())
                    throw ResolveDeadline(late);
            });
        }
    }

    LabelOracle oracle(anno_graph);
    Timer timer;
    SupportProfile profile;
    std::optional<SeedSelection> selection;
    std::string not_made;
    try {
        profile = resolve_support(oracle, req.sequence, req.options);
        if (client_gone && client_gone())
            abandon();
        if (req.select) {
            req.policy.release_id = release;
            // under a stop the selection is the prefix's, made from it unless an
            // explicit seed reaches beyond what the prefix can tell
            if (profile.stop && req.policy.policy == SelectionPolicy::EXPLICIT)
                not_made = explicit_selection_blocked(profile, req.policy,
                                                      req.options.discover);
            if (not_made.empty())
                selection = select_seeds(profile, req.sequence, req.policy, oracle.regime() != Regime::BASIC);
        }
    } catch (const std::invalid_argument &e) {
        throw InvalidRequest(e.what());
    }
    // the selection is not polled: read once it is done
    if (req.options.finish_check)
        req.options.finish_check();
    Json::Value out = profile_to_json(profile, selection ? &*selection : nullptr, req.run_format,
                                      req.options.finish_check);
    if (!not_made.empty())
        out["selection"] = Json::Value();
    out["release"] = release;
    out["capabilities"] = capabilities_to_json(oracle, release, identity);
    if (deadline) {
        // the budget the work ran under, stated whether or not it stopped the work
        Json::Value l;
        l["time_budget_ms"] = number_json(deadline->time_budget_ms());
        l["finalize_reserve_ms"] = number_json(deadline->finalize_reserve_ms());
        l["clamped"] = std::move(clamped);
        out["limits"] = std::move(l);
        out["stop"] = profile.stop
            ? resolve_stop_json(*profile.stop, profile.k, deadline->time_budget_ms(),
                                deadline->finalize_reserve_ms(), not_made)
            : Json::Value();
    }
    Json::Value t;
    t["elapsed_ms"] = timer.elapsed() * 1000;
    out["timing"] = std::move(t);
    if (req.options.finish_check)
        req.options.finish_check();
    return out;
}


// ---------------------------------------------------------------- CLI

int traverse_graph(Config *config) {
    assert(config);
    assert(config->infbase_annotators.size() == 1);

    if (config->traverse_index_inventory) {
        // what a manifest of this index must cover (index_manifest.py cross-checks its own
        // list against it): nothing is loaded
        Json::StreamWriterBuilder builder;
        builder["indentation"] = config->output_json ? "" : "  ";
        std::cout << Json::writeString(builder, index_inventory_json(
                             config->infbase, config->infbase_annotators[0],
                             !config->no_coord_mapping)) << std::endl;
        return 0;
    }

    auto loaded = load_graph_with_async_annotation(*config);
    auto graph = loaded.first.get();
    auto anno_graph = loaded.second.get();
    assert(anno_graph);
    graph.reset();

    // the reverse index of the sequence headers, built once before the requests (as the
    // server builds it while the index loads): a first header name does not pay it inside
    // its seed's time budget
    if (const auto *coord_to_header = anno_graph->get_coord_to_header())
        coord_to_header->build_header_index();

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
        // A request with attempt_id states its usage on an error after it was read too, as the
        // server's 400 does; a malformed id is refused
        // without (nothing to reconcile)
        std::unique_ptr<Attempt> attempt;
        try {
            if (!config->traverse_resolve) {
                const AttemptIds ids = attempt_ids(json);
                // a request meant for another server_instance is not started, as the server
                // refuses it (409): this process's instance is its own, random, so any one named
                // is another's
                if (instance_mismatch(ids, local_instance())) {
                    const Json::Value body = instance_mismatch_json(ids, local_instance());
                    logger->error("Request in {} not started: {}", file, body["error"].asString());
                    std::cout << Json::writeString(builder, body) << std::endl;
                    status = 1;
                    continue;
                }
                // a request whose not_after_ms has passed is not started, as the server
                // refuses it (409): its body, and the exit status of a request error
                if (const uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now().time_since_epoch()).count();
                        not_after_passed(ids, now)) {
                    const Json::Value body = expired_json(ids, now, local_instance());
                    logger->error("Request in {} not started: {}", file, body["error"].asString());
                    std::cout << Json::writeString(builder, body) << std::endl;
                    status = 1;
                    continue;
                }
                if (!ids.attempt_id.empty())
                    attempt = local_attempt(ids);
            }
            // operator-run: no caps; the reads are chunked under a deadline as on the server
            TraverseLimits limits;
            limits.chunk_target_ms = static_cast<double>(config->traverse_chunk_target_ms);
            limits.path_cache_bytes = config->traverse_path_cache_mb << 20;
            // a /resolve's deadline (bounds.time_budget_ms, no cap here) holds for the
            // writing of its answer too, as on the server (else its 503 body, exit 1)
            ResolveDelivery delivery;
            Json::Value out = config->traverse_resolve
                ? process_resolve_request(json, *anno_graph, config->index_release, 0, &identity,
                                          nullptr, {}, &delivery)
                : process_traverse_request(json, *anno_graph, config->index_release, limits,
                                           &identity, attempt.get());
            const std::string text = Json::writeString(builder, out);
            delivery.check();
            std::cout << text << std::endl;
        } catch (const ResolveDeadline &e) {
            logger->error("Request in {} not answered: {}", file, e.what());
            std::cout << Json::writeString(builder, resolve_deadline_body(e)) << std::endl;
            status = 1;
        } catch (const InvalidRequest &e) {
            logger->error("Invalid request in {}: {}", file, e.what());
            Json::Value err;
            err["error"] = e.what();
            if (attempt)
                err["usage"] = attempt->usage_json("error");
            std::cout << Json::writeString(builder, err) << std::endl;
            status = 1;
        }
    }
    return status;
}

} // namespace cli
} // namespace mtg
