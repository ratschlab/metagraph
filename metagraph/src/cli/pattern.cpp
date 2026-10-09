#include "pattern.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <sstream>
#include <string_view>
#include <vector>

#include "common/logger.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "config/config.hpp"
#include "json_helpers.hpp"
#include "load/load_annotated_graph.hpp"
#include "pattern_predicate.hpp"
#include "traverse.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph::pattern;
using graph::AnnotatedDBG;
using graph::CanonicalDBG;
using graph::DBGSuccinct;
using graph::DeBruijnGraph;

namespace {

/**
 * The projection of an omitted output.labels: "none" (the design's default is "all", §7.1), so
 * that a request that names no projection reads no annotation. It is stated in every answer
 * (`output`) and in the capabilities (`default_projection`), so a client never has to assume
 * it.
 */
constexpr const char *kDefaultProjection = "none";

// Request fields that are not served (§7.1): refused by name, any value (null included), rather
// than reported as unknown, so that the answer says what to wait for
const char *const kLaterIncrementFields[] = {
    "graphs", "budget_split",
};

// long_search: "anchors", the default, answers a pattern longer than k by its anchors; "paths"
// extends them into paths. A pattern of at most k bases is answered alike under both
constexpr const char kLongSearchAnchors[] = "anchors";
constexpr const char kLongSearchPaths[] = "paths";
// require_support: every label carrying a path is returned with its support
// ("label_intersection", the default), or only the record-verified ones
constexpr const char kSupportIntersection[] = "label_intersection";
constexpr const char kSupportVerified[] = "record_verified";

// The residues a protein pattern may hold (DESIGN §6): the 20 amino acids, the ambiguity codes
// X, B, Z, J and the stop '*' (a stop codon of the request's genetic code), as the engine's
// Pattern::parse reads them; listed in the capabilities
constexpr const char kProteinResidues[] = "ACDEFGHIKLMNPQRSTVWYXBZJ*";

// The two prose fields of the capabilities, references to the SPEC (pattern_capabilities_json):
// §4.5 classifies every cap (a request field's maximum or the server's policy), §12.2 states
// how a peptide is read
constexpr const char kCapsRule[] = "SPEC-pattern-search.md section 4.5";
constexpr const char kProteinRule[] = "SPEC-pattern-search.md section 12.2";

// The route of the full pattern block, which the block of GET /traverse/capabilities names in
// `details` (SPEC §23)
constexpr const char kPatternDetails[] = "GET /pattern/capabilities";

// The keys of the pattern block a client of GET /traverse/capabilities alone gates on or parses
// (SPEC §23, the gate block; what the search service reads, its interface inventory of
// 2026-10-08 §3.10): the contract version and availability, the values a request is checked
// against (modes, projections, kinds, residues, genetic codes, strands, scopes, k, the long
// patterns' search), the budget floor and the defaults its parser reads, what the annotation
// gives (support, placement, annotation), the mask and counting. Each keeps the full block's
// value on that route, whatever else the route leaves to GET /pattern/capabilities
const char *const kPatternGateKeys[] = {
    "pattern_contract_version", "available", "unavailable_reason",
    "modes", "projections", "kinds", "protein_residues", "genetic_codes", "strands",
    "scopes", "scopes_by_graph_mode", "graph_mode", "k",
    "long_patterns", "long_search", "default_long_search",
    "finalize_reserve_ms", "default_time_budget_ms", "default_genetic_code",
    "default_occurrences",
    "support", "placement", "annotation", "mask", "counting",
    "caps",
};
// the caps of the gate block: the ceilings a request is checked against and the chunk size
// (max_patterns), and the two its parser reads (min_information_bits, max_anchors)
const char *const kPatternGateCaps[] = {
    "max_patterns", "max_contexts", "max_anchors", "max_paths", "max_steps", "time_budget_ms",
    "min_information_bits", "max_memory_mb", "max_labels_per_anchor",
};

// The note of an entry answered on a graph without its dummy-edge mask (counting "upper_bound")
// where a count carries an estimate: each such count is the bounds [lower, upper], upper the
// BOSS entries of its ranges (source dummies included), and its estimate is upper x the graph's
// sampled dummy fraction (index.dummy_fraction), not a bound
constexpr const char kNoteEstimate[] = "estimate_sampled_dummy_fraction";

// counting (capabilities and, on a graph without its mask, the answer's index): "exact" with
// the dummy-edge mask, "upper_bound" without it
constexpr const char kCountingExact[] = "exact";
constexpr const char kCountingUpperBound[] = "upper_bound";

// The note of an entry whose request named the labels (output.labels "all", or an
// annotation field) but whose answer reads none: mode count reads no annotation (§5.2), and
// neither does output.labels "none"; the fields were checked and had no effect
constexpr const char kNoteAnnotationNotRead[] = "annotation_not_read";

// The JSON objects built between two readings of the answer's deadline (§5.3:
// "serialisation every 4,096 objects")
constexpr uint64_t kDeliveryStride = 4096;

// Predicates (SPEC §19): output.labels of a predicate request that returns the predicate's
// labels on each selected context; predicate_strands' values ("either", the default, reads a
// context's k-mer and its reverse complement on a BASIC graph, never mixing the two
// orientations: a label is present when it annotates the one k-mer or the other); the scope of
// a predicate's claim (§19.5: per context, per index); and the notes of a predicate entry, in
// this order after the others (§19.10)
constexpr const char kLabelsPredicateOnly[] = "predicate_only";
constexpr const char kPredicateStrandsEither[] = "either";
constexpr const char kPredicateStrandsContext[] = "context";
constexpr const char kPredicateScope[] = "shard_context";
constexpr const char kNotePredicateConstant[] = "predicate_constant";
constexpr const char kNoteProjectionNotRead[] = "projection_not_read";
// the operators of the predicate language, in the order the capabilities list them
const char *const kPredicateOperators[] = {
    "any", "all", "none", "at_least", "and", "or", "not",
};

// The deepest nesting of arrays and objects a request body may have (jsoncpp's stackLimit,
// its default): a deeper one is refused invalid_request
constexpr int kMaxBodyNesting = 1000;

constexpr InvalidPatternRequest invalid{};

using Fields = StrictObject<InvalidPatternRequest>;

PatternRefusal later(const std::string &message) {
    return PatternRefusal(400, "later_increment", message);
}

// One pattern of the request: parsed, or refused in its slot (bad_alphabet)
struct PatternSpec {
    Json::Value id;                      // the string given, or null
    PatternKind kind = PatternKind::DNA;
    // a peptide's text, parsed once the request's genetic code is known
    std::string protein;
    std::optional<Pattern> pattern;
    std::optional<std::pair<std::string, std::string>> error;  // code, message
};

struct ParsedRequest {
    std::vector<PatternSpec> patterns;
    Request request;
    uint64_t max_steps = 0;
    double time_budget_ms = 0;
    // {field, requested, effective} per request value lowered to its cap
    Json::Value clamped = Json::Value(Json::arrayValue);
    // output.labels "all"
    bool labels_all = false;
    // the request named the labels: output.labels "all" or an annotation field
    bool annotation_named = false;
    // the effective annotation limits (used with labels_all in a retrieval mode)
    RetrievalLimits retrieval;
    // long_search "paths": patterns longer than k are extended into paths
    // (Request::extend_paths); false for "anchors" or when omitted
    bool long_paths = false;
    // require_support "record_verified": the labels of paths that one record verifies only
    bool require_verified = false;
    // the genetic code of the request's peptides: genetic_code, default 1
    const GeneticCode *genetic_code = &GeneticCode::standard();
    // the request's predicate as parsed (SPEC §19; bound to the index once, before the first
    // pattern); none without one
    std::optional<predicate::Predicate> predicate;
    // output.labels "predicate_only": the predicate's labels on each selected context
    bool labels_predicate_only = false;
    // a field only a projection reads was named: output.labels "all" or "predicate_only",
    // max_labels_per_anchor, max_annotation_work, max_labels, max_occurrences_per_label,
    // require_support (a predicate answer that builds no projection says so,
    // projection_not_read; max_memory_mb and allow_unbudgeted_annotation act on the selection
    // too)
    bool projection_named = false;
    // max_predicate_contexts (effective): the raw contexts a pattern's selection may test
    uint64_t max_predicate_contexts = 0;
    // max_predicate_work (effective) and predicate_strands
    SelectionLimits selection;
};

// A cap of the request (max_contexts, max_anchors, max_steps): the server's cap when omitted,
// else the request's integer >= |min|, lowered to the cap and listed when above it
uint64_t capped_integer(Fields &f, const char *key, uint64_t cap, uint64_t min,
                        Json::Value *clamped) {
    if (!f.has(key))
        return cap;
    const Json::Value &v = f.raw(key);
    if (!v.isIntegral() || (v.isInt64() && v.asInt64() < 0))
        throw invalid(f.path(key) + ": expected a non-negative integer");
    const uint64_t x = v.asUInt64();
    if (x < min)
        throw invalid(f.path(key) + ": expected an integer >= " + std::to_string(min));
    if (x > cap) {
        note_clamped(clamped, key, uint_json(x), uint_json(cap));
        return cap;
    }
    return x;
}

/**
 * The request as the route serves it (§7.1; SPEC §5): the first error wins, in this order — a
 * field that is not served or resident-only (named, whatever its value), the patterns, mode,
 * output, scope, strands, stop_at_threshold, the caps, the time budget, the annotation caps,
 * long_search, max_paths, require_support, genetic_code (then the peptides are read in it: a
 * slot error, never a refusal), the predicate (its form, then predicate_too_large),
 * max_predicate_contexts, max_predicate_work, predicate_strands, then "predicate_only" without
 * a predicate and a predicate with long_search "paths" — and a field nothing read is refused
 * last, as /traverse's Strict refuses it.
 */
ParsedRequest parse_request(const Json::Value &json, const PatternLimits &limits) {
    if (!json.isObject())
        throw invalid("request: expected an object");

    for (const std::string &name : json.getMemberNames()) {
        if (name == "in_ram") {
            throw PatternRefusal(400, "resident_only",
                                 "request.in_ram: /pattern serves the resident index only and "
                                 "never loads one inside a request (a load cannot be stopped at "
                                 "the deadline)");
        }
        for (const char *field : kLaterIncrementFields) {
            if (name == field)
                throw later("request." + name + ": '" + name + "' is not served by this build");
        }
    }

    Fields f(json, "request");
    ParsedRequest req;

    if (!f.has("patterns"))
        throw invalid("request.patterns: required: a list of 1 to "
                      + std::to_string(limits.max_patterns) + " pattern objects");
    const Json::Value &patterns = f.raw("patterns");
    if (!patterns.isArray() || patterns.empty() || patterns.size() > limits.max_patterns) {
        throw invalid("request.patterns: expected a list of 1 to "
                      + std::to_string(limits.max_patterns) + " pattern objects"
                      + (patterns.isArray() && patterns.size() > limits.max_patterns
                            ? " (got " + std::to_string(patterns.size())
                                + "; the server's --pattern-max-patterns)"
                            : std::string()));
    }
    for (Json::ArrayIndex i = 0; i < patterns.size(); ++i) {
        const std::string path = "request.patterns[" + std::to_string(i) + "]";
        Fields p(patterns[i], path);
        PatternSpec spec;
        if (p.has("id"))
            spec.id = p.str("id", "");
        // the kinds (§4.2): dna, iupac and protein (a peptide)
        const bool dna = p.has("dna");
        const bool iupac = p.has("iupac");
        const bool protein = p.has("protein");
        if (dna + iupac + protein != 1) {
            // (a request that names no peptide is not told about 'protein')
            throw invalid(path + (protein ? ": expected exactly one of 'dna', 'iupac', 'protein'"
                                          : ": expected exactly one of 'dna', 'iupac'"));
        }
        const std::string text = p.str(dna ? "dna" : iupac ? "iupac" : "protein", "");
        p.finish();
        spec.kind = dna ? PatternKind::DNA : iupac ? PatternKind::IUPAC : PatternKind::PROTEIN;
        if (protein) {
            // parsed with the request's genetic code, once that is read (below)
            spec.protein = text;
        } else {
            try {
                spec.pattern = Pattern::parse(spec.kind, text);
            } catch (const PatternError &e) {
                // the pattern's own error (§7.2): the other patterns are still answered
                spec.error = std::make_pair(e.code(), std::string(e.what()));
            }
        }
        req.patterns.push_back(std::move(spec));
    }

    const std::string mode = f.str("mode", to_string(Mode::ALL_OR_COUNT));
    if (auto m = parse_mode(mode)) {
        req.request.mode = *m;
    } else {
        throw invalid(f.path("mode") + ": expected one of count|all_or_count|partial");
    }

    if (f.has("output")) {
        Fields o(f.raw("output"), f.path("output"));
        if (o.has("labels")) {
            const std::string labels = o.str("labels", "");
            if (labels != "none" && labels != "all" && labels != kLabelsPredicateOnly)
                throw invalid(o.path("labels") + ": expected one of none|all|predicate_only");
            req.labels_all = labels == "all";
            // the predicate's labels on each selected context (§19.10); a request without a
            // predicate is refused below, once the predicate is read
            req.labels_predicate_only = labels == kLabelsPredicateOnly;
            req.annotation_named |= req.labels_all || req.labels_predicate_only;
            req.projection_named |= req.labels_all || req.labels_predicate_only;
        }
        if (o.has("occurrences")) {
            const bool occurrences = o.boolean("occurrences", false);
            // placed occurrences are the labels' (§4.3): without labels there is nothing
            // to place (the predicate's labels, with "predicate_only", are placed alike)
            if (occurrences && !req.labels_all && !req.labels_predicate_only) {
                throw invalid(o.path("occurrences") + ": true needs output.labels \"all\""
                              + std::string(json.isMember("predicate")
                                                ? " or \"predicate_only\"" : "")
                              + " (occurrences are placed per label)");
            }
            req.retrieval.occurrences = occurrences;
        }
        // accepted with either value and changes nothing, since a path result always carries
        // its node path (nodes, rows), as a context its node and row
        o.boolean("paths", false);
        o.finish();
    }

    const std::string scope = f.str("scope", to_string(Scope::ANY_OFFSET));
    if (auto s = parse_scope(scope)) {
        req.request.scope = *s;
    } else {
        throw invalid(f.path("scope") + ": expected one of suffix|any_offset (a pattern longer "
                      "than k is searched as 'long' whatever is named)");
    }
    const std::string strands = f.str("strands", to_string(Strands::BOTH));
    if (auto s = parse_strands(strands)) {
        req.request.strands = *s;
    } else {
        throw invalid(f.path("strands") + ": expected one of both|forward|reverse");
    }
    req.request.stop_at_threshold = f.boolean("stop_at_threshold",
                                              req.request.stop_at_threshold);

    req.request.max_contexts = capped_integer(f, "max_contexts", limits.max_contexts, 0,
                                              &req.clamped);
    req.request.max_anchors = capped_integer(f, "max_anchors", limits.max_anchors, 0,
                                             &req.clamped);
    req.max_steps = capped_integer(f, "max_steps", limits.max_steps, 1, &req.clamped);
    req.request.min_information_bits = limits.min_information_bits;
    // the server's policy, read on a graph without its mask only
    req.request.max_checked_entries = limits.max_checked_entries;

    req.time_budget_ms = limits.default_time_ms;
    if (f.has("time_budget_ms")) {
        const Json::Value &v = f.raw("time_budget_ms");
        if (!v.isNumeric())
            throw invalid(f.path("time_budget_ms") + ": expected a number (ms)");
        const double t = v.asDouble();
        // the reserve is inside the budget: a budget not above it leaves no time to work
        if (!(t > limits.finalize_ms)) {
            throw invalid(f.path("time_budget_ms") + ": expected more than the finalisation "
                          "reserve of " + ms_text(limits.finalize_ms) + " ms");
        }
        if (t > limits.max_time_ms) {
            // lowered to the cap, not to the default, which lies below it
            note_clamped(&req.clamped, "time_budget_ms", v, number_json(limits.max_time_ms));
            req.time_budget_ms = limits.max_time_ms;
        } else {
            req.time_budget_ms = t;
        }
    }

    // the annotation caps (§5.3): accepted with any projection and mode, used by output.labels
    // "all" in a retrieval mode (elsewhere stated: annotation_not_read)
    for (const char *key : { "max_labels_per_anchor", "max_annotation_work", "max_memory_mb",
                             "max_labels", "max_occurrences_per_label",
                             "allow_unbudgeted_annotation" }) {
        req.annotation_named |= f.has(key);
    }
    // (with a predicate the account and its access act on the selection in every mode: only
    // these four are a projection's alone)
    for (const char *key : { "max_labels_per_anchor", "max_annotation_work", "max_labels",
                             "max_occurrences_per_label" }) {
        req.projection_named |= f.has(key);
    }
    RetrievalLimits &r = req.retrieval;
    r.max_labels_per_anchor = capped_integer(f, "max_labels_per_anchor",
                                             limits.max_labels_per_anchor, 1, &req.clamped);
    r.max_annotation_work = capped_integer(f, "max_annotation_work", limits.max_annotation_work,
                                           1, &req.clamped);
    r.max_memory_bytes = capped_integer(f, "max_memory_mb", limits.max_memory_mb, 1,
                                        &req.clamped) << 20;
    r.max_labels = capped_integer(f, "max_labels", limits.max_labels, 0, &req.clamped);
    r.max_occurrences_per_label = capped_integer(f, "max_occurrences_per_label",
                                                 limits.max_occurrences_per_label, 0,
                                                 &req.clamped);
    r.allow_unbudgeted = f.boolean("allow_unbudgeted_annotation", r.allow_unbudgeted);
    r.chunk_target_ms = limits.chunk_target_ms;

    // the paths of a pattern longer than k, opt-in
    const std::string long_search = f.str("long_search", kLongSearchAnchors);
    if (long_search != kLongSearchAnchors && long_search != kLongSearchPaths)
        throw invalid(f.path("long_search") + ": expected one of anchors|paths");
    req.long_paths = long_search == kLongSearchPaths;
    req.request.extend_paths = req.long_paths;
    // accepted with any request; it bounds the paths of long_search "paths" only
    req.request.max_paths = capped_integer(f, "max_paths", limits.max_paths, 0, &req.clamped);
    // an annotation field (its effect is on the labels of paths): named in a request that
    // reads no labels, it is stated as not read (annotation_not_read)
    req.annotation_named |= f.has("require_support");
    req.projection_named |= f.has("require_support");
    const std::string support = f.str("require_support", kSupportIntersection);
    if (support != kSupportIntersection && support != kSupportVerified) {
        throw invalid(f.path("require_support") + ": expected one of label_intersection|"
                      "record_verified");
    }
    req.require_verified = support == kSupportVerified;
    if (req.require_verified && !r.occurrences) {
        // the verification reads the labels' coordinates, which output.occurrences false
        // declines: the two contradict each other
        throw invalid(f.path("require_support") + ": \"record_verified\" needs the labels' "
                      "coordinates, which output.occurrences false does not read");
    }

    // the genetic code of the request's peptides, an NCBI translation table id; accepted with
    // any request, it acts on protein patterns only
    if (f.has("genetic_code")) {
        const Json::Value &v = f.raw("genetic_code");
        if (!v.isIntegral()) {
            throw invalid(f.path("genetic_code") + ": expected an integer (an NCBI translation "
                          "table id: " + GeneticCode::ids_text() + ")");
        }
        const GeneticCode *code = v.isInt() ? GeneticCode::find(v.asInt()) : nullptr;
        if (!code) {
            const std::string given = v.isInt64() ? std::to_string(v.asInt64())
                                                  : std::to_string(v.asUInt64());
            throw PatternRefusal(400, "genetic_code_unknown",
                                 f.path("genetic_code") + ": " + given + " is not an NCBI "
                                 "translation table id; the tables served are "
                                 + GeneticCode::ids_text() + " (capabilities genetic_codes), "
                                 "1, the standard code, the default");
        }
        req.genetic_code = code;
    }
    for (PatternSpec &spec : req.patterns) {
        if (spec.kind != PatternKind::PROTEIN)
            continue;
        try {
            spec.pattern = Pattern::parse(PatternKind::PROTEIN, spec.protein, *req.genetic_code);
        } catch (const PatternError &e) {
            // bad_alphabet: the slot's error (§8.9); the stop '*' is a residue
            spec.error = std::make_pair(e.code(), std::string(e.what()));
        }
    }

    // the predicate (SPEC §19.2, §19.3; its form, then its size: 400 predicate_too_large above
    // the server's max_predicate_labels), then the selection's two caps, accepted with any
    // request and acting only with a predicate (lowered and listed after max_paths), then
    // predicate_strands
    if (f.has("predicate")) {
        req.predicate = predicate::Predicate::parse(f.raw("predicate"),
                                                    limits.max_predicate_labels,
                                                    f.path("predicate"));
    }
    req.max_predicate_contexts = capped_integer(f, "max_predicate_contexts",
                                                limits.max_predicate_contexts, 0, &req.clamped);
    req.selection.max_predicate_work = capped_integer(f, "max_predicate_work",
                                                      limits.max_predicate_work, 1,
                                                      &req.clamped);
    const std::string predicate_strands = f.str("predicate_strands",
                                                       kPredicateStrandsEither);
    if (predicate_strands != kPredicateStrandsEither
            && predicate_strands != kPredicateStrandsContext) {
        throw invalid(f.path("predicate_strands") + ": expected one of either|context");
    }
    req.selection.either = predicate_strands == kPredicateStrandsEither;
    // the combinations: the predicate's labels need a predicate; a predicate selects among
    // supported paths for L > k, which long_search "paths" (every graph walk) does not search
    if (req.labels_predicate_only && !req.predicate) {
        throw invalid("request.output.labels: \"predicate_only\" returns the labels a "
                      "predicate names: it needs a predicate (request.predicate)");
    }
    if (req.predicate && req.long_paths) {
        throw invalid("request.long_search: \"paths\" with a predicate: a predicate selects "
                      "among the supported paths of a pattern longer than k (long_search "
                      "\"supported_paths\", not served by this build), never among every graph "
                      "walk; send long_search \"anchors\" (a pattern longer than k is then "
                      "answered by its anchors, its selection not_started)");
    }

    f.finish();
    return req;
}

std::string support_message(const GraphSupport &support) {
    // (no reason is mask_required: a graph without its mask is served, with upper bounds)
    if (support.reason == "representation_unsupported") {
        return "pattern: the graph is not a succinct graph (a DBGSuccinct, or a PRIMARY one "
               "wrapped in CanonicalDBG): the pattern lookup narrows BOSS ranges";
    }
    if (support.reason == "primary_unwrapped") {
        return "pattern: a PRIMARY graph is served only wrapped in CanonicalDBG";
    }
    if (support.reason == "alphabet_unsupported") {
        return "pattern: the graph's alphabet '" + support.alphabet + "' is not $ACGT (a DNA4 "
               "build; $ACGTN, a DNA5 build, is not served yet either)";
    }
    if (support.reason == "alphabet_untested") {
        return "pattern: the graph's alphabet is $ACGTN (a DNA5 build): the pattern search is "
               "not served on the $ACGTN alphabet until a DNA5 build passes the pattern tests "
               "(its DNA5 code has not been run in a test yet); serve a $ACGT (DNA4) graph";
    }
    if (support.reason == "mask_invalid") {
        return "pattern: the graph's dummy-edge mask (.edgemask) marks dummy edges with W = $ "
               "valid (as an older build's `metagraph extend` writes it on a masked graph, or a "
               "stale mask): counts on it would take dummy edges for k-mers and could be claimed "
               "exact while too large; rebuild the mask with `metagraph transform --mask-dummy "
               "--force <graph>.dbg` (writes the .edgemask beside the graph; node ids and "
               "annotation unchanged) and restart the server";
    }
    return "pattern: the graph is not supported (" + support.reason + ")";
}

/**
 * The estimate of a count of a graph served without its dummy-edge mask: its upper bound times
 * the graph's dummy fraction, rounded, and kept inside the bounds (lower, a true lower bound,
 * can exceed the product). Not a bound: what the count would be if the source dummies among the
 * upper bound's entries were as frequent as among the graph's.
 */
uint64_t estimate(const Count &c, const DummyFraction &fraction) {
    assert(c.relation == Relation::BOUNDS);
    const double x = std::round(static_cast<double>(c.upper) * fraction.value);
    const uint64_t e = x <= 0 ? 0 : x >= static_cast<double>(c.upper) ? c.upper
                                                                    : static_cast<uint64_t>(x);
    return std::max(c.lower, e);
}

/**
 * The JSON of one count (§7.4). |fraction|: the graph's dummy fraction when it is served
 * without its dummy-edge mask (counting "upper_bound"), null with the mask: a count with
 * relation bounds then also carries `estimate`, and |estimated| is set. With the mask nothing
 * is added.
 */
Json::Value count_json(const Count &c, const DummyFraction *fraction = nullptr,
                       bool *estimated = nullptr) {
    Json::Value v;
    // the conservative number: the count, or the lower bound; null when nothing is known
    v["value"] = c.relation == Relation::UNKNOWN ? Json::Value() : uint_json(c.value);
    v["relation"] = to_string(c.relation);
    v["unit"] = to_string(c.unit);
    if (c.relation == Relation::BOUNDS) {
        v["lower"] = uint_json(c.lower);
        v["upper"] = uint_json(c.upper);
        if (fraction) {
            v["estimate"] = uint_json(estimate(c, *fraction));
            if (estimated)
                *estimated = true;
        }
    }
    return v;
}

// by_strand (keys +, -, both) on a BASIC graph, by_orientation (forward, reverse,
// palindromic) on the others, where no strand is known (§3, "Strand")
void put_orientations(Json::Value *count, const std::map<Orientation, Count> &by,
                      bool strand_stated, const DummyFraction *fraction = nullptr,
                      bool *estimated = nullptr) {
    Json::Value o(Json::objectValue);
    for (const auto &[orientation, c] : by) {
        o[strand_stated ? strand_key(orientation) : orientation_key(orientation)]
                = count_json(c, fraction, estimated);
    }
    (*count)[strand_stated ? "by_strand" : "by_orientation"] = std::move(o);
}

Json::Value stop_json(const std::optional<Stop> &stop) {
    if (!stop)
        return Json::Value();
    Json::Value s;
    s["phase"] = to_string(stop->phase);
    s["reason"] = to_string(stop->reason);
    return s;
}

Json::Value error_json(const std::string &code, const std::string &message) {
    Json::Value e;
    e["code"] = code;
    e["message"] = message;
    return e;
}

/**
 * The answer's entry for one pattern (§7.2). |results| holds the JSON
 * of the contexts enumerate() released (empty after count()); it is published only when the
 * engine withheld nothing, and moved into the entry. |released| is the number of contexts the
 * engine released: results.size(), or more with output.labels "all" when the memory account
 * admitted only the first of them (their objects not built; apply_labels states the cut).
 */
Json::Value entry_json(const PatternSpec &spec, const Result *result, Mode mode,
                       bool strand_stated, Json::Value results, uint64_t released,
                       bool long_paths, const DummyFraction *fraction) {
    Json::Value e;
    e["id"] = spec.id;
    e["kind"] = to_string(spec.kind);
    if (spec.error) {
        e["error"] = error_json(spec.error->first, spec.error->second);
        return e;
    }
    assert(spec.pattern && result);
    e["pattern"] = spec.pattern->text();
    // L, in bases, whatever the kind (a peptide's 3m)
    e["length"] = uint_json(spec.pattern->length());
    if (spec.kind == PatternKind::PROTEIN) {
        // a peptide's residues (m) and the genetic code it was read in
        e["residues"] = uint_json(spec.pattern->text().size());
        e["genetic_code"] = spec.pattern->genetic_code();
    }
    e["information_bits"] = result->information_bits;
    e["anchor_information_bits"] = result->anchor_information_bits
            ? Json::Value(*result->anchor_information_bits) : Json::Value();
    // the least informative searched anchor window, the floor's operand for L > k (an addition
    // under contract version 1)
    e["min_anchor_information_bits"] = result->min_anchor_information_bits
            ? Json::Value(*result->min_anchor_information_bits) : Json::Value();
    if (result->refusal) {
        e["error"] = error_json(result->refusal->code, result->refusal->message);
        return e;
    }

    e["mode"] = to_string(mode);
    e["scope"] = to_string(result->scope);
    Json::Value searched(Json::arrayValue);
    for (Orientation o : result->searched) {
        searched.append(strand_stated ? strand_symbol(o) : orientation_key(o));
    }
    e["strands"] = std::move(searched);
    e["palindromic"] = result->palindromic;

    // a graph without its mask: every count with relation bounds carries its estimate, and the
    // entry then says what the estimate rests on (kNoteEstimate)
    bool estimated = false;
    Json::Value counts;
    if (result->contexts) {
        Json::Value c = count_json(result->contexts->total, fraction, &estimated);
        c["suffix"] = count_json(result->contexts->suffix, fraction, &estimated);
        Json::Value by_offset(Json::objectValue);
        for (const auto &[offset, count] : result->contexts->by_offset) {
            by_offset[std::to_string(offset)] = count_json(count, fraction, &estimated);
        }
        c["by_offset"] = std::move(by_offset);
        put_orientations(&c, result->contexts->by_orientation, strand_stated, fraction,
                         &estimated);
        counts["contexts"] = std::move(c);
    } else if (result->anchors) {
        Json::Value a = count_json(result->anchors->total, fraction, &estimated);
        put_orientations(&a, result->anchors->by_orientation, strand_stated, fraction,
                         &estimated);
        counts["anchors"] = std::move(a);
        Json::Value paths = count_json(result->anchors->paths, fraction, &estimated);
        if (long_paths) {
            // long_search "paths" (§4.2): the paths counted by the extension, per orientation,
            // beside the branches it entered (work, not a count of the pattern) and what it did
            // (no_anchors, not_started, not_admitted, stopped, completed)
            put_orientations(&paths, result->anchors->paths_by_orientation, strand_stated,
                             fraction, &estimated);
            paths["candidates_examined"] = uint_json(result->anchors->candidates_examined);
            paths["extension"] = to_string(result->anchors->extension);
        }
        counts["paths"] = std::move(paths);
    } else {
        throw std::logic_error("pattern: the engine answered a pattern without counts");
    }
    // the engine reads no annotation (§7.2): the label counts are stated unknown, never omitted
    counts["labels"] = count_json(Count::unknown(Unit::LABELS));
    counts["occurrences"] = count_json(Count::unknown(Unit::PLACED_OCCURRENCES));
    e["counts"] = std::move(counts);

    Json::Value work;
    work["ranges_visited"] = uint_json(result->work.ranges_visited);
    work["mask_scans"] = uint_json(result->work.mask_scans);
    work["steps"] = uint_json(result->work.steps);
    if (long_paths && result->anchors) {
        // the outgoing edges the extension examined, one step each (part of steps)
        work["extension_edges"] = uint_json(result->work.extension_edges);
        // (not steps) the anchors whose extension began, each spelled once (k - 1 BOSS steps no
        // step charges), and the nodes the extension expanded with two or more k-mers allowed
        // at the next position
        work["extension_anchors"] = uint_json(result->work.extension_anchors);
        work["extension_branches"] = uint_json(result->work.extension_branches);
    }
    e["work"] = std::move(work);
    e["stop"] = stop_json(result->stop);

    if (mode == Mode::COUNT) {
        // a count names no context and licenses no absence (§5.1)
        e["retrieval_complete"] = false;
    } else {
        if (!result->extraction)
            throw std::logic_error("pattern: the engine answered a retrieval without extraction");
        const Extraction &x = *result->extraction;
        if (x.withheld) {
            results = Json::Value(Json::arrayValue);
        } else if (released != x.returned || results.size() > released) {
            throw std::logic_error("pattern: the engine released "
                                   + std::to_string(released) + " contexts and stated "
                                   + std::to_string(x.returned));
        }
        // every context of the pattern is in results: an absence claim over graph contexts
        // in this scope, never over labels (none were read)
        e["retrieval_complete"] = x.complete && !x.withheld;
        e["withheld"] = x.withheld ? reason_json(to_string(*x.withheld)) : Json::Value();
        e["returned"] = uint_json(results.size());
        // all_or_count returns all or nothing: it has no cut
        e["cut"] = mode == Mode::PARTIAL && x.cut && !x.withheld && !x.complete
                ? reason_json(to_string(*x.cut)) : Json::Value();
        e["results"] = std::move(results);
    }

    e["absence_scope"] = absence_scope(result->scope);
    e["determinism"] = result->time_limited ? "time_limited" : "full";
    Json::Value notes(Json::arrayValue);
    for (const std::string &note : result->notes) {
        notes.append(note);
    }
    // (the engine's notes first, no_stop_codon among them; then the route's note of the
    // estimates)
    if (estimated)
        notes.append(kNoteEstimate);
    e["notes"] = std::move(notes);
    Json::Value timing;
    timing["elapsed_ms"] = result->elapsed_ms;
    if (long_paths && result->anchors)
        timing["extension_ms"] = result->extension_ms;
    e["timing"] = std::move(timing);
    return e;
}

/**
 * The labelled retrieval's counters of one entry, which apply_labels does not merge: in every
 * entry with labels (0 where none were read, e.g. a withheld count), the distinct rows read
 * beside the reads (work.annotation_rows counts both steps' reads); for the paths of a long
 * pattern also the verification's units of work and the time of the paths' label lists and of
 * their verification. Additive fields of work and timing; nothing for a refused slot.
 */
void put_retrieval_counters(Json::Value *entry, const RetrievalCounters &c, bool paths) {
    Json::Value &e = *entry;
    if (e.isMember("error"))
        return;
    e["work"]["annotation_rows_distinct"] = uint_json(c.rows_distinct);
    if (!paths)
        return;
    e["work"]["verification_steps"] = uint_json(c.verification_steps);
    e["timing"]["label_intersection_ms"] = c.label_intersection_ms;
    e["timing"]["verification_ms"] = c.verification_ms;
}

/**
 * What the route made of one pattern's selection (SPEC §19.6-§19.10), beside the
 * SelectionAnswer of PatternRetrieval::select (or of constant_selection,
 * selection_without_pass): how the entry is composed from the engine's answer and the pass's.
 */
struct SelectionEntry {
    SelectionAnswer answer;
    // the engine ran with the selection's raw threshold (max_contexts = max_predicate_contexts):
    // its withheld count_above_threshold, its stop and cut max_contexts are the predicate's
    // (predicate_above_threshold, max_predicate_contexts), and the results are the selected
    // contexts the route built (a pass that ran, or could not start); false for a constant
    // normal form and a pattern longer than k (the engine's own thresholds and results)
    bool pass_path = false;
    // selection.support and selection.access
    const char *support = "kmer";
    const char *access = "rows";
    // the results of the selected contexts (built after the pass): a stop while building them
    // (output: time, or max_memory without a projection), and what it does to the list
    std::optional<std::pair<std::string, std::string>> output_stop;
    std::optional<std::string> output_withheld;
    std::optional<std::string> output_cut;
    // the annotation was read without the budget-aware decode (stated by a note when no
    // projection states it)
    bool unbudgeted_read = false;
    // the request's memory account at its peak once the pattern's work was done
    uint64_t memory_bytes = 0;
};

/**
 * Merges one pattern's selection into its |entry| (entry_json's, of the engine's answer; the
 * results the route built for the selected contexts in it, for a pass): counts.tested and
 * counts.selected, selection, absence_filter, work.predicate_*, timing.selection_ms, the stop
 * (the engine's first, then the pass's, then the results' output), determinism, and on the pass
 * path in a retrieval mode withheld (the engine's, then the pass's, then the output's), cut
 * (partial, §19.8: the engine's stop, the pass's stop, the raw release's cut, the selected
 * list's cut; a cut of the output replaces them, as it does without a predicate) and
 * retrieval_complete (every selected context returned). The projection's labels are merged
 * after it (apply_labels); the notes after those.
 */
void put_selection(Json::Value *entry, const SelectionEntry &s, Mode mode) {
    Json::Value &e = *entry;
    if (e.isMember("error"))
        return;
    const SelectionAnswer &sel = s.answer;
    e["counts"]["tested"] = count_json(sel.tested);
    e["counts"]["selected"] = count_json(sel.selected);
    Json::Value selection;
    selection["pass"] = to_string(sel.pass);
    selection["support"] = s.support;
    selection["access"] = s.access;
    e["selection"] = std::move(selection);
    // the absence claims of this entry are the predicate's (§19.11): absence_scope keeps its
    // closed values
    e["absence_filter"] = "predicate";
    e["work"]["predicate_rows"] = uint_json(sel.rows);
    e["work"]["predicate_units"] = uint_json(sel.units);
    e["work"]["predicate_lookups"] = uint_json(sel.lookups);
    e["work"]["memory_bytes"] = uint_json(s.memory_bytes);
    e["timing"]["selection_ms"] = sel.ms;

    // (read without adding the member: an entry of mode count has no withheld and no cut)
    auto reason_of = [&e](const char *key) {
        return e.isMember(key) && e[key].isObject() ? e[key]["reason"].asString()
                                                    : std::string();
    };
    if (s.pass_path) {
        // the engine counted the raw contexts against max_predicate_contexts
        if (reason_of("stop") == "max_contexts")
            e["stop"]["reason"] = "max_predicate_contexts";
        if (reason_of("withheld") == "count_above_threshold")
            e["withheld"] = reason_json("predicate_above_threshold");
    }
    for (const auto &stop : { sel.stop, s.output_stop }) {
        if (stop && e["stop"].isNull()) {
            Json::Value v;
            v["phase"] = stop->first;
            v["reason"] = stop->second;
            e["stop"] = std::move(v);
        }
    }
    if (sel.time_limited || (s.output_stop && s.output_stop->second == "time"))
        e["determinism"] = "time_limited";
    if (mode == Mode::COUNT || !s.pass_path)
        return;

    if (!e["withheld"].isObject()) {
        const std::optional<std::string> withheld = sel.withheld ? sel.withheld
                                                                 : s.output_withheld;
        if (withheld) {
            e["withheld"] = reason_json(withheld->c_str());
            e["retrieval_complete"] = false;
            e["returned"] = 0;
            e["results"] = Json::Value(Json::arrayValue);
        }
    }
    if (e["withheld"].isObject()) {
        e["cut"] = Json::Value();
        return;
    }
    // every selected context of the pattern is in results: the engine released every raw
    // context, the pass decided each, and every selected one was built
    const uint64_t returned = e["results"].size();
    e["returned"] = uint_json(returned);
    e["retrieval_complete"] = e["retrieval_complete"].asBool()
                                && sel.pass == SelectionPass::COMPLETED
                                && sel.selected.relation == Relation::EXACT
                                && returned == sel.selected.value;
    if (mode != Mode::PARTIAL)
        return;
    // the cut (§19.8): the first that applies
    std::optional<std::string> cut;
    const std::string engine_cut = reason_of("cut");
    if (!engine_cut.empty() && engine_cut != "max_contexts") {
        cut = engine_cut;                   // the engine's stop (max_steps, time)
    } else if (sel.cut) {
        cut = *sel.cut;                     // the pass's stop, refused rows, descriptors
    } else if (engine_cut == "max_contexts") {
        cut = "max_predicate_contexts";     // the raw release held the first of them only
    } else if (sel.list_cut) {
        cut = "max_contexts";               // more selected than the list holds
    }
    if (s.output_cut)
        cut = *s.output_cut;                // the output's: the list's length
    e["cut"] = cut && !e["retrieval_complete"].asBool() ? reason_json(cut->c_str())
                                                        : Json::Value();
}

} // namespace


std::string alphabet_refusal(const std::string &alphabet) {
    if (alphabet == "$ACGT")
        return "";
    // the engine counts on $ACGTN, but no DNA5 build has run its tests; the route serves it
    // once one passes
    if (alphabet == "$ACGTN")
        return "alphabet_untested";
    return "alphabet_unsupported";
}

GraphSupport route_support(const DeBruijnGraph &graph) {
    GraphSupport support = PatternSearch::support(graph);
    // (a graph without its mask is served, counting upper bounds)
    if (support.supported) {
        const std::string refusal = alphabet_refusal(support.alphabet);
        if (!refusal.empty()) {
            support.supported = false;
            support.reason = refusal;
            return support;
        }
    }
    // a mask that marks a W = $ edge valid, found once at load (check_mask_at_load): a loaded
    // mask is trusted only once checked; a graph without a mask counts upper bounds and needs
    // no such check
    if (support.supported && support.mask_present && mask_invalid_at_load(graph)) {
        support.supported = false;
        support.reason = "mask_invalid";
    }
    return support;
}

namespace {

// The dummy fractions sampled so far, one per succinct graph served without its mask (weak
// references: a graph freed and another allocated at its address must not pass for it, as in
// load_annotated_graph.cpp's mask registries)
std::mutex dummy_fractions_mutex;
std::vector<std::pair<std::weak_ptr<const DBGSuccinct>, DummyFraction>> dummy_fractions;

// the succinct graph |graph| is, or the PRIMARY one its CanonicalDBG wraps; null otherwise
std::shared_ptr<const DBGSuccinct> succinct_of(std::shared_ptr<const DeBruijnGraph> graph) {
    if (auto canonical = std::dynamic_pointer_cast<const CanonicalDBG>(graph))
        graph = canonical->get_graph_ptr();
    return std::dynamic_pointer_cast<const DBGSuccinct>(graph);
}

// the kept fraction of |dbg_succ|, sampled now when there is none; the caller holds
// dummy_fractions_mutex
const DummyFraction& fraction_of(const std::shared_ptr<const DBGSuccinct> &dbg_succ,
                                 bool *sampled = nullptr) {
    for (auto it = dummy_fractions.begin(); it != dummy_fractions.end(); ) {
        if (auto held = it->first.lock()) {
            if (held == dbg_succ)
                return it->second;
            ++it;
        } else {
            it = dummy_fractions.erase(it);
        }
    }
    if (sampled)
        *sampled = true;
    dummy_fractions.emplace_back(dbg_succ, sample_real_fraction(*dbg_succ));
    return dummy_fractions.back().second;
}

} // namespace

std::optional<DummyFraction> dummy_fraction(const AnnotatedDBG &anno_graph) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());
    auto dbg_succ = graph ? succinct_of(graph) : nullptr;
    if (!dbg_succ || dbg_succ->get_mask())
        return std::nullopt;
    std::lock_guard<std::mutex> lock(dummy_fractions_mutex);
    return fraction_of(dbg_succ);
}

void sample_dummy_fraction_at_load(const std::shared_ptr<DeBruijnGraph> &graph,
                                   bool stdout_reserved) {
    auto dbg_succ = succinct_of(graph);
    if (!dbg_succ || dbg_succ->get_mask())
        return;
    bool sampled = false;
    DummyFraction f;
    const auto start = std::chrono::steady_clock::now();
    {
        std::lock_guard<std::mutex> lock(dummy_fractions_mutex);
        f = fraction_of(dbg_succ, &sampled);
    }
    if (!sampled)
        return;
    logger->log(stdout_reserved ? spdlog::level::trace : spdlog::level::info,
                "Dummy fraction sampled for the pattern search in {:.3f} s (the graph has no "
                "dummy-edge mask: counts are upper bounds with estimates): {} real k-mers "
                "among {} entries with W != $ drawn of {} edges ({} with W = $), f = {:.6f}, "
                "95% interval [{:.6f}, {:.6f}], seed {}",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(),
                f.real, f.samples, f.edges, f.sentinel_edges, f.value, f.lower, f.upper,
                f.seed);
}

Json::Value dummy_fraction_json(const DummyFraction &f) {
    Json::Value v;
    v["value"] = f.value;
    Json::Value interval(Json::arrayValue);
    interval.append(f.lower);
    interval.append(f.upper);
    v["interval"] = std::move(interval);
    v["samples"] = uint_json(f.samples);
    v["source"] = "sampled";
    return v;
}

Json::Value PatternRefusal::body() const {
    Json::Value b;
    b["error"] = what();
    b["code"] = code_;
    return b;
}

void PatternDelivery::check() const {
    if (abort_ && abort_())
        throw Aborted();
    if (deadline_ && deadline_->respond_expired()) {
        throw PatternRefusal(503, "deadline",
                             "pattern: the answer could not be written within time_budget_ms ("
                             + ms_text(deadline_->time_budget_ms())
                             + " ms, the finalisation reserve of "
                             + ms_text(deadline_->finalize_reserve_ms())
                             + " ms included): nothing partial is sent");
    }
}

PatternLimits pattern_limits(const Config &config) {
    PatternLimits limits;
    limits.max_contexts = config.pattern_max_contexts;
    limits.max_anchors = config.pattern_max_anchors;
    limits.max_paths = config.pattern_max_paths;
    limits.max_steps = config.pattern_max_steps;
    limits.default_time_ms = static_cast<double>(config.pattern_default_time_ms);
    limits.max_time_ms = static_cast<double>(config.pattern_max_time_ms);
    limits.finalize_ms = static_cast<double>(config.pattern_finalize_ms);
    limits.min_information_bits = config.pattern_min_information_bits;
    limits.max_checked_entries = config.pattern_max_checked_entries;
    limits.max_patterns = config.pattern_max_patterns;
    limits.max_labels_per_anchor = config.pattern_max_labels_per_anchor;
    limits.max_annotation_work = config.pattern_max_annotation_work;
    limits.max_memory_mb = config.pattern_max_memory_mb;
    limits.max_labels = config.pattern_max_labels;
    limits.max_occurrences_per_label = config.pattern_max_occurrences;
    limits.max_predicate_contexts = config.pattern_max_predicate_contexts;
    limits.max_predicate_work = config.pattern_max_predicate_work;
    limits.max_predicate_labels = config.pattern_max_predicate_labels;
    limits.chunk_target_ms = static_cast<double>(config.traverse_chunk_target_ms);
    limits.delivery_build_mbps = config.pattern_delivery_build_mbps;
    limits.delivery_compress_mbps = config.pattern_delivery_compress_mbps;
    return limits;
}

Json::Value parse_pattern_body(const std::string &content) {
    // one RFC 8259 JSON text, its members' names unique: no comments, no trailing comma,
    // nothing after the value, and no duplicated member name, whose earlier value jsoncpp would
    // drop without a word (a budget given twice would be replaced silently). Not
    // CharReaderBuilder::strictMode(): its strictRoot would refuse a root that is not an object
    // here, before the graph check (§5 step 4), where step 5 does. The route's own reader:
    // /search, /align and /traverse keep theirs
    // jsoncpp skips a comment before an object's member name whatever allowComments says: a
    // '/' outside a string is no JSON at all, so the text is refused before it is parsed
    bool in_string = false;
    for (size_t i = 0; i < content.size(); ++i) {
        if (in_string) {
            if (content[i] == '\\')
                ++i;
            else if (content[i] == '"')
                in_string = false;
        } else if (content[i] == '"') {
            in_string = true;
        } else if (content[i] == '/') {
            throw invalid("request: not JSON: a comment (or a '/' outside a string) at byte "
                          + std::to_string(i));
        }
    }
    Json::CharReaderBuilder builder;
    builder["allowComments"] = false;
    builder["allowTrailingCommas"] = false;
    builder["failIfExtra"] = true;
    builder["rejectDupKeys"] = true;
    builder["stackLimit"] = kMaxBodyNesting;
    std::unique_ptr<Json::CharReader> reader { builder.newCharReader() };
    Json::Value json;
    std::string errors;
    try {
        if (!reader->parse(content.data(), content.data() + content.size(), &json, &errors))
            throw invalid("request: not JSON: " + errors);
    } catch (const Json::Exception &e) {
        // jsoncpp throws rather than returns past its nesting limit: the client's to fix, not a
        // server bug (a 400 without a code)
        throw invalid("request: not JSON: " + std::string(e.what()) + " (more than "
                      + std::to_string(kMaxBodyNesting) + " nested arrays and objects)");
    }
    return json;
}

Json::Value process_pattern_request(const Json::Value &json,
                                    const AnnotatedDBG &anno_graph,
                                    const Config &config,
                                    const IndexIdentity *identity,
                                    PatternDelivery *delivery) {
    return process_pattern_request(json, anno_graph, pattern_limits(config),
                                   config.index_release, identity, delivery);
}

Json::Value process_pattern_request(
        const Json::Value &json,
        const AnnotatedDBG &anno_graph,
        const PatternLimits &limits,
        const std::string &release,
        const IndexIdentity *identity,
        PatternDelivery *delivery,
        const std::function<Deadline::Clock::time_point()> &clock,
        const RetrievalHooks *hooks) {
    // the request's one deadline starts here, its body parsed (§5.3)
    const std::function<Deadline::Clock::time_point()> now
            = clock ? clock : std::function<Deadline::Clock::time_point()>(&Deadline::Clock::now);
    const Deadline::Clock::time_point start = now();
    const DeBruijnGraph &graph = anno_graph.get_graph();

    // the graph first: no request is answerable on a graph the engine
    // cannot count on (or the route does not serve), whatever it asks
    const GraphSupport support = route_support(graph);
    if (!support.supported)
        throw PatternRefusal(400, support.reason, support_message(support));

    ParsedRequest req = parse_request(json, limits);
    const Mode mode = req.request.mode;

    const Deadline deadline(start, req.time_budget_ms, limits.finalize_ms, now);
    if (delivery)
        delivery->set_deadline(deadline);
    /**
     * The work's deadline: the request's, but read through a clock that runs ahead of |now| by
     * the estimated time to write what the answer holds so far (AnswerVolume), so that the work
     * stops finalize_ms plus that estimate before the deadline, and a request that buffered
     * many results still answers within its budget (stopped by time, its counts kept) rather
     * than 503 with all of its work lost. Every reading of the work time — the engine's
     * (discovery, the release), the annotation reads' and the output of their labels — reads
     * it; the answer's own deadline (|deadline|: the delivery check, timing.elapsed_ms) reads
     * |now|.
     */
    AnswerVolume volume(limits.delivery_build_mbps, limits.delivery_compress_mbps,
                        limits.delivery_text_scale);
    auto work_clock = [&now, &volume]() {
        // (bounded far above any budget, so that absurd rates cannot overflow the clock)
        return now() + std::chrono::duration_cast<Deadline::Clock::duration>(
                std::chrono::duration<double, std::milli>(std::min(volume.finalize_ms(), 1e12)));
    };
    // one budget for the request, spent by the patterns in request order (§5.3, §5.5)
    Budget budget(req.max_steps, Deadline(start, req.time_budget_ms, limits.finalize_ms,
                                          work_clock));
    // a caller that left is not answered: the work ends at its next clock reading
    if (delivery && delivery->abort())
        budget.set_abort(delivery->abort());
    const PatternSearch search(graph);

    const size_t k = graph.get_k();
    const bool strand_stated = support.strand_stated;
    const uint64_t num_rows = anno_graph.get_annotator().num_objects();
    // without the dummy-edge mask: the counts are upper bounds, each with its estimate from the
    // graph's dummy fraction (sampled once per graph)
    const std::optional<DummyFraction> fraction = support.mask_present
            ? std::nullopt : dummy_fraction(anno_graph);
    if (!support.mask_present && !fraction)
        throw std::logic_error("pattern: no dummy fraction for a graph without its mask");

    // output.labels "all" in a retrieval mode reads the annotation (§4.3): on the budget-aware
    // path, or unbudgeted by the request's explicit opt-in
    const bool read_labels = req.labels_all && mode != Mode::COUNT;
    // a predicate's selection (SPEC §19) reads the annotation in every mode, and its projection
    // is the selected contexts' labels: none, the predicate's ("predicate_only") or all of
    // them, in a retrieval mode
    const bool has_predicate = req.predicate.has_value();
    const Projection projection = mode == Mode::COUNT ? Projection::NONE
                                : req.labels_all ? Projection::ALL
                                : req.labels_predicate_only ? Projection::PREDICATE_ONLY
                                                            : Projection::NONE;
    std::optional<PatternRetrieval> retrieval;
    if (read_labels || has_predicate) {
        retrieval.emplace(anno_graph, support.mode, req.retrieval, budget, hooks, &volume);
        if (!retrieval->description().budgeted && !req.retrieval.allow_unbudgeted) {
            if (has_predicate) {
                throw PatternRefusal(400, "annotation_unbudgeted",
                                     "pattern: a predicate reads the annotation in every mode "
                                     "(its selection reads the rows of the contexts it tests), "
                                     "and this index's annotation has no budget-aware decode "
                                     "(only the row-diff family has one): its reads would run "
                                     "without a memory bound. Set allow_unbudgeted_annotation: "
                                     "true to read it anyway (selection.access then says how: "
                                     "\"rows\", or \"columns\" for at most 16 labels), or ask "
                                     "without the predicate");
            }
            throw PatternRefusal(400, "annotation_unbudgeted",
                                 "pattern: output.labels \"all\" reads the annotation, and this "
                                 "index's annotation has no budget-aware decode (only the "
                                 "row-diff family has one): its reads would run without a "
                                 "memory bound. Set allow_unbudgeted_annotation: true to read "
                                 "it anyway (the answer then says annotation: unbudgeted), or "
                                 "ask for output.labels \"none\" or mode count");
        }
        // require_support "record_verified" (DESIGN §4.3) keeps the labels one record verifies,
        // which needs a BASIC index with coordinates and its record mapping; an index that
        // cannot verify refuses it rather than answer in the weaker mode
        if (read_labels && req.long_paths && req.require_verified
                && std::string(retrieval->description().support) != kSupportVerified) {
            throw PatternRefusal(400, "support_unavailable",
                                 "pattern: require_support \"record_verified\" needs a BASIC "
                                 "index with coordinates and its record mapping (.seqs), to "
                                 "check that one record holds the whole path; this index's "
                                 "best support is \""
                                 + std::string(retrieval->description().support)
                                 + "\" (placement \""
                                 + std::string(retrieval->description().placement)
                                 + "\"). Ask without require_support: each label of a path "
                                 "then states its support");
        }
    }
    // the graph's name in by_label (the index's --index-name; null without one)
    const Json::Value graph_name = identity && !identity->name.empty()
            ? Json::Value(identity->name) : Json::Value();

    /**
     * One released context as its JSON result (§7.2, output.labels none), built in the engine's
     * callback: spelling the k-mer is graph work, done while the engine still reads the clock
     * (every 64 contexts handed over, in partial's release and in all_or_count's delivery of
     * its buffered release), so that the finalisation reserve is left to assembling and writing
     * the answer. The row is the annotation row the context's k-mer is annotated in, named
     * without reading it: the stored k-mer's (base_node) on BASIC and wrapped PRIMARY graphs;
     * on a native CANONICAL graph the canonical k-mer's (the annotation key of every route,
     * LabelOracle::key_of), since the other orientation's row carries no labels.
     */
    auto context_json = [&](const Context &c, size_t length,
                            RetrievalContext *collected = nullptr) {
        Json::Value r;
        std::string kmer = graph.get_node_sequence(c.node);
        if (kmer.size() != k || c.offset + length > k) {
            throw std::logic_error("pattern: the engine released node " + std::to_string(c.node)
                                   + " at offset " + std::to_string(c.offset)
                                   + ", outside its k-mer");
        }
        DeBruijnGraph::node_index key = c.base_node;
        if (support.mode == GraphMode::CANONICAL) {
            key = DeBruijnGraph::npos;
            graph.map_to_nodes(kmer, [&](DeBruijnGraph::node_index n) { key = n; });
        }
        r["instance"] = kmer.substr(c.offset, length);
        if (collected) {
            // what the labelled retrieval reads: the row's key, the k-mer naming it
            collected->orientation = c.orientation;
            collected->offset = c.offset;
            collected->kmer = kmer;
            collected->key = key != DeBruijnGraph::npos
                                && AnnotatedDBG::graph_to_anno_index(key) < num_rows
                    ? key : DeBruijnGraph::npos;
        }
        r["kmer"] = std::move(kmer);
        r["offset"] = c.offset;
        if (strand_stated) {
            r["strand"] = strand_symbol(c.orientation);
        } else {
            r["orientation"] = orientation_key(c.orientation);
        }
        r["node"] = uint_json(c.node);
        // null: the k-mer has no row of this annotation (never expected on a compatible
        // index; stated rather than guessed)
        r["row"] = key != DeBruijnGraph::npos
                        && AnnotatedDBG::graph_to_anno_index(key) < num_rows
                ? uint_json(AnnotatedDBG::graph_to_anno_index(key)) : Json::Value();
        // in the answer from now on: the work time is read with it counted
        volume.add(compact_json_bytes(r));
        return r;
    };

    /**
     * One released path of a pattern longer than k as its JSON result (long_search "paths"):
     * the fields sequence (the L bases it spells), anchor_kmer (its anchor's k bases, as the
     * graph spells them) and the node path with the row of each k-mer (nodes, rows), never
     * kmer, which keeps its meaning (a context's k-mer); instance is the sequence, offset 0,
     * the strand or orientation as for contexts. The rows are named without reading them, as a
     * context's row (§7.10).
     */
    auto path_json = [&](const Context &c, size_t length, RetrievalPath *collected = nullptr) {
        const size_t n = length - k + 1;
        if (length <= k || c.path.size() != n || c.sequence.size() != length
                || c.path.front() != c.node || c.offset != 0) {
            throw std::logic_error("pattern: the engine released a path of "
                                   + std::to_string(c.path.size()) + " k-mers and "
                                   + std::to_string(c.sequence.size()) + " bases for a pattern "
                                   "of " + std::to_string(length) + " bases");
        }
        std::string anchor = graph.get_node_sequence(c.node);
        if (anchor.size() != k || c.sequence.compare(0, k, anchor) != 0) {
            throw std::logic_error("pattern: the engine released a path whose sequence does "
                                   "not start with its anchor's k-mer");
        }
        // the annotation key of each k-mer: the stored k-mer's (base_node) on BASIC and
        // wrapped PRIMARY graphs, the canonical k-mer's on a native CANONICAL graph
        std::vector<DeBruijnGraph::node_index> keys(n, DeBruijnGraph::npos);
        if (support.mode == GraphMode::CANONICAL) {
            size_t j = 0;
            graph.map_to_nodes(c.sequence, [&](DeBruijnGraph::node_index x) {
                if (j < n)
                    keys[j] = x;
                ++j;
            });
        } else {
            for (size_t j = 0; j < n; ++j) {
                keys[j] = search.base_node(c.path[j]);
            }
        }
        Json::Value r;
        Json::Value nodes(Json::arrayValue), rows(Json::arrayValue);
        for (size_t j = 0; j < n; ++j) {
            nodes.append(uint_json(c.path[j]));
            const bool valid = keys[j] != DeBruijnGraph::npos
                                && AnnotatedDBG::graph_to_anno_index(keys[j]) < num_rows;
            rows.append(valid ? uint_json(AnnotatedDBG::graph_to_anno_index(keys[j]))
                              : Json::Value());
            if (!valid)
                keys[j] = DeBruijnGraph::npos;
        }
        if (collected) {
            collected->orientation = c.orientation;
            collected->sequence = c.sequence;
            collected->keys.assign(keys.begin(), keys.end());
        }
        r["sequence"] = c.sequence;
        r["anchor_kmer"] = std::move(anchor);
        r["instance"] = c.sequence;
        r["offset"] = 0;
        if (strand_stated) {
            r["strand"] = strand_symbol(c.orientation);
        } else {
            r["orientation"] = orientation_key(c.orientation);
        }
        r["nodes"] = std::move(nodes);
        r["rows"] = std::move(rows);
        // in the answer from now on: the work time is read with it counted
        volume.add(compact_json_bytes(r));
        return r;
    };

    Json::Value entries(Json::arrayValue);
    struct Answered {
        std::optional<Result> result;
        Json::Value results;
        // the contexts the engine released (results.size(), or more when the memory account
        // of output.labels "all" admitted only the first of them; with a predicate's pass, the
        // raw contexts released into it)
        uint64_t released = 0;
        // output.labels "all" (or "predicate_only"): the pattern's labels, read in the work
        // phase
        std::optional<LabelsAnswer> labels;
        // the labels are a long pattern's paths' (retrieve_paths)
        bool label_paths = false;
        // the pattern's selection (a request with a predicate)
        std::optional<SelectionEntry> selection;
        // its normal form is a constant (no pass, no read: note predicate_constant)
        bool constant = false;
    };

    /**
     * One pattern answered by the engine's own thresholds (max_contexts, max_paths), as without
     * a predicate: its contexts (paths) released as results, with their labels when |labels|
     * (output.labels "all", or a projection's empty release of a long pattern's anchors).
     */
    auto answer_plain = [&](const Pattern &pattern, Answered &a, bool labels) {
        // long_search "paths": a pattern longer than k is answered by its paths
        const bool paths = req.long_paths && pattern.length() > k;
        if (mode == Mode::COUNT) {
            a.result = search.count(pattern, req.request, budget);
        } else if (!labels) {
            const size_t length = pattern.length();
            a.result = search.enumerate(pattern, req.request, budget, [&](const Context &c) {
                a.results.append(paths ? path_json(c, length) : context_json(c, length));
            });
            a.released = a.results.size();
        } else if (paths) {
            const size_t length = pattern.length();
            std::vector<RetrievalPath> collected;
            retrieval->begin_release(mode);
            a.result = search.enumerate(pattern, req.request, budget, [&](const Context &c) {
                ++a.released;
                // the path's descriptor, sequence and arrays are charged to the memory
                // account before its object is built; after the first that does not fit
                // none is built
                if (!retrieval->admit_path(length))
                    return;
                collected.emplace_back();
                a.results.append(path_json(c, length, &collected.back()));
            });
            if (!a.result->refusal && a.result->extraction) {
                a.labels = retrieval->retrieve_paths(collected, a.released, length, mode,
                                                     *a.result->extraction, graph_name,
                                                     req.require_verified);
                a.label_paths = true;
            }
        } else {
            const size_t length = pattern.length();
            std::vector<RetrievalContext> collected;
            retrieval->begin_release(mode);
            a.result = search.enumerate(pattern, req.request, budget, [&](const Context &c) {
                ++a.released;
                // the context's descriptor is charged to the memory account before its
                // object is built; after the first that does not fit none is built
                if (!retrieval->admit_context())
                    return;
                collected.emplace_back();
                a.results.append(context_json(c, length, &collected.back()));
            });
            // the labels of what was released (nothing when it was withheld), work
            // still: the reads end where the deadline's work time does (§5.3)
            if (!a.result->refusal && a.result->extraction) {
                a.labels = retrieval->retrieve(collected, a.released, length, mode,
                                               *a.result->extraction, graph_name);
            }
        }
    };

    // the request's predicate bound to the index's columns, once, before the first pattern (its
    // bytes in the memory account; the work time read every 4,096 names); its echo (the normal
    // form, the unknown names) is text the answer will hold. A binding the time or the account
    // stopped leaves no bound predicate: every pattern's selection is then not_started, stopped
    // as the binding was (stop {selection, time | max_memory})
    const predicate::Bound *bound = nullptr;
    if (has_predicate) {
        retrieval->bind(*req.predicate, req.selection);
        bound = retrieval->bound();
        if (bound)
            volume.add_pending(bound->text_bytes());
    }
    // the raw search of a pattern's selection (§19.6 steps 1 and 2): every raw context released
    // into the pass up to max_predicate_contexts, all or none in all_or_count (and in count,
    // which reads them too), the first ones in partial; stop_at_threshold on that threshold
    Request selection_request = req.request;
    selection_request.max_contexts = req.max_predicate_contexts;
    if (mode == Mode::COUNT)
        selection_request.mode = Mode::ALL_OR_COUNT;
    SelectionRequest pass_request;
    pass_request.mode = mode;
    pass_request.max_contexts = req.request.max_contexts;
    pass_request.stop_at_threshold = req.request.stop_at_threshold;
    pass_request.projection = projection;

    /**
     * One pattern of L <= k with a predicate (SPEC §19.6): a constant normal form without a
     * pass (false: discovery as mode count, nothing selected; true: the unfiltered request),
     * else the engine's release into the pass (each raw context's 64-byte descriptor admitted
     * before it is kept), the pass (PatternRetrieval::select), the results of the chosen
     * contexts (their objects admitted, 512 + 2k, before each is built; the clock every 64),
     * then the projection: none, the predicate's labels on the chosen rows (retrieve_given) or
     * every label of them (retrieve). A pass that cannot start (the predicate not bound, or the
     * request's selection work spent by an earlier pattern: sticky) leaves the raw search as
     * mode count does it, nothing retained.
     */
    auto answer_selection = [&](const Pattern &pattern, Answered &a) {
        const size_t length = pattern.length();
        SelectionEntry s;
        s.access = retrieval->selection_access();
        s.unbudgeted_read = !retrieval->description().budgeted;
        auto finish = [&]() {
            s.memory_bytes = retrieval->memory_peak();
            a.selection = std::move(s);
        };
        // an empty release (nothing selected) as the projection reads it: nothing read, every
        // count exact 0
        auto nothing_selected = [&](const Extraction &x) {
            if (projection == Projection::NONE)
                return;
            retrieval->begin_release(mode);
            a.labels = retrieval->retrieve({}, 0, length, mode, x, graph_name);
        };

        if (bound && bound->constant()) {
            a.constant = true;
            if (!*bound->constant()) {
                // constant false: discovery as mode count does it, nothing retained, no read;
                // nothing can pass, also after a discovery stop
                a.result = search.count(pattern, req.request, budget);
                if (a.result->refusal)
                    return;
                s.answer = constant_selection(false, a.result->contexts->total);
                if (mode != Mode::COUNT) {
                    // every selected context (none) is returned
                    a.result->extraction = Extraction();
                    a.result->extraction->complete = true;
                    nothing_selected(*a.result->extraction);
                }
                finish();
                return;
            }
            // true: the unfiltered request, the engine's thresholds; output.labels "all" reads
            // the labels as without a predicate, "predicate_only" has none to read (the
            // predicate names no column of the index)
            answer_plain(pattern, a, projection == Projection::ALL);
            if (a.result->refusal)
                return;
            s.answer = constant_selection(true, a.result->contexts->total);
            if (mode != Mode::COUNT && projection != Projection::NONE) {
                for (Json::Value &r : a.results) {
                    r["selection_labels"] = Json::Value(Json::arrayValue);
                    r["selection_strands"] = Json::Value(Json::arrayValue);
                }
                volume.add(a.results.size()
                           * std::string(",\"selection_labels\":[],\"selection_strands\":[]")
                                     .size());
            }
            if (projection == Projection::PREDICATE_ONLY && a.result->extraction) {
                Extraction none;
                none.complete = true;
                if (a.result->extraction->withheld)
                    none.withheld = a.result->extraction->withheld;
                nothing_selected(none);
                Json::Value f;
                f["support"] = "kmer";
                f["labels_status"] = "complete";
                f["labels_total"] = 0;
                f["labels"] = Json::Value(Json::arrayValue);
                a.labels->result_fields.assign(a.results.size(), f);
                volume.add(a.results.size() * compact_json_bytes(f));
            }
            finish();
            return;
        }

        s.pass_path = true;
        if (!bound || retrieval->predicate_units() >= req.selection.max_predicate_work) {
            // the pass cannot start: the raw discovery as mode count does it (against the
            // selection's threshold), nothing released; the pass states why it did not start
            a.result = search.count(pattern, selection_request, budget);
            if (a.result->refusal)
                return;
            Extraction x;
            std::vector<TestedContext> none;
            retrieval->begin_selection(mode);
            s.answer = retrieval->select(none, 0, a.result->contexts->total, x, pass_request);
            retrieval->end_selection();
            if (mode != Mode::COUNT) {
                a.result->extraction = x;
                nothing_selected(x);
            }
            finish();
            return;
        }

        // ---- the raw contexts released into the pass
        std::vector<TestedContext> tested;
        retrieval->begin_selection(mode);
        a.result = search.enumerate(pattern, selection_request, budget, [&](const Context &c) {
            ++a.released;
            // its descriptor is charged before it is kept; after the first that does not fit
            // none is kept (the pass's set ends there, stated)
            if (!retrieval->admit_tested())
                return;
            tested.push_back(retrieval->tested_context(c));
        });
        if (a.result->refusal || !a.result->extraction) {
            retrieval->end_selection();
            return;
        }
        const Extraction &x = *a.result->extraction;
        s.answer = retrieval->select(tested, a.released, a.result->contexts->total, x,
                                     pass_request);
        const SelectionAnswer &sel = s.answer;
        if (mode == Mode::COUNT) {
            retrieval->end_selection();
            finish();
            return;
        }

        // ---- the results of the chosen contexts, in answer order
        std::vector<RetrievalContext> collected;
        retrieval->begin_release(mode);
        bool output_memory = false;
        if (!x.withheld && !sel.withheld) {
            for (size_t j = 0; j < sel.chosen.size(); ++j) {
                // the clock every 64 result objects (spelling a k-mer is graph work)
                if (j % Budget::kReleaseClockStride == 0 && !budget.check_time()) {
                    s.output_stop = std::make_pair(std::string("output"), std::string("time"));
                    if (mode == Mode::PARTIAL) {
                        s.output_cut = "time";
                    } else {
                        s.output_withheld = "deadline";
                    }
                    break;
                }
                // its object charged before it is built (512 + 2k); after the first that does
                // not fit none is built (a projection states it: its released is the chosen)
                if (!retrieval->admit_context()) {
                    output_memory = true;
                    break;
                }
                const TestedContext &t = tested[sel.chosen[j]];
                Context c;
                c.orientation = t.orientation;
                c.offset = t.offset;
                c.node = t.node;
                c.base_node = t.base_node;
                collected.emplace_back();
                Json::Value r = context_json(c, length, &collected.back());
                if (projection != Projection::NONE) {
                    // why it was selected: the predicate's labels in the set it was evaluated
                    // on, and per label the orientation whose row carries it, their bytes
                    // charged by the pass
                    r["selection_labels"] = retrieval->selection_labels_json(sel, j);
                    r["selection_strands"] = retrieval->selection_strands_json(sel, j);
                    volume.add(compact_json_bytes(r["selection_labels"])
                               + std::string(",\"selection_labels\":").size()
                               + compact_json_bytes(r["selection_strands"])
                               + std::string(",\"selection_strands\":").size());
                }
                a.results.append(std::move(r));
            }
        }
        if (output_memory && projection == Projection::NONE) {
            // the account could not hold the next result object: the list ends there
            s.output_stop = std::make_pair(std::string("output"), std::string("max_memory"));
            if (mode == Mode::PARTIAL) {
                s.output_cut = "max_memory";
            } else {
                s.output_withheld = "output_budget";
            }
        }
        if (s.output_withheld)
            a.results = Json::Value(Json::arrayValue);

        // ---- the projection of what was built (nothing is read for what is withheld)
        if (projection != Projection::NONE) {
            Extraction given;
            given.returned = collected.size();
            given.complete = sel.selected.relation == Relation::EXACT
                                && sel.chosen.size() == sel.selected.value;
            if (x.withheld || sel.withheld || s.output_withheld) {
                given.withheld = x.withheld ? *x.withheld : Withheld::DEADLINE;
                collected.clear();
            }
            // a time stop of the output ended the list (cut: time): the projection's release
            // is what was built; a memory refusal is the projection's own cut (max_memory)
            const uint64_t released = s.output_stop || given.withheld ? collected.size()
                                                                       : sel.chosen.size();
            a.labels = projection == Projection::PREDICATE_ONLY
                    ? retrieval->retrieve_given(collected, released, length, mode, given,
                                                graph_name)
                    : retrieval->retrieve(collected, released, length, mode, given,
                                          graph_name);
        }
        retrieval->end_selection();
        finish();
    };

    std::vector<Answered> answered;
    answered.reserve(req.patterns.size());
    for (const PatternSpec &spec : req.patterns) {
        Answered a;
        a.results = Json::Value(Json::arrayValue);
        if (spec.pattern) {
            if (!has_predicate) {
                answer_plain(*spec.pattern, a, read_labels);
            } else if (spec.pattern->length() <= k) {
                answer_selection(*spec.pattern, a);
            } else {
                // a pattern longer than k under long_search "anchors" (with "paths" a predicate
                // is refused): its anchors as without a predicate, a projection reading nothing
                // (no context is released); its selection did not start (§19.5)
                answer_plain(*spec.pattern, a, projection != Projection::NONE);
                if (!a.result->refusal) {
                    SelectionEntry s;
                    s.answer = selection_without_pass(SelectionPass::NOT_STARTED,
                                                      a.result->anchors->total);
                    // what a long pattern's selection counts: its supported walks (§19.7)
                    s.answer.tested.unit = Unit::PATHS;
                    s.answer.selected.unit = Unit::PATHS;
                    s.support = retrieval->description().support;
                    s.access = retrieval->selection_access();
                    s.memory_bytes = retrieval->memory_peak();
                    a.selection = std::move(s);
                }
            }
        }
        answered.push_back(std::move(a));
    }
    if (hooks && hooks->work_done_hook)
        hooks->work_done_hook();
    const double elapsed_ms = deadline.elapsed_ms();

    // the finalisation: assembling the answer, read against the deadline every
    // kDeliveryStride objects (the writing and the compression then check it too)
    uint64_t objects = 0;
    for (size_t i = 0; i < req.patterns.size(); ++i) {
        Answered &a = answered[i];
        objects += 1 + a.results.size();
        if (a.labels)
            objects += a.labels->result_fields.size();
        if (delivery && objects >= kDeliveryStride) {
            delivery->check();
            objects = 0;
        }
        Json::Value entry = entry_json(req.patterns[i], a.result ? &*a.result : nullptr, mode,
                                       strand_stated, std::move(a.results), a.released,
                                       req.long_paths, fraction ? &*fraction : nullptr);
        // the selection (before the projection's labels, whose withheld, cut and stop come
        // after the pass's)
        if (a.selection)
            put_selection(&entry, *a.selection, mode);
        const bool projected = a.labels.has_value();
        if (a.labels) {
            const RetrievalCounters counters = a.labels->counters;
            apply_labels(&entry, std::move(*a.labels), mode);
            put_retrieval_counters(&entry, counters, a.label_paths);
        }
        if (a.selection && !entry.isMember("error")) {
            // the rows the account refused: the selection's (phase "selection") first, then the
            // projection's (§14.4), each once per phase
            Json::Value refused = std::move(a.selection->answer.rows_refused);
            if (entry["rows_refused"].isArray()) {
                for (Json::Value &row : entry["rows_refused"]) {
                    refused.append(std::move(row));
                }
            }
            entry["rows_refused"] = std::move(refused);
        }
        if (!projected && req.annotation_named && !read_labels && !has_predicate
                && entry.isMember("notes")) {
            // the labels were asked for (or bounded) and none are read here: said, not
            // ignored (mode count, or output.labels "none")
            entry["notes"].append(kNoteAnnotationNotRead);
        }
        if (a.selection && entry.isMember("notes")) {
            // a predicate's notes, after the others (§19.10): the pass read without the
            // budget-aware decode (when no projection said so), a constant normal form, a
            // projection field named and no projection built (mode count, or labels "none"; a
            // predicate answer never says annotation_not_read: its pass read rows)
            Json::Value &notes = entry["notes"];
            bool unbudgeted_stated = false;
            for (const Json::Value &note : notes) {
                unbudgeted_stated |= note.asString() == "annotation_unbudgeted";
            }
            if (a.selection->unbudgeted_read && a.selection->answer.rows && !unbudgeted_stated)
                notes.append("annotation_unbudgeted");
            if (a.constant)
                notes.append(kNotePredicateConstant);
            if (req.projection_named && !projected)
                notes.append(kNoteProjectionNotRead);
        }
        entries.append(std::move(entry));
    }

    Json::Value out;
    out["pattern_contract_version"] = kPatternContractVersion;
    out["mode"] = to_string(mode);
    if (mode == Mode::COUNT) {
        out["output"] = Json::Value();
    } else if (projection == Projection::PREDICATE_ONLY) {
        out["output"]["labels"] = kLabelsPredicateOnly;
        out["output"]["occurrences"] = req.retrieval.occurrences;
    } else if (!read_labels) {
        out["output"]["labels"] = "none";
    } else {
        out["output"]["labels"] = "all";
        out["output"]["occurrences"] = req.retrieval.occurrences;
    }
    if (has_predicate) {
        // the request's predicate as bound to this index (§19.10): its normal form, the names
        // of its lists and of them the index's columns, the others (a typo shows here, never as
        // an absence), whether it holds on a context carrying none of its labels; null where
        // the binding stopped (every entry says how: stop {selection, time | max_memory}).
        // Written under the answer's delivery check
        auto check = [delivery]() {
            if (delivery)
                delivery->check();
        };
        Json::Value p;
        if (bound) {
            p["normal_form"] = bound->normal_form_json(check);
            p["names"] = uint_json(bound->num_names());
            p["known"] = uint_json(bound->num_known());
            p["unknown_labels"] = bound->unknown_labels_json(check);
            p["vacuous"] = bound->vacuous();
            volume.settle(bound->text_bytes());
        } else {
            p["normal_form"] = Json::Value();
            p["names"] = uint_json(req.predicate->names.size());
            p["known"] = Json::Value();
            p["unknown_labels"] = Json::Value();
            p["vacuous"] = Json::Value();
        }
        p["scope"] = kPredicateScope;
        p["strands"] = retrieval->selection_strands();
        out["predicate"] = std::move(p);
    }

    Json::Value index;
    index["index_ns"] = identity ? string_or_null(identity->name) : Json::Value();
    index["index_fp"] = identity ? string_or_null(identity->fp) : Json::Value();
    index["release"] = release;
    index["k"] = uint_json(k);
    index["graph_mode"] = to_string(support.mode);
    index["alphabet"] = support.alphabet;
    index["strand_stated"] = support.strand_stated;
    if (fraction) {
        // in the answers on a graph without its mask only (an answer on a masked graph does not
        // state it: its counting is exact): what the counts are and the dummy fraction the
        // estimates rest on
        index["counting"] = kCountingUpperBound;
        index["dummy_fraction"] = dummy_fraction_json(*fraction);
    }
    out["index"] = std::move(index);

    Json::Value l;
    l["max_contexts"] = uint_json(req.request.max_contexts);
    l["max_anchors"] = uint_json(req.request.max_anchors);
    l["max_steps"] = uint_json(req.max_steps);
    l["time_budget_ms"] = number_json(req.time_budget_ms);
    l["finalize_reserve_ms"] = number_json(limits.finalize_ms);
    l["min_information_bits"] = number_json(limits.min_information_bits);
    l["max_patterns"] = uint_json(limits.max_patterns);
    l["stop_at_threshold"] = req.request.stop_at_threshold;
    if (read_labels || has_predicate) {
        // the annotation limits, in the answers that read annotation only; with a predicate in
        // every mode (its selection reads rows under the account)
        const RetrievalLimits &r = req.retrieval;
        l["max_labels_per_anchor"] = uint_json(r.max_labels_per_anchor);
        l["max_annotation_work"] = uint_json(r.max_annotation_work);
        l["max_memory_mb"] = uint_json(r.max_memory_bytes >> 20);
        l["max_labels"] = uint_json(r.max_labels);
        l["max_occurrences_per_label"] = uint_json(r.max_occurrences_per_label);
        l["allow_unbudgeted_annotation"] = r.allow_unbudgeted;
    }
    if (req.long_paths) {
        // in the answers that ask for paths only
        l["long_search"] = kLongSearchPaths;
        l["max_paths"] = uint_json(req.request.max_paths);
        if (read_labels)
            l["require_support"] = req.require_verified ? kSupportVerified : kSupportIntersection;
    }
    if (has_predicate) {
        // in the answers with a predicate only (§19.10): the selection's caps (effective), the
        // server's cap on the names, and predicate_strands as requested
        l["max_predicate_contexts"] = uint_json(req.max_predicate_contexts);
        l["max_predicate_work"] = uint_json(req.selection.max_predicate_work);
        l["max_predicate_labels"] = uint_json(limits.max_predicate_labels);
        l["predicate_strands"] = req.selection.either ? kPredicateStrandsEither
                                                      : kPredicateStrandsContext;
    }
    l["clamped"] = std::move(req.clamped);
    out["limits"] = std::move(l);

    Json::Value timing;
    timing["elapsed_ms"] = elapsed_ms;
    out["timing"] = std::move(timing);
    out["patterns"] = std::move(entries);

    if (delivery)
        delivery->check();
    return out;
}

Json::Value pattern_capabilities_json(const AnnotatedDBG *anno_graph,
                                      const PatternLimits &limits, bool multi_graph) {
    Json::Value p;
    p["pattern_contract_version"] = kPatternContractVersion;
    if (multi_graph) {
        // the fan-out over shards with its barriers is not served (§8)
        p["available"] = false;
        p["unavailable_reason"] = "multi_graph_later_increment";
        return p;
    }

    p["modes"] = strings_json({ "count", "all_or_count", "partial" });
    p["default_mode"] = to_string(Mode::ALL_OR_COUNT);
    // the predicate's labels ("predicate_only", with a predicate) are served; the default is
    // "none" (a client sends "predicate_only" explicitly)
    p["projections"] = strings_json({ "none", "all", kLabelsPredicateOnly });
    p["default_projection"] = kDefaultProjection;
    p["projections_later_increment"] = Json::Value(Json::arrayValue);
    // output.occurrences with output.labels "all": placed where the index can place
    p["default_occurrences"] = true;
    p["kinds"] = strings_json({ "dna", "iupac", "protein" });
    p["kinds_later_increment"] = Json::Value(Json::arrayValue);
    // the residues a protein pattern may hold, the genetic codes (NCBI translation table ids)
    // and the default
    Json::Value residues(Json::arrayValue);
    for (char c : std::string(kProteinResidues)) {
        residues.append(std::string(1, c));
    }
    p["protein_residues"] = std::move(residues);
    Json::Value codes(Json::arrayValue);
    for (int id : GeneticCode::ids()) {
        codes.append(id);
    }
    p["genetic_codes"] = std::move(codes);
    p["default_genetic_code"] = GeneticCode::kStandard;
    // a reference, not the rule: the capabilities document a service's MCP tool returns in one
    // piece has a ceiling of 32 KiB (api/python/metagraph/traverse/mcp_tools.py
    // CAPABILITIES_MAX_BYTES), which the document of the mini index nearly filled; the rule is
    // the SPEC's, and the machine-readable part is in the fields above. ASCII only: the writers
    // escape any other byte as \uXXXX
    p["protein_rule"] = kProteinRule;
    p["default_scope"] = to_string(Scope::ANY_OFFSET);
    Json::Value by_mode;
    by_mode["basic"] = strings_json({ "suffix", "any_offset" });
    by_mode["canonical"] = strings_json({ "suffix", "any_offset" });
    // a virtual suffix of the wrapper is a stored prefix (§4.1): not a BOSS range
    by_mode["primary"] = strings_json({ "any_offset" });
    p["scopes_by_graph_mode"] = std::move(by_mode);
    // what a pattern longer than k gets without the option (paths are opt-in: SPEC §12)
    p["long_patterns"] = "anchors_counted";
    // the long_search values served, and the default
    p["long_search"] = strings_json({ kLongSearchAnchors, kLongSearchPaths });
    p["default_long_search"] = kLongSearchAnchors;
    p["strands"] = strings_json({ "both", "forward", "reverse" });
    p["default_strands"] = to_string(Strands::BOTH);
    p["graph_cleaned"] = "unknown";
    p["records_shorter_than_k"] = "not_indexed";
    p["resident_only"] = true;
    Json::Value caps;
    caps["max_contexts"] = uint_json(limits.max_contexts);
    caps["max_anchors"] = uint_json(limits.max_anchors);
    // long_search "paths"
    caps["max_paths"] = uint_json(limits.max_paths);
    // not a request field (caps_rule)
    caps["max_checked_entries"] = uint_json(limits.max_checked_entries);
    caps["max_steps"] = uint_json(limits.max_steps);
    caps["time_budget_ms"] = number_json(limits.max_time_ms);
    caps["min_information_bits"] = number_json(limits.min_information_bits);
    caps["max_patterns"] = uint_json(limits.max_patterns);
    // the labelled retrieval's (output.labels "all")
    caps["max_labels_per_anchor"] = uint_json(limits.max_labels_per_anchor);
    caps["max_annotation_work"] = uint_json(limits.max_annotation_work);
    caps["max_memory_mb"] = uint_json(limits.max_memory_mb);
    caps["max_labels"] = uint_json(limits.max_labels);
    caps["max_occurrences_per_label"] = uint_json(limits.max_occurrences_per_label);
    // a predicate's selection (SPEC §19.2, §19.3); max_predicate_labels is not a request field
    // (caps_rule)
    caps["max_predicate_contexts"] = uint_json(limits.max_predicate_contexts);
    caps["max_predicate_work"] = uint_json(limits.max_predicate_work);
    caps["max_predicate_labels"] = uint_json(limits.max_predicate_labels);
    p["caps"] = std::move(caps);
    // the predicate (SPEC §19.12): the operators served, the predicate_strands values (on
    // CANONICAL and PRIMARY graphs both answer "either": one row serves a k-mer and its reverse
    // complement) and the access the index gives a selection: "rows" (budget-aware, or an
    // unbudgeted annotation without direct access) or "columns" (an unbudgeted annotation with
    // direct access: single cells for at most 16 labels); null while it is not known
    Json::Value predicate;
    Json::Value operators(Json::arrayValue);
    for (const char *op : kPredicateOperators) {
        operators.append(op);
    }
    predicate["operators"] = std::move(operators);
    predicate["strands"] = strings_json({ kPredicateStrandsEither, kPredicateStrandsContext });
    predicate["access"] = Json::Value();
    p["predicate"] = std::move(predicate);
    // the budget of a request that names none: unlike the other caps, below the maximum
    p["default_time_budget_ms"] = number_json(limits.default_time_ms);
    p["finalize_reserve_ms"] = number_json(limits.finalize_ms);
    // the rates of the time kept back for the answer (SPEC §7.6), MB/s, stated as numbers
    Json::Value delivery;
    delivery["build"] = number_json(limits.delivery_build_mbps);
    delivery["compress"] = number_json(limits.delivery_compress_mbps);
    p["delivery_mbps"] = std::move(delivery);
    // which caps are request fields' maxima and which are the server's policy: a reference to
    // the SPEC's section that classifies every cap (as protein_rule)
    p["caps_rule"] = kCapsRule;

    const char *graph_fields[] = { "graph_mode", "k", "alphabet", "strand_stated", "mask",
                                   "counting", "dummy_fraction", "scopes", "placement",
                                   "support", "annotation" };
    for (const char *field : graph_fields) {
        p[field] = Json::Value();
    }
    if (!anno_graph) {
        // the index is loading: what it supports is not known yet
        p["available"] = Json::Value();
        p["unavailable_reason"] = Json::Value();
        return p;
    }

    const DeBruijnGraph &graph = anno_graph->get_graph();
    const GraphSupport support = route_support(graph);
    p["available"] = support.supported;
    p["unavailable_reason"] = support.supported ? Json::Value() : Json::Value(support.reason);
    p["k"] = uint_json(graph.get_k());
    // the engine recognised the representation (only its alphabet is not served, or its
    // mask: alphabet_unsupported, alphabet_untested, mask_invalid)
    const bool recognised = support.supported || support.reason == "alphabet_unsupported"
                                || support.reason == "alphabet_untested"
                                || support.reason == "mask_invalid";
    if (!recognised)
        return p;
    p["graph_mode"] = to_string(support.mode);
    p["alphabet"] = support.alphabet.empty() ? Json::Value() : Json::Value(support.alphabet);
    p["strand_stated"] = support.strand_stated;
    // §4: file | built_at_load | absent: the .edgemask read beside the .dbg, or the same mask
    // built in memory at load (--pattern-build-mask) when there was no file
    p["mask"] = !support.mask_present ? "absent"
              : mask_built_at_load(graph) ? "built_at_load" : "file";
    if (support.supported) {
        // exact with the mask, upper bounds and estimates without it (the rule is SPEC's; no
        // prose here: see protein_rule's note on the document's ceiling); the dummy fraction
        // the estimates rest on, null with the mask
        p["counting"] = support.mask_present ? kCountingExact : kCountingUpperBound;
        if (!support.mask_present) {
            if (auto fraction = dummy_fraction(*anno_graph))
                p["dummy_fraction"] = dummy_fraction_json(*fraction);
        }
    }
    Json::Value scopes(Json::arrayValue);
    for (Scope s : support.scopes) {
        scopes.append(to_string(s));
    }
    p["scopes"] = std::move(scopes);

    // what the annotation gives the labelled retrieval (output.labels "all"): placement and
    // support need BASIC coordinates (and the .seqs mapping for records), the budgeted reads a
    // row-diff annotation (§4.3); support is the best per-label support a path's labels can
    // have on this index
    try {
        const graph::traversal::LabelOracle oracle(*anno_graph);
        const AnnotationDescription d = describe_annotation(oracle, support.mode);
        p["placement"] = d.placement;
        p["support"] = d.support;
        p["annotation"] = d.budgeted ? "budgeted" : "unbudgeted";
        p["predicate"]["access"] = !d.budgeted && oracle.supports_direct() ? "columns" : "rows";
    } catch (const std::exception &e) {
        // stated as unknown (null) rather than failing the capabilities
        logger->warn("[Server] pattern capabilities: the annotation could not be described: {}",
                     e.what());
    }
    return p;
}

std::vector<std::string> pattern_gate_keys() {
    std::vector<std::string> keys(std::begin(kPatternGateKeys), std::end(kPatternGateKeys));
    for (const char *cap : kPatternGateCaps) {
        keys.push_back(std::string("caps.") + cap);
    }
    return keys;
}

Json::Value pattern_traverse_block(Json::Value full) {
    // every key of the full block is kept until the search service reads the full block from
    // GET /pattern/capabilities (SPEC §23); the gate keys keep their values after that too
    full["details"] = kPatternDetails;
    return full;
}

bool write_pattern_answer(const std::string &content,
                          const AnnotatedDBG &anno_graph,
                          const PatternLimits &limits,
                          const std::string &release,
                          const IndexIdentity *identity,
                          const Json::StreamWriterBuilder &builder,
                          std::ostream &out,
                          const std::string &name,
                          const RetrievalHooks *hooks) {
    PatternDelivery delivery;
    try {
        const Json::Value json = parse_pattern_body(content);
        const std::string text = Json::writeString(
                builder, process_pattern_request(json, anno_graph, limits, release, identity,
                                                 &delivery, nullptr, hooks));
        // written by the deadline, as the server's answer must be (else 503 there)
        delivery.check();
        out << text << std::endl;
        return true;
    } catch (const PatternRefusal &e) {
        logger->error("Request in {} refused ({} {}): {}", name, e.status(), e.code(),
                      e.what());
        out << Json::writeString(builder, e.body()) << std::endl;
    } catch (const std::exception &e) {
        // what the server answers 400 {"error"} without a code (an unexpected failure): the
        // same body, and the next request file is still answered (an exception that escaped
        // would abort the run with the later files unanswered)
        logger->error("Request in {} failed: {}", name, e.what());
        Json::Value body;
        body["error"] = e.what();
        out << Json::writeString(builder, body) << std::endl;
    }
    return false;
}

int pattern_graph(Config *config) {
    assert(config);
    assert(config->infbase_annotators.size() == 1);

    std::unique_ptr<AnnotatedDBG> anno_graph = initialize_annotated_dbg(*config);

    // the index identity every answer states, as the server's (--index-name,
    // --index-manifest)
    IndexIdentity identity;
    try {
        identity = index_identity(*config, *anno_graph);
    } catch (const std::exception &e) {
        logger->error("{}", e.what());
        return 1;
    }

    Json::StreamWriterBuilder builder;
    builder["indentation"] = config->output_json ? "" : "  ";
    PatternLimits limits = pattern_limits(*config);
    // the indented text writes about twice the compact one's bytes (AnswerVolume)
    limits.delivery_text_scale = config->output_json ? 1 : 2;
    int status = 0;
    for (const auto &file : config->fnames) {
        std::ifstream in(file);
        if (!in.good()) {
            logger->error("Cannot open request file {}", file);
            return 1;
        }
        const std::string content((std::istreambuf_iterator<char>(in)),
                                  std::istreambuf_iterator<char>());
        if (!write_pattern_answer(content, *anno_graph, limits, config->index_release,
                                  &identity, builder, std::cout, file)) {
            status = 1;
        }
    }
    return status;
}

} // namespace cli
} // namespace mtg
