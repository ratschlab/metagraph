#include "pattern.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string_view>
#include <vector>

#include "common/logger.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "config/config.hpp"
#include "load/load_annotated_graph.hpp"
#include "traverse.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph::pattern;
using graph::AnnotatedDBG;
using graph::DeBruijnGraph;

namespace {

/**
 * The projection of an omitted output.labels in this increment. The design's default is
 * "all" (§7.1), which needs annotation and is a later increment; the owner's instruction for
 * this milestone is that an omitted projection means "none", the one served. It is stated in
 * every answer (`output`) and in the capabilities (`default_projection`), so a client never
 * has to assume it; a later increment that serves "all" changes it together with the
 * capabilities.
 */
constexpr const char *kDefaultProjection = "none";

// Request fields of later increments (§7.1): refused by name, any value (null included),
// rather than reported as unknown, so that the answer says what to wait for
const char *const kLaterIncrementFields[] = {
    "max_paths", "max_labels_per_anchor", "max_annotation_work", "max_memory_mb",
    "max_labels", "max_occurrences_per_label", "allow_unbudgeted_annotation",
    "require_support", "predicate", "max_predicate_contexts", "max_predicate_work",
    "graphs", "genetic_code", "budget_split",
};

// The JSON objects built between two readings of the answer's deadline (§5.3:
// "serialisation every 4,096 objects")
constexpr uint64_t kDeliveryStride = 4096;

PatternRefusal invalid(const std::string &message) {
    return PatternRefusal(400, "invalid_request", message);
}

PatternRefusal later(const std::string &message) {
    return PatternRefusal(400, "later_increment", message);
}

Json::Value uint_json(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

// a number of milliseconds or bits as JSON: an integer when it is one (the flags are integers
// and a client compares them as written), else the double
Json::Value number_json(double x) {
    if (x >= 0 && x == std::floor(x) && x <= 9007199254740991.0)
        return uint_json(static_cast<uint64_t>(x));
    return Json::Value(x);
}

Json::Value strings_json(std::initializer_list<const char*> values) {
    Json::Value a(Json::arrayValue);
    for (const char *v : values) {
        a.append(v);
    }
    return a;
}

/**
 * Strict access to one JSON object of the request, as /traverse's: every field read is
 * remembered, and finish() refuses the first one nothing read ("unknown field"), after the
 * known ones were checked — the guarantee rule: a field this increment does not know is never
 * ignored.
 */
class Fields {
  public:
    Fields(const Json::Value &value, std::string path) : v_(value), path_(std::move(path)) {
        if (!v_.isObject())
            throw invalid(path_ + ": expected an object");
    }

    bool has(const std::string &key) {
        seen_.insert(key);
        return v_.isMember(key);
    }
    const Json::Value& raw(const std::string &key) {
        seen_.insert(key);
        return v_[key];
    }
    std::string path(const std::string &key) const { return path_ + "." + key; }

    void finish() const {
        for (const std::string &name : v_.getMemberNames()) {
            if (!seen_.count(name))
                throw invalid(path_ + ": unknown field '" + name + "'");
        }
    }

  private:
    const Json::Value &v_;
    std::string path_;
    std::set<std::string> seen_;
};

// One pattern of the request: parsed, or refused in its slot with bad_alphabet
struct PatternSpec {
    Json::Value id;                      // the string given, or null
    PatternKind kind = PatternKind::DNA;
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
};

void note_clamped(Json::Value *clamped, const char *field, Json::Value requested,
                  Json::Value effective) {
    Json::Value c;
    c["field"] = field;
    c["requested"] = std::move(requested);
    c["effective"] = std::move(effective);
    clamped->append(std::move(c));
}

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

std::string string_field(Fields &f, const char *key, const std::string &def) {
    if (!f.has(key))
        return def;
    const Json::Value &v = f.raw(key);
    if (!v.isString())
        throw invalid(f.path(key) + ": expected a string");
    return v.asString();
}

/**
 * The request as this increment serves it (§7.1): the first error wins, in
 * this order — a later-increment or resident-only field (named, whatever its value), the
 * patterns, mode, output, scope, strands, stop_at_threshold, the caps, the time budget — and
 * a field nothing read is refused last, as /traverse's Strict refuses it.
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
                throw later("request." + name + ": '" + name + "' in a later increment");
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
        if (p.has("protein"))
            throw later(p.path("protein") + ": protein patterns in a later increment");
        if (p.has("id")) {
            if (!p.raw("id").isString())
                throw invalid(p.path("id") + ": expected a string");
            spec.id = p.raw("id").asString();
        }
        const bool dna = p.has("dna");
        const bool iupac = p.has("iupac");
        if (dna == iupac)
            throw invalid(path + ": expected exactly one of 'dna', 'iupac'");
        const char *key = dna ? "dna" : "iupac";
        if (!p.raw(key).isString())
            throw invalid(p.path(key) + ": expected a string");
        p.finish();
        spec.kind = dna ? PatternKind::DNA : PatternKind::IUPAC;
        try {
            spec.pattern = Pattern::parse(spec.kind, p.raw(key).asString());
        } catch (const PatternError &e) {
            // the pattern's own error (§7.2): the other patterns are still answered
            spec.error = std::make_pair(e.code(), std::string(e.what()));
        }
        req.patterns.push_back(std::move(spec));
    }

    const std::string mode = string_field(f, "mode", to_string(Mode::ALL_OR_COUNT));
    if (auto m = parse_mode(mode)) {
        req.request.mode = *m;
    } else {
        throw invalid(f.path("mode") + ": expected one of count|all_or_count|partial");
    }

    if (f.has("output")) {
        Fields o(f.raw("output"), f.path("output"));
        if (o.has("labels")) {
            const Json::Value &v = o.raw("labels");
            if (!v.isString())
                throw invalid(o.path("labels") + ": expected a string");
            const std::string labels = v.asString();
            if (labels == "all" || labels == "predicate_only") {
                // even in mode count, where no projection applies: the request names a
                // reading of the annotation this increment cannot do
                throw later(o.path("labels") + ": \"" + labels + "\" (labels read from the "
                            "annotation) in a later increment; this increment serves \"none\"");
            }
            if (labels != "none")
                throw invalid(o.path("labels") + ": expected one of none|all|predicate_only");
        }
        for (const char *key : { "occurrences", "paths" }) {
            if (!o.has(key))
                continue;
            if (!o.raw(key).isBool())
                throw invalid(o.path(key) + ": expected a boolean");
            if (o.raw(key).asBool())
                throw later(o.path(key) + ": true in a later increment");
        }
        o.finish();
    }

    const std::string scope = string_field(f, "scope", to_string(Scope::ANY_OFFSET));
    if (auto s = parse_scope(scope)) {
        req.request.scope = *s;
    } else {
        throw invalid(f.path("scope") + ": expected one of suffix|any_offset (a pattern longer "
                      "than k is searched as 'long' whatever is named)");
    }
    const std::string strands = string_field(f, "strands", to_string(Strands::BOTH));
    if (auto s = parse_strands(strands)) {
        req.request.strands = *s;
    } else {
        throw invalid(f.path("strands") + ": expected one of both|forward|reverse");
    }
    if (f.has("stop_at_threshold")) {
        if (!f.raw("stop_at_threshold").isBool())
            throw invalid(f.path("stop_at_threshold") + ": expected a boolean");
        req.request.stop_at_threshold = f.raw("stop_at_threshold").asBool();
    }

    req.request.max_contexts = capped_integer(f, "max_contexts", limits.max_contexts, 0,
                                              &req.clamped);
    req.request.max_anchors = capped_integer(f, "max_anchors", limits.max_anchors, 0,
                                             &req.clamped);
    req.max_steps = capped_integer(f, "max_steps", limits.max_steps, 1, &req.clamped);
    req.request.min_information_bits = limits.min_information_bits;

    req.time_budget_ms = limits.default_time_ms;
    if (f.has("time_budget_ms")) {
        const Json::Value &v = f.raw("time_budget_ms");
        if (!v.isNumeric())
            throw invalid(f.path("time_budget_ms") + ": expected a number (ms)");
        const double t = v.asDouble();
        // the reserve is inside the budget: a budget not above it leaves no time to work
        if (!(t > limits.finalize_ms)) {
            throw invalid(f.path("time_budget_ms") + ": expected more than the finalisation "
                          "reserve of " + std::to_string(static_cast<uint64_t>(limits.finalize_ms))
                          + " ms");
        }
        if (t > limits.max_time_ms) {
            // lowered to the cap, not to the default, which lies below it
            note_clamped(&req.clamped, "time_budget_ms", v, number_json(limits.max_time_ms));
            req.time_budget_ms = limits.max_time_ms;
        } else {
            req.time_budget_ms = t;
        }
    }

    f.finish();
    return req;
}

std::string support_message(const GraphSupport &support) {
    if (support.reason == "mask_required") {
        return "pattern: the graph was loaded without its dummy-edge mask (.edgemask): without "
               "it every dummy edge would count as a k-mer and no count would be right; give "
               "the graph its mask once with `metagraph transform --mask-dummy <graph>.dbg` "
               "(writes the .edgemask beside the graph; node ids and annotation unchanged), or "
               "pass --pattern-build-mask to server_query or pattern (builds it in memory at "
               "load)";
    }
    if (support.reason == "representation_unsupported") {
        return "pattern: the graph is not a succinct graph (a DBGSuccinct, or a PRIMARY one "
               "wrapped in CanonicalDBG): the pattern lookup narrows BOSS ranges";
    }
    if (support.reason == "primary_unwrapped") {
        return "pattern: a PRIMARY graph is served only wrapped in CanonicalDBG";
    }
    if (support.reason == "alphabet_unsupported") {
        return "pattern: the graph's alphabet '" + support.alphabet + "' is neither $ACGT nor "
               "$ACGTN";
    }
    return "pattern: the graph is not supported (" + support.reason + ")";
}

Json::Value count_json(const Count &c) {
    Json::Value v;
    // the conservative number: the count, or the lower bound; null when nothing is known
    v["value"] = c.relation == Relation::UNKNOWN ? Json::Value() : uint_json(c.value);
    v["relation"] = to_string(c.relation);
    v["unit"] = to_string(c.unit);
    if (c.relation == Relation::BOUNDS) {
        v["lower"] = uint_json(c.lower);
        v["upper"] = uint_json(c.upper);
    }
    return v;
}

// by_strand (keys +, -, both) on a BASIC graph, by_orientation (forward, reverse,
// palindromic) on the others, where no strand is known (§3, "Strand")
void put_orientations(Json::Value *count, const std::map<Orientation, Count> &by,
                      bool strand_stated) {
    Json::Value o(Json::objectValue);
    for (const auto &[orientation, c] : by) {
        o[strand_stated ? strand_key(orientation) : orientation_key(orientation)] = count_json(c);
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

Json::Value reason_json(const char *reason) {
    Json::Value r;
    r["reason"] = reason;
    return r;
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
 * engine withheld nothing, and moved into the entry.
 */
Json::Value entry_json(const PatternSpec &spec, const Result *result, Mode mode,
                       bool strand_stated, Json::Value results) {
    Json::Value e;
    e["id"] = spec.id;
    e["kind"] = to_string(spec.kind);
    if (spec.error) {
        e["error"] = error_json(spec.error->first, spec.error->second);
        return e;
    }
    assert(spec.pattern && result);
    e["pattern"] = spec.pattern->text();
    e["length"] = uint_json(spec.pattern->length());
    e["information_bits"] = result->information_bits;
    e["anchor_information_bits"] = result->anchor_information_bits
            ? Json::Value(*result->anchor_information_bits) : Json::Value();
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

    Json::Value counts;
    if (result->contexts) {
        Json::Value c = count_json(result->contexts->total);
        c["suffix"] = count_json(result->contexts->suffix);
        Json::Value by_offset(Json::objectValue);
        for (const auto &[offset, count] : result->contexts->by_offset) {
            by_offset[std::to_string(offset)] = count_json(count);
        }
        c["by_offset"] = std::move(by_offset);
        put_orientations(&c, result->contexts->by_orientation, strand_stated);
        counts["contexts"] = std::move(c);
    } else if (result->anchors) {
        Json::Value a = count_json(result->anchors->total);
        put_orientations(&a, result->anchors->by_orientation, strand_stated);
        counts["anchors"] = std::move(a);
        counts["paths"] = count_json(result->anchors->paths);
    } else {
        throw std::logic_error("pattern: the engine answered a pattern without counts");
    }
    // nothing of the annotation is read in this increment (§7.2): stated, never omitted
    counts["labels"] = count_json(Count::unknown(Unit::LABELS));
    counts["occurrences"] = count_json(Count::unknown(Unit::PLACED_OCCURRENCES));
    e["counts"] = std::move(counts);

    Json::Value work;
    work["ranges_visited"] = uint_json(result->work.ranges_visited);
    work["mask_scans"] = uint_json(result->work.mask_scans);
    work["steps"] = uint_json(result->work.steps);
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
        } else if (results.size() != x.returned) {
            throw std::logic_error("pattern: the engine released "
                                   + std::to_string(results.size()) + " contexts and stated "
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
    e["notes"] = std::move(notes);
    Json::Value timing;
    timing["elapsed_ms"] = result->elapsed_ms;
    e["timing"] = std::move(timing);
    return e;
}

} // namespace


Json::Value PatternRefusal::body() const {
    Json::Value b;
    b["error"] = what();
    b["code"] = code_;
    return b;
}

void PatternDelivery::check() const {
    if (deadline_ && deadline_->respond_expired()) {
        throw PatternRefusal(503, "deadline",
                             "pattern: the answer could not be written within time_budget_ms ("
                             + std::to_string(static_cast<uint64_t>(deadline_->time_budget_ms()))
                             + " ms, the finalisation reserve of "
                             + std::to_string(static_cast<uint64_t>(
                                     deadline_->finalize_reserve_ms()))
                             + " ms included): nothing partial is sent");
    }
}

PatternLimits pattern_limits(const Config &config) {
    PatternLimits limits;
    limits.max_contexts = config.pattern_max_contexts;
    limits.max_anchors = config.pattern_max_anchors;
    limits.max_steps = config.pattern_max_steps;
    limits.default_time_ms = static_cast<double>(config.pattern_default_time_ms);
    limits.max_time_ms = static_cast<double>(config.pattern_max_time_ms);
    limits.finalize_ms = static_cast<double>(config.pattern_finalize_ms);
    limits.min_information_bits = config.pattern_min_information_bits;
    limits.max_patterns = config.pattern_max_patterns;
    return limits;
}

Json::Value parse_pattern_body(const std::string &content) {
    Json::CharReaderBuilder builder;
    std::unique_ptr<Json::CharReader> reader { builder.newCharReader() };
    Json::Value json;
    std::string errors;
    if (!reader->parse(content.data(), content.data() + content.size(), &json, &errors))
        throw invalid("request: not JSON: " + errors);
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
        const std::function<Deadline::Clock::time_point()> &clock) {
    // the request's one deadline starts here, its body parsed (§5.3)
    const std::function<Deadline::Clock::time_point()> now
            = clock ? clock : std::function<Deadline::Clock::time_point()>(&Deadline::Clock::now);
    const Deadline::Clock::time_point start = now();
    const DeBruijnGraph &graph = anno_graph.get_graph();

    // the graph first: no request is answerable on a graph the engine
    // cannot count on, whatever it asks
    const GraphSupport support = PatternSearch::support(graph);
    if (!support.supported)
        throw PatternRefusal(400, support.reason, support_message(support));

    ParsedRequest req = parse_request(json, limits);
    const Mode mode = req.request.mode;

    const Deadline deadline(start, req.time_budget_ms, limits.finalize_ms, now);
    if (delivery)
        delivery->set_deadline(deadline);
    // one budget for the request, spent by the patterns in request order (§5.3, §5.5)
    Budget budget(req.max_steps, deadline);
    const PatternSearch search(graph);

    const size_t k = graph.get_k();
    const bool strand_stated = support.strand_stated;
    const uint64_t num_rows = anno_graph.get_annotator().num_objects();

    /**
     * One released context as its JSON result (§7.2, output.labels none), built in the engine's
     * callback: spelling the k-mer is graph work, done while the engine still reads the
     * clock between releases, so that the finalisation reserve is left to assembling and
     * writing the answer. The row is the annotation row the context's k-mer is annotated in,
     * named without reading it: the stored k-mer's (base_node) on BASIC and wrapped PRIMARY
     * graphs; on a native CANONICAL graph the canonical k-mer's (the annotation key of every
     * route, LabelOracle::key_of), since the other orientation's row carries no labels.
     */
    auto context_json = [&](const Context &c, size_t length) {
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
        return r;
    };

    Json::Value entries(Json::arrayValue);
    std::vector<std::pair<std::optional<Result>, Json::Value>> answered;
    answered.reserve(req.patterns.size());
    for (const PatternSpec &spec : req.patterns) {
        Json::Value results(Json::arrayValue);
        std::optional<Result> result;
        if (spec.pattern) {
            if (mode == Mode::COUNT) {
                result = search.count(*spec.pattern, req.request, budget);
            } else {
                const size_t length = spec.pattern->length();
                result = search.enumerate(*spec.pattern, req.request, budget,
                                          [&](const Context &c) {
                    results.append(context_json(c, length));
                });
            }
        }
        answered.emplace_back(std::move(result), std::move(results));
    }
    const double elapsed_ms = deadline.elapsed_ms();

    // the finalisation: assembling the answer, read against the deadline every
    // kDeliveryStride objects (the writing and the compression then check it too)
    uint64_t objects = 0;
    for (size_t i = 0; i < req.patterns.size(); ++i) {
        auto &[result, results] = answered[i];
        objects += 1 + results.size();
        if (delivery && objects >= kDeliveryStride) {
            delivery->check();
            objects = 0;
        }
        entries.append(entry_json(req.patterns[i], result ? &*result : nullptr, mode,
                                  strand_stated, std::move(results)));
    }

    Json::Value out;
    out["pattern_contract_version"] = kPatternContractVersion;
    out["mode"] = to_string(mode);
    if (mode == Mode::COUNT) {
        out["output"] = Json::Value();
    } else {
        out["output"]["labels"] = "none";
    }

    Json::Value index;
    index["index_ns"] = identity && !identity->name.empty() ? Json::Value(identity->name)
                                                            : Json::Value();
    index["index_fp"] = identity && !identity->fp.empty() ? Json::Value(identity->fp)
                                                          : Json::Value();
    index["release"] = release;
    index["k"] = uint_json(k);
    index["graph_mode"] = to_string(support.mode);
    index["alphabet"] = support.alphabet;
    index["strand_stated"] = support.strand_stated;
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
        // the fan-out over shards with its barriers is a later increment (§8)
        p["available"] = false;
        p["unavailable_reason"] = "multi_graph_later_increment";
        return p;
    }

    p["modes"] = strings_json({ "count", "all_or_count", "partial" });
    p["default_mode"] = to_string(Mode::ALL_OR_COUNT);
    p["projections"] = strings_json({ "none" });
    p["default_projection"] = kDefaultProjection;
    p["projections_later_increment"] = strings_json({ "all", "predicate_only" });
    p["kinds"] = strings_json({ "dna", "iupac" });
    p["kinds_later_increment"] = strings_json({ "protein" });
    p["default_scope"] = to_string(Scope::ANY_OFFSET);
    Json::Value by_mode;
    by_mode["basic"] = strings_json({ "suffix", "any_offset" });
    by_mode["canonical"] = strings_json({ "suffix", "any_offset" });
    // a virtual suffix of the wrapper is a stored prefix (§4.1): not a BOSS range
    by_mode["primary"] = strings_json({ "any_offset" });
    p["scopes_by_graph_mode"] = std::move(by_mode);
    p["long_patterns"] = "anchors_counted";
    p["strands"] = strings_json({ "both", "forward", "reverse" });
    p["default_strands"] = to_string(Strands::BOTH);
    p["graph_cleaned"] = "unknown";
    p["records_shorter_than_k"] = "not_indexed";
    p["resident_only"] = true;
    Json::Value caps;
    caps["max_contexts"] = uint_json(limits.max_contexts);
    caps["max_anchors"] = uint_json(limits.max_anchors);
    caps["max_steps"] = uint_json(limits.max_steps);
    caps["time_budget_ms"] = number_json(limits.max_time_ms);
    caps["min_information_bits"] = number_json(limits.min_information_bits);
    caps["max_patterns"] = uint_json(limits.max_patterns);
    p["caps"] = std::move(caps);
    // the budget of a request that names none: unlike the other caps, below the maximum
    p["default_time_budget_ms"] = number_json(limits.default_time_ms);
    p["finalize_reserve_ms"] = number_json(limits.finalize_ms);
    p["caps_rule"] = "each cap is the maximum of its request field: a larger request value is "
                     "lowered to it and listed in limits.clamped; each is also the field's "
                     "default, except time_budget_ms, whose default is default_time_budget_ms; "
                     "max_patterns is not lowered: a longer list is refused";

    const char *graph_fields[] = { "graph_mode", "k", "alphabet", "strand_stated", "mask",
                                   "scopes", "placement", "support", "annotation" };
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
    const GraphSupport support = PatternSearch::support(graph);
    p["available"] = support.supported;
    p["unavailable_reason"] = support.supported ? Json::Value() : Json::Value(support.reason);
    p["k"] = uint_json(graph.get_k());
    // the engine recognised the representation (only its mask or alphabet is missing)
    const bool recognised = support.supported || support.reason == "mask_required"
                                || support.reason == "alphabet_unsupported";
    if (!recognised)
        return p;
    p["graph_mode"] = to_string(support.mode);
    p["alphabet"] = support.alphabet.empty() ? Json::Value() : Json::Value(support.alphabet);
    p["strand_stated"] = support.strand_stated;
    // §4: file | built_at_load | absent: the .edgemask read beside the .dbg, or the same mask
    // built in memory at load (--pattern-build-mask) when there was no file
    p["mask"] = !support.mask_present ? "absent"
              : mask_built_at_load(graph) ? "built_at_load" : "file";
    Json::Value scopes(Json::arrayValue);
    for (Scope s : support.scopes) {
        scopes.append(to_string(s));
    }
    p["scopes"] = std::move(scopes);

    // what the annotation could give a later increment's labelled retrieval (none of it is
    // served now): placement and support need BASIC coordinates (and the .seqs mapping for
    // records), the budgeted reads a row-diff annotation (§4.3)
    try {
        const graph::traversal::LabelOracle oracle(*anno_graph);
        const bool basic = support.mode == GraphMode::BASIC;
        const bool coords = oracle.has_coordinates();
        const bool records = coords && oracle.coord_to_header();
        p["placement"] = !basic ? "none_canonical" : records ? "record" : coords ? "global"
                                                                                 : "none";
        p["support"] = basic && records ? "record_verified" : "label_intersection";
        p["annotation"] = oracle.decode_charged() ? "budgeted" : "unbudgeted";
    } catch (const std::exception &e) {
        // stated as unknown (null) rather than failing the capabilities
        logger->warn("[Server] pattern capabilities: the annotation could not be described: {}",
                     e.what());
    }
    return p;
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
    int status = 0;
    for (const auto &file : config->fnames) {
        std::ifstream in(file);
        if (!in.good()) {
            logger->error("Cannot open request file {}", file);
            return 1;
        }
        const std::string content((std::istreambuf_iterator<char>(in)),
                                  std::istreambuf_iterator<char>());
        PatternDelivery delivery;
        try {
            const Json::Value json = parse_pattern_body(content);
            const std::string text = Json::writeString(
                    builder, process_pattern_request(json, *anno_graph, *config, &identity,
                                                     &delivery));
            // written by the deadline, as the server's answer must be (else 503 there)
            delivery.check();
            std::cout << text << std::endl;
        } catch (const PatternRefusal &e) {
            logger->error("Request in {} refused ({} {}): {}", file, e.status(), e.code(),
                          e.what());
            std::cout << Json::writeString(builder, e.body()) << std::endl;
            status = 1;
        }
    }
    return status;
}

} // namespace cli
} // namespace mtg
