#ifndef __METAGRAPH_CLI_PATTERN_HPP__
#define __METAGRAPH_CLI_PATTERN_HPP__

/**
 * POST /pattern and `metagraph pattern`: count and extract the graph contexts of short motifs
 * and IUPAC patterns, read their labels (docs/DESIGN-pattern-search.md), and find the paths of
 * patterns longer than k (opt-in). The engine is graph::pattern::PatternSearch
 * (src/graph/alignment/pattern_search.hpp), the labelled retrieval PatternRetrieval
 * (pattern_retrieval.hpp); this file turns their results into the JSON of the route's contract
 * (pattern_contract_version 1):
 *  - modes count, all_or_count and partial; the two retrieval modes with output.labels "none"
 *    (the label-free path, §4.3: contexts with k-mer, instance, offset, strand, node and row
 *    ids, no annotation row read) or "all" (each context's labels, placed where the index can
 *    place them, under the annotation budgets);
 *  - a pattern longer than k: its anchors counted and nothing extracted (long_search "anchors",
 *    the default), or with long_search "paths" the anchors extended into paths (§4.2):
 *    counts.paths, path results with the fields sequence, anchor_kmer, nodes and rows (never
 *    kmer), and with output.labels "all" each path's labels, each with its support
 *    (label_intersection, record_verified) and require_support "record_verified" listing the
 *    verified ones only; or with long_search "supported_paths" (SPEC §20) into the walks some
 *    label supports along their whole length, the annotation read during the extension
 *    (pattern_support.hpp: the trackers and the list of supported paths), counts.supported_paths
 *    beside a plain counts.paths;
 *  - a predicate (SPEC §19): the request's predicate bound once to the index's columns; for a
 *    pattern of L <= k each of its raw contexts (at most max_predicate_contexts) tested by the
 *    selection pass of PatternRetrieval (their rows, with predicate_strands "either" on a BASIC
 *    graph also their reverse complements', under max_predicate_work), the selected ones
 *    returned with output.labels "none", "predicate_only" or "all"; counts.tested and
 *    counts.selected, selection, absence_filter and the top-level predicate block. A pattern
 *    longer than k: its supported paths selected under long_search "supported_paths" (SPEC
 *    §20.9, pattern_supported.hpp); under "anchors" it keeps its anchors' answer, its selection
 *    not_started; with "paths" a predicate is refused (it selects among supported paths, never
 *    among every graph walk). With predicate_scope "motif" (SPEC §25) the predicate is also
 *    asked once of each pattern of L <= k, on the union of its contexts' labels (the entry's
 *    motif);
 *  - one graph per call: a multi-graph server calls it once per (graph, annotation) pair a
 *    request selects (`graphs`, as /search's) and concatenates the answers (server.cpp);
 *  - graphs with their dummy-edge mask (counting "exact") and without it (counting
 *    "upper_bound"): a count is then the bounds [lower, U], U the BOSS entries of its ranges
 *    (source dummies among them), with the additive estimate U x f (f the graph's sampled dummy
 *    fraction, DummyFraction), while every list stays exact (the engine drops the dummies it
 *    releases); a pattern with at most max_checked_entries unchecked candidates has them tested
 *    and its counts exact;
 *  - peptides with the stop '*': a stop codon of the request's genetic code; a table without an
 *    unconditional stop codon matches nothing there, stated in a note.
 * A request field that is not served is refused (400 "later_increment"), never ignored, so that
 * nothing is weakened silently.
 */

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <json/json.h>

#include "graph/alignment/pattern_search.hpp"
#include "pattern_retrieval.hpp"


namespace mtg {

namespace graph {
class AnnotatedDBG;
}

namespace cli {

class Config;
struct IndexIdentity;

/**
 * A refusal of a whole /pattern request: the HTTP status and the body {"error", "code"}
 * (400 invalid_request, later_increment and the graph support reasons of
 * route_support, alphabet_untested and mask_invalid among them; 503 deadline). The server
 * answers it as is (HttpError); the CLI writes the same body and exits 1. A refusal of one
 * pattern is not this: it is the `error` of that pattern's slot in a 200 answer.
 */
class PatternRefusal : public std::runtime_error {
  public:
    PatternRefusal(int status, std::string code, const std::string &message)
          : std::runtime_error(message), status_(status), code_(std::move(code)) {}

    int status() const { return status_; }
    const std::string& code() const { return code_; }
    Json::Value body() const;

  private:
    int status_;
    std::string code_;
};

// The refusal of a request that cannot be read (400 invalid_request), as the strict readers of
// the request's parsers throw it (StrictObject<InvalidPatternRequest>)
struct InvalidPatternRequest {
    PatternRefusal operator()(const std::string &message) const {
        return PatternRefusal(400, "invalid_request", message);
    }
};

/**
 * The server's caps of a /pattern request (the --pattern-* flags; the CLI applies the same
 * ones, so that both answer alike). Each is a request field's maximum: a larger request value
 * is lowered to it and listed in limits.clamped (the /traverse convention), never refused,
 * never silently kept; and each is the field's default, except the time budget, whose default
 * (60 s) lies under its maximum (600 s) (§5.3), and max_anchors and max_labels, whose defaults
 * (default_max_anchors, default_max_labels) are at most their maxima (equal by default).
 */
struct PatternLimits {
    uint64_t max_contexts = 10'000;
    uint64_t max_anchors = 1'000;
    // max_anchors of a request that names none (--pattern-default-max-anchors): at most
    // max_anchors, the server refuses a larger one at start-up
    uint64_t default_max_anchors = 1'000;
    // long_search "paths": the retrieval threshold on the completed paths of a pattern longer
    // than k (all_or_count) and partial's cap on them
    uint64_t max_paths = 1'000;
    uint64_t max_steps = 100'000'000;
    double default_time_ms = 60'000;
    double max_time_ms = 600'000;
    // the finalisation reserve inside the time budget (§5.3): work stops at least this long
    // before the request's deadline so that the answer can still be written by it (longer by
    // the estimated finalisation of what the answer buffers: delivery_* below)
    double finalize_ms = 250;
    double min_information_bits = 24;
    // a graph without its dummy-edge mask only: a pattern whose unchecked candidates number at
    // most this has each of them tested (k - 1 steps each), its counts then exact
    // (graph::pattern::Request::max_checked_entries); 0 tests none. Not a request field: the
    // server's policy (--pattern-max-checked-entries), stated in caps
    uint64_t max_checked_entries = graph::pattern::kDefaultMaxCheckedEntries;
    // patterns per request: above it the request is refused (a list is not cut)
    uint64_t max_patterns = 16;
    // the labelled retrieval (output.labels "all"; §4.3, §5.3): the labels kept per row, the
    // annotation work (the oracle's units), the request's memory account (MiB), and partial's
    // lists of labels per pattern and of occurrences per label
    uint64_t max_labels_per_anchor = 64;
    uint64_t max_annotation_work = 100'000'000;
    uint64_t max_memory_mb = 256;
    uint64_t max_labels = 1'000;
    // max_labels of a request that names none (--pattern-default-max-labels): at most
    // max_labels, likewise
    uint64_t default_max_labels = 1'000;
    uint64_t max_occurrences_per_label = 16;
    // long_search "supported_paths" (SPEC §20.3): the ceiling of a pattern's row cache (MiB),
    // which gets a quarter of what the request's account has left when the pattern's
    // extension begins, at most this. Not a request field: the server's policy
    // (--pattern-row-cache-mb), stated in caps
    uint64_t row_cache_mb = 64;
    // a predicate's selection (SPEC §19.2): the raw contexts a pattern's selection may test
    // (its compute admission) and the selection's work per request (the oracle's units),
    // request fields' maxima like the others; and the names a predicate may list, the server's
    // policy (a larger predicate is refused, predicate_too_large)
    uint64_t max_predicate_contexts = 100'000;
    uint64_t max_predicate_work = 100'000'000;
    uint64_t max_predicate_labels = 10'000;
    // not a cap: the annotation reads under the deadline are decoded in chunks of about this
    // many ms (the server's --traverse-chunk-target-ms, as /traverse's reads); 0: one piece
    double chunk_target_ms = 50;
    // not caps: the finalisation of what the answer buffers (AnswerVolume): the rates (MB/s) at
    // which its JSON text is assumed to be built and written, and compressed
    // (--pattern-delivery-build-mbps, --pattern-delivery-compress-mbps), and the text written
    // per byte of compact JSON (1 on the server; the CLI's indented text 2). The work stops
    // finalize_ms plus that estimate before the deadline. 0 or infinity: the reserve alone
    double delivery_build_mbps = 10;
    double delivery_compress_mbps = 50;
    double delivery_text_scale = 1;
};

PatternLimits pattern_limits(const Config &config);

/**
 * The alphabet half of the route's support decision, a pure function of the graph's BOSS
 * alphabet: "" for "$ACGT" (served), "alphabet_untested" for "$ACGTN" (the engine supports it,
 * but the route does not serve it until a DNA5 build passes the pattern tests),
 * "alphabet_unsupported" for any other.
 */
std::string alphabet_refusal(const std::string &alphabet);

/**
 * Whether /pattern and `metagraph pattern` serve |graph|: the engine's PatternSearch::support,
 * narrowed by the route's own reasons, "alphabet_untested" (alphabet_refusal) and
 * "mask_invalid" (a loaded mask that marks an edge with W = $ valid, found once at load:
 * check_mask_at_load). Its reason is the 400 refusal's code and the capabilities'
 * unavailable_reason. The engine itself keeps serving $ACGTN (its tests). A graph without a
 * mask is served with counting "upper_bound"; no configuration answers mask_required.
 */
graph::pattern::GraphSupport route_support(const graph::DeBruijnGraph &graph);

/**
 * f, the fraction of real k-mers among the entries of a succinct graph that a pattern can
 * count: on a graph served without its dummy-edge mask a count is an upper bound U (the BOSS
 * entries of its ranges, the source dummies among them), stated with the additive estimate U x
 * f. Sampled by the engine (graph::pattern::sample_real_fraction: 10,000 entries with W != $
 * drawn with a seed fixed by the graph's number of edges, the same value in every process;
 * Wilson's 95% interval), once per graph by the route (dummy_fraction).
 */
using DummyFraction = graph::pattern::RealFraction;

/**
 * The dummy fraction of the graph |anno_graph| serves when the pattern search counts on it
 * without a dummy-edge mask (counting "upper_bound"); nullopt when it has its mask (counting
 * "exact"), or is not a succinct graph. Sampled once per graph and kept: in the loading thread
 * (sample_dummy_fraction_at_load; every graph a server serves is loaded that way, so no request
 * samples), or at the first call for a graph not loaded that way. Thread-safe: the sample is
 * drawn outside the registry's lock (two first calls on one graph may both draw it, the same
 * draws, and one is kept), so a reader of another graph's fraction never waits for it.
 */
std::optional<DummyFraction> dummy_fraction(const graph::AnnotatedDBG &anno_graph);

/**
 * In the loading thread of a graph the pattern search serves (async_load_critical_dbg): when
 * |graph| is a succinct graph without its dummy-edge mask, samples its dummy fraction and keeps
 * it for dummy_fraction(), logging f, its interval and the time; nothing otherwise.
 * |stdout_reserved| as for build_mask_at_load (the CLI logs at trace level).
 */
void sample_dummy_fraction_at_load(const std::shared_ptr<graph::DeBruijnGraph> &graph,
                                   bool stdout_reserved = false);

// The JSON of a dummy fraction: {value, interval: [lower, upper], samples, source: "sampled"}
// ("counted" for an exact one, which the route never states)
Json::Value dummy_fraction_json(const DummyFraction &fraction);

/**
 * What the transport of one /pattern answer needs from its processing: the request's
 * deadline, known once the request is parsed. The answer is written (serialised, compressed)
 * under check(), which refuses with 503 "deadline" once time_budget_ms has passed, so that
 * nothing partial is sent as if whole (§5.3). Used by one request's thread only.
 */
class PatternDelivery {
  public:
    void set_deadline(const graph::pattern::Deadline &deadline) { deadline_ = deadline; }
    // throws PatternRefusal(503, "deadline") once the deadline's respond time passed; does
    // nothing before a deadline was set (the request was refused before it was parsed); and
    // graph::pattern::Aborted when the abort predicate answers true
    void check() const;

    /**
     * The server's "the client left, or the server stops", so that /pattern does not run to its
     * deadline for a caller that is gone: process_pattern_request gives it to the request's
     * Budget (Budget::set_abort, read at every clock reading of the work), and check() asks it
     * while the answer is assembled, written and compressed. Its answer true throws
     * graph::pattern::Aborted; the server then writes nothing. Unset (the CLI, tests): never
     * asked.
     */
    void set_abort(std::function<bool()> aborted) { abort_ = std::move(aborted); }
    const std::function<bool()>& abort() const { return abort_; }

  private:
    std::optional<graph::pattern::Deadline> deadline_;
    std::function<bool()> abort_;
};

// The request body as JSON: one RFC 8259 JSON text with unique member names and nothing after
// it, nested at most 1,000 deep; anything else is refused as invalid_request (the other routes'
// generic parse error carries no code, and their parser is more lenient)
Json::Value parse_pattern_body(const std::string &content);

/**
 * The checks of process_pattern_request that need no graph (SPEC §5 steps 5 to 9: the body
 * is an object, its fields, the patterns, the predicate), thrown as process_pattern_request
 * throws them (400 invalid_request, later_increment, genetic_code_unknown,
 * predicate_too_large); nothing else is done. A multi-graph server asks it once before its
 * pairs (SPEC §24.1), so that such a refusal is the request's and not every pair's; each pair
 * then starts from its graph's support (step 4) and refuses on its own.
 */
void validate_pattern_request(const Json::Value &json, const PatternLimits &limits);

/**
 * Answers one /pattern request (the JSON contract of pattern_contract_version 1) on the
 * single graph of |anno_graph| under |limits|, stating |release|. The deadline starts on
 * entry (the body was parsed), read from |clock| (the steady clock when null; injectable so
 * that tests stop at a chosen instant); |delivery|, when given, receives it for the writing
 * of the answer. The work stops finalize_ms before the deadline, earlier by the estimated
 * time to write what the answer holds (AnswerVolume, limits.delivery_*). |identity| (may be
 * null) states the index in the answer's `index`. Throws PatternRefusal for a
 * whole-request refusal, 503 "deadline" included when the answer could not be assembled by
 * the deadline; refused patterns are answered in their slots. Reads
 * annotation rows only for output.labels "all" or "predicate_only" in a retrieval mode, and
 * for a predicate's selection in every mode (PatternRetrieval);
 * |hooks| (tests): a record mapping instead of the index's, a hook on every read.
 */
Json::Value process_pattern_request(
        const Json::Value &json,
        const graph::AnnotatedDBG &anno_graph,
        const PatternLimits &limits,
        const std::string &release,
        const IndexIdentity *identity = nullptr,
        PatternDelivery *delivery = nullptr,
        const std::function<graph::pattern::Deadline::Clock::time_point()> &clock = nullptr,
        const RetrievalHooks *hooks = nullptr);

// The same with the caps (pattern_limits) and the release (--index-release) of |config|
Json::Value process_pattern_request(const Json::Value &json,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const Config &config,
                                    const IndexIdentity *identity = nullptr,
                                    PatternDelivery *delivery = nullptr);

/**
 * The full `pattern` block, of GET /pattern/capabilities and GET /capabilities (SPEC §10.2,
 * §23; GET /traverse/capabilities carries pattern_traverse_block of it): the contract
 * version, whether this server can answer /pattern (available: true | false | null while the single index loads, with the
 * reason when false), the modes, projections, kinds, scopes and strands, the caps and the
 * finalisation reserve with the delivery rates (delivery_mbps; caps_rule and protein_rule are
 * references to the SPEC), and what the graph is (mode, k, alphabet, mask, counting: exact with
 * the mask, upper_bound without it, and then its dummy_fraction) and what its annotation gives
 * the labelled retrieval (placement, support, annotation: budgeted or unbudgeted) and a
 * predicate's selection (projections with "predicate_only", the caps max_predicate_contexts,
 * max_predicate_work and max_predicate_labels, and the predicate object: operators, strands,
 * access). |anno_graph| is null while the index loads; |multi_graph| servers answer only that
 * they are not served yet.
 */
Json::Value pattern_capabilities_json(const graph::AnnotatedDBG *anno_graph,
                                      const PatternLimits &limits, bool multi_graph);

/**
 * The keys of the pattern block that GET /traverse/capabilities keeps with the full block's
 * values for a client that reads that document alone (SPEC-pattern-search.md §23, the gate
 * block): every key such a client gates on or parses, and counting; a cap as "caps.NAME". The
 * other keys are the full block's (GET /pattern/capabilities, GET /capabilities).
 */
std::vector<std::string> pattern_gate_keys();

/**
 * The pattern block of GET /traverse/capabilities from |full| (pattern_capabilities_json):
 * |full| with `details`, the route of the full block ("GET /pattern/capabilities"). Every key
 * of |full| is kept until the search service reads that route (SPEC §23).
 */
Json::Value pattern_traverse_block(Json::Value full);

/**
 * One request of `metagraph pattern` (the text of a request file, |name| in the log) as the
 * CLI answers it on |out|, one JSON text written with |builder|: the answer (true), or the body
 * of a refusal or, for any other failure, the body the server answers 400 without a code
 * ({"error"}) (false). Never throws for a request, so that the next file is still answered.
 * |hooks|: as process_pattern_request's (tests).
 */
bool write_pattern_answer(const std::string &content,
                          const graph::AnnotatedDBG &anno_graph,
                          const PatternLimits &limits,
                          const std::string &release,
                          const IndexIdentity *identity,
                          const Json::StreamWriterBuilder &builder,
                          std::ostream &out,
                          const std::string &name = "the request",
                          const RetrievalHooks *hooks = nullptr);

// `metagraph pattern -i GRAPH -a ANNOTATION REQUEST.json ...`: one JSON answer per request
// file on stdout, the server's (a refusal's body, exit 1); for tests and offline use
int pattern_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_HPP__
