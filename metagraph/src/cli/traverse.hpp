#ifndef __METAGRAPH_CLI_TRAVERSE_HPP__
#define __METAGRAPH_CLI_TRAVERSE_HPP__

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <json/json.h>

#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"


namespace mtg {

namespace graph {
class AnnotatedDBG;
}

namespace cli {

class Config;
class Attempt;

// A request error: reported as HTTP 400 by the server, exit code 1 by the CLI.
class InvalidRequest : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

// Label-change cost as requested, by label NAME (converted per seed to LabelChangeCost
// over the seed's label dictionary).
struct CostSpec {
    enum Model { FORBID, CONSTANT, TABLE };
    Model model = FORBID;
    double value = 0;
    double default_cost = graph::traversal::kInfiniteLoss;
    std::vector<std::tuple<std::string, std::string, double>> entries;
};

struct TraverseRequest {
    std::string release;
    std::string graph;        // multi-graph mode: index name (resolved by the server)
    std::string graph_path;   // optional disambiguation when a name spans several shards
    // A ledger-managed attempt (DESIGN-traverse-graphlet.md §14.1: the frozen wire contract's
    // request fields; traverse_attempts.hpp): its id, unique per server process while it runs
    // or is retained, and the budget and locus it is charged to, echoed in `usage` only
    std::string attempt_id;
    std::string budget_id;
    std::string locus_id;
    // not started after this instant (Unix epoch ms, the server's clock): enforced by the
    // caller at handler start (the server: 409 state expired; the CLI: the same body, exit 1),
    // echoed in usage
    std::optional<uint64_t> not_after_ms;
    std::vector<graph::traversal::Seed> seeds;
    graph::traversal::Strategy strategy;
    CostSpec cost;
    // summary | tree | full | graphlet: the last is the per-seed JSON summary of
    // DESIGN-traverse-graphlet.md §3 with the lossless MGT text embedded (spec §7.5)
    std::string detail = "full";
    bool timing = true;
};

// The MGT format version written by `detail: graphlet` (H record, capabilities)
constexpr int kGraphletFormatVersion = 1;
// The traversal algorithm every /traverse response names (`algorithm_version`), also stated
// by both capabilities routes, so that a client can tell before asking which walk it gets
constexpr const char *kTraverseAlgorithmVersion = "traverse-0.2";
// What the server offers beyond the base contract (capabilities.feature_level), monotonic:
// each pass that adds capabilities fields or routes bumps it by one (SPEC §10.3 maps every
// level). 2: attempts and the client-gone stop; 3: not_after_ms, per-graph identity and
// GET /traverse/capabilities?graph=, GET /capabilities, algorithm_version in the
// capabilities, attempts.hard_cap_ms (allowance_ms an integer), deadline_check and the
// chunked deadlines, compression_level and the delivery reserve; 4 (the efficiency pass): the
// server's maxima of the budgets (the probe's max_memory_mb / max_work_units, clamped like the
// time cap), the row-diff path cache (the probe's decode_cache) and the delivery reserve's
// calibrated starting estimates (delivery_reserve.calibration)
constexpr int kTraverseFeatureLevel = 4;

/**
 * What identifies the index a response was computed on (DESIGN-traverse-graphlet.md
 * §3.1). Labels are joined across retrievals by (kind, column, seq_id), which is only
 * meaningful on the same index: k, regime, alphabet and release do not identify one
 * (two annotations of one graph agree on all four).
 *   name     --index-name, for humans and routing, NOT identity ("" = unset)
 *   fp       sha256 over the canonical file list of the index bundle's manifest
 *            (--index-manifest): the identity ("" = no manifest: joins are unverifiable)
 *   meta_fp  FNV-1a-64 over k, regime, alphabet, rows and the ordered column names: a
 *            negative check only (a mismatch proves different indexes, equality proves
 *            nothing: swapped column memberships keep it)
 */
struct IndexIdentity {
    std::string name;
    std::string fp;
    std::string meta_fp;
};

// [A-Za-z0-9._-]+ (an MGT token and safe in file names)
bool valid_index_name(const std::string &name);
// index_meta_fp of the index behind |oracle| (cost: one pass over the column names)
std::string index_meta_fingerprint(const graph::traversal::LabelOracle &oracle);
/**
 * The files an index of |graph| and |annotation| is loaded from: the two, and the sidecars the
 * loader reads beside them when they exist — the row-diff anchors and fork successors
 * (<graph>.anchors, <graph>.rd_succ), the succinct graph's dummy-edge mask (<graph without
 * .dbg>.edgemask), the coordinates of a column annotation (<annotation>.coords) and the
 * sequence headers of a coordinate annotation (<annotation without .<type>.annodbg>.seqs) —
 * as scripts/traversal/index_manifest.py collects them.
 */
std::vector<std::string> index_bundle_files(const std::string &graph,
                                            const std::string &annotation);
/**
 * Read an index manifest, validate it and return its fingerprint. The manifest is a JSON
 * object with `files: [{path, size, sha256}, ...]` (paths relative to the manifest,
 * '/'-separated, unique; sha256 64 lowercase hex); other keys (builder, inputs) are
 * metadata and not part of the identity. The fingerprint is the lowercase hex sha256 of
 * the lines "<path>\t<size>\t<sha256>\n" in ascending byte order of path; an `index_fp`
 * the manifest states must equal it. Every file of |loaded| (index_bundle_files) that exists
 * on disk must be listed (by base name) with its size, and the manifest must list no graph
 * (*dbg) or annotation (*.annodbg) file other than those of |loaded|: a manifest of another
 * bundle, or of a directory holding several, must not lend its identity. Throws
 * std::runtime_error naming the problem. |stated_name|, when given, receives the manifest's
 * own `index_ns` (metadata; "" when it states none): the server's name for the index comes
 * from its configuration, and a different one in the manifest is only logged.
 */
std::string index_manifest_fingerprint(const std::string &manifest_path,
                                       const std::vector<std::string> &loaded,
                                       std::string *stated_name = nullptr);
// lowercase hex SHA-256 (FIPS 180-4) of |data|
std::string sha256_hex(std::string_view data);
// the identity of the index |config| names (--index-name, --index-manifest against the
// loaded files -i / -a) with |anno_graph|'s meta fingerprint; throws std::runtime_error
IndexIdentity index_identity(const Config &config, const graph::AnnotatedDBG &anno_graph);

struct ResolveRequest {
    std::string graph;
    std::string graph_path;
    std::string sequence;
    graph::traversal::ResolveOptions options;
    bool select = false;
    graph::traversal::SelectionPolicy policy;
    std::string run_format = "intervals"; // intervals | rle
    uint64_t max_query_bp = 1'000'000;
};

// Strict parsing: unknown keys, wrong types and out-of-range values throw InvalidRequest
// naming the offending field. Defaults are filled in.
TraverseRequest parse_traverse_request(const Json::Value &json);
ResolveRequest parse_resolve_request(const Json::Value &json);

// Normalized strategy echo (what was actually used), resubmittable verbatim:
// output.detail and output.timing come from the request, not from the Strategy
Json::Value strategy_to_json(const graph::traversal::Strategy &strategy, const CostSpec &cost,
                             const std::string &detail = "full", bool timing = true);
// One entry of `results`. With |detail| "graphlet" this is the summary of
// DESIGN-traverse-graphlet.md §3 WITHOUT the graphlet itself: the body needs the seed,
// the index identity and the final limitations (process_traverse_request adds them).
Json::Value seed_result_to_json(const graph::traversal::SeedResult &result,
                                const graph::traversal::Strategy &strategy,
                                const std::string &detail, bool timing);
Json::Value profile_to_json(const graph::traversal::SupportProfile &profile,
                            const graph::traversal::SeedSelection *selection,
                            const std::string &run_format);
// |identity| null: no name, no manifest, meta_fp computed here
Json::Value capabilities_to_json(const graph::traversal::LabelOracle &oracle,
                                 const std::string &release,
                                 const IndexIdentity *identity = nullptr);

// The width at which the delivery model prices an MGT float when nothing makes it wider.
constexpr uint64_t kMgtFloatWidth = 24;
// The widest an MGT float of a result walked under |st| and |cost| can be written (canonical
// MGT writes floats positionally: 1e-300 takes 302 characters), from the costs, the loss
// budget and the time budgets (|requested_time_ms|: before a server clamp, |max_time_ms|: the
// server's cap; 0: none): at least kMgtFloatWidth.
uint64_t mgt_float_width(const graph::traversal::Strategy &st,
                         const graph::traversal::LabelChangeCost &cost,
                         double requested_time_ms = 0, double max_time_ms = 0);
// What one object of a result costs to deliver in |detail| (per object, bytes; upper
// bounds): the output part of the memory budget's model (DESIGN-traverse-graphlet.md §14),
// which process_traverse_request sets as Strategy::delivery, with every MGT float priced at
// |float_width| characters (mgt_float_width)
graph::traversal::DeliveryCosts delivery_costs(const std::string &detail, bool sequences,
                                               uint64_t float_width = kMgtFloatWidth);
// the length the JSON writers here give |s| inside a JSON string, quotes excluded: what
// delivery_costs prices a name's text by
uint64_t json_escaped_size(std::string_view s);

// What the H record states besides the result: the index and the seed's position.
struct GraphletContext {
    uint64_t k = 0;
    std::string regime;
    std::string alphabet;
    size_t seed_index = 0;
    IndexIdentity identity;
};

/**
 * The MGT v1 text of one seed (DESIGN-traverse-graphlet.md §2, spec §7.5), written by a
 * streaming writer straight from |result|. |seed| is the request's seed (its validated,
 * case-mapped sequence is the S record), |result_json| the seed's entry of `results` as
 * returned to the client (its `outcome` and `limitations`, seed and per arm, become the
 * O and K records, so the two never diverge). Every field with a derivation rule is
 * written `*` exactly when the rule reproduces the walker's value; a value the grammar
 * cannot express throws std::logic_error rather than producing a wrong document.
 * |lines|, when given, receives the line count (the Z record's value).
 */
std::string graphlet_text(const graph::traversal::SeedResult &result,
                          const graph::traversal::Seed &seed,
                          const graph::traversal::Strategy &strategy,
                          const GraphletContext &context,
                          const Json::Value &result_json,
                          size_t *lines = nullptr);

/**
 * The MGT v1 token codec (DESIGN-traverse-graphlet.md §2.1-§2.2, the golden vectors of
 * api/python/tests/data/traverse/codec_vectors.tsv). Every value has exactly one valid
 * spelling; the decoders reject every other one (FormatError), so that two readers never
 * disagree on which documents are valid. Exposed for the golden-vector test.
 */
namespace mgt {

class FormatError : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

// x >= 0 (or -0) or +inf: the shortest round-trip digits, expanded positionally
std::string encode_float(double x);
double decode_float(std::string_view token);
// strictly ascending ids; every maximal run of >= 2 written a-b; "." = empty
std::string encode_ranges(const std::vector<uint64_t> &ids);
// |bound|: every id must be below it, checked BEFORE a run is expanded, so that a few
// bytes of a corrupt or hostile token ("0-3000000000") cannot allocate in proportion to
// the range (a reader passes the number of L records)
std::vector<uint64_t> decode_ranges(std::string_view token,
                                    std::optional<uint64_t> bound = std::nullopt);
// against |base| (null: none, explicit only); the shorter of explicit and delta, a tie
// explicit
std::string encode_setexpr(const std::vector<uint64_t> &ids, const std::vector<uint64_t> *base);
std::vector<uint64_t> decode_setexpr(std::string_view token, const std::vector<uint64_t> *base,
                                     std::optional<uint64_t> bound = std::nullopt);
// the free-text last field: % LF CR as %25 %0A %0D, nothing else
std::string pct_escape(std::string_view raw);
std::string pct_unescape(std::string_view token);
// "<prefix_len> <suffix>" of an L record after |previous| (byte prefix cut back to a
// code-point boundary, suffix pct-escaped); decoding returns the full name
std::string front_code(std::string_view previous, std::string_view name);
std::string front_decode(std::string_view previous, std::string_view token);
// strings in MGT (and in the JSON beside it) are UTF-8 (§2.1): well-formed UTF-8 (no
// overlong forms, no surrogates, at most U+10FFFF). A label name that is not is never
// replaced: the seed is refused (spec §6.1 step 4, cause unrepresentable_label_name)
bool valid_utf8(std::string_view s);

// a typed K value: i:<integer> | f:<float> | s:<percent-encoded string> | u (unlimited)
struct KValue {
    enum Type { INTEGER, FLOAT, STRING, UNLIMITED };
    Type type = UNLIMITED;
    uint64_t integer = 0;
    double real = 0;
    std::string string;
};
std::string encode_kvalue(const KValue &value);
KValue decode_kvalue(std::string_view token);

} // namespace mgt

// Server-side caps on POST /traverse, 0 = unlimited. A request that names no labels
// makes the SERVER choose the permitted set, so the request's own bounds no longer
// describe the work it asks for; these do. |max_time_ms| also bounds the derivation
// itself, and a request asking for less is left alone (a clamp only ever lowers a
// budget, never raises one) — see Config::traverse_max_*.
struct TraverseLimits {
    double max_time_ms = 0;        // cap on strategy.bounds.time_budget_ms
    size_t max_seeds = 0;          // cap on |request.seeds|
    uint64_t max_seed_bp = 0;      // cap on the length of one seed
    size_t max_seed_labels = 0;    // cap on strategy.labels.max_seed_labels
    // the maxima of the request's budgets (R16; 0 = off): a larger budget is lowered to it,
    // an omitted one set to it, both echoed in strategy.clamped
    uint64_t max_memory_mb = 0;    // of strategy.bounds.max_memory_mb
    uint64_t max_work_units = 0;   // of strategy.bounds.max_work_units
    // Not a cap: the chunked deadlines (spec §6.8, Config::traverse_chunk_target_ms) — an
    // annotation read a deadline may fall into is decoded in chunks of about this many ms, the
    // deadline checked between them; 0: one piece per read, as before
    double chunk_target_ms = 0;
    // Not a cap: the bound of the request's row-diff path cache (LabelOracle::path_cache,
    // Config::traverse_path_cache_mb); 0: off, every read decodes its rows' whole paths
    uint64_t path_cache_bytes = 0;
};

// The results of a /traverse response written as text, one per seed, as each was built
// (the server: the JSON tree of a seed's result is freed at once, and an attempt knows the
// exact bytes it has to deliver). |active|: the response's `results` are these texts, in
// order; the envelope process_traverse_request returns then has no "results", and
// assemble_traverse_response (server_checks.hpp) writes the whole response from both, byte
// for byte what json_text of the tree with its results would write.
struct ResultTexts {
    std::vector<std::string> texts;
    bool active = false;
};

// Shared by the CLI and the server. |release| is the configured index release id ("" if
// none); |identity| as in capabilities_to_json. |texts|: each seed's result is written as
// compact text once built (ResultTexts), under the attempt's delivery check; null (the CLI,
// tests): the results stay in the returned tree. |attempt|: the server's control of the
// request (traverse_attempts.hpp) — polled between seeds, at the walker's checkpoints and
// while the results are built, and, for a request with attempt_id, the usage it records and
// states and the bound it enforces; null (the CLI, tests): none, and a request with
// attempt_id states its usage with the bound not enforced. Throws
// graph::traversal::AttemptAborted when the client is gone, AttemptAtBound when the attempt
// reached its bound while the results were built.
Json::Value process_traverse_request(const Json::Value &json,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const std::string &release,
                                     const TraverseLimits &limits = {},
                                     const IndexIdentity *identity = nullptr,
                                     Attempt *attempt = nullptr,
                                     ResultTexts *texts = nullptr);
// |client_gone|: polled between the phases of the request (the discovery read, the support
// fetch, the selection); true abandons it (graph::traversal::AttemptAborted)
Json::Value process_resolve_request(const Json::Value &json,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const std::string &release,
                                    uint64_t max_query_bp = 0,
                                    const IndexIdentity *identity = nullptr,
                                    const std::function<bool()> &client_gone = nullptr);

// CLI entry point: `metagraph traverse [--resolve] -i GRAPH -a ANNOTATION request.json ...`
int traverse_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_TRAVERSE_HPP__
