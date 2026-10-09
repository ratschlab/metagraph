#ifndef __METAGRAPH_CLI_TRAVERSE_HPP__
#define __METAGRAPH_CLI_TRAVERSE_HPP__

#include <chrono>
#include <functional>
#include <optional>
#include <tuple>
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
// What the server offers beyond the base contract (capabilities.feature_level). Monotonic: a
// change that adds capabilities fields or routes raises it by one; SPEC §10.3 states what each
// level adds.
constexpr int kTraverseFeatureLevel = 6;

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
 * The loader dependency inventory of an index listed as |graph| and |annotation|: every file
 * of the index's IDENTITY that the server's loaders open
 * for the pair, its path derived from the LISTED spelling exactly as the loaders derive it (a
 * sidecar next to a symlink, not next to its target), in loading order —
 *   graph               the graph file (required)
 *   annotation          the annotation file (required)
 *   row_diff_anchors    <graph>.anchors and
 *   row_diff_fork_succ  <graph>.rd_succ, for a .row_diff.annodbg (required: build_annotated_dbg)
 *   coord_to_header     <annotation without .<type>.annodbg>.seqs, for a coordinate
 *                       annotation when it exists and |coord_mapping| (not --no-coord-mapping)
 * A required file is listed whether it exists or not (the loader fails without it), an
 * optional one only when the loader would read it. Deliberately not listed: the graph's
 * DERIVED data, which the loader reads too (index_derived_files: the dummy-edge mask and the
 * Bloom filter), what the server does not open — a
 * column annotation's .coords (merge_load reads the columns only), the graph's .weights, any
 * leftover <graph>.anchors beside another annotation type (theirs are inside the annotation
 * file) — and the header index of the coordinate mapping, which is built in memory from the
 * .seqs. scripts/traversal/index_manifest.py mirrors this list (load_inventory; an integration
 * test compares the two on every sidecar kind through `traverse --index-inventory`).
 */
struct IndexFile {
    std::string path;
    const char *role;
    bool required;
};
std::vector<IndexFile> index_load_inventory(const std::string &graph,
                                            const std::string &annotation,
                                            bool coord_mapping = true);
/**
 * The derived data of the graph listed as |graph|:
 * files computed from the graph alone that DBGSuccinct::load reads beside it, and that are
 * NOT part of the index identity (index_fp) — an exact answer is the same with and without
 * them, and a manifest must not list them (index_manifest_fingerprint refuses one that does),
 * so that adding, removing or rebuilding them leaves an index's index_fp unchanged:
 *   graph_mask   <graph without .dbg>.edgemask: loaded when it opens
 *   graph_bloom  <graph without .dbg>.bloom: loaded when it exists and the mask was loaded
 * Both candidates of a .dbg graph are listed, whether they exist or not, with |exists| and
 * whether the loader reads them (|loaded|), so that an operator sees what is loaded; none for
 * another graph type. Paths derived from the LISTED spelling, as the loader derives them.
 * scripts/traversal/index_manifest.py mirrors it (derived_files).
 */
struct IndexDerivedFile {
    std::string path;
    const char *role;
    bool exists;
    bool loaded;
};
std::vector<IndexDerivedFile> index_derived_files(const std::string &graph);
// The rule of index_derived_files, stated by `traverse --index-inventory` and the refusal of
// a manifest that lists derived data
extern const char *const kIndexDerivedDataRule;
// the role of a file whose base name marks it as the derived data of a graph (it ends in
// .edgemask or .bloom), nullptr for any other: a manifest lists none (by extension, whichever
// graph it belongs to)
const char* index_derived_role(const std::string &path);
// the annotation extensions the loader knows (parse_annotation_type's, in its order) and what
// it reads beside each: the row-diff anchors beside the graph, the sequence headers
struct IndexAnnotationKind {
    std::string extension;
    bool row_diff_anchors;
    bool coordinates;
};
const std::vector<IndexAnnotationKind>& index_annotation_kinds();
// the paths of index_load_inventory (what a manifest is checked against; never a derived file)
std::vector<std::string> index_bundle_files(const std::string &graph,
                                            const std::string &annotation,
                                            bool coord_mapping = true);
// the optional identity files the inventory derives for the pair — a coordinate annotation's
// .seqs — that index_load_inventory leaves out (missing beside the listed spelling,
// |coord_mapping| false): a manifest must not list them (by base name). The graph's mask and
// Bloom filter are not among them: a manifest lists those never (index_derived_files)
std::vector<std::string> index_unloaded_optional_files(const std::string &graph,
                                                       const std::string &annotation,
                                                       bool coord_mapping = true);
// `traverse --index-inventory`: {graph, annotation, coord_mapping, files: [{path, role,
// required, exists}] (the identity: what a manifest lists), derived: [{path, role, exists,
// loaded}] (index_derived_files: loaded or not, never in a manifest), derived_rule (the text of
// kIndexDerivedDataRule), annotation_kinds: [{extension, row_diff_anchors, coordinates}]}
Json::Value index_inventory_json(const std::string &graph, const std::string &annotation,
                                 bool coord_mapping = true);
/**
 * Read an index manifest, validate it and return its fingerprint. The manifest is a JSON
 * object with `files: [{path, size, sha256}, ...]` (paths relative to the manifest,
 * '/'-separated, unique; sha256 64 lowercase hex); other keys (builder, inputs) are
 * metadata and not part of the identity. The fingerprint is the lowercase hex sha256 of
 * the lines "<path>\t<size>\t<sha256>\n" in ascending byte order of path; an `index_fp`
 * the manifest states must equal it. Its entries' base names must be distinct (loaded files
 * are matched by base name). Every file of |loaded| (index_bundle_files: the loader
 * dependency inventory) that exists on disk must be listed (by base name) with its size —
 * sizes, not digests: the server does not re-hash the bundle, so a replacement of one file by
 * another of the same size is not detected (index_manifest.py --verify re-hashes) — no file
 * of |not_loaded| (index_unloaded_optional_files) may be listed, no derived data of a graph
 * (an entry whose base name ends in .edgemask or .bloom, index_derived_role: not part of
 * index_fp, kIndexDerivedDataRule) may be listed, and the manifest must list
 * no graph (*dbg) or annotation (*.annodbg) file other than those of |loaded|: a manifest of
 * another bundle, or of a directory holding several, must not lend its identity, and
 * index_fp identifies the loaded identity files. Throws
 * std::runtime_error naming the problem. |stated_name|, when given, receives the manifest's
 * own `index_ns` (metadata; "" when it states none): the server's name for the index comes
 * from its configuration, and a different one in the manifest is only logged.
 */
std::string index_manifest_fingerprint(const std::string &manifest_path,
                                       const std::vector<std::string> &loaded,
                                       std::string *stated_name = nullptr,
                                       const std::vector<std::string> &not_loaded = {});
// lowercase hex SHA-256 (FIPS 180-4) of |data|
std::string sha256_hex(std::string_view data);
// the identity of the index |config| names (--index-name, --index-manifest against the
// loaded files -i / -a) with |anno_graph|'s meta fingerprint; throws std::runtime_error
IndexIdentity index_identity(const Config &config, const graph::AnnotatedDBG &anno_graph);

struct ResolveRequest {
    std::string sequence;
    graph::traversal::ResolveOptions options;
    bool select = false;
    graph::traversal::SelectionPolicy policy;
    std::string run_format = "intervals"; // intervals | rle
    uint64_t max_query_bp = 1'000'000;
    // bounds.time_budget_ms as given (opt-in: none, no deadline); the route checks it against
    // the finalisation reserve and the server's cap, which parsing does not know
    std::optional<double> time_budget_ms;
    Json::Value time_budget_given;          // the value as sent (a clamp states it)
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
// |check|: read every graph::traversal::kResolveCheckLabels objects built (a /resolve's
// deadline: it throws to abandon the answer)
Json::Value profile_to_json(const graph::traversal::SupportProfile &profile,
                            const graph::traversal::SeedSelection *selection,
                            const std::string &run_format,
                            const std::function<void()> &check = nullptr);
// The `coordinates` block of GET /traverse/capabilities (feature level 6; not in the
// per-request capabilities, which change only in their feature_level): whether this index
// reports record coordinates (supports_trace), the knobs and the cap's default, the block
// kinds it can report, the limitation and the action, and what bounds the block's size and the
// rule (positions, the interval rule, the column record-end numbering, the assumptions) as
// references to the SPEC section that states them
Json::Value coordinates_capabilities_json(const graph::traversal::LabelOracle &oracle);
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
// What a result's record coordinates (Strategy::coordinates, DESIGN-traverse-graphlet.md §18)
// add to its output: nothing (not asked for), `coordinates: null` with its reason (asked for
// and ruled out by the index or the support), or the block, where they can be recorded (or
// the null form, which it bounds)
enum class CoordinatesOutput { NONE, REASON, BLOCK };
// What one object of a result costs to deliver in |detail| (per object, bytes; upper
// bounds): the output part of the memory budget's model (DESIGN-traverse-graphlet.md §14),
// which process_traverse_request sets as Strategy::delivery, with every MGT float priced at
// |float_width| characters (mgt_float_width); |coordinates|: the coordinates' share (nothing
// without them; the other prices do not depend on it)
graph::traversal::DeliveryCosts delivery_costs(const std::string &detail, bool sequences,
                                               uint64_t float_width = kMgtFloatWidth,
                                               CoordinatesOutput coordinates
                                                   = CoordinatesOutput::NONE);
// the length the JSON writers here give |s| inside a JSON string, quotes excluded: what
// delivery_costs prices a name's text by
uint64_t json_escaped_size(std::string_view s);
// The length of json_text(|value|, true) — the compact text the server writes — counted from
// digits and fixed punctuation, without writing it (a string's or a real's own text is
// written by the writer itself, so that its escaping and digits are the writer's)
uint64_t compact_json_size(const Json::Value &value);
// The bytes of a seed's compact result text (one entry of `results`, as the server writes
// it) that are there because record coordinates were asked for: the `coordinates` member
// (the block, or null with `coordinates_reason`), a cut list's `coordinates` limitation, the
// `drop_coordinates` action of a resource stop and, in a graphlet, the coordinates K record,
// the Q record's drop_coordinates token and the digits these add to the Z record,
// graphlet_lines and graphlet_bytes. Counted exactly (compact_json_size), so that the text of
// the same result without them — the opt-out request's, when the walk is the same — is the
// result's text less this: what the server's delivery reserve leaves out of its ratio
// samples beside the account's coordinate share. 0 without them. |check|:
// the attempt's delivery check, called every 4096 values counted (a large block is counted
// under the attempt's bound, as it is written)
uint64_t coordinates_text_bytes(const Json::Value &result,
                                const std::function<void()> &check = {});

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
 * The MGT v1 token encoders (DESIGN-traverse-graphlet.md §2.1-§2.2, the golden vectors of
 * api/python/tests/data/traverse/codec_vectors.tsv). Every value has exactly one valid
 * spelling, which these write; a reader rejects every other one, so that two readers never
 * disagree on which documents are valid (the conformance tests' reader:
 * tests/cli/mgt_reader_for_tests.hpp). Exposed for the golden-vector test.
 */
namespace mgt {

// x >= 0 (or -0) or +inf: the shortest round-trip digits, expanded positionally
std::string encode_float(double x);
// the free-text last field: % LF CR as %25 %0A %0D, nothing else
std::string pct_escape(std::string_view raw);
// "<prefix_len> <suffix>" of an L record after |previous| (byte prefix cut back to a
// code-point boundary, suffix pct-escaped)
std::string front_code(std::string_view previous, std::string_view name);
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

} // namespace mgt

// Server-side caps on POST /traverse, 0 = unlimited. A request that names no labels
// makes the SERVER choose the permitted set, so the request's own bounds do not
// describe the work it asks for; these do. |max_time_ms| also bounds the derivation
// itself, and a request asking for less is left alone (a clamp only ever lowers a
// budget, never raises one) — see Config::traverse_max_*.
struct TraverseLimits {
    double max_time_ms = 0;        // cap on strategy.bounds.time_budget_ms
    size_t max_seeds = 0;          // cap on |request.seeds|
    uint64_t max_seed_bp = 0;      // cap on the length of one seed
    size_t max_seed_labels = 0;    // cap on strategy.labels.max_seed_labels
    // the maxima of the request's budgets (SPEC §10.3; 0 = off): a larger budget is lowered
    // to it, an omitted one set to it, both echoed in strategy.clamped
    uint64_t max_memory_mb = 0;    // of strategy.bounds.max_memory_mb
    uint64_t max_work_units = 0;   // of strategy.bounds.max_work_units
    // Not a cap: the chunked deadlines (spec §6.8, Config::traverse_chunk_target_ms) — an
    // annotation read a deadline may fall into is decoded in chunks of about this many ms, the
    // deadline checked between them; 0: one piece per read
    double chunk_target_ms = 0;
    // Not a cap: the bound of the request's row-diff path cache (LabelOracle::path_cache,
    // Config::traverse_path_cache_mb); 0: off, every read decodes its rows' whole paths
    uint64_t path_cache_bytes = 0;
    // Tests: the cache's retention rule (RowDiffCache::keeps: checkpoint, successors,
    // narrow_bytes) instead of its defaults; what it keeps never changes a response
    std::optional<std::tuple<uint32_t, uint32_t, uint64_t>> path_cache_retention;
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
// The finalisation reserve of a /resolve with bounds.time_budget_ms: its work stops this long
// before the deadline, so that the answer can still be built, written and compressed by it
// (as /pattern's, DESIGN-pattern-search.md §5.3; the same 250 ms as --pattern-finalize-ms's
// default). Fixed in this build: no flag sets it (stated in the capabilities)
constexpr double kResolveFinalizeMs = 250;

// The time limits of a /resolve: the cap of bounds.time_budget_ms (the server's
// --traverse-max-time-ms, as /traverse's; 0: none — the CLI, which an operator runs) and the
// finalisation reserve inside every budget. A request without the field has no deadline,
// whatever the cap: it is answered without a time limit
struct ResolveTimeLimits {
    double max_time_ms = 0;
    double finalize_ms = kResolveFinalizeMs;
};

// A /resolve answer that could not be built and written within its bounds.time_budget_ms (the
// finalisation overran the reserve): nothing partial is sent. The server answers 503 with
// resolve_deadline_body(), the CLI writes that body and exits 1
class ResolveDeadline : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};
// {"error": the message, "code": "deadline"} (the body of /pattern's 503 deadline)
Json::Value resolve_deadline_body(const ResolveDeadline &e);

// What the transport of one /resolve answer needs from its processing: the request's
// deadline, known once the request is parsed (none without bounds.time_budget_ms). The answer
// is written and compressed under check(), which throws ResolveDeadline once the budget has
// passed. Used by one request's thread only
class ResolveDelivery {
  public:
    void set_check(std::function<void()> check) { check_ = std::move(check); }
    void check() const {
        if (check_)
            check_();
    }

  private:
    std::function<void()> check_;
};

// The `resolve` block of both capabilities routes (GET /capabilities, GET
// /traverse/capabilities): bounds.time_budget_ms accepted, its cap and the reserve, and the
// rule (where the deadline is read, what a stop answers, what is not polled) as a reference to
// the SPEC section that states it
Json::Value resolve_capabilities_json(const ResolveTimeLimits &limits);

// |client_gone|: polled between the phases of the request (the discovery read, the support
// fetch, the selection); true abandons it (graph::traversal::AttemptAborted).
// |time|: the cap of bounds.time_budget_ms and the finalisation reserve; the deadline of a
// request with the field starts on entry (its body parsed), read from |clock| (the steady
// clock when null; injectable so that tests stop at a chosen instant). Its work stops at the
// budget less the reserve — the answer then is exactly the resolve of the query's first
// stop.resolved_kmers k-mers, with the `stop` block —; its finalisation is read against the
// whole budget, and |delivery|, when given, receives the check for the writing of the answer.
// Throws ResolveDeadline when the answer could not be built by the deadline.
Json::Value process_resolve_request(
        const Json::Value &json,
        const graph::AnnotatedDBG &anno_graph,
        const std::string &release,
        uint64_t max_query_bp = 0,
        const IndexIdentity *identity = nullptr,
        const std::function<bool()> &client_gone = nullptr,
        const ResolveTimeLimits &time = {},
        ResolveDelivery *delivery = nullptr,
        const std::function<std::chrono::steady_clock::time_point()> &clock = nullptr);

// CLI entry point: `metagraph traverse [--resolve] -i GRAPH -a ANNOTATION request.json ...`
int traverse_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_TRAVERSE_HPP__
