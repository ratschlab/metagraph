#ifndef __METAGRAPH_CLI_PATTERN_HPP__
#define __METAGRAPH_CLI_PATTERN_HPP__

/**
 * POST /pattern and `metagraph pattern`: count and extract the graph contexts of short motifs
 * and IUPAC patterns, and read their labels (docs/DESIGN-pattern-search.md, increments 0-3).
 * The engine is graph::pattern::PatternSearch (src/graph/alignment/pattern_search.hpp), the
 * labelled retrieval PatternRetrieval (pattern_retrieval.hpp); this file turns their results
 * into the JSON of the route's contract (pattern_contract_version 1):
 *  - modes count, all_or_count and partial; the two retrieval modes with output.labels "none"
 *    (the label-free path, §4.3: contexts with k-mer, instance, offset, strand, node and row
 *    ids, no annotation row read) or "all" (increment 3: each context's labels, placed where
 *    the index can place them, under the annotation budgets);
 *  - single-graph servers only; no predicate, no extension beyond k (a pattern longer than k
 *    has its anchors counted and nothing extracted).
 * Everything a later increment adds is refused (400 "later_increment"), never ignored: the
 * owner's guarantee rule, nothing weakened silently.
 */

#include <cstdint>
#include <functional>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>

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
 * (400 invalid_request, later_increment, resident_only, mask_required and the other graph
 * support reasons of route_support, alphabet_untested and mask_invalid among them; 503
 * deadline). The server answers it as is (HttpError); the CLI writes the
 * same body and exits 1. A refusal of one pattern is not this: it is the `error` of that
 * pattern's slot in a 200 answer.
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

/**
 * The server's caps of a /pattern request (the --pattern-* flags; the CLI applies the same
 * ones, so that both answer alike). Each is a request field's maximum: a larger request value
 * is lowered to it and listed in limits.clamped (the /traverse convention), never refused,
 * never silently kept; and each is the field's default, except the time budget, whose
 * default (60 s) lies under its maximum (600 s), the owner's decision of 2026-10-07 (§5.3).
 */
struct PatternLimits {
    uint64_t max_contexts = 10'000;
    uint64_t max_anchors = 1'000;
    uint64_t max_steps = 100'000'000;
    double default_time_ms = 60'000;
    double max_time_ms = 600'000;
    // the finalisation reserve inside the time budget (§5.3): work stops at least this long
    // before the request's deadline so that the answer can still be written by it (longer by
    // the estimated finalisation of what the answer buffers: delivery_* below)
    double finalize_ms = 250;
    double min_information_bits = 24;
    // patterns per request: above it the request is refused (a list is not cut)
    uint64_t max_patterns = 16;
    // the labelled retrieval (output.labels "all", increment 3; §4.3, §5.3): the labels kept
    // per row, the annotation work (the oracle's units), the request's memory account (MiB),
    // and partial's lists of labels per pattern and of occurrences per label
    uint64_t max_labels_per_anchor = 64;
    uint64_t max_annotation_work = 100'000'000;
    uint64_t max_memory_mb = 256;
    uint64_t max_labels = 1'000;
    uint64_t max_occurrences_per_label = 16;
    // not a cap: the annotation reads under the deadline are decoded in chunks of about this
    // many ms (the server's --traverse-chunk-target-ms, as /traverse's reads); 0: one piece
    double chunk_target_ms = 50;
    // not caps: the finalisation of what the answer buffers (AnswerVolume; review of
    // 2026-10-07, X-EFFICIENCY-04): the rates (MB/s) at which its JSON text is assumed to be
    // built and written, and compressed (--pattern-delivery-build-mbps,
    // --pattern-delivery-compress-mbps), and the text written per byte of compact JSON (1 on
    // the server; the CLI's indented text 2). The work stops finalize_ms plus that estimate
    // before the deadline. 0 or infinity: the reserve alone
    double delivery_build_mbps = 10;
    double delivery_compress_mbps = 50;
    double delivery_text_scale = 1;
};

PatternLimits pattern_limits(const Config &config);

/**
 * The alphabet half of the route's support decision, a pure function of the graph's BOSS
 * alphabet (owner decision #4 of 2026-10-07, review I26): "" for "$ACGT" (served),
 * "alphabet_untested" for "$ACGTN" (the engine supports it, but the route does not serve it
 * until a DNA5 build passes the pattern tests), "alphabet_unsupported" for any other.
 */
std::string alphabet_refusal(const std::string &alphabet);

/**
 * Whether /pattern and `metagraph pattern` serve |graph|: the engine's PatternSearch::support,
 * narrowed by the route's own reasons, "alphabet_untested" (alphabet_refusal; it takes
 * precedence over mask_required, as alphabet_unsupported does) and "mask_invalid" (a mask
 * that marks an edge with W = $ valid, found once at load: check_mask_at_load; review of
 * 2026-10-07, I17, owner decision #6). Its reason is the 400 refusal's code and the
 * capabilities' unavailable_reason. The engine itself keeps serving $ACGTN (its tests).
 */
graph::pattern::GraphSupport route_support(const graph::DeBruijnGraph &graph);

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
     * The server's "the client left, or the server stops" (review of 2026-10-07,
     * X-CONCURRENCY-01, R2-02: /pattern ran to its deadline for a caller that was gone):
     * process_pattern_request gives it to the request's Budget (Budget::set_abort, read at
     * every clock reading of the work), and check() asks it while the answer is assembled,
     * written and compressed. Its answer true throws graph::pattern::Aborted; the server
     * then writes nothing. Unset (the CLI, tests): never asked.
     */
    void set_abort(std::function<bool()> aborted) { abort_ = std::move(aborted); }
    const std::function<bool()>& abort() const { return abort_; }

  private:
    std::optional<graph::pattern::Deadline> deadline_;
    std::function<bool()> abort_;
};

// The request body as JSON: one RFC 8259 JSON text with unique member names and nothing after
// it, nested at most 1,000 deep; anything else is refused as invalid_request (the generic
// parse error of the other routes carries no code; their parser stays as lenient as it was)
Json::Value parse_pattern_body(const std::string &content);

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
 * annotation rows only for output.labels "all" in a retrieval mode (PatternRetrieval);
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
 * The `pattern` block of GET /capabilities (§7.3): the contract version, whether this server
 * can answer /pattern (available: true | false | null while the single index loads, with the
 * reason when false), the modes, projections, kinds, scopes and strands, the caps and the
 * finalisation reserve, and what the graph is (mode, k, alphabet, mask) and what its
 * annotation gives the labelled retrieval (placement, support, annotation: budgeted or
 * unbudgeted). |anno_graph| is null while the index loads; |multi_graph| servers answer only
 * that they are not served yet.
 */
Json::Value pattern_capabilities_json(const graph::AnnotatedDBG *anno_graph,
                                      const PatternLimits &limits, bool multi_graph);

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
