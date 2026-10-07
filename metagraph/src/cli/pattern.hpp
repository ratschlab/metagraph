#ifndef __METAGRAPH_CLI_PATTERN_HPP__
#define __METAGRAPH_CLI_PATTERN_HPP__

/**
 * POST /pattern and `metagraph pattern`: count, and extract without reading any annotation,
 * the graph contexts of short motifs and IUPAC patterns (docs/DESIGN-pattern-search.md,
 * increments 0-2). The engine is graph::pattern::PatternSearch
 * (src/graph/alignment/pattern_search.hpp); this file turns its results into the JSON of the
 * route's contract (pattern_contract_version 1):
 *  - modes count, all_or_count and partial; the two retrieval modes with output.labels "none"
 *    only (the label-free path, §4.3): contexts with k-mer, instance, offset, strand, node and
 *    row ids, and no annotation row read anywhere;
 *  - single-graph servers only; no placement, no predicate, no extension beyond k (a pattern
 *    longer than k has its anchors counted and nothing extracted).
 * Everything a later increment adds is refused (400 "later_increment"), never ignored: the
 * owner's guarantee rule, nothing weakened silently.
 */

#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>

#include <json/json.h>

#include "graph/alignment/pattern_search.hpp"


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
 * support reasons; 503 deadline). The server answers it as is (HttpError); the CLI writes the
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
    // the finalisation reserve inside the time budget (§5.3): work stops this long before the
    // request's deadline so that the answer can still be written by it
    double finalize_ms = 250;
    double min_information_bits = 24;
    // patterns per request: above it the request is refused (a list is not cut)
    uint64_t max_patterns = 16;
};

PatternLimits pattern_limits(const Config &config);

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
    // nothing before a deadline was set (the request was refused before it was parsed)
    void check() const;

  private:
    std::optional<graph::pattern::Deadline> deadline_;
};

// The request body as JSON; a body that is not JSON is refused as invalid_request (the
// generic parse error of the other routes carries no code)
Json::Value parse_pattern_body(const std::string &content);

/**
 * Answers one /pattern request (the JSON contract of pattern_contract_version 1) on the
 * single graph of |anno_graph| under |limits|, stating |release|. The deadline starts on
 * entry (the body was parsed), read from |clock| (the steady clock when null; injectable so
 * that tests stop at a chosen instant); |delivery|, when given, receives it for the writing
 * of the answer. |identity| (may be null) states the index in the answer's `index`. Throws
 * PatternRefusal for a whole-request refusal, 503 "deadline" included when the answer could
 * not be assembled by the deadline; refused patterns are answered in their slots. Never reads
 * an annotation row.
 */
Json::Value process_pattern_request(
        const Json::Value &json,
        const graph::AnnotatedDBG &anno_graph,
        const PatternLimits &limits,
        const std::string &release,
        const IndexIdentity *identity = nullptr,
        PatternDelivery *delivery = nullptr,
        const std::function<graph::pattern::Deadline::Clock::time_point()> &clock = nullptr);

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
 * annotation could give a later increment (placement, support, annotation). |anno_graph| is
 * null while the index loads; |multi_graph| servers answer only that they are not served yet.
 */
Json::Value pattern_capabilities_json(const graph::AnnotatedDBG *anno_graph,
                                      const PatternLimits &limits, bool multi_graph);

// `metagraph pattern -i GRAPH -a ANNOTATION REQUEST.json ...`: one JSON answer per request
// file on stdout, the server's (a refusal's body, exit 1); for tests and offline use
int pattern_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_PATTERN_HPP__
