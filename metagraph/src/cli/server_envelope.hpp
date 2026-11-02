#ifndef __METAGRAPH_SERVER_ENVELOPE_HPP__
#define __METAGRAPH_SERVER_ENVELOPE_HPP__

#include <chrono>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include <json/json.h>

#include "server_selection.hpp"


namespace mtg {
namespace cli {

// The envelope a multi-graph server answers /resolve, /traverse and /pattern with, and the
// capabilities routes (SPEC-multigraph-envelope.md section 1): one entry per selected pair,
// each the pair's single-graph answer or its refusal, tagged with the pair

// The envelope's version (envelope_version)
constexpr int kEnvelopeVersion = 1;

// The pair as an entry and the lists not_reached and unfinished name it: {graph, graph_path,
// annotation_path}
Json::Value envelope_pair_ref(const SelectedPair &pair);

/**
 * The entry of a pair that answered: its single-graph answer |answer| with the tags merged
 * in, graph (the name it was selected by), graph_path, annotation_path, index_fp (null when
 * |index_fp| is empty: no manifest in the graph list) and outcome "answered"
 */
Json::Value envelope_entry_answered(const SelectedPair &pair, const std::string &index_fp,
                                    Json::Value answer);

/**
 * The entry of a refused pair: the same tags, outcome "refused", and `refusal`: the members
 * of |refusal_body| ({error, code?, reason?, usage?}: what the pair alone would have been
 * answered, or its not_reached / unfinished outcome of section 4) with |http_status|
 */
Json::Value envelope_entry_refused(const SelectedPair &pair, const std::string &index_fp,
                                   int http_status, const Json::Value &refusal_body);

/**
 * The request-level members of a POST envelope that its route states beside the entries:
 * |request_id| (the client's request_id, or /traverse's attempt_id; none: null),
 * |server_instance|, |deadline| (deadline_json; null: every member null) and |usage| (the
 * attempt's usage on a ledger-managed /traverse; null: the member is omitted)
 */
struct EnvelopeContext {
    std::optional<std::string> request_id;
    std::string server_instance;
    Json::Value deadline;
    Json::Value usage;
};

/**
 * The envelope's `deadline` (section 4): {not_after_ms, finish_by_ms} as the request gave
 * them (integers, ms since the epoch), {effective_deadline_at, cutoff_at} as ISO-8601 UTC
 * with milliseconds (iso_utc, as the attempt ledger writes its instants); null where not
 * given or not applicable
 */
Json::Value deadline_json(std::optional<uint64_t> not_after_ms,
                          std::optional<uint64_t> finish_by_ms,
                          std::optional<std::chrono::system_clock::time_point> effective,
                          std::optional<std::chrono::system_clock::time_point> cutoff);

/**
 * The envelope of a POST route (|route|: "/resolve", "/traverse" or "/pattern"):
 * envelope_version, route, pattern_contract_version (/pattern only), graphs (|names|, the
 * selected names in byte order, as given), answered and refused (the entries by outcome),
 * not_reached and unfinished (the pairs of the refused entries whose refusal.code says so),
 * request_id, server_instance, deadline, usage (when |context|.usage is not null), answers
 * (|entries| in order), timing.elapsed_ms and work.memory_bytes (the answered entries'
 * work.memory_bytes summed).
 */
Json::Value post_envelope(const std::string &route,
                          const std::vector<std::string> &names,
                          std::vector<Json::Value> entries,
                          const EnvelopeContext &context,
                          double elapsed_ms);

/**
 * The envelope of a GET capabilities route (|route|: "/traverse/capabilities" or
 * "/pattern/capabilities"): envelope_version, route, graphs (the entries' names, in their
 * order, each once), answered and refused, server_instance and answers (|entries| in order)
 */
Json::Value get_envelope(const std::string &route,
                         std::vector<Json::Value> entries,
                         const std::string &server_instance);

// The request registry's route (section 5), as GET /capabilities names it (routes.requests
// and multigraph.requests_route)
constexpr char kRequestsRoute[] = "GET /requests/{request_id}";

// How composite_fp is computed, stated beside it in GET /capabilities (composite_rule)
constexpr char kCompositeRule[] = "sha256 over the pairs' index_fp sorted as bytes, joined by "
                                  "\"\\n\" without a trailing newline; null when any pair has "
                                  "no index_fp; each loaded (graph_path, annotation_path) pair "
                                  "counted once, whatever names list it";

/**
 * GET /capabilities' composite_fp on a multi-graph server (section 8, kCompositeRule): the
 * lowercase hex sha256 over |index_fps|, the index_fp of each of the server's pairs, sorted
 * as bytes and joined by "\n" without a trailing newline; null when any of them is empty (a
 * pair listed without a manifest) or there is none. The caller passes each pair once
 */
Json::Value composite_fp_json(std::vector<std::string> index_fps);

/**
 * GET /capabilities' `multigraph` block on a multi-graph server (section 8): what this
 * release's multi-graph routes do, for clients to switch on (never on dates):
 * envelope_version, the routes that answer in the envelope (envelope_routes), the selectors
 * of each route (section 6), traverse_pairs_per_request, the routes that read not_after_ms,
 * finish_by_ms and request_id, the registry's route (requests_route), the flags
 * unknown_labels (section 2), search_cuts (section 3) and resolve_boundary_ties (section 7),
 * discovery_trace_scale, and |max_request_memory_mb| (the server's request-wide memory bound of
 * /traverse, --traverse-max-request-memory-mb; 0 = off), which every /traverse capabilities
 * probe states too: a graph-less probe sees it without reading a pair
 */
Json::Value multigraph_capabilities_json(uint64_t max_request_memory_mb);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_ENVELOPE_HPP__
