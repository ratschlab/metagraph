#ifndef __METAGRAPH_SERVER_SELECTION_HPP__
#define __METAGRAPH_SERVER_SELECTION_HPP__

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>
#include <tsl/hopscotch_map.h>

#include "server_checks.hpp"


namespace mtg {
namespace cli {

// Which (graph, annotation) pairs of a multi-graph server a request reads
// (SPEC-multigraph-envelope.md section 6), apart from the server and without the HTTP
// library, so that the selection is tested on its own and shared by the routes

// A physical pair of the graph list: (graph path, annotation path)
using GraphPair = std::pair<std::string, std::string>;
// The graph list: each name with its pairs in the list's order (a pair may be listed twice
// under a name)
using GraphIndexes = tsl::hopscotch_map<std::string, std::vector<GraphPair>>;

// The routes that select pairs, each with its own selectors (the table of section 6)
enum class SelectionRoute { SEARCH, RESOLVE, PATTERN, TRAVERSE, CAPABILITIES };

// One selected pair and the name it was selected by: the tags of its envelope entry, its
// /search rows and its cut
struct SelectedPair {
    std::string name;
    GraphPair pair;
};
bool operator==(const SelectedPair &a, const SelectedPair &b);
inline bool operator!=(const SelectedPair &a, const SelectedPair &b) { return !(a == b); }

// The pair as a request or an answer names it: {graph, graph_path, annotation_path}
Json::Value pair_ref_json(const SelectedPair &pair);

// The pairs |name| holds, each once, in the graph list's order, as the candidates of a
// refusal: [{graph, graph_path, annotation_path}]; an empty array for an unknown name
Json::Value selection_candidates(const std::string &name, const GraphIndexes &indexes);

/**
 * A selection refused before any pair starts: HTTP 400 {error, code, candidates?} (section
 * 6), |code| one of selection_conflict, pair_not_found, pair_ambiguous, traverse_single_pair,
 * graph_unknown, pairs_unsupported, and |candidates| the pairs to choose from where a choice
 * exists (null otherwise: the member is then absent). An HttpError, so that a route that does
 * not catch it answers it as it is.
 */
class SelectionRefusal : public HttpError {
  public:
    SelectionRefusal(const std::string &error, const std::string &code,
                     Json::Value candidates = Json::Value());
    const std::string& code() const { return code_; }
    const Json::Value& candidates() const { return candidates_; }

  private:
    static Json::Value body_of(const std::string &error, const std::string &code,
                               const Json::Value &candidates);

    std::string code_;
    Json::Value candidates_;
};

/**
 * Address exactly one physical (graph, annotation) pair for a traversal request or the probe.
 * A traversal must stay inside one graph, so a name covering several shards is rejected
 * unless the request disambiguates it with |graph_path|; a graph_path naming one graph with
 * several annotations under the name is rejected too — a traversal reads ONE annotation, and
 * answering from the first of them would leave the others unread without saying so. Throws
 * InvalidRequest (traverse.hpp) naming the problem.
 */
const GraphPair& select_traverse_pair(const std::string &name, const std::string *graph_path,
                                      const GraphIndexes &indexes);

/**
 * The names a POST /search request selects: its `graphs` (an array; each name known, else
 * std::invalid_argument "Request with an uninitialized graph <name>"), deduplicated and
 * sorted; without `graphs` (or with a `graphs` that is no array) every name of the list, in
 * the map's order, when there are at most |max_names_without_filtering| of them (the
 * server's --max-graphs-without-selection), std::invalid_argument naming the threshold
 * otherwise. |request_id| is not read.
 */
std::vector<std::string> filter_graphs_from_list(const GraphIndexes &indexes,
                                                 const Json::Value &content_json,
                                                 size_t request_id,
                                                 size_t max_names_without_filtering = 10);

// The pairs of a multi-graph server in the order its answers list them: the names in byte
// order, each name's pairs in the graph list's order, a pair listed twice under a name once
std::vector<std::pair<std::string, GraphPair>>
ordered_pairs(const GraphIndexes &indexes, const std::vector<std::string> &names);

// Every name of the graph list, in byte order
std::vector<std::string> sorted_names(const GraphIndexes &indexes);

// The parameters of a GET probe of one pair on a multi-graph server (?graph=<name>
// [&graph_path=<path>]), refused as GET /traverse/capabilities refuses them (|route| names the
// route in the messages): InvalidRequest for an unknown or a repeated parameter, or a missing
// graph
std::pair<std::string, std::optional<std::string>>
probe_parameters(const std::vector<std::pair<std::string, std::string>> &query,
                 const std::string &route);
template <class Query>  // SimpleWeb::CaseInsensitiveMultimap, or any range of (key, value)
std::pair<std::string, std::optional<std::string>>
probe_parameters(const Query &query, const std::string &route) {
    return probe_parameters(std::vector<std::pair<std::string, std::string>>(query.begin(),
                                                                             query.end()),
                            route);
}

/**
 * The pairs a request on |route| selects on a multi-graph server (|indexes|), in the
 * envelope's order (ordered_pairs), before any pair starts. These are TODAY's rules of each
 * route, so that a handler moved onto this function answers as before, with the refusals'
 * messages and exception types it states today; lane E12 (section 6) replaces the bodies with
 * the `pairs` selector and the coded refusals:
 *  - every route: a `pairs` member is refused, SelectionRefusal pairs_unsupported;
 *  - SEARCH: the names of `graphs` as POST /search reads them (filter_graphs_from_list, its
 *    refusals and the names cap |max_names_without_selection|), each name's pairs. POST
 *    /search itself still queries a pair listed twice under a name twice;
 *  - PATTERN: `graphs` as POST /pattern reads it (pattern_graph_names, its refusals and the
 *    names cap), each name's pairs;
 *  - RESOLVE (the envelope, any number of pairs): `graph` (or `graphs` naming one) with
 *    `graph_path` its one pair, `graph` alone each pair of its name, else `graphs` and no
 *    selection as PATTERN reads them;
 *  - TRAVERSE: the one pair of today's /traverse: `graph`, or
 *    `graphs` ([name]), and `graph_path` where the name spans several graphs
 *    (traverse_graph_selection, select_traverse_pair, their refusals); no names cap, since
 *    the name is required;
 *  - CAPABILITIES: the members `graph` and `graph_path` read as select_pairs_from_parameters
 *    reads the GET probes' parameters.
 * Throws std::invalid_argument (InvalidRequest for the pair rules) or SelectionRefusal.
 */
std::vector<SelectedPair> select_pairs(SelectionRoute route,
                                       const Json::Value &request,
                                       const GraphIndexes &indexes,
                                       size_t max_names_without_selection);

/**
 * The pair a GET capabilities probe selects on a multi-graph server, by today's rule: one pair
 * per probe, ?graph=<name>[&graph_path=<path>] (probe_parameters, then select_traverse_pair,
 * and their refusals; a probe without ?graph= is refused). |route| names the route in the
 * messages: TRAVERSE for GET /traverse/capabilities, PATTERN for GET /pattern/capabilities
 * (CAPABILITIES names both). |max_names_without_selection| is not read yet: it bounds the
 * envelope over every pair that answers a probe without ?graph= (section 1), which replaces
 * today's refusal. Throws InvalidRequest.
 */
std::vector<SelectedPair>
select_pairs_from_parameters(SelectionRoute route,
                             const std::vector<std::pair<std::string, std::string>> &query,
                             const GraphIndexes &indexes,
                             size_t max_names_without_selection);
template <class Query>  // SimpleWeb::CaseInsensitiveMultimap, or any range of (key, value)
std::vector<SelectedPair> select_pairs_from_query(SelectionRoute route,
                                                  const Query &query,
                                                  const GraphIndexes &indexes,
                                                  size_t max_names_without_selection) {
    return select_pairs_from_parameters(
            route, std::vector<std::pair<std::string, std::string>>(query.begin(), query.end()),
            indexes, max_names_without_selection);
}

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_SELECTION_HPP__
