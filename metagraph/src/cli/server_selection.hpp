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
 * The release-9 choice of one pair by a name and a graph_path (/resolve, /traverse and the GET
 * probes of 8fa32e553): |pair| when the name holds exactly one pair on |graph_path| (when
 * given), else |outcome| and |error|, the text the base answered for it, byte for byte (a
 * multi-graph client matches these texts whole; SPEC-multigraph-envelope.md section 6 adds a
 * code beside them, never a new wording), and |matching| the name's pairs the choice could
 * not decide between (SEVERAL_ANNOTATIONS: the pairs on that one graph; SPANS_GRAPHS: every
 * pair of the name), each once, in the graph list's order.
 */
struct TraversePairChoice {
    enum Outcome { ONE, UNKNOWN_GRAPH, NOT_PART, SEVERAL_ANNOTATIONS, SPANS_GRAPHS };
    Outcome outcome = ONE;
    const GraphPair *pair = nullptr;
    std::string error;
    std::vector<GraphPair> matching;
};
TraversePairChoice choose_traverse_pair(const std::string &name, const std::string *graph_path,
                                        const GraphIndexes &indexes);

/**
 * Address exactly one physical (graph, annotation) pair for a traversal request or the probe.
 * A traversal must stay inside one graph, so a name covering several shards is rejected
 * unless the request disambiguates it with |graph_path|; a graph_path naming one graph with
 * several annotations under the name is rejected too — a traversal reads ONE annotation, and
 * answering from the first of them would leave the others unread without saying so. Throws
 * InvalidRequest (traverse.hpp) naming the problem (choose_traverse_pair's text).
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
// graph. The release-9 reading: the GET routes select through select_pairs_from_query, which
// also reads annotation_path and codes its refusals
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
 * envelope's order (ordered_pairs: the names in byte order, each name's pairs in the graph
 * list's order, a pair once), before any pair starts: the selectors of
 * SPEC-multigraph-envelope.md section 6.
 *  - `pairs` ([{graph, graph_path?, annotation_path?}], RESOLVE, PATTERN, TRAVERSE): each
 *    object names exactly one pair of its name (graph_path may be omitted when the name has
 *    one graph; annotation_path is needed where the name's pairs share a graph), none
 *    (pair_not_found) or several (pair_ambiguous) refuse the request with the candidates; a
 *    pair named twice is selected once. `pairs` beside `graphs` (or, on RESOLVE and TRAVERSE,
 *    beside `graph` or `graph_path`) is a selection_conflict. SEARCH refuses it
 *    (pairs_unsupported), CAPABILITIES too (a probe selects by graph, graph_path and
 *    annotation_path).
 *  - `graphs` (a non-empty array of names; an unknown one graph_unknown, the first in the
 *    list's order) selects every pair of its names; RESOLVE and TRAVERSE also take `graph`
 *    (one name; with `graphs`, a selection_conflict) and `graph_path` beside one name, which
 *    keeps the name's pairs on that graph (pair_not_found when it holds none).
 *  - without a selector, every pair of every name, when there are at most
 *    |max_names_without_selection| names (else the names cap's 400, a std::invalid_argument
 *    without a code, as before; on RESOLVE, new in release 10, in /pattern's words).
 *  - TRAVERSE runs on exactly one pair: a selection of several is refused with
 *    traverse_single_pair and the selected pairs as the candidates.
 *  - Every request the release-9 server (8fa32e553) refused is refused in its words, byte for
 *    byte (a client matches them whole), the code added: on RESOLVE and TRAVERSE without
 *    `pairs` the texts of traverse_graph_selection and choose_traverse_pair ("Bad request:
 *    unknown graph 'n'", "... spans several graphs ...", "... with several annotations ...",
 *    "... is not part of index ...", "'graphs' names the one graph a traversal reads: expected
 *    [name]", "'graph' (index name) is required in multi-graph mode", "give the graph as
 *    'graph' or as 'graphs', not both"); on PATTERN pattern_graph_names' texts.
 *  - SEARCH reads `graphs` as POST /search always has (filter_graphs_from_list: a `graphs`
 *    that is no array selects nothing, every name is then read), an unknown name refused
 *    with its message and graph_unknown.
 *  - CAPABILITIES reads the members `graph`, `graph_path` and `annotation_path` as
 *    select_pairs_from_parameters reads the GET probes' parameters.
 * Throws SelectionRefusal (HTTP 400 {error, code, candidates?}), InvalidRequest for a
 * malformed selector, std::invalid_argument for the names cap.
 */
std::vector<SelectedPair> select_pairs(SelectionRoute route,
                                       const Json::Value &request,
                                       const GraphIndexes &indexes,
                                       size_t max_names_without_selection);

/**
 * What select_pairs refuses on RESOLVE and TRAVERSE before the route reads `in_ram`: the
 * release-9 server (8fa32e553) checked the selection's shape (traverse_graph_selection:
 * `graph` beside `graphs`, a `graphs` that is not a list of one name, a missing or non-string
 * `graph`, `graph_path` without a name) before in_ram and resolved the name after it, so that
 * a request it refused for two faults keeps its words (section 6): the routes call this, then
 * in_ram_field, then select_pairs. Also `pairs` beside another selector (a
 * selection_conflict), the names cap of a /resolve without a selector, and what the base
 * refused because it named no single graph (several names; no name on /traverse). Nothing on
 * the other routes (/pattern read its graphs before in_ram, as select_pairs does).
 */
void check_selection_shape(SelectionRoute route,
                           const Json::Value &request,
                           const GraphIndexes &indexes,
                           size_t max_names_without_selection);

/**
 * The pairs a GET capabilities probe selects on a multi-graph server, from its query
 * parameters (section 6): ?graph=<name>[&graph_path=<path>][&annotation_path=<path>] names
 * exactly one pair of the name, as an object of `pairs` does (graph_unknown, pair_not_found,
 * pair_ambiguous, with the candidates; without annotation_path in the words of
 * choose_traverse_pair, as before release 10); without parameters, every pair in the
 * envelope's order, when the server hosts at most |max_names_without_selection| names: the
 * envelope over every pair answers such a probe (section 1); beyond the cap, and with
 * graph_path or annotation_path but no ?graph=, the probe is refused in the words a probe
 * without ?graph= always had, byte for byte ("Bad request: in multi-graph mode GET <route>
 * needs ?graph=<name> (and graph_path=<path> when the name spans several graphs); GET
 * /capabilities lists the graphs"). |route| names the route in the messages: TRAVERSE for GET
 * /traverse/capabilities, PATTERN for GET /pattern/capabilities (CAPABILITIES names both). An
 * unknown or a repeated parameter (in the base's words), or a missing ?graph=, is an
 * InvalidRequest.
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
