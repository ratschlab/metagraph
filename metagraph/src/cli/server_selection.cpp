#include "server_selection.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

#include <spdlog/fmt/fmt.h>

#include "traverse.hpp"


namespace mtg {
namespace cli {

bool operator==(const SelectedPair &a, const SelectedPair &b) {
    return a.name == b.name && a.pair == b.pair;
}

Json::Value pair_ref_json(const SelectedPair &pair) {
    Json::Value ref(Json::objectValue);
    ref["graph"] = pair.name;
    ref["graph_path"] = pair.pair.first;
    ref["annotation_path"] = pair.pair.second;
    return ref;
}

Json::Value selection_candidates(const std::string &name, const GraphIndexes &indexes) {
    Json::Value candidates(Json::arrayValue);
    auto it = indexes.find(name);
    if (it == indexes.end())
        return candidates;
    for (const auto &[_, pair] : ordered_pairs(indexes, { name })) {
        candidates.append(pair_ref_json({ name, pair }));
    }
    return candidates;
}

SelectionRefusal::SelectionRefusal(const std::string &error, const std::string &code,
                                   Json::Value candidates)
      : HttpError(400, body_of(error, code, candidates)), code_(code),
        candidates_(std::move(candidates)) {}

Json::Value SelectionRefusal::body_of(const std::string &error, const std::string &code,
                                      const Json::Value &candidates) {
    Json::Value body(Json::objectValue);
    body["error"] = error;
    body["code"] = code;
    // where a choice exists only: a refusal that names no candidates has no member
    if (!candidates.isNull())
        body["candidates"] = candidates;
    return body;
}

const GraphPair& select_traverse_pair(const std::string &name, const std::string *graph_path,
                                      const GraphIndexes &indexes) {
    auto it = indexes.find(name);
    if (it == indexes.end())
        throw InvalidRequest("Bad request: unknown graph '" + name + "'");

    // the same pair listed twice under the name is one index (/search keeps its list as is)
    std::vector<const GraphPair*> pairs;
    for (const auto &pair : it->second) {
        if (std::none_of(pairs.begin(), pairs.end(),
                         [&](const GraphPair *p) { return *p == pair; })) {
            pairs.push_back(&pair);
        }
    }
    if (pairs.size() == 1)
        return *pairs[0];

    // a name whose pairs share one graph has several annotations, which graph_path cannot
    // choose between
    auto several_annotations = [&](const std::string &graph) {
        std::string annotations;
        for (const GraphPair *pair : pairs) {
            if (pair->first == graph)
                annotations += (annotations.empty() ? "" : ", ") + pair->second;
        }
        return InvalidRequest("Bad request: " + (graph_path
                                  ? "graph_path '" + graph + "' names a graph"
                                  : "index '" + name + "' lists one graph (" + graph + ")")
                              + " with several annotations" + (graph_path ? " under index '"
                                  + name + "'" : std::string()) + " (" + annotations
                              + "): a traversal reads one annotation, which graph_path cannot "
                              "choose; list each annotation under a name of its own");
    };
    if (graph_path) {
        const GraphPair *found = nullptr;
        size_t matching = 0;
        for (const GraphPair *pair : pairs) {
            if (pair->first != *graph_path)
                continue;
            matching++;
            if (!found)
                found = pair;
        }
        if (!found) {
            throw InvalidRequest("Bad request: graph_path '" + *graph_path + "' is not part of "
                                 "index '" + name + "'");
        }
        if (matching > 1)
            throw several_annotations(*graph_path);
        return *found;
    }
    std::vector<std::string> graphs;
    for (const GraphPair *pair : pairs) {
        if (std::find(graphs.begin(), graphs.end(), pair->first) == graphs.end())
            graphs.push_back(pair->first);
    }
    if (graphs.size() == 1)
        throw several_annotations(graphs[0]);
    std::string candidates;
    for (const std::string &graph : graphs) {
        candidates += (candidates.empty() ? "" : ", ") + graph;
    }
    throw InvalidRequest("Bad request: index '" + name + "' spans several graphs (" + candidates
                         + "); pass 'graph_path' to pick one");
}

std::vector<std::string> filter_graphs_from_list(const GraphIndexes &indexes,
                                                 const Json::Value &content_json,
                                                 size_t /* request_id */,
                                                 size_t max_names_without_filtering) {
    std::vector<std::string> graphs_to_query;

    if (content_json.isMember("graphs") && content_json["graphs"].isArray()) {
        for (const auto &item : content_json["graphs"]) {
            graphs_to_query.push_back(item.asString());
            if (!indexes.count(graphs_to_query.back()))
                throw std::invalid_argument("Request with an uninitialized graph " + graphs_to_query.back());
        }
        // deduplicate
        std::sort(graphs_to_query.begin(), graphs_to_query.end());
        graphs_to_query.erase(std::unique(graphs_to_query.begin(), graphs_to_query.end()), graphs_to_query.end());
    } else {
        if (indexes.size() > max_names_without_filtering) {
            throw std::invalid_argument(
                fmt::format("Bad request: requests without names (no \"graphs\" field) are "
                            "only supported for small indexes (<={} names)",
                            max_names_without_filtering));
        }
        // query all graphs from list `config->fnames`
        for (const auto &[name, _] : indexes) {
            graphs_to_query.push_back(name);
        }
    }
    return graphs_to_query;
}

std::vector<std::pair<std::string, GraphPair>>
ordered_pairs(const GraphIndexes &indexes, const std::vector<std::string> &names) {
    std::vector<std::pair<std::string, GraphPair>> pairs;
    for (const std::string &name : names) {
        const size_t first = pairs.size();
        for (const GraphPair &pair : indexes.at(name)) {
            if (std::none_of(pairs.begin() + first, pairs.end(),
                             [&](const auto &p) { return p.second == pair; })) {
                pairs.emplace_back(name, pair);
            }
        }
    }
    return pairs;
}

std::vector<std::string> sorted_names(const GraphIndexes &indexes) {
    std::vector<std::string> names;
    for (const auto &[name, _] : indexes) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

std::pair<std::string, std::optional<std::string>>
probe_parameters(const std::vector<std::pair<std::string, std::string>> &query,
                 const std::string &route) {
    std::optional<std::string> name, graph_path;
    for (const auto &[key, value] : query) {
        std::optional<std::string> *slot = key == "graph" ? &name
                                         : key == "graph_path" ? &graph_path : nullptr;
        if (!slot) {
            throw InvalidRequest("Bad request: unknown parameter '" + key + "' of GET "
                                 + route + " (graph, graph_path)");
        }
        if (*slot)
            throw InvalidRequest("Bad request: the parameter '" + key + "' is given twice");
        *slot = value;
    }
    if (!name) {
        throw InvalidRequest("Bad request: in multi-graph mode GET " + route + " "
                             "needs ?graph=<name> (and graph_path=<path> when the name "
                             "spans several graphs); GET /capabilities lists the graphs");
    }
    return { *name, graph_path };
}

namespace {

// The GET route a probe on |route| answers, as its refusals name it: the two probes take the
// same parameters, so the route is all that tells their messages apart
std::string probe_route(SelectionRoute route) {
    switch (route) {
        case SelectionRoute::TRAVERSE:
            return "/traverse/capabilities";
        case SelectionRoute::PATTERN:
            return "/pattern/capabilities";
        default:
            return "/traverse/capabilities or /pattern/capabilities";
    }
}

// |request| when it is an object, else an empty one: the routes read a request that is no
// object as one without selectors (its parser refuses it). A reference, never a copy of a
// request that may carry megabytes of sequence
const Json::Value& object_or_empty(const Json::Value &request) {
    static const Json::Value empty(Json::objectValue);
    return request.isObject() ? request : empty;
}

// Each pair of |names| in the envelope's order
std::vector<SelectedPair> pairs_of(const std::vector<std::string> &names,
                                   const GraphIndexes &indexes) {
    std::vector<SelectedPair> selected;
    for (auto &[name, pair] : ordered_pairs(indexes, names)) {
        selected.push_back(SelectedPair { std::move(name), std::move(pair) });
    }
    return selected;
}

Json::Value candidates_of(const std::vector<SelectedPair> &pairs) {
    Json::Value candidates(Json::arrayValue);
    for (const SelectedPair &pair : pairs) {
        candidates.append(pair_ref_json(pair));
    }
    return candidates;
}

std::string in_quotes(const std::string &text) { return "'" + text + "'"; }

// The names cap of a request without a selector: its existing 400 (no code), as POST /pattern
// states it; |selector| names what the request would add
void check_names_cap(const GraphIndexes &indexes, size_t max_names,
                     const std::string &selector) {
    if (indexes.size() <= max_names)
        return;
    throw std::invalid_argument(selector + ": required on this server, which hosts "
                                + std::to_string(indexes.size()) + " graph names (more than "
                                + std::to_string(max_names) + "; GET /capabilities lists them)");
}

// What names one pair of |name|: an object of `pairs`, or a probe's parameters
struct PairSpec {
    std::string name;
    std::optional<std::string> graph_path;
    std::optional<std::string> annotation_path;
};

// "graph_path 'g'", "annotation_path 'a'", or both
std::string paths_text(const PairSpec &spec) {
    std::string text;
    if (spec.graph_path)
        text = "graph_path " + in_quotes(*spec.graph_path);
    if (spec.annotation_path) {
        text += (text.empty() ? "" : " and ") + std::string("annotation_path ")
              + in_quotes(*spec.annotation_path);
    }
    return text;
}

// The pairs of |spec|'s name its paths keep, in the list's order, a pair once; the name is
// known
std::vector<SelectedPair> matching_pairs(const PairSpec &spec, const GraphIndexes &indexes) {
    std::vector<SelectedPair> matching;
    for (SelectedPair &pair : pairs_of({ spec.name }, indexes)) {
        if ((!spec.graph_path || pair.pair.first == *spec.graph_path)
                && (!spec.annotation_path || pair.pair.second == *spec.annotation_path)) {
            matching.push_back(std::move(pair));
        }
    }
    return matching;
}

// graph_unknown for |name| unless the list has it; |where| prefixes the message ("" or
// "request.<field>: ")
void check_known(const std::string &name, const GraphIndexes &indexes,
                 const std::string &where) {
    if (!indexes.count(name)) {
        throw SelectionRefusal(where + "unknown graph " + in_quotes(name)
                               + " (GET /capabilities lists the graphs)", "graph_unknown");
    }
}

/**
 * The one pair |spec| names (section 6): its name known (graph_unknown), and exactly one of
 * the name's pairs on its graph_path and annotation_path (each when given): none is
 * pair_not_found with the name's pairs as the candidates, several pair_ambiguous with them
 * (graph_path to add where they span graphs, annotation_path where they share one).
 * |where| prefixes the messages, |subject| names the name in them ("'n'" in an object of
 * `pairs`, "graph 'n'" in a probe).
 */
SelectedPair one_pair_of(const PairSpec &spec, const GraphIndexes &indexes,
                         const std::string &where, const std::string &subject) {
    check_known(spec.name, indexes, where);
    std::vector<SelectedPair> matching = matching_pairs(spec, indexes);
    if (matching.size() == 1)
        return std::move(matching.front());
    if (matching.empty()) {
        throw SelectionRefusal(where + "no pair of " + subject + " has " + paths_text(spec),
                               "pair_not_found", selection_candidates(spec.name, indexes));
    }
    std::vector<std::string> graphs;
    for (const SelectedPair &pair : matching) {
        if (std::find(graphs.begin(), graphs.end(), pair.pair.first) == graphs.end())
            graphs.push_back(pair.pair.first);
    }
    const std::string count = std::to_string(matching.size()) + " pairs";
    std::string error;
    if (graphs.size() == 1) {
        // one graph with several annotations (sra-human): only annotation_path tells them apart
        error = subject + " has " + count + " on graph_path " + in_quotes(graphs.front())
              + ": add annotation_path";
    } else {
        // no graph_path was given (one would leave a single graph): the pairs span graphs,
        // and where one of them carries several annotations annotation_path is needed too
        const bool shared = matching.size() > graphs.size();
        error = subject + " has " + count
              + (spec.annotation_path ? " with " + paths_text(spec) : std::string()) + " on "
              + std::to_string(graphs.size()) + " graphs: add graph_path"
              + (shared && !spec.annotation_path ? " and annotation_path" : "");
    }
    throw SelectionRefusal(where + error, "pair_ambiguous", candidates_of(matching));
}

// A string member of |object| (|where| names it in the refusal): none when absent
std::optional<std::string> string_member(const Json::Value &object, const char *key,
                                         const std::string &where, const char *expected) {
    if (!object.isMember(key))
        return std::nullopt;
    const Json::Value &value = object[key];
    if (!value.isString())
        throw InvalidRequest(where + "." + key + ": expected " + expected);
    return value.asString();
}

// `pairs`: each object's one pair, a pair named twice once, in the envelope's order
std::vector<SelectedPair> select_listed_pairs(const Json::Value &pairs,
                                              const GraphIndexes &indexes) {
    if (!pairs.isArray() || pairs.empty()) {
        throw InvalidRequest("request.pairs: expected a non-empty array of pair objects "
                             "{graph, graph_path?, annotation_path?} (GET /capabilities lists "
                             "the pairs)");
    }
    std::set<std::pair<std::string, GraphPair>> named;
    for (Json::ArrayIndex i = 0; i < pairs.size(); ++i) {
        const std::string where = "request.pairs[" + std::to_string(i) + "]";
        const Json::Value &object = pairs[i];
        if (!object.isObject()) {
            throw InvalidRequest(where + ": expected an object {graph, graph_path?, "
                                 "annotation_path?}");
        }
        for (const std::string &key : object.getMemberNames()) {
            if (key != "graph" && key != "graph_path" && key != "annotation_path") {
                throw InvalidRequest(where + ": unknown member " + in_quotes(key)
                                     + " (graph, graph_path, annotation_path)");
            }
        }
        PairSpec spec;
        const std::optional<std::string> name
                = string_member(object, "graph", where, "a graph name");
        if (!name)
            throw InvalidRequest(where + ".graph: required (the pair's graph name)");
        spec.name = *name;
        spec.graph_path = string_member(object, "graph_path", where, "a string");
        spec.annotation_path = string_member(object, "annotation_path", where, "a string");
        SelectedPair pair = one_pair_of(spec, indexes, where + ": ", in_quotes(spec.name));
        named.emplace(std::move(pair.name), std::move(pair.pair));
    }
    // the envelope's order: the names in byte order (the set's), each name's pairs in the
    // list's order
    std::vector<std::string> names;
    for (const auto &[name, _] : named) {
        if (names.empty() || names.back() != name)
            names.push_back(name);
    }
    std::vector<SelectedPair> selected;
    for (SelectedPair &pair : pairs_of(names, indexes)) {
        if (named.count({ pair.name, pair.pair }))
            selected.push_back(std::move(pair));
    }
    return selected;
}

// `graphs` of RESOLVE, PATTERN and TRAVERSE: a non-empty array of known names, each once, in
// byte order; the first unknown name in the list's order is refused (graph_unknown). Each
// name is checked once, in the list's map: the check runs on the HTTP thread before any pair
// is queued, and a list may repeat a name many times over
std::vector<std::string> listed_names(const Json::Value &graphs, const GraphIndexes &indexes) {
    if (!graphs.isArray() || graphs.empty()) {
        throw InvalidRequest("request.graphs: expected a non-empty array of graph names "
                             "(GET /capabilities lists them)");
    }
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (const Json::Value &name : graphs) {
        if (!name.isString()) {
            throw InvalidRequest("request.graphs: expected a non-empty array of graph names "
                                 "(GET /capabilities lists them)");
        }
        std::string value = name.asString();
        if (!seen.insert(value).second)
            continue;
        check_known(value, indexes, "request.graphs: ");
        names.push_back(std::move(value));
    }
    std::sort(names.begin(), names.end());
    return names;
}

// TRAVERSE runs on one pair: a selection of several is refused with them as the candidates;
// |selection| says what selected them, |lever| how to name one
void check_single_pair(const std::vector<SelectedPair> &selected, const std::string &selection,
                       const std::string &lever) {
    if (selected.size() == 1)
        return;
    throw SelectionRefusal("request: " + selection + " " + std::to_string(selected.size())
                           + " pairs; a traverse runs on exactly one pair: " + lever,
                           "traverse_single_pair", candidates_of(selected));
}

const char kNamePairLever[] = "name it with pairs (graph, graph_path, annotation_path)";

// RESOLVE, PATTERN and TRAVERSE: `pairs`, or the names of `graphs` (`graph` beside them on
// RESOLVE and TRAVERSE, with `graph_path`), or every name up to the cap
std::vector<SelectedPair> select_posted(SelectionRoute route, const Json::Value &request,
                                        const GraphIndexes &indexes, size_t max_names) {
    const bool by_graph = route != SelectionRoute::PATTERN;
    const bool traverse = route == SelectionRoute::TRAVERSE;
    if (request.isMember("pairs")) {
        for (const char *selector : { "graphs", "graph", "graph_path" }) {
            if (!request.isMember(selector) || (!by_graph && std::string(selector) != "graphs"))
                continue;
            const std::string name = selector;
            throw SelectionRefusal(
                    "request: " + name + " and pairs are exclusive: "
                    + (name == "graphs" ? "name graphs, or list pairs, not both"
                       : name == "graph" ? "name a graph, or list pairs, not both"
                                         : "give each pair's graph_path in its object, not both"),
                    "selection_conflict");
        }
        std::vector<SelectedPair> selected = select_listed_pairs(request["pairs"], indexes);
        if (traverse)
            check_single_pair(selected, "pairs names", "list one");
        return selected;
    }

    std::optional<std::string> graph_path;
    std::vector<std::string> names;
    std::string selection;
    if (by_graph && request.isMember("graph")) {
        if (request.isMember("graphs")) {
            throw SelectionRefusal("request: graph and graphs are exclusive: name one graph "
                                   "as graph, or the graphs as graphs, not both",
                                   "selection_conflict");
        }
        const std::optional<std::string> name
                = string_member(request, "graph", "request", "a graph name");
        check_known(*name, indexes, "request.graph: ");
        names = { *name };
        selection = "graph " + in_quotes(*name) + " holds";
    } else if (request.isMember("graphs")) {
        names = listed_names(request["graphs"], indexes);
        selection = names.size() == 1 ? "graph " + in_quotes(names.front()) + " holds"
                                      : "graphs name";
    }
    if (by_graph)
        graph_path = string_member(request, "graph_path", "request", "a string");
    if (graph_path) {
        // graph_path picks among one name's pairs: the alias of a pair object without
        // annotation_path
        if (names.empty()) {
            throw InvalidRequest("request.graph_path: picks a graph of the name in graph: add "
                                 "graph, or list pairs");
        }
        if (names.size() > 1) {
            throw SelectionRefusal("request: graph_path picks a graph of one name, and graphs "
                                   "names " + std::to_string(names.size()) + ": name one, or "
                                   "list pairs", "selection_conflict");
        }
        PairSpec spec;
        spec.name = names.front();
        spec.graph_path = graph_path;
        std::vector<SelectedPair> selected = matching_pairs(spec, indexes);
        if (selected.empty()) {
            throw SelectionRefusal("request: no pair of graph " + in_quotes(spec.name) + " has "
                                   + paths_text(spec), "pair_not_found",
                                   selection_candidates(spec.name, indexes));
        }
        if (traverse) {
            check_single_pair(selected, "graph " + in_quotes(spec.name) + " holds",
                              kNamePairLever);
        }
        return selected;
    }
    if (names.empty()) {
        check_names_cap(indexes, max_names, "request.graphs");
        names = sorted_names(indexes);
        selection = "this server holds";
    }
    std::vector<SelectedPair> selected = pairs_of(names, indexes);
    if (traverse) {
        check_single_pair(selected, selection, selection == "this server holds"
                ? "name it with graph, or with pairs (graph, graph_path, annotation_path)"
                : kNamePairLever);
    }
    return selected;
}

} // namespace

std::vector<SelectedPair> select_pairs(SelectionRoute route,
                                       const Json::Value &request,
                                       const GraphIndexes &indexes,
                                       size_t max_names_without_selection) {
    const Json::Value &object = object_or_empty(request);
    switch (route) {
        case SelectionRoute::SEARCH: {
            if (object.isMember("pairs")) {
                throw SelectionRefusal("request: /search selects by graphs only (it has no "
                                       "re-ask): pairs is not accepted", "pairs_unsupported");
            }
            // its own message for an unknown name, as before, now with the code
            if (object.isMember("graphs") && object["graphs"].isArray()) {
                for (const Json::Value &name : object["graphs"]) {
                    if (!indexes.count(name.asString())) {
                        throw SelectionRefusal("Request with an uninitialized graph "
                                               + name.asString(), "graph_unknown");
                    }
                }
            }
            std::vector<std::string> names = filter_graphs_from_list(
                    indexes, object, 0, max_names_without_selection);
            // Without `graphs` the names come in the hash map's order, but the rows and cuts
            // of a multi-graph answer follow the index's pair order (section 3), which starts
            // from the names in byte order
            std::sort(names.begin(), names.end());
            return pairs_of(names, indexes);
        }
        case SelectionRoute::RESOLVE:
        case SelectionRoute::PATTERN:
        case SelectionRoute::TRAVERSE:
            return select_posted(route, object, indexes, max_names_without_selection);
        case SelectionRoute::CAPABILITIES: {
            if (object.isMember("pairs")) {
                throw SelectionRefusal("request: a capabilities probe selects one pair by "
                                       "graph, graph_path and annotation_path: pairs is not "
                                       "accepted", "pairs_unsupported");
            }
            // the members that spell the probe's parameters; any other is the request's
            std::vector<std::pair<std::string, std::string>> parameters;
            for (const char *key : { "graph", "graph_path", "annotation_path" }) {
                if (!object.isMember(key))
                    continue;
                if (!object[key].isString()) {
                    throw InvalidRequest(std::string("Bad request: '") + key
                                         + "': expected a string");
                }
                parameters.emplace_back(key, object[key].asString());
            }
            return select_pairs_from_parameters(route, parameters, indexes,
                                                max_names_without_selection);
        }
    }
    throw std::logic_error("select_pairs: unknown route");
}

std::vector<SelectedPair>
select_pairs_from_parameters(SelectionRoute route,
                             const std::vector<std::pair<std::string, std::string>> &query,
                             const GraphIndexes &indexes,
                             size_t max_names_without_selection) {
    const std::string get = "GET " + probe_route(route);
    PairSpec spec;
    std::optional<std::string> name;
    for (const auto &[key, value] : query) {
        std::optional<std::string> *slot = key == "graph" ? &name
                                         : key == "graph_path" ? &spec.graph_path
                                         : key == "annotation_path" ? &spec.annotation_path
                                                                    : nullptr;
        if (!slot) {
            throw InvalidRequest("Bad request: unknown parameter " + in_quotes(key) + " of " + get
                                 + " (graph, graph_path, annotation_path)");
        }
        if (*slot)
            throw InvalidRequest("Bad request: the parameter " + in_quotes(key) + " is given twice");
        *slot = value;
    }
    if (!name) {
        if (spec.graph_path || spec.annotation_path) {
            throw InvalidRequest("Bad request: " + get + ": graph_path and annotation_path pick "
                                 "a pair of the name ?graph= names: add ?graph=<name>");
        }
        // every pair: the envelope over them answers the probe (section 1)
        if (indexes.size() > max_names_without_selection) {
            throw std::invalid_argument(
                    "Bad request: " + get + " needs ?graph=<name> on this server, which hosts "
                    + std::to_string(indexes.size()) + " graph names (more than "
                    + std::to_string(max_names_without_selection)
                    + "; GET /capabilities lists them)");
        }
        return pairs_of(sorted_names(indexes), indexes);
    }
    spec.name = *name;
    return { one_pair_of(spec, indexes, "", "graph " + in_quotes(spec.name)) };
}

} // namespace cli
} // namespace mtg
