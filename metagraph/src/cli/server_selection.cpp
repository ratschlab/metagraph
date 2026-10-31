#include "server_selection.hpp"

#include <algorithm>
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

// The one pair |name| (and |graph_path|) addresses, as selected
std::vector<SelectedPair> one_pair(const std::string &name, const std::string *graph_path,
                                   const GraphIndexes &indexes) {
    return { SelectedPair { name, select_traverse_pair(name, graph_path, indexes) } };
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

// Section 6's `pairs` selector is not served yet (E12); /search never takes it
void refuse_pairs(SelectionRoute route, const Json::Value &request) {
    if (!request.isObject() || !request.isMember("pairs"))
        return;
    if (route == SelectionRoute::SEARCH) {
        throw SelectionRefusal("request: /search selects by graphs only (it has no re-ask): "
                               "pairs is not accepted", "pairs_unsupported");
    }
    throw SelectionRefusal("request: pairs is not served by this build: select by graphs",
                           "pairs_unsupported");
}

} // namespace

std::vector<SelectedPair> select_pairs(SelectionRoute route,
                                       const Json::Value &request,
                                       const GraphIndexes &indexes,
                                       size_t max_names_without_selection) {
    refuse_pairs(route, request);
    switch (route) {
        case SelectionRoute::SEARCH: {
            std::vector<std::string> names = filter_graphs_from_list(
                    indexes, object_or_empty(request), 0, max_names_without_selection);
            // Without `graphs` the names come in the hash map's order, but the rows and cuts
            // of a multi-graph answer follow the index's pair order (section 3), which starts
            // from the names in byte order. POST /search itself still walks its names in the
            // map's order and queries a pair listed twice under a name twice.
            std::sort(names.begin(), names.end());
            return pairs_of(names, indexes);
        }
        case SelectionRoute::PATTERN: {
            return pairs_of(pattern_graph_names(object_or_empty(request), sorted_names(indexes),
                                                max_names_without_selection), indexes);
        }
        case SelectionRoute::RESOLVE:
        case SelectionRoute::TRAVERSE: {
            // today /resolve reads one pair as /traverse does; its fan-out over the names of
            // `graphs` (section 6) comes with its envelope
            const GraphSelection selection = traverse_graph_selection(request, true);
            return one_pair(selection.name,
                            selection.graph_path ? &*selection.graph_path : nullptr, indexes);
        }
        case SelectionRoute::CAPABILITIES: {
            // the members that spell the probe's parameters; any other is the request's
            std::vector<std::pair<std::string, std::string>> parameters;
            if (request.isObject()) {
                for (const char *key : { "graph", "graph_path" }) {
                    if (!request.isMember(key))
                        continue;
                    if (!request[key].isString()) {
                        throw InvalidRequest(std::string("Bad request: '") + key
                                             + "': expected a string");
                    }
                    parameters.emplace_back(key, request[key].asString());
                }
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
                             size_t /* max_names_without_selection */) {
    // Today a probe describes one pair and needs ?graph=. The names cap is not read until a
    // probe without ?graph= answers the envelope over every pair (section 1).
    const auto [name, graph_path] = probe_parameters(query, probe_route(route));
    return one_pair(name, graph_path ? &*graph_path : nullptr, indexes);
}

} // namespace cli
} // namespace mtg
