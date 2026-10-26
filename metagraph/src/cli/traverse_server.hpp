#pragma once

#include <atomic>
#include <functional>
#include <span>
#include <string>
#include <utility>

#include "server_utils.hpp"
#include "traverse.hpp"

namespace mtg::cli {

struct GraphSelection {
    std::string name;
    std::optional<std::string> graph_path;
};

GraphSelection traverse_graph_selection(const Json::Value &request, bool multi_graph);

using GraphPair = std::pair<std::string, std::string>;
const GraphPair& select_traverse_pair(const std::string &name,
                                      const std::optional<std::string> &graph_path,
                                      std::span<const GraphPair> pairs);

// References to an index already owned by the server. They outlive its handlers.
struct TraversalIndex {
    const graph::AnnotatedDBG &graph;
    const IndexIdentity &identity;
    double prepare_ms = 0;
};

using TraversalLookup = std::function<TraversalIndex(const Json::Value &)>;

void add_traversal_routes(HttpServer &server, const Config &config,
                          TraversalLookup lookup, std::function<Json::Value()> catalog,
                          std::atomic<size_t> &request_ids,
                          const std::atomic<bool> &stopping);

} // namespace mtg::cli
