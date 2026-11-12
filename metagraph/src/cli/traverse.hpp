#ifndef __CLI_TRAVERSE_HPP__
#define __CLI_TRAVERSE_HPP__

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <json/json.h>

#include "strict_json.hpp"


class ThreadPool;

namespace mtg {

namespace graph {
class AnnotatedDBG;
} // namespace graph

namespace cli {

class Config;

// The server's bounds of a traversal request (--traverse-max-*)
struct TraverseLimits {
    double max_time_ms = 30'000;     // a seed's time budget is clamped to it
    size_t max_seeds = 10;
};
TraverseLimits traverse_limits(const Config &config);

// The longest seed of a /traverse request
constexpr uint64_t kMaxSeedBp = 10'000;
// The largest labels.max_seed_labels of a /traverse request
constexpr uint64_t kMaxSeedLabels = 100'000;

// The request and response schemas are documented in the manual (docs/source/traverse.rst)
Json::Value process_traverse_request(const Json::Value &request,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const TraverseLimits &limits);

// One index of a server started with a list of graphs
struct TraverseIndex {
    std::string graph;  // its name in the list
    const graph::AnnotatedDBG *anno_graph;
};

// The same request on several indexes, walked in |pool|. Each seed gets one walk per index
// that holds the seed and at least one of its labels.
Json::Value process_traverse_request(const Json::Value &request,
                                     const std::vector<TraverseIndex> &indexes,
                                     const TraverseLimits &limits,
                                     ThreadPool &pool);

} // namespace cli
} // namespace mtg

#endif // __CLI_TRAVERSE_HPP__
