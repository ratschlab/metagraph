#ifndef __CLI_TRAVERSE_HPP__
#define __CLI_TRAVERSE_HPP__

#include <cstdint>
#include <stdexcept>

#include <json/json.h>


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
    uint64_t max_resolve_query_bp = 1'000'000;   // the longest /resolve query
};
TraverseLimits traverse_limits(const Config &config);

// The longest seed of a /traverse request
constexpr uint64_t kMaxSeedBp = 10'000;
// The largest labels.max_seed_labels of a /traverse request
constexpr uint64_t kMaxSeedLabels = 100'000;

// A request that is not valid under the v1 schema (an unknown field included)
class InvalidRequest : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

// The v1 request and response schemas are documented in traverse.md
Json::Value process_traverse_request(const Json::Value &request,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const TraverseLimits &limits);
Json::Value process_resolve_request(const Json::Value &request,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const TraverseLimits &limits);
} // namespace cli
} // namespace mtg

#endif // __CLI_TRAVERSE_HPP__
