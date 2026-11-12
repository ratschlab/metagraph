#ifndef __CLI_PATTERN_HPP__
#define __CLI_PATTERN_HPP__

#include <cstdint>
#include <string>
#include <vector>

#include <json/json.h>


namespace mtg {

namespace graph {
class DeBruijnGraph;
} // namespace graph

namespace cli {

// The most patterns of a /pattern request
constexpr uint64_t kMaxPatterns = 16;
// The largest max_contexts of a /pattern request
constexpr uint64_t kMaxPatternContexts = 10'000;
// The most search steps of a /pattern request
constexpr uint64_t kMaxPatternSteps = 10'000'000;

// The request and response schemas are documented in the manual (docs/source/pattern.rst).
// The request's time budget is clamped to |max_time_ms| (--pattern-max-time-ms).
Json::Value process_pattern_request(const Json::Value &request,
                                    const graph::DeBruijnGraph &graph,
                                    double max_time_ms);

// One graph of a server started with a list of graphs
struct PatternGraph {
    std::string name;  // its name in the list
    const graph::DeBruijnGraph *graph;
};

// The same request on several graphs, searched one after another
Json::Value process_pattern_request(const Json::Value &request,
                                    const std::vector<PatternGraph> &graphs,
                                    double max_time_ms);

} // namespace cli
} // namespace mtg

#endif // __CLI_PATTERN_HPP__
