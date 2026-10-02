#ifndef __METAGRAPH_CLI_TRAVERSE_HPP__
#define __METAGRAPH_CLI_TRAVERSE_HPP__

#include <string>
#include <vector>

#include <json/json.h>

#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"


namespace mtg {

namespace graph {
class AnnotatedDBG;
}

namespace cli {

class Config;

// A request error: reported as HTTP 400 by the server, exit code 1 by the CLI.
class InvalidRequest : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

// Label-change cost as requested, by label NAME (converted per seed to LabelChangeCost
// over the seed's label dictionary).
struct CostSpec {
    enum Model { FORBID, CONSTANT, TABLE };
    Model model = FORBID;
    double value = 0;
    double default_cost = graph::traversal::kInfiniteLoss;
    std::vector<std::tuple<std::string, std::string, double>> entries;
};

struct TraverseRequest {
    std::string release;
    std::string graph;        // multi-graph mode: index name (resolved by the server)
    std::string graph_path;   // optional disambiguation when a name spans several shards
    std::vector<graph::traversal::Seed> seeds;
    graph::traversal::Strategy strategy;
    CostSpec cost;
    std::string detail = "full";          // summary | tree | full
    bool timing = true;
};

struct ResolveRequest {
    std::string graph;
    std::string graph_path;
    std::string sequence;
    graph::traversal::ResolveOptions options;
    bool select = false;
    graph::traversal::SelectionPolicy policy;
    std::string run_format = "intervals"; // intervals | rle
    uint64_t max_query_bp = 1'000'000;
};

// Strict parsing: unknown keys, wrong types and out-of-range values throw InvalidRequest
// naming the offending field. Defaults are filled in.
TraverseRequest parse_traverse_request(const Json::Value &json);
ResolveRequest parse_resolve_request(const Json::Value &json);

// Normalized strategy echo (what was actually used)
Json::Value strategy_to_json(const graph::traversal::Strategy &strategy, const CostSpec &cost);
Json::Value seed_result_to_json(const graph::traversal::SeedResult &result,
                                const graph::traversal::Strategy &strategy,
                                const std::string &detail, bool timing);
Json::Value profile_to_json(const graph::traversal::SupportProfile &profile,
                            const graph::traversal::SeedSelection *selection,
                            const std::string &run_format);
Json::Value capabilities_to_json(const graph::traversal::LabelOracle &oracle,
                                 const std::string &release);

// Server-side caps on POST /traverse, 0 = unlimited. A request that names no labels
// makes the SERVER choose the permitted set, so the request's own bounds no longer
// describe the work it asks for; these do. |max_time_ms| also bounds the derivation
// itself, and a request asking for less is left alone (a clamp only ever lowers a
// budget, never raises one) — see Config::traverse_max_*.
struct TraverseLimits {
    double max_time_ms = 0;        // cap on strategy.bounds.time_budget_ms
    size_t max_seeds = 0;          // cap on |request.seeds|
    uint64_t max_seed_bp = 0;      // cap on the length of one seed
    size_t max_seed_labels = 0;    // cap on strategy.labels.max_seed_labels
};

// Shared by the CLI and the server. |release| is the configured index release id ("" if none).
Json::Value process_traverse_request(const Json::Value &json,
                                     const graph::AnnotatedDBG &anno_graph,
                                     const std::string &release,
                                     const TraverseLimits &limits = {});
Json::Value process_resolve_request(const Json::Value &json,
                                    const graph::AnnotatedDBG &anno_graph,
                                    const std::string &release,
                                    uint64_t max_query_bp = 0);

// CLI entry point: `metagraph traverse [--resolve] -i GRAPH -a ANNOTATION request.json ...`
int traverse_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_TRAVERSE_HPP__
