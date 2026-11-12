#include "pattern.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "config/config.hpp"
#include "strict_json.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;

namespace {

struct PatternRequest {
    std::vector<std::pair<Json::Value /* id */, std::string>> patterns;
    uint64_t max_contexts;
    std::chrono::milliseconds time_budget;
};

PatternRequest parse_pattern_request(Strict &root, double max_time_ms) {
    PatternRequest request;
    const Json::Value &patterns = root.raw("patterns");
    if (!patterns.isArray() || patterns.empty() || patterns.size() > kMaxPatterns) {
        throw InvalidRequest(root.child("patterns") + ": expected an array of 1 to "
                             + std::to_string(kMaxPatterns) + " patterns");
    }
    for (Json::ArrayIndex i = 0; i < patterns.size(); ++i) {
        Strict pattern(patterns[i], root.child("patterns") + "[" + std::to_string(i) + "]");
        if (!pattern.has("sequence"))
            throw InvalidRequest(pattern.child("sequence") + ": missing");
        Json::Value id = pattern.has("id") ? Json::Value(pattern.str("id", "")) : Json::Value(i);
        request.patterns.emplace_back(std::move(id), pattern.str("sequence", ""));
        pattern.done();
    }
    request.max_contexts = root.uint("max_contexts", 1000, 0, kMaxPatternContexts);
    request.time_budget = std::chrono::milliseconds(static_cast<int64_t>(std::min(
        root.number("time_budget_ms", max_time_ms, 0, std::numeric_limits<double>::max()),
        max_time_ms
    )));
    return request;
}

const DBGSuccinct& get_dbg_succ(const DeBruijnGraph &graph) {
    const DeBruijnGraph *base = &graph;
    // a primary graph is loaded under a wrapper: the stored graph is searched
    if (const auto *wrapper = dynamic_cast<const DBGWrapper<>*>(base))
        base = &wrapper->get_graph();

    if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base))
        return *dbg_succ;

    throw InvalidRequest("pattern search needs a succinct graph");
}

struct Context {
    bool reverse;
    size_t offset;
    DeBruijnGraph::node_index node;
    std::string kmer;
};

// the answer for one pattern on one graph: count, complete, stop and contexts
Json::Value search_pattern(const DBGSuccinct &graph,
                           const std::string &sequence,
                           uint64_t max_contexts,
                           boss::SearchBudget *budget) {
    if (sequence.empty() || sequence.size() > graph.get_k()) {
        throw InvalidRequest("pattern '" + sequence + "': expected 1 to "
                             + std::to_string(graph.get_k()) + " symbols");
    }
    // the pattern, and its reverse complement if the graph does not hold both strands
    std::vector<boss::Pattern> strands { boss::parse_pattern(graph.get_boss(), sequence) };
#if ! _PROTEIN_GRAPH
    if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        boss::Pattern reverse = boss::reverse_complement_pattern(strands[0]);
        if (reverse != strands[0])
            strands.push_back(std::move(reverse));
    }
#endif

    uint64_t min = 0;
    uint64_t max = 0;
    bool complete = true;
    for (size_t i = 0; i < strands.size() && complete; ++i) {
        uint64_t strand_min, strand_max;
        complete = graph.count_pattern(strands[i], budget, &strand_min, &strand_max);
        min += strand_min;
        max += strand_max;
    }

    Json::Value result;
    if (complete && max <= max_contexts) {
        std::vector<Context> contexts;
        for (size_t i = 0; i < strands.size() && complete; ++i) {
            complete = graph.call_pattern_contexts(strands[i], budget,
                [&](DeBruijnGraph::node_index node, const std::string &kmer, size_t offset) {
#if ! _PROTEIN_GRAPH
                    if (i) {
                        // a k-mer equal to its reverse complement was listed on the first strand
                        std::string reverse = kmer;
                        reverse_complement(reverse.begin(), reverse.end());
                        if (reverse == kmer)
                            return;
                    }
#endif
                    contexts.push_back({ i > 0, offset, node, kmer });
                }
            );
        }
        if (complete) {
            std::sort(contexts.begin(), contexts.end(), [](const Context &a, const Context &b) {
                return std::tie(a.reverse, a.offset, a.node) < std::tie(b.reverse, b.offset, b.node);
            });
            min = max = contexts.size();
            result["contexts"] = Json::arrayValue;
            for (const Context &context : contexts) {
                Json::Value &c = result["contexts"].append(Json::objectValue);
                c["kmer"] = context.kmer;
                c["offset"] = static_cast<Json::UInt64>(context.offset);
                c["strand"] = context.reverse ? "-" : "+";
            }
        }
    }
    result["count"]["min"] = static_cast<Json::UInt64>(min);
    // null if the search was stopped
    result["count"]["max"] = complete ? Json::Value(static_cast<Json::UInt64>(max)) : Json::Value();
    result["complete"] = complete;
    if (!complete)
        result["stop"] = budget->stop();

    return result;
}

} // namespace

Json::Value process_pattern_request(const Json::Value &request,
                                    const DeBruijnGraph &graph,
                                    double max_time_ms) {
    Strict root(request, "request");
    const PatternRequest parsed = parse_pattern_request(root, max_time_ms);
    root.done();

    const DBGSuccinct &dbg_succ = get_dbg_succ(graph);
    // one budget for all patterns of the request
    boss::SearchBudget budget(kMaxPatternSteps, parsed.time_budget);

    Json::Value out;
    out["k"] = static_cast<Json::UInt64>(dbg_succ.get_k());
    out["mode"] = Config::graphmode_to_string(dbg_succ.get_mode());
    out["mask"] = dbg_succ.get_mask() != nullptr;
    out["results"] = Json::arrayValue;
    for (const auto &[id, sequence] : parsed.patterns) {
        Json::Value result = search_pattern(dbg_succ, sequence, parsed.max_contexts, &budget);
        result["id"] = id;
        result["length"] = static_cast<Json::UInt64>(sequence.size());
        out["results"].append(std::move(result));
    }
    return out;
}

Json::Value process_pattern_request(const Json::Value &request,
                                    const std::vector<PatternGraph> &graphs,
                                    double max_time_ms) {
    Strict root(request, "request");
    const PatternRequest parsed = parse_pattern_request(root, max_time_ms);
    root.has("graphs");  // the server selected |graphs| by it
    root.done();

    // one budget for all patterns and graphs of the request
    boss::SearchBudget budget(kMaxPatternSteps, parsed.time_budget);

    Json::Value out;
    out["results"] = Json::arrayValue;
    for (const auto &[id, sequence] : parsed.patterns) {
        Json::Value &result = out["results"].append(Json::objectValue);
        result["id"] = id;
        result["length"] = static_cast<Json::UInt64>(sequence.size());
        result["graphs"] = Json::arrayValue;
        for (const PatternGraph &graph : graphs) {
            Json::Value found = search_pattern(get_dbg_succ(*graph.graph), sequence,
                                               parsed.max_contexts, &budget);
            found["graph"] = graph.name;
            result["graphs"].append(std::move(found));
        }
    }
    return out;
}

} // namespace cli
} // namespace mtg
