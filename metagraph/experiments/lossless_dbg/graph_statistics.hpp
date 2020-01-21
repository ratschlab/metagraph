//
// Created by studenyj on 6/6/19.
//

#ifndef __GRAPH_STATISTICS_HPP__
#define __GRAPH_STATISTICS_HPP__

#define STATS_INCOMING_HISTOGRAM 1u
#define STATS_OUTGOING_HISTOGRAM 2u
#include <nlohmann/json.hpp>
#include "graph_patch.hpp"
#include "utilities.hpp"
#include "unix_tools.hpp"

using namespace std;
using namespace nlohmann;
using namespace std::string_literals;

json get_statistics(DBGSuccinct& graph,int64_t verbosity=~0) {
    VerboseTimer graph_statistics_timer ("computation of graph statistics");
    PRINT_VAR(verbosity,verbosity & STATS_INCOMING_HISTOGRAM,verbosity & STATS_OUTGOING_HISTOGRAM)(std::string());
    int64_t joins = 0;
    int64_t splits = 0;
    int64_t incoming_histogram[6] = {0,0,0,0,0,0};
    int64_t outgoing_histogram[6] = {0,0,0,0,0,0};
#pragma omp parallel for reduction(+:joins,splits,incoming_histogram[:6],outgoing_histogram[:6])
    for(uint64_t node = 1; node <= graph.num_nodes();node++) {
        int64_t indegree = graph.indegree(node);
        int64_t outdegree = graph.outdegree(node);
        if (verbosity & STATS_INCOMING_HISTOGRAM) {
            incoming_histogram[indegree]++;
        }
        if (verbosity & STATS_OUTGOING_HISTOGRAM) {
            outgoing_histogram[outdegree]++;
        }
        if (indegree > 1) {
            joins++;
        }
        if (outdegree > 1) {
            splits++;
        }
    }
    json result{{"joins",joins},
                {"splits",splits},
                {"incoming_histogram", incoming_histogram},
                {"outgoing_histogram", outgoing_histogram},
                {"num_of_nodes", graph.num_nodes()}
    };
    graph_statistics_timer.finished();
    return result;
}

#endif // __GRAPH_STATISTICS_HPP__
