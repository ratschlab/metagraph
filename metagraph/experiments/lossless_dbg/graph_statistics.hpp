//
// Created by studenyj on 6/6/19.
//

#ifndef METAGRAPH_GRAPH_STATISTICS_HPP
#define METAGRAPH_GRAPH_STATISTICS_HPP

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
    Timer timer;
    cerr << "Starting computation of graph statistics" << endl;
    PRINT_VAR(verbosity,verbosity & STATS_INCOMING_HISTOGRAM,verbosity & STATS_OUTGOING_HISTOGRAM);
    int64_t joins = 0;
    int64_t splits = 0;
    int64_t incoming_histogram[6] = {0,0,0,0,0,0};
    int64_t outgoing_histogram[6] = {0,0,0,0,0,0};
#pragma omp parallel for reduction(+:joins,splits,incoming_histogram[:6],outgoing_histogram[:6])
    for (int64_t node = 1; node <= graph.num_nodes();node++) {
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
    cerr << "Computation of statistics finished in " << timer.elapsed() << " sec." << endl;
    return result;
}

#endif //METAGRAPH_GRAPH_STATISTICS_HPP
