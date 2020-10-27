#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <cstdlib>

namespace mtg {

namespace graph {
    class DeBruijnGraph;
    namespace boss {
        class BOSS;
    } // namespace boss
} // namespace graph

namespace cli {

class Config;


void print_stats(const graph::DeBruijnGraph &graph);

void print_boss_stats(const graph::boss::BOSS &boss_graph,
                      bool count_dummy = false,
                      size_t num_threads = 0);

int print_stats(Config *config);

int compare(Config *config);

} // namespace cli
} // namespace mtg

#endif // __STATS_HPP__
