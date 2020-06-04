#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <cstdlib>

class DeBruijnGraph;
class BOSS;


namespace mtg {
namespace cli {

class Config;


void print_stats(const DeBruijnGraph &graph);

void print_boss_stats(const BOSS &boss_graph,
                      bool count_dummy = false,
                      size_t num_threads = 0,
                      bool verbose = false);

int print_stats(Config *config);

int compare(Config *config);

} // namespace cli
} // namespace mtg

#endif // __STATS_HPP__
