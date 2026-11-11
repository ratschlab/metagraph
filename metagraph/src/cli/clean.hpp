#ifndef __CLEAN_GRAPH_HPP__
#define __CLEAN_GRAPH_HPP__


#include <memory>
#include <utility>

#include <sdsl/int_vector.hpp>


namespace mtg {

namespace graph {
class DeBruijnGraph;
} // namespace graph

namespace cli {

class Config;

int clean_graph(Config *config);

// Same, for a graph that is already loaded from config->fnames[0]
// With |succ_pred| (basic succinct graphs only), also fill the succ_is_next and
// pred_is_prev bitmaps of the SSHash graph built from the written sequences.
int clean_graph(Config *config, std::shared_ptr<graph::DeBruijnGraph> graph,
                std::pair<sdsl::bit_vector, sdsl::bit_vector> *succ_pred = nullptr);

} // namespace cli
} // namespace mtg

#endif // __CLEAN_GRAPH_HPP__
