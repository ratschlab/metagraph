#ifndef __ANNOTATE_GRAPH_HPP__
#define __ANNOTATE_GRAPH_HPP__


namespace mtg {
namespace cli {

class Config;

int annotate_graph(Config *config);

int annotate_graph_with_genome_coordinates(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ANNOTATE_GRAPH_HPP__
