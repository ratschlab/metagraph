#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <string>

class AnnotatedDBG;
class IDBGAligner;
class Config;

void execute_query(const std::string &seq_name,
                   const std::string &sequence,
                   bool count_labels,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream,
                   IDBGAligner *aligner = nullptr);

int query_graph(Config *config);

#endif // __QUERY_GRAPH_HPP__
