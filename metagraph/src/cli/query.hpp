#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>

class AnnotatedDBG;
class IDBGAligner;
class Config;

void execute_query(const std::string &seq_name,
                   const std::string &sequence,
                   bool count_labels,
                   bool print_signature,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream);


using StringGenerator = std::function<void(std::function<void(const std::string&)>)>;

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads);

int query_graph(Config *config);

#endif // __QUERY_GRAPH_HPP__
