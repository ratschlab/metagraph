#ifndef __BENCHMARK_GRAPH_HELPERS_HPP__
#define __BENCHMARK_GRAPH_HELPERS_HPP__

#include <memory>
#include <string>

class AnnotatedDBG;

namespace annotate {
template <typename Label>
class ColumnCompressed;
};


template <class Annotation = annotate::ColumnCompressed<std::string>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename);

std::unique_ptr<AnnotatedDBG> build_query_graph(const AnnotatedDBG &anno_graph,
                                                const std::string &query_filename);

#endif // __BENCHMARK_GRAPH_HELPERS_HPP__
