#ifndef __BENCHMARK_GRAPH_HELPERS_HPP__
#define __BENCHMARK_GRAPH_HELPERS_HPP__

#include <memory>
#include <string>

class AnnotatedDBG;


std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename);

#endif // __BENCHMARK_GRAPH_HELPERS_HPP__
