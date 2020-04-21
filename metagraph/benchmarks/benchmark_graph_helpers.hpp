#ifndef __BENCHMARK_GRAPH_HELPERS_HPP__
#define __BENCHMARK_GRAPH_HELPERS_HPP__

#include <memory>
#include <string>


// TODO: move these into mg namespace later
class DeBruijnGraph;
class AnnotatedDBG;

namespace annotate {
template <typename Label>
class ColumnCompressed;
}


namespace mg {
namespace bm {

std::shared_ptr<DeBruijnGraph> build_graph(const std::string &filename);

template <class Annotation = annotate::ColumnCompressed<std::string>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename);

std::unique_ptr<AnnotatedDBG> build_query_graph(const AnnotatedDBG &anno_graph,
                                                const std::string &query_filename);

} //namespace bm
} //namespace mg

#endif // __BENCHMARK_GRAPH_HELPERS_HPP__
