#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "annotated_dbg.hpp"
#include "bitmap.hpp"
#include "dbg_aligner.hpp"
#include "threading.hpp"
#include "aligner_helper.hpp"
#include "aligner_methods.hpp"


typedef std::function<uint64_t()> UInt64Callback;

namespace annotated_graph_algorithm {

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    std::function<bool(uint64_t, uint64_t)> keep_node,
                    double lazy_evaluation_density_cutoff = 0.05);

// Allows for lazy evaluation of the counts passed to keep_node
std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    std::function<bool(UInt64Callback, UInt64Callback)> keep_node,
                    double lazy_evaluation_density_cutoff = 0.05);

template <class Index, typename... Args>
using VariantCallback = std::function<void(Alignment<Index>&&,
                                           const std::string&, // query sequence
                                           Args&&...)>;

typedef Alignment<DeBruijnGraph::node_index> DBGAlignment;
typedef VariantCallback<DeBruijnGraph::node_index,
                        std::vector<AnnotatedDBG::Annotator::Label>&&> VariantLabelCallback;


// These functions will callback nothing if graph is equal to the graph stored
// in anno_graph

void call_bubbles(const DeBruijnGraph &graph,
                  const AnnotatedDBG &anno_graph,
                  const VariantLabelCallback &callback,
                  ThreadPool *thread_pool = nullptr,
                  const std::function<bool()> &terminate = []() { return false; });

void call_breakpoints(const DeBruijnGraph &graph,
                      const AnnotatedDBG &anno_graph,
                      const VariantLabelCallback &callback,
                      ThreadPool *thread_pool = nullptr,
                      const std::function<bool()> &terminate = []() { return false; });

// Given a seed, construct a masked graph that only includes the labels of the seed node,
// then extend on that.
Extender<DeBruijnGraph::node_index>
build_masked_graph_extender(const AnnotatedDBG &anno_graph,
                            double seed_label_discovery_fraction = 1.0,
                            Extender<DeBruijnGraph::node_index>&& extender
                                = default_extender<DeBruijnGraph::node_index>);



} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
