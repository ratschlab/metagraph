#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "annotated_dbg.hpp"
#include "bitmap.hpp"
#include "masked_graph.hpp"
#include "threading.hpp"


typedef std::function<uint64_t()> UInt64Callback;

namespace annotated_graph_algorithm {

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(uint64_t, uint64_t)> &keep_node);

// Allows for lazy evaluation of the counts passed to keep_node
std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(UInt64Callback, UInt64Callback)> &keep_node);

template <class Index, class VLabels>
using IndexRefVarLabelCallback = std::function<void(const Index&,
                                                    const std::string&,
                                                    const std::string&,
                                                    VLabels&&)>;

typedef IndexRefVarLabelCallback<MaskedDeBruijnGraph::node_index,
                                 std::vector<AnnotatedDBG::Annotator::Label>>
    AnnotatedDBGIndexRefVarLabelsCallback;

void call_bubbles(const MaskedDeBruijnGraph &masked_graph,
                  const AnnotatedDBG &anno_graph,
                  const AnnotatedDBGIndexRefVarLabelsCallback &callback,
                  ThreadPool *thread_pool = nullptr,
                  const std::function<bool()> &terminate = []() { return false; });

void call_breakpoints(const MaskedDeBruijnGraph &masked_graph,
                      const AnnotatedDBG &anno_graph,
                      const AnnotatedDBGIndexRefVarLabelsCallback &callback,
                      ThreadPool *thread_pool = nullptr,
                      const std::function<bool()> &terminate = []() { return false; });

} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
