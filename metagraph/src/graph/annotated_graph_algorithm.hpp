#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "annotated_dbg.hpp"
#include "bitmap.hpp"
#include "masked_graph.hpp"


namespace annotated_graph_algorithm {

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(uint64_t, uint64_t)> &keep_node);


} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
