#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <vector>
#include <functional>

#include "annotated_dbg.hpp"
#include "bit_vector.hpp"
#include "masked_graph.hpp"


namespace annotated_graph_algorithm {

MaskedDeBruijnGraph
mask_insignificant_nodes(std::shared_ptr<DeBruijnGraph> graph,
                         const AnnotatedDBG &anno_graph,
                         const std::vector<AnnotatedDBG::Annotator::Label> &ingroup,
                         const std::vector<AnnotatedDBG::Annotator::Label> &outgroup,
                         const std::function<bool(uint64_t, uint64_t)> &is_significant,
                         const std::vector<const bit_vector*> &ingroup_masks = {},
                         const std::vector<const bit_vector*> &outgroup_masks = {});

} // namespace annotated_graph_algorithm

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
