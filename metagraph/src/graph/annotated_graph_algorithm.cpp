#include "annotated_graph_algorithm.hpp"

#include "utils.hpp"


namespace annotated_graph_algorithm {

MaskedDeBruijnGraph
mask_insignificant_nodes(std::shared_ptr<DeBruijnGraph> graph,
                         const AnnotatedDBG &anno_graph,
                         const std::vector<AnnotatedDBG::Annotator::Label> &ingroup,
                         const std::vector<AnnotatedDBG::Annotator::Label> &outgroup,
                         const std::function<bool(uint64_t, uint64_t)> &is_significant,
                         const std::vector<const bit_vector*> &ingroup_masks,
                         const std::vector<const bit_vector*> &outgroup_masks) {
    assert(graph.get() == dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph()));

    if (!graph->num_nodes())
        return MaskedDeBruijnGraph(graph);

    std::unique_ptr<bit_vector_stat> mask(
        new bit_vector_stat(graph->num_nodes() + 1, false)
    );

    if (ingroup.empty() && ingroup_masks.empty())
        return MaskedDeBruijnGraph(graph, mask.release());

    // store incount and outcount interleaved
    sdsl::int_vector<> counts(
        (graph->num_nodes() + 1) << 1, 0,
        utils::code_length(std::max(ingroup.size() + ingroup_masks.size(),
                                    outgroup.size() + outgroup_masks.size()))
    );

    for (const auto *ingroup_mask : ingroup_masks) {
        if (ingroup_mask)
            ingroup_mask->call_ones([&](const auto &i) { counts[i << 1]++; });
    }

    for (const auto &inlabel : ingroup) {
        anno_graph.call_indices(inlabel,
                                [&](const auto &i) { counts[i << 1]++; });
    }

    for (const auto *outgroup_mask : outgroup_masks) {
        if (outgroup_mask)
            outgroup_mask->call_ones([&](const auto &i) { counts[(i << 1) | 0x1]++; });
    }

    if (ingroup.size()) {
        anno_graph.call_indices(
            ingroup[0],
            [&](const auto &graph_index) {
                assert(((graph_index << 1) + 1) < counts.size());
                auto incount = counts[graph_index << 1];
                auto outcount = counts[(graph_index << 1) + 1]
                    + anno_graph.count_labels(graph_index, outgroup);
                if (is_significant(incount, outcount))
                    mask->set(graph_index, true);
            }
        );
    } else {
        assert(ingroup_masks.size());
        ingroup_masks[0]->call_ones(
            [&](const auto &graph_index) {
                auto incount = counts[graph_index << 1];
                auto outcount = counts[(graph_index << 1) + 1]
                    + anno_graph.count_labels(graph_index, outgroup);
                if (is_significant(incount, outcount))
                    mask->set(graph_index, true);
            }
        );
    }

    return MaskedDeBruijnGraph(graph, mask.release());
}

} // namespace annotated_graph_algorithm
