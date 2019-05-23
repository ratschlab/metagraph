#include "annotated_graph_algorithm.hpp"

#include "annotate_column_compressed.hpp"
#include "utils.hpp"


namespace annotated_graph_algorithm {

uint64_t count_node_labels(const AnnotatedDBG &anno_graph,
                           const DeBruijnGraph::node_index &index,
                           const std::vector<AnnotatedDBG::Annotator::Label> &labels_to_check) {
    uint64_t count = 0;
    for (const auto &label : labels_to_check) {
        if (anno_graph.has_label(index, label))
            count++;
    }
    return count;
}

// TODO: optimize this
constexpr double density_cutoff = 0.1;

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(uint64_t, uint64_t)> &keep_node) {
    if (!anno_graph.get_graph().num_nodes())
        return {};

    std::unique_ptr<bitmap> mask(
        new bitmap_vector(anno_graph.get_graph().num_nodes() + 1, false)
    );

    if (mask_in.empty())
        return mask;

    // store counts interleaved
    size_t width = utils::code_length(std::max(mask_in.size(), mask_out.size()));
    size_t int_mask = (size_t(1) << width) - 1;
    sdsl::int_vector<> counts(mask->size(), 0, width * 2);

    const auto *columns = dynamic_cast<const annotate::ColumnCompressed<>*>(
        &anno_graph.get_annotation()
    );

    if (columns) {
        // Pick how to query annotation based on column density

        std::vector<AnnotatedDBG::Annotator::Label> mask_in_dense, mask_out_dense;
        size_t density_cutoff_count = mask->size() * density_cutoff;

        for (const auto &label_in : mask_in) {
            if (columns->get_column(label_in).num_set_bits() < density_cutoff_count) {
                anno_graph.call_annotated_nodes(
                    label_in,
                    [&](const auto &i) { counts[i]++; }
                );
            } else {
                mask_in_dense.push_back(label_in);
            }
        }

        counts.width(width);
        assert(counts.size() == mask->size() * 2);

        for (const auto &label_out : mask_out) {
            if (columns->get_column(label_out).num_set_bits() < density_cutoff_count) {
                anno_graph.call_annotated_nodes(
                    label_out,
                    [&](const auto &i) { counts[(i << 1) + 1]++; }
                );
            } else {
                mask_out_dense.push_back(label_out);
            }
        }

        counts.width(width * 2);
        assert(counts.size() == mask->size());

        for (size_t i = 1; i < counts.size(); ++i) {
            size_t count = counts[i];

            size_t count_in = (count & int_mask)
                + count_node_labels(anno_graph, i, mask_in_dense);

            size_t count_out = (count >> width)
                + count_node_labels(anno_graph, i, mask_out_dense);

            if (!count_in && !count_out)
                continue;

            if (keep_node(count_in, count_out))
                mask->set(i, true);
        }
    } else {
        // TODO: make this more efficient for row-major annotations
        for (const auto &label_in : mask_in) {
            anno_graph.call_annotated_nodes(
                label_in,
                [&](const auto &i) { counts[i]++; }
            );
        }

        counts.width(width);
        assert(counts.size() == mask->size() * 2);

        for (const auto &label_out : mask_out) {
            anno_graph.call_annotated_nodes(
                label_out,
                [&](const auto &i) { counts[(i << 1) + 1]++; }
            );
        }

        counts.width(width * 2);
        assert(counts.size() == mask->size());

        for (size_t i = 1; i < counts.size(); ++i) {
            size_t count = counts[i];
            if (!count)
                continue;

            if (keep_node(count & int_mask, count >> width))
                mask->set(i, true);
        }
    }

    return std::unique_ptr<bitmap>(mask.release());
}

} // namespace annotated_graph_algorithm
