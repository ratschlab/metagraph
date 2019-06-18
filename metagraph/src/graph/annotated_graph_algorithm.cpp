#include "annotated_graph_algorithm.hpp"

#include "annotate_column_compressed.hpp"
#include "utils.hpp"

namespace annotated_graph_algorithm {

typedef DeBruijnGraph::node_index NodeIndex;

uint64_t count_node_labels(const AnnotatedDBG &anno_graph,
                           const NodeIndex &index,
                           const std::vector<AnnotatedDBG::Annotator::Label> &labels_to_check) {
    uint64_t count = 0;
    for (const auto &label : labels_to_check) {
        if (anno_graph.has_label(index, label))
            count++;
    }
    return count;
}

// TODO: optimize this
constexpr double density_cutoff = 0.05;

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(uint64_t, uint64_t)> &keep_node) {
    return mask_nodes_by_label(
        anno_graph,
        mask_in,
        mask_out,
        [&](UInt64Callback call_in, UInt64Callback call_out) {
            return keep_node(call_in(), call_out());
        }
    );
}

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(UInt64Callback, UInt64Callback)> &keep_node) {
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

            if (keep_node(
                    [&]() {
                        return (count & int_mask)
                            + count_node_labels(anno_graph, i, mask_in_dense);
                    },
                    [&]() {
                        return (count >> width)
                            + count_node_labels(anno_graph, i, mask_out_dense);
                    }
                ))
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
            if (keep_node([&]() { return count & int_mask; },
                          [&]() { return count >> width; }))
                mask->set(i, true);
        }
    }

    return std::unique_ptr<bitmap>(mask.release());
}

void
call_paths_from_branch(const MaskedDeBruijnGraph &masked_graph,
                       std::function<void(NodeIndex, NodeIndex, std::string&&)> callback,
                       std::function<bool(NodeIndex, NodeIndex, const std::string&)> stop_path,
                       std::function<bool()> terminate = []() { return false; }) {
    bit_vector_stat visited(masked_graph.num_nodes() + 1, false);

    masked_graph.call_nodes(
        [&](auto start) {
            if (visited[start] || masked_graph.unmasked_outdegree(start) <= 1)
                return;

            visited.set(start, true);

            auto node_seq = masked_graph.get_node_sequence(start);

            if (stop_path(start, start, node_seq)) {
                callback(start, start, std::move(node_seq));
                return;
            }

            std::stack<std::tuple<NodeIndex, NodeIndex, std::string>> paths;
            paths.emplace(start, start, std::move(node_seq));
            while (paths.size()) {
                auto path = std::move(paths.top());
                paths.pop();

                if (!visited[std::get<1>(path)]
                        && masked_graph.unmasked_outdegree(std::get<1>(path)) > 1) {
                    visited.set(std::get<1>(path), true);
                    paths.emplace(std::get<1>(path),
                                  std::get<1>(path),
                                  std::string(std::get<2>(path).end() - masked_graph.get_k(),
                                              std::get<2>(path).end()));
                }

                while (!stop_path(std::get<0>(path), std::get<1>(path), std::get<2>(path))) {
                    size_t cur_size = std::get<2>(path).size();

                    masked_graph.call_outgoing_kmers(
                        std::get<1>(path),
                        [&](const auto &index, char out_char) {
                            if (std::get<2>(path).size() == cur_size) {
                                std::get<2>(path).push_back(out_char);
                                std::get<1>(path) = index;
                            } else {
                                paths.emplace(std::get<0>(path),
                                              index,
                                              std::get<2>(path).substr(0, std::get<2>(path).length() - 1)
                                                  + out_char);
                            }
                        }
                    );

                    if (std::get<2>(path).size() == cur_size)
                        break;
                }

                callback(std::get<0>(path),
                         std::get<1>(path),
                         std::move(std::get<2>(path)));
                if (terminate())
                    break;
            }
        },
        terminate
    );
}


void call_breakpoints(const MaskedDeBruijnGraph &masked_graph,
                      const AnnotatedDBG &anno_graph,
                      AnnotatedDBGIndexRefVarLabelsCallback callback,
                      ThreadPool *thread_pool,
                      std::function<bool()> terminate) {
    assert(&masked_graph.get_graph()
        == dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph()));

    thread_pool = nullptr;

    call_paths_from_branch(masked_graph,
        [&](auto first, auto, auto&& sequence) {
            MaskedDeBruijnGraph background(
                masked_graph.get_graph_ptr(),
                [&](const auto &i) {
                    return i != DeBruijnGraph::npos
                        && (i == first || !masked_graph.in_graph(i));
                }
            );

            background.call_outgoing_kmers(
                first,
                [&](const auto &next_index, char c) {
                    callback(first,
                             std::move(sequence), std::string(1, c),
                             anno_graph.get_labels(next_index));
                }
            );
        },
        [&](auto, auto, const auto&) { return true; },
        terminate
    );

    if (thread_pool)
        thread_pool->join();
}

template <class Index, class VLabels>
using IndexSeqLabelCallback = std::function<void(const Index&,
                                                 const std::string&,
                                                 VLabels&&)>;

typedef IndexSeqLabelCallback<NodeIndex, std::vector<AnnotatedDBG::Annotator::Label>>
    AnnotatedDBGIndexSeqLabelsCallback;

void call_bubbles_from_path(const MaskedDeBruijnGraph &foreground,
                            const MaskedDeBruijnGraph &background,
                            const AnnotatedDBG &anno_graph,
                            NodeIndex first,
                            std::string seq,
                            AnnotatedDBGIndexSeqLabelsCallback callback,
                            std::function<bool()> &terminate) {
    if (foreground.get_graph_ptr() != background.get_graph_ptr())
        throw std::runtime_error("ERROR: bubble calling in matching graphs not implemented");

    assert(seq.length() == foreground.get_k() * 2 + 1);

    char ref = seq[foreground.get_k()];
    background.call_outgoing_kmers(
        first,
        [&](auto, char c) {
            if (terminate())
                return;

            if (c == ref)
                return;

            // For each outgoing edge, traverse path and intersect labels
            seq[foreground.get_k()] = c;

            bool in_foreground = true;
            bool in_background = true;
            anno_graph.get_graph().map_to_nodes_sequentially(
                seq.begin() + 1, seq.end(),
                [&](const auto &i) {
                    in_foreground &= foreground.in_graph(i);
                    in_background &= background.in_graph(i);
                },
                [&]() { return !in_background; }
            );

            if (in_foreground || !in_background)
                return;

            auto count_labels = anno_graph.get_top_labels(
                seq,
                anno_graph.get_annotation().num_labels(),
                1.0
            );

            if (count_labels.size()) {
                std::vector<std::string> labels;
                labels.reserve(count_labels.size());

                std::transform(std::make_move_iterator(count_labels.begin()),
                               std::make_move_iterator(count_labels.end()),
                               std::back_inserter(labels),
                               [&](auto&& pair) { return std::move(pair.first); });

                callback(first, seq, std::move(labels));
            }
        }
    );
}

void call_bubbles(const MaskedDeBruijnGraph &masked_graph,
                  const AnnotatedDBG &anno_graph,
                  AnnotatedDBGIndexRefVarLabelsCallback callback,
                  ThreadPool *thread_pool,
                  std::function<bool()> terminate) {
    assert(&masked_graph.get_graph()
        == dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph()));

    size_t path_length = masked_graph.get_k() * 2 + 1;

    call_paths_from_branch(masked_graph,
        [&](auto first, auto last, auto&& sequence) {
            if (sequence.size() != path_length
                    || masked_graph.unmasked_indegree(last) <= 1)
                return;

            auto process_path = [&, first, last, sequence]() {
                call_bubbles_from_path(
                    masked_graph,
                    MaskedDeBruijnGraph(masked_graph.get_graph_ptr()),
                    anno_graph,
                    first,
                    sequence,
                    [&](const auto &index, const auto &variant, auto&& labels) {
                        callback(index, sequence, variant, std::move(labels));
                    },
                    terminate
                );
            };

            if (thread_pool) {
                thread_pool->enqueue(process_path);
            } else {
                process_path();
            }
        },
        [&](auto, auto, auto&& sequence) { return sequence.size() >= path_length; },
        terminate
    );

    if (thread_pool)
        thread_pool->join();
}


} // namespace annotated_graph_algorithm
