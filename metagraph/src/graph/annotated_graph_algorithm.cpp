#include "annotated_graph_algorithm.hpp"

#include "annotate_column_compressed.hpp"
#include "utils.hpp"

#include <mutex>

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

void call_outgoing_nmers(const DeBruijnGraph &graph,
                         DeBruijnGraph::node_index node,
                         const std::function<void(const std::string&)> &callback,
                         size_t n) {
    if (!n)
        return;

    size_t path_length = n + graph.get_k() - 1;
    std::stack<std::pair<DeBruijnGraph::node_index, std::string>> paths;
    paths.emplace(node, graph.get_node_sequence(node));
    while (paths.size()) {
        auto path = std::move(paths.top());
        paths.pop();
        while (path.second.size() < path_length) {
            size_t cur_size = path.second.size();
            graph.call_outgoing_kmers(
                path.first,
                [&](const auto &index, char out_char) {
                    if (path.second.size() == cur_size) {
                        path.second.push_back(out_char);
                        path.first = index;
                    } else {
                        paths.emplace(index,
                                      path.second.substr(0, path.second.length() - 1)
                                          + out_char);
                    }
                }
            );
            if (path.second.size() == cur_size)
                break;
        }
        if (path.second.size() == path_length)
            callback(path.second);
    }
}

template <class Index, class VLabels>
using IndexSeqLabelCallback = std::function<void(const Index&,
                                                 const std::string&,
                                                 VLabels&&)>;

typedef IndexSeqLabelCallback<MaskedDeBruijnGraph::node_index,
                              std::vector<AnnotatedDBG::Annotator::Label>>
    AnnotatedDBGIndexSeqLabelsCallback;

void call_bubbles_from_path(const MaskedDeBruijnGraph &foreground,
                            const MaskedDeBruijnGraph &background,
                            const AnnotatedDBG &anno_graph,
                            MaskedDeBruijnGraph::node_index path_start,
                            std::string seq,
                            const AnnotatedDBGIndexSeqLabelsCallback &callback) {
    if (foreground.get_graph_ptr() != background.get_graph_ptr())
        throw std::runtime_error("ERROR: bubble calling in matching graphs not implemented");

    // Don't check if sequence is too short for bubble to exist
    if (seq.length() < foreground.get_k() * 2 + 1)
        return;

    auto candidates = anno_graph.get_labels(path_start);
    if (candidates.empty())
        return;

    std::sort(candidates.begin(), candidates.end());

    char ref = seq[foreground.get_k()];
    background.call_outgoing_kmers(
        path_start,
        [&](auto, char c) {
            if (c == ref)
                return;

            // For each outgoing edge, traverse path and intersect labels
            seq[foreground.get_k()] = c;

            bool in_foreground = true;
            bool in_background = true;
            anno_graph.get_graph().map_to_nodes(
                seq.substr(1),
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

                callback(path_start, seq, std::move(labels));
            }
        }
    );
}

void call_bubbles(const MaskedDeBruijnGraph &masked_graph,
                  const AnnotatedDBG &anno_graph,
                  const AnnotatedDBGIndexRefVarLabelsCallback &callback,
                  ThreadPool *thread_pool,
                  const std::function<bool()> &terminate) {
    assert(&masked_graph.get_graph()
        == dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph()));

    std::mutex visited_mutex;
    bit_vector_stat visited(masked_graph.num_nodes() + 1, false);
    uint64_t path_length = masked_graph.get_k() + 2;

    masked_graph.call_nodes(
        [&](const auto &index) {
            if (visited[index])
                return;

            if (masked_graph.unmasked_outdegree(index) <= 1)
                return;

            auto process_path =
                [&](auto index, auto sequence) {
                    if (thread_pool) {
                        std::lock_guard<std::mutex> lock(visited_mutex);
                        visited.set(index, true);
                    } else {
                        visited.set(index, true);
                    }
                    call_bubbles_from_path(
                        masked_graph,
                        MaskedDeBruijnGraph(masked_graph.get_graph_ptr()),
                        anno_graph,
                        index,
                        sequence,
                        [&](const auto &index, const auto &variant, auto&& labels) {
                            callback(index, sequence, variant, std::move(labels));
                        }
                    );
                };

            if (thread_pool) {
                call_outgoing_nmers(
                    masked_graph,
                    index,
                    [&](const auto &sequence) {
                        thread_pool->enqueue(process_path, index, sequence);
                    },
                    path_length
                );
            } else {
                call_outgoing_nmers(masked_graph,
                                    index,
                                    [&](const auto &sequence) {
                                        process_path(index, sequence);
                                    },
                                    path_length);
            }
        }, terminate
    );

    if (thread_pool)
        thread_pool->join();
}

void call_breakpoints(const MaskedDeBruijnGraph &masked_graph,
                      const AnnotatedDBG &anno_graph,
                      const AnnotatedDBGIndexRefVarLabelsCallback &callback,
                      ThreadPool *thread_pool,
                      const std::function<bool()> &terminate) {
    assert(&masked_graph.get_graph()
        == dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph()));

    thread_pool = nullptr;

    masked_graph.call_nodes(
        [&](const auto &index) {
            if (masked_graph.unmasked_outdegree(index) <= 1)
                return;

            auto seq = masked_graph.get_node_sequence(index);

            MaskedDeBruijnGraph background(
                masked_graph.get_graph_ptr(),
                [&](const auto &i) {
                    return i != DeBruijnGraph::npos
                        && (i == index || !masked_graph.in_graph(i));
                }
            );

            background.call_outgoing_kmers(
                index,
                [&](const auto &next_index, char c) {
                    callback(index,
                             seq, std::string(1, c),
                             anno_graph.get_labels(next_index));
                }
            );
        }, terminate
    );

    if (thread_pool)
        thread_pool->join();
}


} // namespace annotated_graph_algorithm
