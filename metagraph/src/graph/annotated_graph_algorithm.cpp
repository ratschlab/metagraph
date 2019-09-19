#include "annotated_graph_algorithm.hpp"

#include "annotate_column_compressed.hpp"
#include "utils.hpp"
#include "int_vector.hpp"
#include "masked_graph.hpp"

namespace annotated_graph_algorithm {

typedef DeBruijnGraph::node_index node_index;
typedef AnnotatedDBG::Annotator::Label Label;


sdsl::int_vector<>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<Label> &labels_in,
                  const std::vector<Label> &labels_out) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    size_t width = utils::code_length(std::max(labels_in.size(), labels_out.size()));
    sdsl::int_vector<> counts(anno_graph.get_graph().num_nodes() + 1, 0, width << 1);

    for (const auto &label_in : labels_in) {
        anno_graph.call_annotated_nodes(label_in,
                                        [&](const auto &i) { counts[i]++; });
    }

    // correct the width of counts, making it single-width
    counts.width(width);

    for (const auto &label_out : labels_out) {
        anno_graph.call_annotated_nodes(label_out,
                                        [&](const auto &i) { counts[(i << 1) + 1]++; });
    }

    // set the width to be double again
    counts.width(width << 1);

    return counts;
}

std::function<uint64_t(SequenceGraph::node_index i)>
build_label_counter(const AnnotatedDBG &anno_graph,
                    const std::vector<Label> &labels_to_check) {
    return [&anno_graph, labels_to_check](auto i) {
        return i == SequenceGraph::npos ? 0 : std::count_if(
            labels_to_check.begin(),
            labels_to_check.end(),
            [&](const auto &label) { return anno_graph.has_label(i, label); }
        );
    };
}

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<Label> &labels_in,
                    const std::vector<Label> &labels_out,
                    std::function<bool(LabelCountCallback /* get_num_labels_in */,
                                       LabelCountCallback /* get_num_labels_out */)> is_node_in_mask,
                    double lazy_evaluation_label_frequency_cutoff) {
    if (!anno_graph.get_graph().num_nodes())
        return std::make_unique<bitmap_lazy>([](auto) { return false; },
                                             anno_graph.get_graph().num_nodes() + 1,
                                             0);

    const auto *columns = dynamic_cast<const annotate::ColumnCompressed<>*>(
        &anno_graph.get_annotation()
    );

    if (columns) {
        size_t frequent_column_label_min_count = (anno_graph.get_graph().num_nodes() + 1)
                                                * lazy_evaluation_label_frequency_cutoff;

        // Partition labels into frequent and infrequent sets
        std::vector<Label> labels_in_infrequent,
                           labels_in_frequent,
                           labels_out_infrequent,
                           labels_out_frequent;

        for (const auto &label_in : labels_in) {
            if (columns->get_column(label_in).num_set_bits()
                    >= frequent_column_label_min_count) {
                labels_in_frequent.push_back(label_in);
            } else {
                labels_in_infrequent.push_back(label_in);
            }
        }

        for (const auto &label_out : labels_out) {
            if (columns->get_column(label_out).num_set_bits()
                    >= frequent_column_label_min_count) {
                labels_out_frequent.push_back(label_out);
            } else {
                labels_out_infrequent.push_back(label_out);
            }
        }

        // If at least one infrequent label exists, construct a count vector to
        // reduce calls to the annotator
        if (labels_in_infrequent.size() || labels_out_infrequent.size()) {
            auto counts = fill_count_vector(anno_graph,
                                            labels_in_infrequent,
                                            labels_out_infrequent);

            // the width of counts is double, since it's both in and out counts interleaved
            auto width = counts.width() >> 1;
            size_t int_mask = (size_t(1) << width) - 1;

            // Flatten count vector to bitmap if all labels were infrequent
            if (labels_in_frequent.empty() && labels_out_frequent.empty()) {
                auto mask_vector = std::make_unique<bitmap_vector>(
                    anno_graph.get_graph().num_nodes() + 1, false
                );

                call_nonzeros(counts,
                    [&](auto i, auto count) {
                        if (is_node_in_mask([&]() { return count & int_mask; },
                                            [&]() { return count >> width; }))
                            mask_vector->set(i, true);
                    }
                );

                return mask_vector;
            }

            // If at least one of the labels was frequent, construct a lazy bitmap
            // which references both the Annotator and the count vector
            auto count_frequent_in_labels = build_label_counter(anno_graph,
                                                                labels_in_frequent);
            auto count_frequent_out_labels = build_label_counter(anno_graph,
                                                                 labels_out_frequent);
            return std::make_unique<bitmap_lazy>([=](auto i) {
                auto count = counts[i];
                return is_node_in_mask(
                    [&]() { return (count & int_mask) + count_frequent_in_labels(i); },
                    [&]() { return (count >> width) + count_frequent_out_labels(i); }
                );
            }, counts.size());
        }
    }

    // If all of the labels were frequent, or if the Annotator is not ColumnCompressed,
    // construct a lazy bitmap
    auto count_frequent_in_labels = build_label_counter(anno_graph, labels_in);
    auto count_frequent_out_labels = build_label_counter(anno_graph, labels_out);
    return std::make_unique<bitmap_lazy>([=](auto i) {
        return is_node_in_mask([&]() { return count_frequent_in_labels(i); },
                               [&]() { return count_frequent_out_labels(i); });
    }, anno_graph.get_graph().num_nodes() + 1);
}

void
call_paths_from_branch(const DeBruijnGraph &graph,
                       const DeBruijnGraph &full_graph,
                       const std::function<void(node_index, node_index, std::string&&)> &callback,
                       const std::function<bool(node_index, node_index, const std::string&)> &stop_path,
                       const std::function<bool()> &terminate = []() { return false; }) {
    assert(&graph != &full_graph);

    sdsl::bit_vector visited(graph.num_nodes() + 1, false);

    bool stop = false;
    graph.call_nodes(
        [&](auto start) {
            if (visited[start] || full_graph.outdegree(start) <= 1)
                return;

            visited[start] = true;

            auto node_seq = graph.get_node_sequence(start);

            if (stop_path(start, start, node_seq)) {
                if (!(stop = terminate()))
                    callback(start, start, std::move(node_seq));

                return;
            }

            std::stack<std::tuple<node_index, node_index, std::string>> paths;
            paths.emplace(start, start, std::move(node_seq));
            while (paths.size()) {
                auto path = std::move(paths.top());
                paths.pop();

                if (!visited[std::get<1>(path)]
                        && full_graph.outdegree(std::get<1>(path)) > 1) {
                    visited[std::get<1>(path)] = true;
                    paths.emplace(std::get<1>(path),
                                  std::get<1>(path),
                                  std::string(std::get<2>(path).end() - graph.get_k(),
                                              std::get<2>(path).end()));
                }

                while (!std::apply(stop_path, path)) {
                    size_t cur_size = std::get<2>(path).size();

                    graph.call_outgoing_kmers(
                        std::get<1>(path),
                        [&](const auto &index, char out_char) {
                            if (std::get<2>(path).size() == cur_size) {
                                std::get<2>(path).push_back(out_char);
                                std::get<1>(path) = index;
                            } else {
                                paths.emplace(
                                    std::get<0>(path),
                                    index,
                                    std::get<2>(path).substr(
                                        0, std::get<2>(path).length() - 1
                                    ) + out_char
                                );
                            }
                        }
                    );

                    if (std::get<2>(path).size() == cur_size)
                        break;
                }

                if ((stop = terminate()))
                    break;

                std::apply(callback, std::move(path));
            }
        },
        [&]() { return stop; }
    );
}

void call_breakpoints(const DeBruijnGraph &graph,
                      const AnnotatedDBG &anno_graph,
                      const VariantLabelCallback &callback,
                      ThreadPool *thread_pool,
                      const std::function<bool()> &terminate) {
    const auto dbg_succ = dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    assert(dbg_succ.get());

    if (&graph == dbg_succ.get())
        return;

    thread_pool = nullptr;

    const DBGAlignerConfig variant_config(
        DBGAlignerConfig::unit_scoring_matrix(1, graph.alphabet()),
        -1, -1
    );

    bool stop = false;
    call_paths_from_branch(
        graph,
        *dbg_succ,
        [&](auto first, auto, auto&& sequence) {
            assert(sequence.size() == graph.get_k());

            std::vector<std::pair<DeBruijnGraph::node_index, char>> outgoing;
            graph.call_outgoing_kmers(
                first,
                [&](const auto &next_index, char c) {
                    outgoing.emplace_back(next_index, c);
                }
            );

            MaskedDeBruijnGraph background_masked(
                dbg_succ,
                [&](const auto &i) { return i == first || !graph.in_graph(i); }
            );

            // if outgoing is empty, we don't have to check for it to be excluded
            // from the mask
            const DeBruijnGraph& background = outgoing.size()
                ? background_masked : *dbg_succ;

            background.call_outgoing_kmers(
                first,
                [&](const auto &next_index, char c) {
                    if ((stop = terminate()))
                        return;

                    auto var = sequence + c;
                    auto ref = const_cast<const std::string&>(sequence);

                    if (outgoing.size())
                        ref += outgoing.front().second;

                    auto score = variant_config.score_sequences(ref.begin(),
                                                                ref.end(),
                                                                var.begin());

                    if (ref.size() != var.size()) {
                        auto difference = ref.size() > var.size()
                            ? ref.size() - var.size() - 1
                            : var.size() - ref.size() - 1;

                        score += variant_config.gap_opening_penalty
                            + difference * variant_config.gap_extension_penalty;
                    }

                    DBGAlignment breakpoint(ref.c_str(),
                                            &*ref.end(),
                                            { first, next_index },
                                            std::move(var),
                                            score);

                    assert(breakpoint.is_valid(background));

                    callback(std::move(breakpoint),
                             ref,
                             anno_graph.get_labels(next_index));
                }
            );
        },
        [&](const auto&...) { return true; },
        [&]() { return stop; }
    );

    if (thread_pool)
        thread_pool->join();
}

void call_bubbles_from_path(const DeBruijnGraph &foreground,
                            const DeBruijnGraph &background,
                            const AnnotatedDBG &anno_graph,
                            node_index first,
                            const std::string &ref,
                            const VariantLabelCallback &callback,
                            const std::function<bool()> &terminate,
                            const DBGAlignerConfig &variant_config) {
    assert(&foreground != &background);
    assert(ref.length() == foreground.get_k() * 2 + 1);

    const char ref_char = ref[foreground.get_k()];
    bool stop = false;
    background.call_outgoing_kmers(
        first,
        [&](auto, char c) {
            if ((stop = terminate()))
                return;

            if (c == ref_char)
                return;

            // For each outgoing edge, traverse path and intersect labels
            auto var = ref;
            var[foreground.get_k()] = c;

            bool in_foreground = true;
            bool in_background = true;
            std::vector<DeBruijnGraph::node_index> nodes { first };
            anno_graph.get_graph().map_to_nodes_sequentially(
                var.begin() + 1,
                var.end(),
                [&](const auto &i) {
                    nodes.emplace_back(i);
                    in_foreground &= foreground.in_graph(i);
                    in_background &= background.in_graph(i);
                },
                [&]() { return !in_background; }
            );

            if (in_foreground || !in_background)
                return;

            auto labels = anno_graph.get_labels(var, 1.0);
            if (labels.size()) {
                auto score = variant_config.score_sequences(ref.begin(),
                                                            ref.end(),
                                                            var.begin());

                DBGAlignment bubble(ref.c_str(),
                                    &*ref.end(),
                                    std::move(nodes),
                                    std::move(var),
                                    score);
                assert(bubble.is_valid(background));

                callback(std::move(bubble), ref, std::move(labels));
            }
        }
    );
}

void call_bubbles(const DeBruijnGraph &graph,
                  const AnnotatedDBG &anno_graph,
                  const VariantLabelCallback &callback,
                  ThreadPool *thread_pool,
                  const std::function<bool()> &terminate) {
    const auto dbg_succ = dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    assert(dbg_succ.get());

    if (&graph == dbg_succ.get())
        return;

    const DBGAlignerConfig variant_config(
        DBGAlignerConfig::unit_scoring_matrix(1, graph.alphabet()),
        -1, -1
    );

    size_t path_length = graph.get_k() * 2 + 1;

    call_paths_from_branch(
        graph,
        *dbg_succ,
        [&](auto first, auto last, auto&& sequence) {
            if (sequence.size() != path_length || dbg_succ->indegree(last) <= 1)
                return;

            auto process_path = [&, first, sequence]() {
                call_bubbles_from_path(
                    graph,
                    *dbg_succ,
                    anno_graph,
                    first,
                    sequence,
                    callback,
                    terminate,
                    variant_config
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


Extender<DeBruijnGraph::node_index>
build_masked_graph_extender(const AnnotatedDBG &anno_graph,
                            double seed_label_discovery_fraction,
                            Extender<DeBruijnGraph::node_index>&& extender) {
    assert(dynamic_cast<const DeBruijnGraph*>(anno_graph.get_graph_ptr().get()));

    return [=,&anno_graph,
            extender{std::move(extender)}](const DeBruijnGraph &graph,
                                           const auto &path,
                                           std::vector<auto>* next_paths,
                                           const char* sequence_end,
                                           const DBGAlignerConfig &config,
                                           bool orientation,
                                           auto min_path_score) {
        auto labels = anno_graph.get_labels(path.get_sequence(),
                                            seed_label_discovery_fraction);

        if (labels.empty())
            return;

        extender(MaskedDeBruijnGraph(
                     std::dynamic_pointer_cast<const DeBruijnGraph>(
                         anno_graph.get_graph_ptr()
                     ),
                     [&, labels{std::move(labels)}](const auto &node) {
                         assert(node != DeBruijnGraph::npos);
                         return graph.in_graph(node)
                             && std::all_of(labels.begin(), labels.end(),
                                            [&](const auto &label) {
                                                return anno_graph.has_label(node, label);
                                            });
                     }
                 ),
                 path,
                 next_paths,
                 sequence_end,
                 config,
                 orientation,
                 min_path_score);
    };
}


Extender<DeBruijnGraph::node_index>
build_background_graph_extender(const DeBruijnGraph &background,
                                Extender<DeBruijnGraph::node_index>&& extender) {
    return [&,base_extender{std::move(extender)}](const DeBruijnGraph &foreground,
                                                  const auto &path,
                                                  std::vector<auto>* next_paths,
                                                  const char *sequence_end,
                                                  const DBGAlignerConfig &config,
                                                  bool orientation,
                                                  auto min_path_score) {
        if (path.empty())
            return;

        // TODO: fix this hax for generating a shared_ptr from a const reference
        base_extender(
            MaskedDeBruijnGraph({ std::shared_ptr<const DeBruijnGraph>{}, &background },
                                [&, seed_node = path.back()](const auto &node) {
                                    // check if node is not in the foreground, or
                                    // if it's a breakpoint
                                    return node == seed_node
                                        || !foreground.in_graph(node)
                                        || (!background.has_single_incoming(node)
                                               && foreground.has_single_incoming(node));
                                }),
            path,
            next_paths,
            sequence_end,
            config,
            orientation,
            min_path_score
        );
    };
}

} // namespace annotated_graph_algorithm
