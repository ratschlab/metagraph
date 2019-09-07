#include "annotated_graph_algorithm.hpp"

#include "annotate_column_compressed.hpp"
#include "utils.hpp"
#include "progress_bar.hpp"
#include "int_vector.hpp"
#include "masked_graph.hpp"

namespace annotated_graph_algorithm {

typedef DeBruijnGraph::node_index node_index;

uint64_t count_node_labels(const AnnotatedDBG &anno_graph,
                           const node_index &index,
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

void fill_count_vector(const AnnotatedDBG &anno_graph,
                       const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                       const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                       sdsl::int_vector<>* counts,
                       std::vector<AnnotatedDBG::Annotator::Label>* mask_in_dense,
                       std::vector<AnnotatedDBG::Annotator::Label>* mask_out_dense,
                       const std::function<bool(const std::string&)> &pick_label
                           = [](const std::string&) { return true; }) {
    // in the beginning, counts is the correct size, but double width
    assert(!(counts->width() & 1));
    uint64_t width = counts->width() >> 1;

    for (const auto &label_in : mask_in) {
        if (pick_label(label_in)) {
            anno_graph.call_annotated_nodes(
                label_in,
                [&](const auto &i) { (*counts)[i]++; }
            );
        } else {
            mask_in_dense->push_back(label_in);
        }
    }

    // correct the width of counts, making it double-length
    counts->width(width);

    for (const auto &label_out : mask_out) {
        if (pick_label(label_out)) {
            anno_graph.call_annotated_nodes(
                label_out,
                [&](const auto &i) { (*counts)[(i << 1) + 1]++; }
            );
        } else {
            mask_out_dense->push_back(label_out);
        }
    }

    // set the width to be double again
    counts->width(width * 2);
}

std::unique_ptr<bitmap>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_in,
                    const std::vector<AnnotatedDBG::Annotator::Label> &mask_out,
                    const std::function<bool(const UInt64Callback&,
                                             const UInt64Callback&)> &keep_node) {
    if (!anno_graph.get_graph().num_nodes())
        return {};

    std::unique_ptr<bitmap> mask(
        new bitmap_vector(anno_graph.get_graph().num_nodes() + 1, false)
    );

    if (mask_in.empty())
        return mask;

    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    size_t width = utils::code_length(std::max(mask_in.size(), mask_out.size()));
    size_t int_mask = (size_t(1) << width) - 1;

    sdsl::int_vector<> counts(mask->size(), 0, width * 2);
    std::vector<AnnotatedDBG::Annotator::Label> mask_in_dense, mask_out_dense;

    const auto *columns = dynamic_cast<const annotate::ColumnCompressed<>*>(
        &anno_graph.get_annotation()
    );

    if (columns) {
        fill_count_vector(
            anno_graph,
            mask_in, mask_out,
            &counts,
            &mask_in_dense, &mask_out_dense,
            [&, density_cutoff_count = mask->size() * density_cutoff](const std::string &label) {
                return columns->get_column(label).num_set_bits() < density_cutoff_count;
            }
        );
    } else {
        // TODO: make this more efficient for row-major annotations
        fill_count_vector(anno_graph,
                          mask_in, mask_out,
                          &counts,
                          &mask_in_dense, &mask_out_dense);
    }

    if (utils::get_verbose())
        std::cerr << "Generating mask" << std::endl;

    if (mask_in_dense.empty() && mask_out_dense.empty()) {
        call_nonzeros(counts,
            [&](auto i, auto count) {
                if (keep_node([&]() { return count & int_mask; },
                              [&]() { return count >> width; }))
                    mask->set(i, true);
            }
        );
    } else {
        ProgressBar progress_bar(counts.size(),
                                 "Generating mask",
                                 std::cerr,
                                 !utils::get_verbose());
        progress_bar.SetFrequencyUpdate(100000);

        for (size_t i = 1; i < counts.size(); ++i) {
            auto count = counts[i];
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

            ++progress_bar;
        }
    }

    return mask;
}

void
call_paths_from_branch(const DeBruijnGraph &graph,
                       const DeBruijnGraph &full_graph,
                       const std::function<void(node_index, node_index, std::string&&)> &callback,
                       const std::function<bool(node_index, node_index, const std::string&)> &stop_path,
                       const std::function<bool()> &terminate = []() { return false; }) {
    assert(&graph != &full_graph);

    bit_vector_stat visited(graph.num_nodes() + 1, false);

    bool stop = false;
    graph.call_nodes(
        [&](auto start) {
            if (visited[start] || full_graph.outdegree(start) <= 1)
                return;

            visited.set(start, true);

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
                    visited.set(std::get<1>(path), true);
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

                    auto breakpoint = MaskedAlignment(ref.c_str(),
                                                      &*ref.end(),
                                                      { first, next_index },
                                                      std::move(var),
                                                      score);

                    assert(breakpoint.is_valid(background));

                    callback(std::move(breakpoint),
                             std::move(ref),
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

                auto ref_copy = ref;
                auto bubble = MaskedAlignment(ref_copy.c_str(),
                                              &*ref_copy.end(),
                                              std::move(nodes),
                                              std::move(var),
                                              score);
                assert(bubble.is_valid(background));

                callback(std::move(bubble), std::move(ref_copy), std::move(labels));
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


} // namespace annotated_graph_algorithm
