#include "annotated_graph_algorithm.hpp"

#include <unordered_set>

#include "common/vectors/vector_algorithm.hpp"
#include "common/vector_map.hpp"
#include "kmer/alphabets.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/representation/masked_graph.hpp"


namespace mtg {
namespace graph {

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::row_index row_index;
typedef AnnotatedDBG::Annotator::Label Label;
typedef Alignment<DeBruijnGraph::node_index> DBGAlignment;


std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig(const DeBruijnGraph &graph,
                     const KeepUnitigPath &keep_unitig,
                     size_t num_threads) {
    sdsl::bit_vector unitig_mask(graph.max_index() + 1, false);

    const bool atomic = num_threads > 1;

    graph.call_unitigs([&](const std::string &unitig, const auto &path) {
        if (keep_unitig(unitig, path)) {
            for (DeBruijnGraph::node_index node : path) {
                set_bit(unitig_mask.data(), node, atomic);
            }
        }
    }, num_threads);

    return std::make_unique<bitmap_vector>(std::move(unitig_mask));
}

std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig_labels(const AnnotatedDBG &anno_graph,
                            const std::vector<Label> &labels_in,
                            const std::vector<Label> &labels_out,
                            size_t num_threads,
                            double label_mask_in_fraction,
                            double label_mask_out_fraction,
                            double label_other_fraction) {
    const auto &dbg = anno_graph.get_graph();
    const auto &annotation = anno_graph.get_annotation();
    const auto &label_encoder = annotation.get_label_encoder();

    assert(labels_in.size() + labels_out.size() <= annotation.num_labels());

    const double label_in_factor = label_mask_in_fraction * labels_in.size();
    const double label_out_factor = label_mask_out_fraction * labels_out.size();

    std::unordered_set<uint64_t> labels_in_enc;
    for (const auto &label_in : labels_in) {
        labels_in_enc.emplace(label_encoder.encode(label_in));
    }

    std::unordered_set<uint64_t> labels_out_enc;
    for (const auto &label_out : labels_out) {
        labels_out_enc.emplace(label_encoder.encode(label_out));
    }

    return mask_nodes_by_unitig(dbg, [&](const auto &, const auto &path) {
        VectorMap<row_index, size_t> index_counts;
        for (node_index i : path) {
            index_counts[anno_graph.graph_to_anno_index(i)]++;
        }

        size_t other_count = 0;
        size_t in_count = 0;
        size_t out_count = 0;
        const size_t out_count_cutoff = label_out_factor * path.size();

        for (const auto &pair : annotation.get_matrix().sum_rows(index_counts.values_container())) {
            if (labels_in_enc.find(pair.first) != labels_in_enc.end()) {
                in_count += pair.second;

            } else if (labels_out_enc.find(pair.first) != labels_out_enc.end()) {
                // early cutoff
                if ((out_count += pair.second) > out_count_cutoff)
                    return false;

            } else {
                other_count += pair.second;
            }
        }

        return (in_count >= label_in_factor * path.size())
            && (other_count <= label_other_fraction * (in_count + out_count + other_count));
    }, num_threads);
}

sdsl::int_vector<>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const std::vector<Label> &labels_in,
                                  const std::vector<Label> &labels_out,
                                  size_t num_threads) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    auto counts = aligned_int_vector(anno_graph.get_graph().max_index() + 1, 0, width << 1);

    std::mutex backup_mutex;
    constexpr std::memory_order memorder = std::memory_order_relaxed;

    std::atomic_thread_fence(std::memory_order_release);
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < labels_in.size(); ++j) {
        const std::string &label_in = labels_in[j];
        anno_graph.call_annotated_nodes(label_in, [&](DeBruijnGraph::node_index i) {
            atomic_fetch_and_add(counts, i, 1, backup_mutex, memorder);
        });
    }
    std::atomic_thread_fence(std::memory_order_acquire);

    // correct the width of counts, making it single-width
    counts.width(width);

    std::atomic_thread_fence(std::memory_order_release);
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < labels_out.size(); ++j) {
        const std::string &label_out = labels_out[j];
        anno_graph.call_annotated_nodes(label_out, [&](DeBruijnGraph::node_index i) {
            atomic_fetch_and_add(counts, (i << 1) + 1, 1, backup_mutex, memorder);
        });
    }
    std::atomic_thread_fence(std::memory_order_acquire);

    // set the width to be double again
    counts.width(width << 1);

    return counts;
}

std::function<uint64_t(SequenceGraph::node_index)>
build_label_counter(const AnnotatedDBG &anno_graph,
                    const std::vector<Label> &labels_to_check) {
    return [&anno_graph, labels_to_check](SequenceGraph::node_index i) {
        assert(i != SequenceGraph::npos);
        return std::count_if(labels_to_check.begin(), labels_to_check.end(),
            [&](const auto &label) { return anno_graph.has_label(i, label); }
        );
    };
}

std::unique_ptr<bitmap>
mask_nodes_by_node_label(const AnnotatedDBG &anno_graph,
                         const std::vector<Label> &labels_in,
                         const std::vector<Label> &labels_out,
                         const std::function<bool(DeBruijnGraph::node_index,
                                                  const LabelCountCallback & /* get_num_labels_in */,
                                                  const LabelCountCallback & /* get_num_labels_out */)> &is_node_in_mask,
                         size_t num_threads,
                         double min_frequency_for_frequent_label) {
    if (!anno_graph.get_graph().num_nodes())
        return std::make_unique<bitmap_lazy>([](uint64_t) { return false; },
                                             anno_graph.get_graph().max_index() + 1,
                                             0);

    const auto *columns = dynamic_cast<const annot::ColumnCompressed<>*>(
        &anno_graph.get_annotation()
    );

    //TODO: why do we have num_nodes + 1 here and not just num_nodes?
    if (columns) {
        size_t frequent_column_label_min_count
            = (anno_graph.get_graph().num_nodes() + 1)
                * min_frequency_for_frequent_label;

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
            sdsl::int_vector<> counts = construct_diff_label_count_vector(
                anno_graph, labels_in_infrequent, labels_out_infrequent, num_threads
            );

            // the width of counts is double, since it's both in and out counts interleaved
            uint32_t width = counts.width() >> 1;
            uint64_t int_mask = (uint64_t(1) << width) - 1;

            // Flatten count vector to bitmap if all labels were infrequent
            if (labels_in_frequent.empty() && labels_out_frequent.empty()) {
                sdsl::bit_vector mask(anno_graph.get_graph().max_index() + 1, false);

                call_nonzeros(counts, [&](DeBruijnGraph::node_index i, size_t count) {
                    if (i != DeBruijnGraph::npos
                            && is_node_in_mask(i, [&]() { return count & int_mask; },
                                                  [&]() { return count >> width; }))
                        mask[i] = true;
                });

                return std::make_unique<bitmap_vector>(std::move(mask));
            }

            // If at least one of the labels was frequent, construct a lazy bitmap
            // which references both the Annotator and the count vector
            auto count_frequent_in_labels = build_label_counter(anno_graph,
                                                                labels_in_frequent);
            auto count_frequent_out_labels = build_label_counter(anno_graph,
                                                                 labels_out_frequent);
            return std::make_unique<bitmap_lazy>([=](uint64_t i) {
                auto count = counts[i];
                return i != DeBruijnGraph::npos
                    && is_node_in_mask(i,
                            [&]() { return (count & int_mask) + count_frequent_in_labels(i); },
                            [&]() { return (count >> width) + count_frequent_out_labels(i); });
            }, counts.size());
        }
    }

    // If all of the labels were frequent, or if the Annotator is not ColumnCompressed,
    // construct a lazy bitmap
    auto count_frequent_in_labels = build_label_counter(anno_graph, labels_in);
    auto count_frequent_out_labels = build_label_counter(anno_graph, labels_out);

    return std::make_unique<bitmap_lazy>([=](DeBruijnGraph::node_index i) {
        return i != DeBruijnGraph::npos
            && is_node_in_mask(i, [&]() { return count_frequent_in_labels(i); },
                                  [&]() { return count_frequent_out_labels(i); });
    }, anno_graph.get_graph().max_index() + 1);
}

} // namespace graph
} // namespace mtg
