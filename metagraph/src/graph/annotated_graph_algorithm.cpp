#include "annotated_graph_algorithm.hpp"

#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vector_map.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"

using mtg::common::logger;

namespace mtg {
namespace graph {

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::row_index row_index;


std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig(const DeBruijnGraph &graph,
                     const KeepUnitigPath &keep_unitig) {
    sdsl::bit_vector unitig_mask(graph.max_index() + 1, false);

    std::atomic<size_t> kept_unitigs(0);
    std::atomic<size_t> total_unitigs(0);
    std::atomic<size_t> num_kept_nodes(0);
    bool parallel = get_num_threads() > 1;
    auto memorder = std::memory_order_relaxed;

    std::atomic_thread_fence(std::memory_order_release);
    graph.call_unitigs([&](const std::string &unitig, const auto &path) {
        if (keep_unitig(unitig, path)) {
            kept_unitigs.fetch_add(1, memorder);
            for (node_index node : path) {
                num_kept_nodes.fetch_add(
                    !fetch_and_set_bit(unitig_mask.data(), node, parallel),
                    memorder
                );
            }
        }
        total_unitigs.fetch_add(1, memorder);
    }, get_num_threads());
    std::atomic_thread_fence(std::memory_order_acquire);

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (graph.get_k() - 1))
                      / kept_unitigs);

    return std::make_unique<bitmap_vector>(std::move(unitig_mask));
}

MaskedDeBruijnGraph
make_masked_graph_by_unitig_labels(const AnnotatedDBG &anno_graph,
                                   const std::vector<Label> &labels_in,
                                   const std::vector<Label> &labels_out,
                                   double label_mask_in_fraction,
                                   double label_mask_out_fraction,
                                   double label_other_fraction,
                                   bool add_complement) {
    const auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    logger->trace("Generating initial mask");
    // Construct initial masked graph from union of labels in labels_in
    sdsl::bit_vector union_mask(graph_ptr->max_index() + 1, false);

    bool parallel = get_num_threads() > 1;
    int memorder = __ATOMIC_RELAXED;

    __atomic_thread_fence(__ATOMIC_RELEASE);
    #pragma omp parallel num_threads(get_num_threads())
    #pragma omp single
    {
        #pragma omp taskloop default(shared)
        for (size_t i = 0; i < labels_in.size(); ++i) {
            const std::string &label = labels_in[i];
            logger->trace("Adding label: {}", label);
            anno_graph.call_annotated_nodes(label, [&](node_index node) {
                set_bit(union_mask.data(), node, parallel, memorder);
            });
        }
        #pragma omp taskloop default(shared)
        for (size_t i = 0; i < labels_out.size(); ++i) {
            const std::string &label = labels_out[i];
            logger->trace("Adding label: {}", label);
            anno_graph.call_annotated_nodes(label, [&](node_index node) {
                set_bit(union_mask.data(), node, parallel, memorder);
            });
        }
    }
    __atomic_thread_fence(__ATOMIC_ACQUIRE);

    if (graph_ptr->is_canonical_mode()) {
        __atomic_thread_fence(__ATOMIC_RELEASE);
        MaskedDeBruijnGraph forward_masked(
            graph_ptr,
            std::make_unique<bit_vector_stat>(sdsl::bit_vector(union_mask)),
            true // only_valid_nodes_in_mask
        );

        logger->trace("Starting with {} nodes", forward_masked.num_nodes());
        logger->trace("Adding reverse complements");

        forward_masked.call_unitigs([&](const std::string &unitig, const auto &) {
            auto rev = unitig;
            reverse_complement(rev.begin(), rev.end());
            graph_ptr->map_to_nodes_sequentially(rev, [&](node_index node) {
                assert(node);
                set_bit(union_mask.data(), node, parallel, memorder);
            });
        }, get_num_threads());
        __atomic_thread_fence(__ATOMIC_ACQUIRE);
    }

    auto masked_graph = std::make_unique<MaskedDeBruijnGraph>(
        graph_ptr,
        std::make_unique<bit_vector_stat>(std::move(union_mask)),
        true, // only_valid_nodes_in_mask
        graph_ptr->is_canonical_mode()
    );
    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");

    tsl::hopscotch_set<std::string> label_in_set(labels_in.begin(), labels_in.end());
    tsl::hopscotch_set<std::string> label_out_set(labels_out.begin(), labels_out.end());

    size_t num_labels = anno_graph.get_annotation().num_labels();
    double min_frac = std::min(std::min(label_mask_in_fraction, label_mask_out_fraction), label_other_fraction);

    auto check = [&](size_t path_size, const std::string &label, size_t count, size_t *in_label_counter) {
        if (label_in_set.count(label)) {
            if (count < label_mask_in_fraction * path_size) {
                return false;
            } else {
                ++(*in_label_counter);
            }
        } else if (label_out_set.count(label)) {
            if (count > label_mask_out_fraction * path_size)
                return false;

        } else if (count > label_other_fraction * path_size) {
            return false;
        }

        return true;
    };

    if (add_complement && !graph_ptr->is_canonical_mode()) {
        std::mutex add_mutex;
        BOSSConstructor constructor(graph_ptr->get_k() - 1, true);
        constructor.add_sequences([&](const CallString &callback) {
            masked_graph->call_sequences([&](const std::string &seq, const auto &) {
                std::lock_guard<std::mutex> lock(add_mutex);
                callback(seq);
            }, get_num_threads(), true);
        });

        auto inlabel_graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), true);

        return MaskedDeBruijnGraph(
            inlabel_graph,
            mask_nodes_by_unitig(*inlabel_graph, [&](const auto &unitig, const auto &path) {
                VectorMap<std::string, std::pair<sdsl::bit_vector, sdsl::bit_vector>> label_signatures;
                for (auto&& [label, signature] : anno_graph.get_top_label_signatures(unitig, num_labels)) {
                    label_signatures.emplace(std::move(label),
                                             std::make_pair(std::move(signature),
                                                            sdsl::bit_vector()));
                }

                std::string rev = unitig;
                reverse_complement(rev.begin(), rev.end());
                for (auto&& [label, signature] : anno_graph.get_top_label_signatures(rev, num_labels)) {
                    if (label_signatures.count(label)) {
                        label_signatures[label].second = std::move(signature);
                    } else {
                        label_signatures.emplace(std::move(label),
                                                 std::make_pair(sdsl::bit_vector(),
                                                                std::move(signature)));
                    }
                }

                size_t in_label_counter = 0;
                for (const auto &[label, sig_pair] : label_signatures) {
                    const auto &[sig, rev_sig] = sig_pair;
                    size_t count = 0;
                    if (sig.size() && rev_sig.size()) {
                        auto it = sig.begin();
                        auto jt = rev_sig.end();
                        for (--jt; it != sig.end(); ++it, --jt) {
                            count += (*it) || (*jt);
                        }
                    } else if (sig.size()) {
                        count = sdsl::util::cnt_one_bits(sig);
                    } else {
                        count = sdsl::util::cnt_one_bits(rev_sig);
                    }

                    if (!check(path.size(), label, count, &in_label_counter))
                        return false;
                }

                return in_label_counter == labels_in.size();
            }),
            true
        );
    }

    return MaskedDeBruijnGraph(
        graph_ptr,
        mask_nodes_by_unitig(*masked_graph, [&](const auto &unitig, const auto &path) {
            size_t in_label_counter = 0;
            for (const auto &[label, count] : anno_graph.get_top_labels(unitig, num_labels, min_frac)) {
                if (!check(path.size(), label, count, &in_label_counter))
                    return false;
            }

            return in_label_counter == labels_in.size();
        }),
        true
    );
}


sdsl::int_vector<> fill_count_vector(const AnnotatedDBG &anno_graph,
                                     const std::vector<Label> &labels_in,
                                     const std::vector<Label> &labels_out) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());
    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    sdsl::int_vector<> counts(anno_graph.get_graph().max_index() + 1, 0, width << 1);
    sdsl::bit_vector indicator(graph->is_canonical_mode() ? counts.size() : 0);

    size_t num_threads = get_num_threads();
    bool async = num_threads > 1;
    constexpr int memorder = __ATOMIC_RELAXED;

    // TODO: replace locked increment operations on int_vector<> with actual
    //       atomic operations when we figure out how to align int_vector<> storage
    std::mutex count_mutex;

    __atomic_thread_fence(__ATOMIC_RELEASE);
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (size_t i = 0; i < labels_in.size(); ++i) {
        const std::string &label_in = labels_in[i];
        if (graph->is_canonical_mode()) {
            anno_graph.call_annotated_nodes(label_in, [&](node_index i) {
                set_bit(indicator.data(), i, async, memorder);
                std::lock_guard<std::mutex> lock(count_mutex);
                ++counts[i];
            });
        } else {
            anno_graph.call_annotated_nodes(label_in, [&](node_index i) {
                std::lock_guard<std::mutex> lock(count_mutex);
                ++counts[i];
            });
        }
    }
    __atomic_thread_fence(__ATOMIC_ACQUIRE);

    // correct the width of counts, making it single-width
    counts.width(width);

    __atomic_thread_fence(__ATOMIC_RELEASE);
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (size_t i = 0; i < labels_out.size(); ++i) {
        const std::string &label_out = labels_out[i];
        if (graph->is_canonical_mode()) {
            anno_graph.call_annotated_nodes(label_out, [&](node_index i) {
                set_bit(indicator.data(), i, async, memorder);
                std::lock_guard<std::mutex> lock(count_mutex);
                ++counts[2 * i + 1];
            });
        } else {
            anno_graph.call_annotated_nodes(label_out, [&](node_index i) {
                std::lock_guard<std::mutex> lock(count_mutex);
                ++counts[2 * i + 1];
            });
        }
    }
    __atomic_thread_fence(__ATOMIC_ACQUIRE);

    // set the width to be double again
    counts.width(width << 1);

    if (graph->is_canonical_mode()) {
        __atomic_thread_fence(__ATOMIC_RELEASE);
        MaskedDeBruijnGraph masked_graph(
            graph,
            std::make_unique<bit_vector_stat>(std::move(indicator)),
            true
        );

        masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
            auto it = path.rbegin();
            auto rev = unitig;
            reverse_complement(rev.begin(), rev.end());
            graph->map_to_nodes_sequentially(rev, [&](node_index i) {
                assert(i != DeBruijnGraph::npos);
                assert(it != path.rend());
                std::lock_guard<std::mutex> lock(count_mutex);
                counts[i] += counts[*it];
                ++it;
            });
        }, num_threads);
        __atomic_thread_fence(__ATOMIC_ACQUIRE);
    }

    return counts;
}

std::function<uint64_t(SequenceGraph::node_index)>
build_label_counter(const AnnotatedDBG &anno_graph,
                    const std::vector<Label> &labels_to_check) {
    if (anno_graph.get_graph().is_canonical_mode()) {
        return [&anno_graph, labels_to_check](SequenceGraph::node_index i) {
            assert(i != SequenceGraph::npos);
            const auto &graph = anno_graph.get_graph();
            graph.map_to_nodes(graph.get_node_sequence(i), [&](auto j) { i = j; });
            return i != SequenceGraph::npos
                ? std::count_if(labels_to_check.begin(), labels_to_check.end(),
                                [&](const auto &label) {
                                    return anno_graph.has_label(i, label);
                                })
                : 0;
        };
    } else {
        return [&anno_graph, labels_to_check](SequenceGraph::node_index i) {
            assert(i != SequenceGraph::npos);
            return std::count_if(labels_to_check.begin(), labels_to_check.end(),
                [&](const auto &label) { return anno_graph.has_label(i, label); }
            );
        };
    }
}

std::unique_ptr<bitmap>
mask_nodes_by_node_label(const AnnotatedDBG &anno_graph,
                         const std::vector<Label> &labels_in,
                         const std::vector<Label> &labels_out,
                         const std::function<bool(DeBruijnGraph::node_index,
                                                  const LabelCountCallback & /* get_num_labels_in */,
                                                  const LabelCountCallback & /* get_num_labels_out */)> &is_node_in_mask,
                         double min_frequency_for_frequent_label) {
    if (!anno_graph.get_graph().num_nodes())
        return std::make_unique<bitmap_lazy>([](uint64_t) { return false; },
                                             anno_graph.get_graph().max_index() + 1,
                                             0);

    const auto *columns = dynamic_cast<const annot::ColumnCompressed<>*>(
        &anno_graph.get_annotation()
    );

    if (columns) {
        size_t frequent_column_label_min_count
            = anno_graph.get_graph().num_nodes()
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
            sdsl::int_vector<> counts = fill_count_vector(anno_graph,
                                                          labels_in_infrequent,
                                                          labels_out_infrequent);

            // the width of counts is double, since it's both in and out counts interleaved
            uint32_t width = counts.width() >> 1;
            uint64_t int_mask = (uint64_t(1) << width) - 1;

            // Flatten count vector to bitmap if all labels were infrequent
            if (labels_in_frequent.empty() && labels_out_frequent.empty()) {
                sdsl::bit_vector mask(anno_graph.get_graph().max_index() + 1, false);

                call_nonzeros(counts, [&](node_index i, size_t count) {
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

    if (min_frequency_for_frequent_label == 1.0) {
        sdsl::int_vector<> counts = fill_count_vector(anno_graph,
                                                      labels_in,
                                                      labels_out);

        // the width of counts is double, since it's both in and out counts interleaved
        uint32_t width = counts.width() >> 1;
        uint64_t int_mask = (uint64_t(1) << width) - 1;

        sdsl::bit_vector mask(anno_graph.get_graph().max_index() + 1, false);

        call_nonzeros(counts, [&](node_index i, size_t count) {
            if (i != DeBruijnGraph::npos
                    && is_node_in_mask(i, [&]() { return count & int_mask; },
                                          [&]() { return count >> width; }))
                mask[i] = true;
        });

        return std::make_unique<bitmap_vector>(std::move(mask));
    }

    // If all of the labels were frequent, or if the Annotator is not ColumnCompressed,
    // construct a lazy bitmap
    auto count_frequent_in_labels = build_label_counter(anno_graph, labels_in);
    auto count_frequent_out_labels = build_label_counter(anno_graph, labels_out);

    return std::make_unique<bitmap_lazy>([=](node_index i) {
        return i != DeBruijnGraph::npos
            && is_node_in_mask(i, [&]() { return count_frequent_in_labels(i); },
                                  [&]() { return count_frequent_out_labels(i); });
    }, anno_graph.get_graph().max_index() + 1);
}

} // namespace graph
} // namespace mtg
