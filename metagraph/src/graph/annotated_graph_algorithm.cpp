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


namespace mtg {
namespace graph {

using mtg::graph::boss::BOSS;
using mtg::graph::boss::BOSSConstructor;
using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::row_index row_index;


std::pair<sdsl::int_vector<>, sdsl::bit_vector>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<Label> &labels_in,
                  const std::vector<Label> &labels_out) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );
    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    sdsl::int_vector<> counts(graph->max_index() + 1, 0, width * 2);
    sdsl::bit_vector indicator(counts.size(), false);

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
        anno_graph.call_annotated_nodes(label_in, [&](node_index i) {
            set_bit(indicator.data(), i, async, memorder);
            std::lock_guard<std::mutex> lock(count_mutex);
            ++counts[i];
        });
    }

    // correct the width of counts, making it single-width
    counts.width(width);

    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (size_t i = 0; i < labels_out.size(); ++i) {
        const std::string &label_out = labels_out[i];
        anno_graph.call_annotated_nodes(label_out, [&](node_index i) {
            set_bit(indicator.data(), i, async, memorder);
            std::lock_guard<std::mutex> lock(count_mutex);
            ++counts[2 * i + 1];
        });
    }
    __atomic_thread_fence(__ATOMIC_ACQUIRE);

    if (graph->is_canonical_mode()) {
        logger->trace("Adding reverse complements");

        __atomic_thread_fence(__ATOMIC_RELEASE);
        MaskedDeBruijnGraph masked_graph(
            graph,
            std::make_unique<bit_vector_stat>(indicator),
            true
        );

        masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
            auto it = path.rbegin();
            auto rev = unitig;
            reverse_complement(rev.begin(), rev.end());
            graph->map_to_nodes_sequentially(rev, [&](node_index i) {
                assert(i != DeBruijnGraph::npos);
                assert(it != path.rend());
                set_bit(indicator.data(), i, async, memorder);
                std::lock_guard<std::mutex> lock(count_mutex);
                counts[i * 2] += counts[*it * 2];
                counts[i * 2 + 1] += counts[*it * 2 + 1];
                ++it;
            });
        }, num_threads);
        __atomic_thread_fence(__ATOMIC_ACQUIRE);
    }

    // set the width to be double again
    counts.width(width * 2);

    return std::make_pair(std::move(counts), std::move(indicator));
}


std::unique_ptr<bitmap_vector>
mask_nodes_by_unitig(const DeBruijnGraph &graph, const KeepUnitigPath &keep_unitig) {
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
            num_kept_nodes.fetch_add(path.size(), memorder);
            for (node_index node : path) {
                set_bit(unitig_mask.data(), node, parallel, memorder);
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
    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    logger->trace("Generating initial mask");

    // Construct initial masked graph from union of labels in labels_in
    auto [counts, union_mask] = fill_count_vector(anno_graph, labels_in, labels_out);
    uint32_t width = counts.width() / 2;
    uint64_t int_mask = (uint64_t(1) << width) - 1;

    std::shared_ptr<DeBruijnGraph> masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        std::make_unique<bit_vector_stat>(std::move(union_mask)),
        true, // only_valid_nodes_in_mask
        graph_ptr->is_canonical_mode()
    );
    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    if (add_complement && !graph_ptr->is_canonical_mode()) {
        // we can't guarantee that the reverse complement is present, so
        // construct a new subgraph
        logger->trace("Constructing canonical DBG from labeled subgraph");
        std::vector<std::pair<std::string, sdsl::int_vector<>>> contigs;
        std::mutex add_mutex;
        BOSSConstructor constructor(graph_ptr->get_k() - 1, true);
        constructor.add_sequences([&](const CallString &callback) {
            masked_graph->call_sequences([&](const std::string &seq, const auto &path) {
                sdsl::int_vector<> path_counts(path.size(), 0, width << 1);
                auto it = path_counts.begin();
                for (node_index i : path) {
                    *it = counts[i];
                    ++it;
                }

                std::lock_guard<std::mutex> lock(add_mutex);
                contigs.emplace_back(seq, std::move(path_counts));
                callback(seq);
            }, get_num_threads(), true);
        });

        masked_graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), true);
        graph_ptr = masked_graph;

        logger->trace("Reconstructing count vector");
        counts = sdsl::int_vector<>(graph_ptr->max_index() + 1, 0, width * 2);
        for (auto &[seq, seq_counts] : contigs) {
            size_t j = 0;
            graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                assert(!counts[i]);
                counts[i] = seq_counts[j++];
            });
            assert(j == seq_counts.size());

            reverse_complement(seq.begin(), seq.end());

            // counts stores counts for in labels and out labels interleaved,
            // so to add values properly, it should be reshaped first
            counts.width(width);
            seq_counts.width(width);
            j = seq_counts.size();
            graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                j -= 2;
                counts[i * 2] += seq_counts[j];
                counts[(i * 2) + 1] += seq_counts[j + 1];
            });
            assert(!j);

            // reshape back
            counts.width(width * 2);
            seq_counts.width(width * 2);
        }

        logger->trace("Constructed canonical DBG with {} nodes", graph_ptr->num_nodes());
    }

    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");

    tsl::hopscotch_set<std::string> masked_labels;

    if (label_other_fraction != 1.0) {
        masked_labels.insert(labels_in.begin(), labels_in.end());
        masked_labels.insert(labels_out.begin(), labels_out.end());
    }

    return MaskedDeBruijnGraph(
        graph_ptr,
        mask_nodes_by_unitig(*masked_graph, [&](const auto &unitig, const auto &path) {
            size_t in_label_counter = 0;
            size_t out_label_counter = 0;
            size_t label_in_cutoff = std::ceil(label_mask_in_fraction * path.size());
            size_t label_out_cutoff = std::floor(label_mask_out_fraction * path.size());
            double other_frac = label_other_fraction + 1.0 / (path.size() + 1);

            for (size_t i = 0; i < path.size(); ++i) {
                node_index node = path[i];
                uint64_t count = counts[node];
                uint64_t in_count = count & int_mask;
                uint64_t out_count = count >> width;

                if (in_count == labels_in.size())
                    ++in_label_counter;

                if ((path.size() - 1 - i + in_label_counter < label_in_cutoff)
                        || (out_count && ++out_label_counter > label_out_cutoff))
                    return false;
            }

            if (label_other_fraction != 1.0) {
                // TODO: extract these beforehand and construct an annotator
                // discard this unitig if other labels are found with too high frequency
                for (const auto &label : anno_graph.get_labels(unitig, other_frac)) {
                    if (!masked_labels.count(label))
                        return false;
                }
            }

            return true;
        }),
        true
    );
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
            = anno_graph.get_graph().num_nodes() * min_frequency_for_frequent_label;

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
            auto [counts, mask] = fill_count_vector(anno_graph,
                                                    labels_in_infrequent,
                                                    labels_out_infrequent);

            // the width of counts is double, since it's both in and out counts interleaved
            uint32_t width = counts.width() >> 1;
            uint64_t int_mask = (uint64_t(1) << width) - 1;

            // Flatten count vector to bitmap if all labels were infrequent
            if (labels_in_frequent.empty() && labels_out_frequent.empty()) {
                mask[DeBruijnGraph::npos] = false;
                call_nonzeros(counts, [&](node_index i, size_t count) {
                    if (!is_node_in_mask(i, [&]() { return count & int_mask; },
                                            [&]() { return count >> width; }))
                        mask[i] = false;
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
        auto [counts, mask] = fill_count_vector(anno_graph, labels_in, labels_out);
        mask[DeBruijnGraph::npos] = false;

        // the width of counts is double, since it's both in and out counts interleaved
        uint32_t width = counts.width() >> 1;
        uint64_t int_mask = (uint64_t(1) << width) - 1;

        call_nonzeros(counts, [&](node_index i, size_t count) {
            if (!is_node_in_mask(i, [&]() { return count & int_mask; },
                                    [&]() { return count >> width; }))
                mask[i] = false;
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
