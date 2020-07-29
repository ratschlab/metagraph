#include "annotated_graph_algorithm.hpp"

#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vector_map.hpp"
#include "common/threads/threading.hpp"
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
typedef AnnotatedDBG::Annotator::Label Label;

typedef std::function<size_t()> LabelCountCallback;

typedef std::function<bool(const std::string&, const std::vector<node_index>&)> KeepUnitig;
typedef std::function<bool(node_index)> KeepNode;


std::pair<sdsl::int_vector<>, sdsl::bit_vector>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<Label> &labels_in,
                  const std::vector<Label> &labels_out,
                  size_t num_threads);

void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const KeepUnitig &keep_unitig,
                                   size_t num_threads);

void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads);

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> &graph_ptr,
                          sdsl::int_vector<> &counts,
                          sdsl::bit_vector&& union_mask,
                          bool add_complement,
                          size_t num_threads);


MaskedDeBruijnGraph mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                                        const std::vector<Label> &labels_in,
                                        const std::vector<Label> &labels_out,
                                        const DifferentialAssemblyConfig &config,
                                        size_t num_threads) {
    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    logger->trace("Generating initial mask");

    // Construct initial masked graph from union of labels in labels_in
    auto [counts, union_mask] = fill_count_vector(anno_graph,
                                                  labels_in, labels_out,
                                                  num_threads);

    // counts is a double-width, interleaved vector where the significant bits
    // represent the out-label count and the least significant bits represent
    // the in-label count
    uint32_t width = counts.width() / 2;
    uint64_t int_mask = (uint64_t(1) << width) - 1;

    auto masked_graph = make_initial_masked_graph(graph_ptr,
                                                  counts, std::move(union_mask),
                                                  config.add_complement, num_threads);

    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");

    tsl::hopscotch_set<std::string> masked_labels;

    if (config.label_mask_other_unitig_fraction != 1.0) {
        masked_labels.insert(labels_in.begin(), labels_in.end());
        masked_labels.insert(labels_out.begin(), labels_out.end());
    } else if (config.label_mask_in_unitig_fraction == 0.0
            && config.label_mask_out_unitig_fraction == 1.0) {
        if (config.label_mask_in_kmer_fraction == 0.0
                && config.label_mask_out_kmer_fraction == 1.0) {
            logger->trace("Bypassing background filtration");
            return std::move(*masked_graph);
        }

        logger->trace("Filtering by node");
        update_masked_graph_by_node(*masked_graph, [&](node_index node) {
            uint64_t count = counts[node];
            uint64_t in_count = count & int_mask;
            uint64_t out_count = count >> width;
            uint64_t sum = in_count + out_count;
            return in_count >= config.label_mask_in_kmer_fraction * sum
                && out_count <= config.label_mask_out_kmer_fraction * sum;
        }, num_threads);

        return std::move(*masked_graph);
    }

    logger->trace("Filtering by unitig");

    update_masked_graph_by_unitig(*masked_graph, [&](const auto &unitig, const auto &path) {
        size_t in_label_counter = 0;
        size_t out_label_counter = 0;
        size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * path.size());
        size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * path.size());
        double other_frac = config.label_mask_other_unitig_fraction + 1.0 / (path.size() + 1);

        for (size_t i = 0; i < path.size(); ++i) {
            node_index node = path[i];
            uint64_t count = counts[node];
            uint64_t in_count = count & int_mask;
            uint64_t out_count = count >> width;

            if (in_count >= config.label_mask_in_kmer_fraction * labels_in.size())
                ++in_label_counter;

            if ((path.size() - 1 - i + in_label_counter < label_in_cutoff)
                    || (out_count > config.label_mask_out_kmer_fraction * labels_out.size()
                        && ++out_label_counter > label_out_cutoff))
                return false;
        }

        if (config.label_mask_other_unitig_fraction != 1.0) {
            // TODO: extract these beforehand and construct an annotator
            // discard this unitig if other labels are found with too high frequency
            for (const auto &label : anno_graph.get_labels(unitig, other_frac)) {
                if (!masked_labels.count(label))
                    return false;
            }
        }

        return true;
    }, num_threads);

    return std::move(*masked_graph);
}


std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> &graph_ptr,
                          sdsl::int_vector<> &counts,
                          sdsl::bit_vector&& union_mask,
                          bool add_complement,
                          size_t num_threads) {
    std::unique_ptr<bitmap> mask;

    // This hack relies on the fact that MaskedDeBruijnGraph::call_sequences
    // makes a copy of the mask when the underlying graph is DBGSuccinct.
    // By making the mask derive from bitmap_vector, we are indicating that
    // it can be safely modified by the callback of call_sequences.
    if (std::dynamic_pointer_cast<const DBGSuccinct>(graph_ptr)) {
        mask = std::make_unique<bitmap_vector>(std::move(union_mask));
    } else {
        mask = std::make_unique<bit_vector_stat>(std::move(union_mask));
    }

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        std::move(mask),
        true, // only_valid_nodes_in_mask
        graph_ptr->is_canonical_mode()
    );
    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    if (add_complement && !graph_ptr->is_canonical_mode()) {
        // we can't guarantee that the reverse complement is present, so
        // construct a new subgraph
        logger->trace("Constructing canonical DBG from labeled subgraph");
        uint32_t width = counts.width() / 2;
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
            }, num_threads, true);
        });

        auto dbg_succ = std::make_shared<DBGSuccinct>(new BOSS(&constructor), true);
        dbg_succ->mask_dummy_kmers(num_threads, false);
        graph_ptr = dbg_succ;
        auto new_indicator = std::make_unique<bitmap_vector>(graph_ptr->max_index() + 1, true);
        new_indicator->set(DeBruijnGraph::npos, false);
        masked_graph = std::make_shared<MaskedDeBruijnGraph>(
            graph_ptr,
            std::move(new_indicator),
            true
        );

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

    return masked_graph;
}


std::pair<sdsl::int_vector<>, sdsl::bit_vector>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<Label> &labels_in,
                  const std::vector<Label> &labels_out,
                  size_t num_threads) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    // in the beginning, it's the correct size, but double width
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );
    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    sdsl::int_vector<> counts(graph->max_index() + 1, 0, width * 2);
    sdsl::bit_vector indicator(counts.size(), false);

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

        std::unique_ptr<bit_vector_stat> mask;
        sdsl::bit_vector *updated_indicator;
        auto dbg_succ = std::dynamic_pointer_cast<const DBGSuccinct>(graph);
        if (dbg_succ) {
            // This hack relies on the fact that call_sequences on a masked
            // DBGSuccinct makes a copy of the mask before traversing, so directly
            // modifying indicator won't have any side effects. This way, we
            // avoid making another bit vector
            mask = std::make_unique<bit_vector_stat>(std::move(indicator));
            updated_indicator = &const_cast<sdsl::bit_vector&>(mask->data());
        } else {
            mask = std::make_unique<bit_vector_stat>(indicator);
            updated_indicator = &indicator;
        }

        MaskedDeBruijnGraph masked_graph(graph, std::move(mask), true);

        __atomic_thread_fence(__ATOMIC_RELEASE);
        masked_graph.call_sequences([&](const std::string &seq, const auto &path) {
            auto it = path.rbegin();
            auto rev = seq;
            reverse_complement(rev.begin(), rev.end());
            graph->map_to_nodes_sequentially(rev, [&](node_index i) {
                assert(i != DeBruijnGraph::npos);
                assert(it != path.rend());
                set_bit(updated_indicator->data(), i, async, memorder);
                std::lock_guard<std::mutex> lock(count_mutex);
                counts[i * 2] += counts[*it * 2];
                counts[i * 2 + 1] += counts[*it * 2 + 1];
                ++it;
            });
        }, num_threads);
        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        if (dbg_succ)
            std::swap(indicator, *updated_indicator);
    }

    // set the width to be double again
    counts.width(width * 2);

    return std::make_pair(std::move(counts), std::move(indicator));
}


void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const KeepUnitig &keep_unitig,
                                   size_t num_threads) {
    std::atomic<size_t> kept_unitigs(0);
    std::atomic<size_t> total_unitigs(0);
    std::atomic<size_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;
    constexpr auto memorder = std::memory_order_relaxed;

    // If the underlying storage derives from bitmap_dyn, then we can safely
    // update it. Otherwise, populate a new bit vector and replace the one in
    // masked_graph after.
    // Also, we require that the underlying mask is bitmap_vector so we can update
    // it with atomic operations.
    auto updateable_mask = dynamic_cast<bitmap_vector*>(masked_graph.get_updateable_mask());
    std::unique_ptr<sdsl::bit_vector> tmp_mask;
    sdsl::bit_vector &unitig_mask = updateable_mask
        ? const_cast<sdsl::bit_vector&>(updateable_mask->data())
        : *(tmp_mask = std::make_unique<sdsl::bit_vector>(masked_graph.max_index() + 1, false));

    if (updateable_mask) {
        assert(!tmp_mask);
        std::atomic_thread_fence(std::memory_order_release);
        masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
            if (keep_unitig(unitig, path)) {
                kept_unitigs.fetch_add(1, memorder);
                num_kept_nodes.fetch_add(path.size(), memorder);
            } else {
                for (node_index node : path) {
                    unset_bit(unitig_mask.data(), node, parallel, memorder);
                }
            }
            total_unitigs.fetch_add(1, memorder);
        }, num_threads);
        std::atomic_thread_fence(std::memory_order_acquire);
    } else {
        assert(tmp_mask);
        std::atomic_thread_fence(std::memory_order_release);
        masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
            if (keep_unitig(unitig, path)) {
                kept_unitigs.fetch_add(1, memorder);
                num_kept_nodes.fetch_add(path.size(), memorder);
                for (node_index node : path) {
                    set_bit(unitig_mask.data(), node, parallel, memorder);
                }
            }
            total_unitigs.fetch_add(1, memorder);
        }, num_threads);
        std::atomic_thread_fence(std::memory_order_acquire);
    }

    if (tmp_mask)
        masked_graph.set_mask(new bit_vector_stat(std::move(*tmp_mask)));

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs);
}

void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads) {
    std::atomic<size_t> kept_nodes(0);
    bool parallel = num_threads > 1;
    constexpr auto memorder = std::memory_order_relaxed;

    // If the underlying storage derives from bitmap_dyn, then we can safely
    // update it. Otherwise, populate a new bit vector and replace the one in
    // masked_graph after.
    // Also, we require that the underlying mask is bitmap_vector so we can update
    // it with atomic operations.
    auto updateable_mask = dynamic_cast<bitmap_vector*>(masked_graph.get_updateable_mask());
    std::unique_ptr<sdsl::bit_vector> tmp_mask;
    sdsl::bit_vector &unitig_mask = updateable_mask
        ? const_cast<sdsl::bit_vector&>(updateable_mask->data())
        : *(tmp_mask = std::make_unique<sdsl::bit_vector>(masked_graph.max_index() + 1, false));

    if (updateable_mask) {
        assert(!tmp_mask);
        std::atomic_thread_fence(std::memory_order_release);
        call_ones(unitig_mask, [&](node_index node) {
            if (!keep_node(node))
                unset_bit(unitig_mask.data(), node, parallel, memorder);
        }, parallel, memorder);
        std::atomic_thread_fence(std::memory_order_acquire);
    } else {
        // TODO: parallel
        assert(tmp_mask);
        masked_graph.call_nodes([&](node_index node) {
            if (keep_node(node))
                unitig_mask[node] = true;
        });
    }

    if (tmp_mask)
        masked_graph.set_mask(new bit_vector_stat(std::move(*tmp_mask)));
}


} // namespace graph
} // namespace mtg
