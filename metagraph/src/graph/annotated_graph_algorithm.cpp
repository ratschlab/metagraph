#include "annotated_graph_algorithm.hpp"

#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
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

constexpr bool MAKE_BOSS = true;


/**
 * Return an int_vector<>, bitmap pair, each of length anno_graph.get_graph().get_base_graph().max_index().
 * For an index i, the int_vector will contain a packed integer representing the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with. The least significant half of each integer represents the count
 * from labels_in, while the most significant half represents the count from
 * labels_out.
 * The returned bitmap is a binarization of the int_vector
 */
std::pair<sdsl::int_vector<>, std::unique_ptr<bitmap>>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const std::vector<Label> &labels_in,
                                  const std::vector<Label> &labels_out,
                                  size_t num_threads);


template <class GetKeptIntervals>
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads);

template <class KeepNode>
void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads);

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector<> &counts,
                          std::unique_ptr<bitmap>&& mask,
                          bool add_complement,
                          bool make_boss,
                          size_t num_threads);


MaskedDeBruijnGraph mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                                        const std::vector<Label> &labels_in,
                                        const std::vector<Label> &labels_out,
                                        const std::vector<Label> &labels_in_post,
                                        const std::vector<Label> &labels_out_post,
                                        const DifferentialAssemblyConfig &config,
                                        size_t num_threads) {
    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    logger->trace("Generating initial mask");

    // Construct initial masked graph from union of labels in labels_in
    auto count_vector = construct_diff_label_count_vector(
        anno_graph, labels_in, labels_out, num_threads
    );
    auto &[counts, union_mask] = count_vector;

    // in and out counts are stored interleaved in the counts vector
    assert(counts.size() == union_mask->size() * 2);

    auto masked_graph = make_initial_masked_graph(graph_ptr,
                                                  counts, std::move(union_mask),
                                                  config.add_complement, MAKE_BOSS,
                                                  num_threads);

    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");

    tsl::hopscotch_set<std::string> masked_labels;
    tsl::hopscotch_set<std::string> labels_in_post_set(labels_in_post.begin(),
                                                       labels_in_post.end());
    tsl::hopscotch_set<std::string> labels_out_post_set(labels_out_post.begin(),
                                                        labels_out_post.end());

    bool check_other = config.label_mask_other_unitig_fraction != 1.0;

    if (check_other) {
        masked_labels.insert(labels_in.begin(), labels_in.end());
        masked_labels.insert(labels_out.begin(), labels_out.end());
        masked_labels.insert(labels_in_post.begin(), labels_in_post.end());
        masked_labels.insert(labels_out_post.begin(), labels_out_post.end());
    } else if (config.label_mask_in_unitig_fraction == 0.0
            && config.label_mask_out_unitig_fraction == 1.0
            && labels_in_post.empty() && labels_out_post.empty()) {
        if (config.label_mask_in_kmer_fraction == 0.0
                && config.label_mask_out_kmer_fraction == 1.0) {
            logger->trace("Bypassing background filtration");
            return std::move(*masked_graph);
        }

        logger->trace("Filtering by node");
        update_masked_graph_by_node(*masked_graph, [&](node_index node) {
            uint64_t in_count = count_vector.first[node * 2];
            uint64_t out_count = count_vector.first[node * 2 + 1];
            uint64_t sum = in_count + out_count;
            return in_count >= config.label_mask_in_kmer_fraction * sum
                && out_count <= config.label_mask_out_kmer_fraction * sum;
        }, num_threads);

        return std::move(*masked_graph);
    }

    logger->trace("Filtering by unitig");

    size_t num_labels = anno_graph.get_annotation().num_labels();

    auto get_kept_intervals = [&](const auto &unitig, const auto &path)
            -> std::vector<std::pair<size_t, size_t>> {
        sdsl::int_vector<> &counts = count_vector.first;
        sdsl::bit_vector in_mask(path.size(), false);
        sdsl::bit_vector out_mask(path.size(), false);
        sdsl::bit_vector other_mask(check_other ? path.size() : 0, false);

        size_t min_label_in_count = config.label_mask_in_kmer_fraction
                                        * (labels_in.size() + labels_in_post.size());
        size_t max_label_out_count = config.label_mask_out_kmer_fraction
                                        * (labels_out.size() + labels_out_post.size());

        size_t in_kmer_count = 0;
        size_t out_kmer_count = 0;

        for (size_t i = 0; i < path.size(); ++i) {
            if (counts[2 * path[i]] >= min_label_in_count) {
                in_mask[i] = true;
                ++in_kmer_count;
            }

            if (counts[2 * path[i] + 1] > max_label_out_count) {
                out_mask[i] = true;
                ++out_kmer_count;
            }
        }

        if (check_other || labels_in_post.size() || labels_out_post.size()) {
            for (auto &[label, sig] : anno_graph.get_top_label_signatures(unitig, num_labels)) {
                if (check_other && !masked_labels.count(label)) {
                    bitmap_vector(std::move(sig)).add_to(&other_mask);
                    continue;
                }

                if (labels_in_post_set.count(label)) {
                    call_ones(sig, [&](size_t i) {
                        if (++counts[2 * path[i]] >= min_label_in_count) {
                            in_mask[i] = true;
                            ++in_kmer_count;
                        }
                    });
                }

                if (labels_out_post_set.count(label)) {
                    call_ones(sig, [&](size_t i) {
                        if (++counts[2 * path[i] + 1] > max_label_out_count) {
                            out_mask[i] = true;
                            ++out_kmer_count;
                        }
                    });
                }
            }
        }

        size_t begin = next_bit(in_mask, 0);

        if (begin == in_mask.size())
            return {};

        size_t end = prev_bit(in_mask, in_mask.size() - 1) + 1;
        assert(end > begin);

        out_kmer_count -= count_ones(out_mask, 0, begin) + count_ones(out_mask, end, out_mask.size());
        size_t other_kmer_count = check_other ? count_ones(other_mask, begin, end) : 0;

        size_t cur_size = end - begin;
        size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * cur_size);
        size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * cur_size);
        size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * cur_size);

        if (in_kmer_count >= label_in_cutoff && out_kmer_count <= label_out_cutoff
                && other_kmer_count <= other_cutoff)
            return { std::make_pair(begin, end) };


        return {};
    };

    update_masked_graph_by_unitig(*masked_graph, get_kept_intervals, num_threads);

    return std::move(*masked_graph);
}

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector<> &counts,
                          std::unique_ptr<bitmap>&& mask,
                          bool add_complement,
                          bool make_boss,
                          size_t num_threads) {
    // counts is a double-length vector storing the in-label and out-label
    // counts interleaved
    assert(counts.size() == mask->size() * 2);

    graph_ptr = std::shared_ptr<const DeBruijnGraph>(
        std::shared_ptr<const DeBruijnGraph>(), &graph_ptr->get_base_graph()
    );

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        std::move(mask),
        true,
        graph_ptr->get_mode()
    );

    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    bool masked_canonical = add_complement
        || graph_ptr->get_mode() == DeBruijnGraph::CANONICAL;

    if (make_boss
            || (add_complement && graph_ptr->get_mode() != DeBruijnGraph::CANONICAL)) {
        // we can't guarantee that the reverse complement is present, so
        // construct a new subgraph
        logger->trace("Constructing BOSS from labeled subgraph");
        std::vector<std::pair<std::string, sdsl::int_vector<>>> contigs;
        std::mutex add_mutex;
        BOSSConstructor constructor(
            graph_ptr->get_k() - 1,
            masked_canonical,
            0 /* count width */, "" /* suffix */, 1 /* num_threads */,
            masked_graph->num_nodes() * 32
        );
        masked_graph->call_sequences([&](const std::string &seq, const auto &path) {
            sdsl::int_vector<> path_counts(path.size() * 2, 0, counts.width());
            auto it = path_counts.begin();
            for (node_index i : path) {
                *it = counts[i * 2];
                ++it;

                *it = counts[i * 2 + 1];
                ++it;
            }

            std::lock_guard<std::mutex> lock(add_mutex);
            contigs.emplace_back(seq, std::move(path_counts));
            constructor.add_sequence(seq);
        }, num_threads, masked_canonical);

        auto dbg_succ = std::make_shared<DBGSuccinct>(
            new BOSS(&constructor),
            masked_canonical ? DeBruijnGraph::CANONICAL : DeBruijnGraph::BASIC
        );

        // instead of keeping multiple masks (one for valid nodes and another
        // for the masked graph), transfer the valid node mask to the masked graph
        dbg_succ->mask_dummy_kmers(num_threads, false);
        std::unique_ptr<bit_vector> dummy_mask(dbg_succ->release_mask());
        graph_ptr = dbg_succ;
        sdsl::bit_vector new_indicator(graph_ptr->max_index() + 1, false);
        dummy_mask->add_to(&new_indicator);
        masked_graph = std::make_shared<MaskedDeBruijnGraph>(
            graph_ptr,
            std::make_unique<bitmap_vector>(std::move(new_indicator)),
            true,
            masked_canonical ? DeBruijnGraph::CANONICAL : DeBruijnGraph::BASIC
        );

        logger->trace("Reconstructing count vector");
        counts = aligned_int_vector((graph_ptr->max_index() + 1) * 2,
                                    0, counts.width(), 16);

        std::atomic_thread_fence(std::memory_order_release);

        std::mutex vector_backup_mutex;
        constexpr std::memory_order memorder = std::memory_order_relaxed;

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l = 0; l < contigs.size(); ++l) {
            auto &[seq, seq_counts] = contigs[l];
            size_t j = 0;
            graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                atomic_exchange(counts, i * 2, contigs[l].second[j],
                                vector_backup_mutex, memorder);
                atomic_exchange(counts, i * 2 + 1, contigs[l].second[j + 1],
                                vector_backup_mutex, memorder);
                j += 2;
            });
            assert(j == seq_counts.size());
        }

        if (masked_canonical) {
            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t l = 0; l < contigs.size(); ++l) {
                auto &[seq, seq_counts] = contigs[l];
                reverse_complement(seq.begin(), seq.end());

                // counts stores counts for in labels and out labels interleaved,
                // so to add values properly, it should be reshaped first
                size_t j = seq_counts.size();
                graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                    uint64_t old_val = atomic_exchange(counts, i * 2 + 1,
                                                       contigs[l].second[--j],
                                                       vector_backup_mutex, memorder);
                    std::ignore = old_val;
                    assert(!old_val);

                    old_val = atomic_exchange(counts, i * 2, contigs[l].second[--j],
                                              vector_backup_mutex, memorder);
                    std::ignore = old_val;
                    assert(!old_val);
                });
                assert(!j);
            }
        }

        std::atomic_thread_fence(std::memory_order_acquire);

        logger->trace("Constructed BOSS with {} nodes", graph_ptr->num_nodes());
    }

    assert(counts.size() == (graph_ptr->max_index() + 1) * 2);

    return masked_graph;
}

/**
 * Return an int_vector<>, bitmap pair, each of length anno_graph.get_graph().get_base_graph().max_index().
 * For an index i, the int_vector will contain a packed integer representing the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with. The least significant half of each integer represents the count
 * from labels_in, while the most significant half represents the count from
 * labels_out.
 * The returned bitmap is a binarization of the int_vector
 */
std::pair<sdsl::int_vector<>, std::unique_ptr<bitmap>>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const std::vector<Label> &labels_in,
                                  const std::vector<Label> &labels_out,
                                  size_t num_threads) {
    const DeBruijnGraph &graph = anno_graph.get_graph().get_base_graph();

    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    sdsl::bit_vector indicator(graph.max_index() + 1, false);

    // the in and out counts are stored interleaved
    sdsl::int_vector<> counts = aligned_int_vector(indicator.size() * 2, 0, width, 16);

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    const auto &binmat = anno_graph.get_annotation().get_matrix();

    std::vector<uint64_t> label_in_codes(labels_in.size());
    std::vector<uint64_t> label_out_codes(labels_out.size());
    for (size_t i = 0; i < labels_in.size(); ++i) {
        label_in_codes[i] = label_encoder.encode(labels_in[i]);
    }
    for (size_t i = 0; i < labels_out.size(); ++i) {
        label_out_codes[i] = label_encoder.encode(labels_out[i]);
    }

    constexpr std::memory_order memorder = std::memory_order_relaxed;

    std::mutex vector_backup_mutex;
    std::atomic_thread_fence(std::memory_order_release);
    bool async = num_threads > 1;
    binmat.call_columns(label_in_codes,
        [&](auto, const bitmap &rows) {
            rows.call_ones([&](auto r) {
                node_index i = AnnotatedDBG::anno_to_graph_index(r);
                set_bit(indicator.data(), i, async, memorder);
                atomic_fetch_and_add(counts, i * 2, 1, vector_backup_mutex, memorder);
            });
        },
        num_threads
    );

    binmat.call_columns(label_out_codes,
        [&](auto, const bitmap &rows) {
            rows.call_ones([&](auto r) {
                node_index i = AnnotatedDBG::anno_to_graph_index(r);
                set_bit(indicator.data(), i, async, memorder);
                atomic_fetch_and_add(counts, i * 2 + 1, 1, vector_backup_mutex, memorder);
            });
        },
        num_threads
    );

    std::atomic_thread_fence(std::memory_order_acquire);

    std::unique_ptr<bitmap> union_mask = std::make_unique<bitmap_vector>(
        std::move(indicator)
    );

    // TODO: this function doesn't work if MAKE_BOSS is false and the complements
    //       of k-mers should also be present
    static_assert(MAKE_BOSS);

    return std::make_pair(std::move(counts), std::move(union_mask));
}


template <class GetKeptIntervals>
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads) {
    std::atomic<uint64_t> kept_unitigs(0);
    std::atomic<uint64_t> total_unitigs(0);
    std::atomic<uint64_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;
    constexpr std::memory_order memorder = std::memory_order_relaxed;

    sdsl::bit_vector mask = dynamic_cast<const bitmap_vector&>(masked_graph.get_mask()).data();

    std::atomic_thread_fence(std::memory_order_release);

    masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
        total_unitigs.fetch_add(1, memorder);

        size_t last = 0;
        for (const auto &pair : get_kept_intervals(unitig, path)) {
            const auto &[begin, end] = pair;
            kept_unitigs.fetch_add(1, memorder);
            num_kept_nodes.fetch_add(end - begin, memorder);
            for ( ; last < begin; ++last) {
                unset_bit(mask.data(), path[last], parallel, memorder);
            }
            last = end;
        }

        for ( ; last < path.size(); ++last) {
            unset_bit(mask.data(), path[last], parallel, memorder);
        }

    }, num_threads);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs);
}

template <class KeepNode>
void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads) {
    bool parallel = num_threads > 1;
    constexpr std::memory_order memorder = std::memory_order_relaxed;

    std::atomic<size_t> kept_nodes(0);
    size_t total_nodes = masked_graph.num_nodes();

    const auto &in_mask = dynamic_cast<const bitmap_vector&>(masked_graph.get_mask()).data();
    sdsl::bit_vector mask(in_mask);

    std::atomic_thread_fence(std::memory_order_release);
    call_ones(in_mask, [&](node_index node) {
        if (keep_node(node)) {
            kept_nodes.fetch_add(1, memorder);
        } else {
            unset_bit(mask.data(), node, parallel, memorder);
        }
    }, parallel, memorder);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);
}

} // namespace graph
} // namespace mtg
