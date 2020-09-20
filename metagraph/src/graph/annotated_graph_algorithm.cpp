#include "annotated_graph_algorithm.hpp"

#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
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

constexpr bool MAKE_BOSS = true;


template <class GetKeptIntervals>
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads,
                                   bool update_in_place);

template <class KeepNode>
void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads,
                                 bool update_in_place);

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> &graph_ptr,
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
                                        size_t num_threads,
                                        const sdsl::int_vector<> *init_counts) {
    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    // TODO: find a way to avoid this hack
    // Since call_unitigs/sequences on a masked DBGSuccinct makes a copy of the
    // mask before traversal, we can safely update the mask in the callback
    bool update_in_place = static_cast<bool>(
        std::dynamic_pointer_cast<const DBGSuccinct>(graph_ptr)
    );

    logger->trace("Generating initial mask");

    // Construct initial masked graph from union of labels in labels_in
    auto [counts, union_mask] = fill_count_vector(anno_graph,
                                                  labels_in, labels_out,
                                                  num_threads,
                                                  update_in_place,
                                                  init_counts);

    // counts is a double-width, interleaved vector where the significant bits
    // represent the out-label count and the least significant bits represent
    // the in-label count
    uint32_t width = counts.width() / 2;
    uint64_t int_mask = (uint64_t(1) << width) - 1;

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

    if (config.label_mask_other_unitig_fraction != 1.0) {
        masked_labels.insert(labels_in.begin(), labels_in.end());
        masked_labels.insert(labels_out.begin(), labels_out.end());
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
            uint64_t count = counts[node];
            uint64_t in_count = count & int_mask;
            uint64_t out_count = count >> width;
            uint64_t sum = in_count + out_count;
            return in_count >= config.label_mask_in_kmer_fraction * sum
                && out_count <= config.label_mask_out_kmer_fraction * sum;
        }, num_threads, update_in_place);

        return std::move(*masked_graph);
    }

    logger->trace("Filtering by unitig");

    size_t num_labels = anno_graph.get_annotation().num_labels();
    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    counts.width(width);

    update_masked_graph_by_unitig(*masked_graph,
                                  [&](const auto &unitig, const auto &path)
                                      -> std::vector<std::pair<size_t, size_t>> {
        sdsl::bit_vector in_mask(path.size(), false);
        sdsl::bit_vector out_mask(path.size(), false);
        sdsl::bit_vector other_mask(check_other ? path.size() : 0, false);

        size_t min_label_in_count = config.label_mask_in_kmer_fraction * labels_in.size();
        size_t max_label_out_count = config.label_mask_out_kmer_fraction * labels_out.size();

        size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * path.size());
        size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * path.size());
        size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * path.size());
        size_t in_kmer_count = 0;
        size_t out_kmer_count = 0;

        for (size_t i = 0; i < path.size(); ++i) {
            if (counts[2 * path[i]] >= min_label_in_count && !in_mask[i]) {
                in_mask[i] = true;
                ++in_kmer_count;
            }

            if (counts[2 * path[i] + 1] > max_label_out_count && !out_mask[i]) {
                out_mask[i] = true;
                ++out_kmer_count;
            }
        }

        if (check_other
                || (in_kmer_count < label_in_cutoff && labels_in_post.size())
                || (out_kmer_count < label_out_cutoff && labels_out_post.size())) {
            for (auto &[label, sig] : anno_graph.get_top_label_signatures(unitig, num_labels)) {
                if (check_other && !masked_labels.count(label))
                    bitmap_vector(std::move(sig)).add_to(&other_mask);

                if (in_kmer_count < label_in_cutoff && labels_in_post_set.count(label)) {
                    for (size_t i = 0; i < path.size(); ++i) {
                        if (!sig[i])
                            continue;

                        if (!in_mask[i] && ++counts[2 * path[i]] >= min_label_in_count) {
                            in_mask[i] = true;
                            ++in_kmer_count;
                        }
                    }
                }

                if (out_kmer_count < label_out_cutoff && labels_out_post_set.count(label)) {
                    for (size_t i = 0; i < path.size(); ++i) {
                        if (!sig[i])
                            continue;

                        if (!out_mask[i] && ++counts[2 * path[i] + 1] >= max_label_out_count) {
                            out_mask[i] = true;
                            ++out_kmer_count;
                        }
                    }
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

        if (out_kmer_count > label_out_cutoff)
            return {};

        if (other_kmer_count > other_cutoff)
            return {};

        return { std::make_pair(begin, end) };

    }, num_threads, update_in_place);

    return std::move(*masked_graph);
}


std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> &graph_ptr,
                          sdsl::int_vector<> &counts,
                          std::unique_ptr<bitmap>&& mask,
                          bool add_complement,
                          bool make_boss,
                          size_t num_threads) {
    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        std::move(mask),
        graph_ptr->is_canonical_mode()
    );
    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());
    bool masked_canonical = add_complement || graph_ptr->is_canonical_mode();

    if (make_boss || (add_complement && !graph_ptr->is_canonical_mode())) {
        // we can't guarantee that the reverse complement is present, so
        // construct a new subgraph
        logger->trace("Constructing BOSS from labeled subgraph");
        uint32_t width = counts.width() / 2;
        std::vector<std::pair<std::string, sdsl::int_vector<>>> contigs;
        std::mutex add_mutex;
        BOSSConstructor constructor(
            graph_ptr->get_k() - 1,
            masked_canonical,
            0 /* count width */, "" /* suffix */, 1 /* num_threads */,
            masked_graph->num_nodes() * 32
        );
        masked_graph->call_sequences([&](const std::string &seq, const auto &path) {
            sdsl::int_vector<> path_counts(path.size(), 0, width * 2);
            auto it = path_counts.begin();
            for (node_index i : path) {
                *it = counts[i];
                ++it;
            }

            std::lock_guard<std::mutex> lock(add_mutex);
            contigs.emplace_back(seq, std::move(path_counts));
            constructor.add_sequence(seq);
        }, num_threads, masked_canonical);

        auto dbg_succ = std::make_shared<DBGSuccinct>(new BOSS(&constructor),
                                                      masked_canonical);

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
            masked_canonical
        );

        logger->trace("Reconstructing count vector");
        counts = aligned_int_vector(graph_ptr->max_index() + 1, 0, width * 2, 16);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l = 0; l < contigs.size(); ++l) {
            auto &[seq, seq_counts] = contigs[l];
            size_t j = 0;
            graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                atomic_set(counts, i, seq_counts[j++]);
            });
            assert(j == seq_counts.size());
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        if (masked_canonical) {
            // reshape to normal width
            counts.width(width);

            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t l = 0; l < contigs.size(); ++l) {
                auto &[seq, seq_counts] = contigs[l];
                seq_counts.width(width);
                reverse_complement(seq.begin(), seq.end());

                // counts stores counts for in labels and out labels interleaved,
                // so to add values properly, it should be reshaped first
                size_t j = seq_counts.size();
                graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                    j -= 2;
                    atomic_set(counts, i * 2, seq_counts[j]);
                    atomic_set(counts, i* 2 + 1, seq_counts[j + 1]);
                });
                assert(!j);
                seq_counts.width(width * 2);
            }

            __atomic_thread_fence(__ATOMIC_ACQUIRE);

            // reshape back
            counts.width(width * 2);
        }

        logger->trace("Constructed BOSS with {} nodes", graph_ptr->num_nodes());
    }

    assert(counts.size() == graph_ptr->max_index() + 1);

    return masked_graph;
}


std::pair<sdsl::int_vector<>, std::unique_ptr<bitmap>>
fill_count_vector(const AnnotatedDBG &anno_graph,
                  const std::vector<Label> &labels_in,
                  const std::vector<Label> &labels_out,
                  size_t num_threads,
                  bool update_in_place,
                  const sdsl::int_vector<> *init_counts) {
    // at this stage, the width of counts is twice what it should be, since
    // the intention is to store the in label and out label counts interleaved
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );
    size_t width = sdsl::bits::hi(std::max(labels_in.size(), labels_out.size())) + 1;
    sdsl::bit_vector indicator(graph->max_index() + 1, false);
    sdsl::int_vector<> counts = aligned_int_vector(graph->max_index() + 1, 0, width * 2, 16);

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    const auto &binmat = anno_graph.get_annotation().get_matrix();

    if (init_counts) {
        assert(init_counts->size() == graph->max_index() + 1);
        call_nonzeros(*init_counts, [&](uint64_t i, uint64_t val) {
            counts[i] = val;
        });
    }

    std::vector<uint64_t> label_in_codes(labels_in.size());
    std::vector<uint64_t> label_out_codes(labels_out.size());
    for (size_t i = 0; i < labels_in.size(); ++i) {
        label_in_codes[i] = label_encoder.encode(labels_in[i]);
    }
    for (size_t i = 0; i < labels_out.size(); ++i) {
        label_out_codes[i] = label_encoder.encode(labels_out[i]);
    }

    #pragma omp parallel num_threads(num_threads)
    #pragma omp single
    {
        bool async = num_threads > 1;
        binmat.slice_columns(label_in_codes, [&](auto, const bitmap &rows) {
            rows.call_ones([&](auto r) {
                node_index i = AnnotatedDBG::anno_to_graph_index(r);
                set_bit(indicator.data(), i, async, __ATOMIC_RELAXED);
                atomic_increment(counts, i);
            });
        });

        #pragma omp taskwait
        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        // correct the width of counts, making it single-width
        counts.width(width);

        binmat.slice_columns(label_out_codes, [&](auto, const bitmap &rows) {
            rows.call_ones([&](auto r) {
                node_index i = AnnotatedDBG::anno_to_graph_index(r);
                set_bit(indicator.data(), i, async, __ATOMIC_RELAXED);
                atomic_increment(counts, i * 2 + 1);
            });
        });

        #pragma omp taskwait
        __atomic_thread_fence(__ATOMIC_ACQUIRE);
    }

    std::unique_ptr<bitmap> union_mask = std::make_unique<bitmap_vector>(
        std::move(indicator)
    );

    if (!MAKE_BOSS && graph->is_canonical_mode()) {
        logger->trace("Adding reverse complements");

        MaskedDeBruijnGraph masked_graph(graph, std::move(union_mask));

        std::mutex count_mutex;
        masked_graph.update_mask([&](const auto &callback) {
            masked_graph.call_sequences([&](const std::string &seq, const auto &path) {
                auto it = path.rbegin();
                auto rev = seq;
                reverse_complement(rev.begin(), rev.end());
                graph->map_to_nodes_sequentially(rev, [&](node_index i) {
                    std::lock_guard<std::mutex> lock(count_mutex);
                    assert(i != DeBruijnGraph::npos);
                    assert(it != path.rend());
                    callback(i, true);
                    counts[i * 2] += counts[*it * 2];
                    counts[i * 2 + 1] += counts[*it * 2 + 1];
                    ++it;
                });
            }, num_threads);
        }, update_in_place);

        union_mask = std::unique_ptr<bitmap>(masked_graph.release_mask());
    }

    // set the width to be double again
    counts.width(width * 2);

    return std::make_pair(std::move(counts), std::move(union_mask));
}


template <class GetKeptIntervals>
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads,
                                   bool update_in_place) {
    std::atomic<size_t> kept_unitigs(0);
    std::atomic<size_t> total_unitigs(0);
    std::atomic<size_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;
    constexpr auto memorder = std::memory_order_relaxed;

    std::atomic_thread_fence(std::memory_order_release);
    masked_graph.update_mask([&](const auto &callback) {
        masked_graph.call_unitigs([&](const std::string &unitig, const auto &path) {
            total_unitigs.fetch_add(1, memorder);

            size_t last = 0;
            for (const auto &[begin, end] : get_kept_intervals(unitig, path)) {
                kept_unitigs.fetch_add(1, memorder);
                num_kept_nodes.fetch_add(end - begin, memorder);
                for ( ; last < begin; ++last) {
                    callback(path[last], false);
                }
                last = end;
            }

            for ( ; last < path.size(); ++last) {
                callback(path[last], false);
            }

        }, num_threads);
    }, update_in_place, parallel, memorder);
    std::atomic_thread_fence(std::memory_order_acquire);

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs);
}

template <class KeepNode>
void update_masked_graph_by_node(MaskedDeBruijnGraph &masked_graph,
                                 const KeepNode &keep_node,
                                 size_t num_threads,
                                 bool update_in_place) {
    // TODO: parallel
    std::ignore = num_threads;
    bool parallel = false;
    constexpr auto memorder = std::memory_order_relaxed;

    size_t kept_nodes = 0;
    size_t total_nodes = masked_graph.num_nodes();

    masked_graph.update_mask([&](const auto &callback) {
        masked_graph.call_nodes([&](node_index node) {
            if (keep_node(node)) {
                ++kept_nodes;
            } else {
                callback(node, false);
            }
        });
    }, update_in_place, parallel, memorder);

    logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);
}


} // namespace graph
} // namespace mtg
