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
        }, num_threads, update_in_place);

        return std::move(*masked_graph);
    }

    logger->trace("Filtering by unitig");

    size_t num_labels = anno_graph.get_annotation().num_labels();

    update_masked_graph_by_unitig(*masked_graph,
                                  [&](const auto &unitig, const auto &path)
                                      -> std::vector<std::pair<size_t, size_t>> {
        sdsl::bit_vector in_mask(path.size(), false);
        sdsl::bit_vector out_mask(path.size(), true);
        sdsl::bit_vector other_mask(path.size(), true);

        size_t min_label_in_count = config.label_mask_in_kmer_fraction * labels_in.size();
        size_t max_label_out_count = config.label_mask_out_kmer_fraction * labels_out.size();

        size_t in_label_counter = 0;
        size_t out_label_counter = 0;
        size_t other_label_counter = 0;

        for (size_t i = 0; i < path.size(); ++i) {
            uint64_t count = counts[path[i]];
            uint64_t in_count = count & int_mask;
            uint64_t out_count = count >> width;

            if (in_count >= min_label_in_count) {
                ++in_label_counter;
                in_mask[i] = true;
            }

            if (out_count > max_label_out_count) {
                ++out_label_counter;
                out_mask[i] = false;
            }
        }

        if (config.label_mask_other_unitig_fraction != 1.0) {
            for (auto &[label, sig] : anno_graph.get_top_label_signatures(unitig, num_labels)) {
                if (masked_labels.count(label))
                    continue;

                for (size_t i = 0; i < sig.size(); ++i) {
                    ++other_label_counter;
                    other_mask[i] = false;
                }
            }
        }

        size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * path.size());
        size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * path.size());
        size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * path.size());

        if (in_label_counter < label_in_cutoff
                || out_label_counter > label_out_cutoff
                || other_label_counter > other_cutoff) {
            return {};
        }

        return { std::make_pair(0, path.size()) };

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
        BOSSConstructor constructor(graph_ptr->get_k() - 1, masked_canonical);
        constructor.add_sequences([&](const CallString &callback) {
            masked_graph->call_sequences([&](const std::string &seq, const auto &path) {
                sdsl::int_vector<> path_counts(path.size(), 0, width * 2);
                auto it = path_counts.begin();
                for (node_index i : path) {
                    *it = counts[i];
                    ++it;
                }

                std::lock_guard<std::mutex> lock(add_mutex);
                contigs.emplace_back(seq, std::move(path_counts));
                callback(seq);
            }, num_threads, masked_canonical);
        });

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
        counts = sdsl::int_vector<>(graph_ptr->max_index() + 1, 0, width * 2);
        for (auto &[seq, seq_counts] : contigs) {
            size_t j = 0;
            graph_ptr->map_to_nodes_sequentially(seq, [&](node_index i) {
                assert(!counts[i]);
                counts[i] = seq_counts[j++];
            });
            assert(j == seq_counts.size());

            if (masked_canonical) {
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
    sdsl::int_vector<> counts(graph->max_index() + 1, 0, width * 2);
    sdsl::bit_vector indicator(counts.size(), false);

    if (init_counts) {
        assert(init_counts->size() == graph->max_index() + 1);
        call_nonzeros(*init_counts, [&](uint64_t i, uint64_t val) {
            counts[i] = val;
        });
    }

    // TODO: replace locked increment operations on int_vector<> with actual
    //       atomic operations when we figure out how to align int_vector<> storage
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads) default(shared)
    for (size_t i = 0; i < labels_in.size(); ++i) {
        anno_graph.call_annotated_nodes(labels_in[i], [&](node_index j) {
            #pragma omp critical
            {
                indicator[j] = true;
                ++counts[j];
            }
        });
    }

    // correct the width of counts, making it single-width
    counts.width(width);

    #pragma omp parallel for schedule(dynamic) num_threads(num_threads) default(shared)
    for (size_t i = 0; i < labels_out.size(); ++i) {
        anno_graph.call_annotated_nodes(labels_out[i], [&](node_index j) {
            #pragma omp critical
            {
                indicator[j] = true;
                ++counts[j * 2 + 1];
            }
        });
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
                for (size_t i = last; i < begin; ++i) {
                    callback(path[i], false);
                }
                last = end;
            }

            for (size_t i = last; i < path.size(); ++i) {
                callback(path[i], false);
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
