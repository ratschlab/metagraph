#include "annotated_graph_algorithm.hpp"

#include <atomic>
#include <type_traits>
#include <utility>
#include <variant>

#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/masked_graph.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator::Label Label;

typedef std::function<size_t()> LabelCountCallback;

// Byte-width counters have real typed elements for atomic_ref; other widths stay packed.
class LabelCounts {
  public:
    LabelCounts(uint64_t size, uint8_t width) {
        assert(width && width <= 64);
        switch (width) {
            case 8: storage_.emplace<std::vector<uint8_t>>(size); break;
            case 16: storage_.emplace<std::vector<uint16_t>>(size); break;
            case 32: storage_.emplace<std::vector<uint32_t>>(size); break;
            case 64: storage_.emplace<std::vector<uint64_t>>(size); break;
            default: storage_.emplace<sdsl::int_vector<>>(size, 0, width); break;
        }
    }

    uint64_t size() const {
        return std::visit([](const auto &v) -> uint64_t { return v.size(); }, storage_);
    }

    uint64_t operator[](uint64_t i) const {
        return std::visit([&](const auto &v) -> uint64_t { return v[i]; }, storage_);
    }

    // Select the storage type once per parallel phase, outside the inner loops.
    template <class Callback>
    void visit(Callback &&callback) {
        std::visit(std::forward<Callback>(callback), storage_);
    }

  private:
    std::variant<sdsl::int_vector<>, std::vector<uint8_t>, std::vector<uint16_t>,
                 std::vector<uint32_t>, std::vector<uint64_t>> storage_;
};

template <class T>
std::atomic_ref<T> counter_atomic(T &value) {
    static_assert(std::atomic_ref<T>::required_alignment <= alignof(T));
    return std::atomic_ref<T>(value);
}

template <class Vector>
uint64_t atomic_fetch(Vector &vector, uint64_t i, std::mutex &backup_mutex) {
    if constexpr (std::is_same_v<Vector, sdsl::int_vector<>>) {
        std::lock_guard<std::mutex> lock(backup_mutex);
        return vector[i];
    } else {
        return counter_atomic(vector[i]).load(std::memory_order_relaxed);
    }
}

template <class Vector>
void atomic_add(Vector &vector, uint64_t i, uint64_t value, std::mutex &backup_mutex) {
    if constexpr (std::is_same_v<Vector, sdsl::int_vector<>>) {
        std::lock_guard<std::mutex> lock(backup_mutex);
        vector[i] += value;
    } else {
        counter_atomic(vector[i]).fetch_add(value, std::memory_order_relaxed);
    }
}


/**
 * Return a count-vector, bit-vector pair, of lengths anno_graph.get_graph().max_index() * 2
 * and anno_graph.get_graph().max_index(), respectively.
 * For an index i, the count-vector at indices 2*i and 2*i + 1 represent the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with, respectively. The count width is computed to be
 * wide enough to contain counts up to num_labels.
 * The returned bit_vector is a k-mer mask indicating those k-mers annotated
 * with at least one in-label or out-label.
 */
std::pair<LabelCounts, sdsl::bit_vector>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  size_t num_labels,
                                  size_t num_threads);

// Regions of a graph mask which should be kept (i.e., masked in)
typedef std::vector<std::pair<size_t, size_t>> Intervals;

// Returns a vector of kept regions given a unitig and its corresponding path
typedef std::function<Intervals(const std::string&, const std::vector<node_index>&)> GetKeptIntervals;

// Assemble unitigs from the masked graph, then use get_kept_intervals to generate
// regions which should be masked in. Update the graph mask accordingly.
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads);

// Given an initial mask and counts, generate a masked graph. If add_complement
// is true, then add the reverse complements of all nodes to the graph as well.
std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          LabelCounts &counts,
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads);


std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads) {
    bool parallel = num_threads > 1;
    size_t num_labels = anno_graph.get_annotator().num_labels();
    size_t num_in_labels = labels_in.size() + labels_in_round2.size();
    size_t num_out_labels = labels_out.size() + labels_out_round2.size();
    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    logger->trace("Generating initial mask");

    // Construct initial masked graph from union of labels in labels_in
    auto count_vector = construct_diff_label_count_vector(
        anno_graph, labels_in, labels_out,
        std::max(num_in_labels, num_out_labels),
        num_threads
    );
    auto &[counts, init_mask] = count_vector;

    // in and out counts are stored interleaved in the counts vector
    assert(counts.size() == init_mask.size() * 2);

    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    sdsl::bit_vector union_mask;

    sdsl::bit_vector other_mask(init_mask.size() * check_other, false);
    auto masked_graph = make_initial_masked_graph(graph_ptr, counts, std::move(init_mask),
                                                  config.add_complement, num_threads);

    if (check_other || labels_in_round2.size() || labels_out_round2.size())
        union_mask = dynamic_cast<const bitmap_vector&>(masked_graph->get_mask()).data();

    // check all other labels and post labels
    if (check_other || labels_in_round2.size() || labels_out_round2.size()) {
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);

        auto mask_or = [&](sdsl::bit_vector &a,
                           const sdsl::bit_vector &b,
                           const std::vector<node_index> &id_map) {
            call_ones(b, [&](size_t i) {
                if (id_map[i])
                    set_bit(a.data(), id_map[i], parallel, std::memory_order_relaxed);
            });
        };

        counts.visit([&](auto &values) {
            auto count_merge = [&](sdsl::bit_vector &a,
                                   const sdsl::bit_vector &b,
                                   const std::vector<node_index> &id_map,
                                   size_t offset = 0) {
                call_ones(b, [&](size_t i) {
                    if (id_map[i]) {
                        set_bit(a.data(), id_map[i], parallel, std::memory_order_relaxed);
                        atomic_add(values, id_map[i] * 2 + offset, 1, vector_backup_mutex);
                    }
                });
            };

            logger->trace("Checking shared and other labels");
            masked_graph->call_sequences([&](const std::string &contig, const std::vector<node_index> &path) {
                for (const auto &[label, count, sig] : anno_graph.get_top_label_signatures(contig, num_labels)) {
                    bool found_in = labels_in.count(label);
                    bool found_out = labels_out.count(label);
                    bool found_in_round2 = labels_in_round2.count(label);
                    bool found_out_round2 = labels_out_round2.count(label);
                    if (!found_in && !found_out
                            && !found_in_round2 && !found_out_round2 && check_other) {
                        mask_or(other_mask, sig, path);
                    }

                    if (found_in_round2)
                        count_merge(union_mask, sig, path);

                    if (found_out_round2)
                        count_merge(union_mask, sig, path, 1);
                }
            }, num_threads);
        });

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(union_mask)));
    }

    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");

    size_t min_label_in_count = std::ceil(config.label_mask_in_kmer_fraction * num_in_labels);
    size_t max_label_out_count = std::floor(config.label_mask_out_kmer_fraction * num_out_labels);

    if (config.label_mask_in_unitig_fraction == 0.0
            && config.label_mask_out_unitig_fraction == 1.0
            && config.label_mask_other_unitig_fraction == 1.0) {
        logger->trace("Filtering by node");
        size_t total_nodes = masked_graph->num_nodes();
        const auto &in_mask = dynamic_cast<const bitmap_vector&>(masked_graph->get_mask()).data();
        sdsl::bit_vector mask = in_mask;

        // TODO: make this part multithreaded
        size_t kept_nodes = 0;
        call_ones(in_mask, [&](node_index node) {
            uint64_t in_count = count_vector.first[node * 2];
            uint64_t out_count = count_vector.first[node * 2 + 1];

            if (in_count >= min_label_in_count && out_count <= max_label_out_count) {
                ++kept_nodes;
            } else {
                mask[node] = false;
            }
        });

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);

        return masked_graph;
    }

    logger->trace("Filtering by unitig");

    update_masked_graph_by_unitig(*masked_graph,
        [&](const std::string &, const std::vector<node_index> &path) -> Intervals {
            // return a set of intervals to keep in the graph
            const auto &counts = count_vector.first;

            size_t in_kmer_count = 0;

            size_t begin = path.size();
            size_t end = 0;
            for (size_t i = 0; i < path.size(); ++i) {
                if (counts[path[i] * 2] >= min_label_in_count) {
                    if (begin == path.size())
                        begin = i;

                    end = std::max(end, i + 1);

                    ++in_kmer_count;
                }
            }

            if (begin >= end)
                return {};

            size_t size = end - begin;
            size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * size);
            if (in_kmer_count < label_in_cutoff)
                return {};

            size_t out_kmer_count = 0;
            size_t other_kmer_count = 0;
            size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * size);
            size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * size);


            for (size_t i = begin; i < end; ++i) {
                if (counts[path[i] * 2 + 1] > max_label_out_count
                        && ++out_kmer_count > label_out_cutoff) {
                    return {};
                }

                if (check_other && other_mask[path[i]]
                        && ++other_kmer_count > other_cutoff) {
                    return {};
                }
            }

            return { std::make_pair(begin, end) };
        },
        num_threads
    );

    return masked_graph;
}


/**
 * Helpers
 */

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          LabelCounts &counts,
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads) {
    // counts is a double-length vector storing the in-label and out-label
    // counts interleaved
    assert(counts.size() == mask.size() * 2);

    add_complement |= graph_ptr->get_mode() == DeBruijnGraph::CANONICAL;
    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        add_complement
            ? std::make_unique<bitmap_vector>(mask)
            : std::make_unique<bitmap_vector>(std::move(mask)),
        true
    );

    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    if (add_complement) {
        logger->trace("Adding reverse complements");
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);
        counts.visit([&](auto &values) {
            masked_graph->call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
                std::string rc_seq = seq;
                std::vector<DeBruijnGraph::node_index> rc_path = path;
                reverse_complement_seq_path(*graph_ptr, rc_seq, rc_path);

                auto it = rc_path.rbegin();
                for (size_t i = 0; i < path.size(); ++i, ++it) {
                    if (*it) {
                        uint64_t in_count = atomic_fetch(values, path[i] * 2, vector_backup_mutex);
                        uint64_t out_count = atomic_fetch(values, path[i] * 2 + 1, vector_backup_mutex);
                        atomic_add(values, *it * 2, in_count, vector_backup_mutex);
                        atomic_add(values, *it * 2 + 1, out_count, vector_backup_mutex);

                        set_bit(mask.data(), *it, num_threads > 1, std::memory_order_relaxed);
                    }
                }
            }, num_threads, true);
        });

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Updated masked graph has {} nodes", masked_graph->num_nodes());
    }

    return masked_graph;
}

std::pair<LabelCounts, sdsl::bit_vector>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  size_t num_labels,
                                  size_t num_threads) {
    size_t width = sdsl::bits::hi(num_labels) + 1;
    sdsl::bit_vector indicator(anno_graph.get_graph().max_index() + 1, false);

    // the in and out counts are stored interleaved
    LabelCounts counts(indicator.size() * 2, width);

    const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
    const auto &binmat = anno_graph.get_annotator().get_matrix();

    tsl::hopscotch_map<uint64_t, uint8_t> code_to_indicator;
    for (const std::string &label_in : labels_in) {
        code_to_indicator[label_encoder.encode(label_in)] = 1;
    }
    for (const std::string &label_out : labels_out) {
        code_to_indicator[label_encoder.encode(label_out)] |= 2;
    }

    std::vector<annot::matrix::BinaryMatrix::Column> label_codes;
    label_codes.reserve(code_to_indicator.size());
    for (const auto &[code, indicator] : code_to_indicator) {
        label_codes.push_back(code);
    }

    std::mutex vector_backup_mutex;
    std::atomic_thread_fence(std::memory_order_release);
    bool parallel = num_threads > 1;

    counts.visit([&](auto &values) {
        binmat.call_columns(label_codes,
            [&](auto col_idx, const bitmap &rows) {
                uint8_t col_indicator = code_to_indicator[label_codes[col_idx]];
                rows.call_ones([&](auto r) {
                    node_index i = AnnotatedDBG::anno_to_graph_index(r);
                    set_bit(indicator.data(), i, parallel, std::memory_order_relaxed);
                    if (col_indicator & 1)
                        atomic_add(values, i * 2, 1, vector_backup_mutex);

                    if (col_indicator & 2)
                        atomic_add(values, i * 2 + 1, 1, vector_backup_mutex);
                });
            },
            num_threads
        );
    });

    std::atomic_thread_fence(std::memory_order_acquire);

    return std::make_pair(std::move(counts), std::move(indicator));
}


void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads) {
    std::atomic<uint64_t> kept_unitigs(0);
    std::atomic<uint64_t> total_unitigs(0);
    std::atomic<uint64_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;

    sdsl::bit_vector mask = dynamic_cast<const bitmap_vector&>(masked_graph.get_mask()).data();

    std::atomic_thread_fence(std::memory_order_release);

    masked_graph.call_unitigs([&](const std::string &unitig, const std::vector<node_index> &path) {
        total_unitigs.fetch_add(1, std::memory_order_relaxed);

        size_t last = 0;
        for (const auto &pair : get_kept_intervals(unitig, path)) {
            const auto &[begin, end] = pair;
            kept_unitigs.fetch_add(1, std::memory_order_relaxed);
            num_kept_nodes.fetch_add(end - begin, std::memory_order_relaxed);
            for ( ; last < begin; ++last) {
                unset_bit(mask.data(), path[last], parallel, std::memory_order_relaxed);
            }
            last = end;
        }

        for ( ; last < path.size(); ++last) {
            unset_bit(mask.data(), path[last], parallel, std::memory_order_relaxed);
        }

    }, num_threads);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs);
}

} // namespace graph
} // namespace mtg
