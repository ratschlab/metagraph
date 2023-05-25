#include "annotated_graph_algorithm.hpp"

#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "common/vector_map.hpp"
#include "common/vector_set.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator::Label Label;

typedef std::function<size_t()> LabelCountCallback;

constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

uint64_t atomic_fetch(const sdsl::int_vector<> &vector,
                      uint64_t i,
                      std::mutex &backup_mutex,
                      int mo) {
    const uint8_t width = vector.width();
    if (width == 64)
        return __atomic_load_n(vector.data() + i, mo);

    if (width == 32)
        return __atomic_load_n((uint32_t*)vector.data() + i, mo);

    if (width == 16)
        return __atomic_load_n((uint16_t*)vector.data() + i, mo);

    if (width == 8)
        return __atomic_load_n((uint8_t*)vector.data() + i, mo);

    std::lock_guard<std::mutex> lock(backup_mutex);
    return vector[i];
}

uint64_t atomic_fetch_and_add(sdsl::int_vector<> &vector,
                              uint64_t i,
                              uint64_t val,
                              std::mutex &backup_mutex,
                              int mo) {
    const uint8_t width = vector.width();
    if (width == 64)
        return __atomic_fetch_add(vector.data() + i, val, mo);

    if (width == 32)
        return __atomic_fetch_add((uint32_t*)vector.data() + i, val, mo);

    if (width == 16)
        return __atomic_fetch_add((uint16_t*)vector.data() + i, val, mo);

    if (width == 8)
        return __atomic_fetch_add((uint8_t*)vector.data() + i, val, mo);

    std::lock_guard<std::mutex> lock(backup_mutex);
    uint64_t old_val = vector[i];
    vector[i] += val;
    return old_val;
}


/**
 * Return an int_vector<>, bit_vector pair, of lengths anno_graph.get_graph().max_index() * 2
 * and anno_graph.get_graph().max_index(), respectively.
 * For an index i, the int_vector at indices 2*i and 2*i + 1 represent the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with, respectively. The width of the int_vector<> is computed to be
 * wide enough to contain counts up to num_labels.
 * The returned bit_vector is a k-mer mask indicating those k-mers annotated
 * with at least one in-label or out-label.
 */
std::pair<sdsl::int_vector<>, sdsl::bit_vector>
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
                          sdsl::int_vector<> &counts,
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
                    set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
            });
        };

        auto count_merge = [&](sdsl::bit_vector &a,
                               const sdsl::bit_vector &b,
                               const std::vector<node_index> &id_map,
                               size_t offset = 0) {
            call_ones(b, [&](size_t i) {
                if (id_map[i]) {
                    set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
                    atomic_fetch_and_add(count_vector.first, id_map[i] * 2 + offset, 1,
                                         vector_backup_mutex, MO_RELAXED);
                }
            });
        };

        logger->trace("Checking shared and other labels");
        masked_graph->call_sequences([&](const std::string &contig, const std::vector<node_index> &path) {
            for (const auto &[label, sig] : anno_graph.get_top_label_signatures(contig, num_labels)) {
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
                          sdsl::int_vector<> &counts,
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
        masked_graph->call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
            std::string rc_seq = seq;
            std::vector<DeBruijnGraph::node_index> rc_path = path;
            reverse_complement_seq_path(*graph_ptr, rc_seq, rc_path);

            auto it = rc_path.rbegin();
            for (size_t i = 0; i < path.size(); ++i, ++it) {
                if (*it) {
                    uint64_t in_count = atomic_fetch(counts, path[i] * 2, vector_backup_mutex, MO_RELAXED);
                    uint64_t out_count = atomic_fetch(counts, path[i] * 2 + 1, vector_backup_mutex, MO_RELAXED);
                    atomic_fetch_and_add(counts, *it * 2, in_count, vector_backup_mutex, MO_RELAXED);
                    atomic_fetch_and_add(counts, *it * 2 + 1, out_count, vector_backup_mutex, MO_RELAXED);

                    set_bit(mask.data(), *it, in_count, MO_RELAXED);
                }
            }
        }, num_threads, true);

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Updated masked graph has {} nodes", masked_graph->num_nodes());
    }

    return masked_graph;
}

std::pair<sdsl::int_vector<>, sdsl::bit_vector>
construct_diff_label_count_vector(const AnnotatedDBG &anno_graph,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  size_t num_labels,
                                  size_t num_threads) {
    size_t width = sdsl::bits::hi(num_labels) + 1;
    sdsl::bit_vector indicator(anno_graph.get_graph().max_index() + 1, false);

    // the in and out counts are stored interleaved
    sdsl::int_vector<> counts = aligned_int_vector(indicator.size() * 2, 0, width, 16);

    const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
    const auto &binmat = anno_graph.get_annotator().get_matrix();

    tsl::hopscotch_map<uint64_t, uint8_t> code_to_indicator;
    for (const std::string &label_in : labels_in) {
        code_to_indicator[label_encoder.encode(label_in)] = 1;
    }
    for (const std::string &label_out : labels_out) {
        code_to_indicator[label_encoder.encode(label_out)] |= 2;
    }

    std::vector<uint64_t> label_codes;
    label_codes.reserve(code_to_indicator.size());
    for (const auto &[code, indicator] : code_to_indicator) {
        label_codes.push_back(code);
    }

    std::mutex vector_backup_mutex;
    std::atomic_thread_fence(std::memory_order_release);
    bool parallel = num_threads > 1;

    binmat.call_columns(label_codes,
        [&](auto col_idx, const bitmap &rows) {
            uint8_t col_indicator = code_to_indicator[label_codes[col_idx]];
            rows.call_ones([&](auto r) {
                node_index i = AnnotatedDBG::anno_to_graph_index(r);
                set_bit(indicator.data(), i, parallel, MO_RELAXED);
                if (col_indicator & 1)
                    atomic_fetch_and_add(counts, i * 2, 1, vector_backup_mutex, MO_RELAXED);

                if (col_indicator & 2)
                    atomic_fetch_and_add(counts, i * 2 + 1, 1, vector_backup_mutex, MO_RELAXED);
            });
        },
        num_threads
    );

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
        total_unitigs.fetch_add(1, MO_RELAXED);

        size_t last = 0;
        for (const auto &pair : get_kept_intervals(unitig, path)) {
            const auto &[begin, end] = pair;
            kept_unitigs.fetch_add(1, MO_RELAXED);
            num_kept_nodes.fetch_add(end - begin, MO_RELAXED);
            for ( ; last < begin; ++last) {
                unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
            }
            last = end;
        }

        for ( ; last < path.size(); ++last) {
            unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
        }

    }, num_threads);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} unitigs with average length {}",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs);
}

void assemble_with_coordinates(size_t k,
                               const SequenceGenerator &input_generator,
                               const SequenceCoordCallback &callback) {
    DBGHashOrdered dbg(k);
    input_generator([&](const std::string &seq) { dbg.add_sequence(seq); });

    // assemble unitigs
    std::vector<std::pair<std::string, std::vector<node_index>>> unitigs;
    tsl::hopscotch_map<node_index, size_t> node_to_unitig;
    std::mutex mu;
    dbg.call_unitigs([&](const std::string &seq, const auto &path) {
        std::lock_guard<std::mutex> lock(mu);
        size_t unitig_id = unitigs.size();
        unitigs.emplace_back(seq, path);
        node_to_unitig[path.front()] = unitig_id;
        node_to_unitig[path.back()] = unitig_id;
    }, get_num_threads());

    size_t num_unitigs = unitigs.size();

    logger->trace("Assembled {} unitigs", num_unitigs);

    // enumerate superbubbles
    sdsl::bit_vector in_superbubble(num_unitigs, false);
    sdsl::bit_vector is_superbubble_start(num_unitigs, false);
    std::vector<std::pair<size_t, VectorMap<size_t, std::pair<std::vector<int64_t>, std::vector<int64_t>>>>> superbubbles(unitigs.size());
    for (size_t i = 0; i < unitigs.size(); ++i) {
        auto &[terminus, seens] = superbubbles[i];
        terminus = std::numeric_limits<size_t>::max();

        tsl::hopscotch_set<size_t> visited;
        std::vector<std::pair<size_t, size_t>> traversal_stack;
        traversal_stack.emplace_back(i, 0);
        seens[i].first.emplace_back(0);
        // std::cerr << "starting at\t" << unitigs[i].first << std::endl;

        while (traversal_stack.size()) {
            auto [unitig_id, dist] = traversal_stack.back();
            traversal_stack.pop_back();

            const auto &[unitig, unipath] = unitigs[unitig_id];
            // std::cerr << "\tvisiting\t" << unitig << std::endl;
            node_index back = unipath.back();
            size_t length = unipath.size();

            bool found_cycle = false;
            size_t num_children = 0;
            dbg.adjacent_outgoing_nodes(back, [&](node_index next) {
                if (found_cycle)
                    return;

                ++num_children;

                size_t next_id = node_to_unitig[next];
                // std::cerr << "\tseen\t" << unitigs[next_id].first << std::endl;
                if (next_id == i || unitig_id == next_id) {
                    // std::cerr << "\tfound cycle" << std::endl;
                    found_cycle = true;
                    traversal_stack.clear();
                    return;
                }

                seens[next_id].first.emplace_back(dist + length);
                bool all_visited = true;
                dbg.adjacent_incoming_nodes(next, [&](node_index sibling) {
                    size_t sibling_id = node_to_unitig[sibling];
                    if (sibling_id == unitig_id)
                        return;

                    if (all_visited && !visited.count(sibling_id))
                        all_visited = false;
                });

                if (all_visited)
                    traversal_stack.emplace_back(next_id, dist + length);
            });

            if (!num_children) {
                // std::cerr << "\tfound dead end" << std::endl;
                seens.clear();
                break;
            }

            visited.emplace(unitig_id);

            if (traversal_stack.size() == 1 && visited.size() + 1 == seens.size()) {
                auto [unitig_id, dist] = traversal_stack.back();
                traversal_stack.pop_back();
                terminus = unitig_id;
            }
        }

        if (terminus == std::numeric_limits<size_t>::max()) {
            // std::cerr << "\tno terminus found" << std::endl;
            seens.clear();
            continue;
        }

        dbg.adjacent_outgoing_nodes(unitigs[terminus].second.back(), [&](node_index next) {
            size_t next_id = node_to_unitig[next];
            if (next_id == i) {
                // std::cerr << "\tfound cycle" << std::endl;
                seens.clear();
            }
        });

        if (seens.empty())
            continue;

        is_superbubble_start[i] = true;

        // std::cerr << "\tfound superbubble" << std::endl;

        assert(seens.size() >= 4);

        for (auto it = seens.begin(); it != seens.end(); ++it) {
            auto &d = it.value().first;
            in_superbubble[it->first] = true;
            std::sort(d.begin(), d.end());
            d.erase(std::unique(d.begin(), d.end()), d.end());
        }

        // logger->info("Found superbubble with {} unitigs with widths {}", seens.size(),
        //     fmt::join(seens[terminus].first, ","));

        if (seens[terminus].first.size() == 1) {
            // easy case
            int64_t sb_size = seens[terminus].first[0] + 1;
            for (auto it = seens.begin(); it != seens.end(); ++it) {
                auto &[d, de] = it.value();
                assert(d.size());
                de.reserve(d.size());
                std::transform(d.rbegin(), d.rend(), std::back_inserter(de),
                               [&](int64_t dd) {
                                   assert(sb_size > dd);
                                   return dd - sb_size;
                               });
                assert(d.size());
            }
        } else {
            // hard case, multiple ways to get to the end
            std::vector<std::tuple<size_t, tsl::hopscotch_set<int64_t>, tsl::hopscotch_set<int64_t>>> traversal;
            traversal.emplace_back(terminus,
                                   tsl::hopscotch_set<int64_t>(seens[terminus].first.begin(),
                                                               seens[terminus].first.end()),
                                   tsl::hopscotch_set<int64_t>{ -1 });
            while (traversal.size()) {
                auto [u_id, d_begin, d_end] = std::move(traversal.back());
                traversal.pop_back();
                assert(d_end.size());
                assert(std::all_of(d_end.begin(), d_end.end(), [&](int64_t d) { return d < 0; }));
                seens[u_id].second = std::vector<int64_t>(d_end.begin(), d_end.end());
                seens[u_id].second.erase(std::unique(seens[u_id].second.begin(),
                                                     seens[u_id].second.end()),
                                         seens[u_id].second.end());

                if (u_id == i)
                    continue;

                node_index front = unitigs[u_id].second.front();
                dbg.adjacent_incoming_nodes(front, [&](node_index prev) {
                    size_t prev_id = node_to_unitig[prev];
                    size_t length = unitigs[prev_id].second.size();
                    auto &[p, d_begin_prev, d_end_prev] = traversal.emplace_back(prev_id, tsl::hopscotch_set<int64_t>{}, tsl::hopscotch_set<int64_t>{});
                    const auto &d = seens[prev_id].first;
                    for (int64_t dd : d) {
                        if (d_begin.count(dd + length)) {
                            d_begin_prev.emplace(dd);
                            for (int64_t dde : d_end) {
                                d_end_prev.emplace(dde - length);
                            }
                        }
                    }

                    if (d_begin_prev.empty())
                        traversal.pop_back();
                });
            }

            assert(std::all_of(seens.begin(), seens.end(), [&](const auto &a) {
                return a.second.first.size() && a.second.second.size()
                    && std::all_of(a.second.first.begin(), a.second.first.end(),
                        [&](int64_t c) { return c >= 0; })
                    && std::all_of(a.second.second.begin(), a.second.second.end(),
                        [&](int64_t c) { return c < 0; });
            }));
        }
    }

    size_t cur_unitig_id = 0;
    logger->info("Outputting {} unitigs not in a superbubble",
                  num_unitigs - sdsl::util::cnt_one_bits(in_superbubble));
    call_zeros(in_superbubble, [&](size_t i) {
        const auto &[unitig, path] = unitigs[i];
        ++cur_unitig_id;
        // logger->info("U:\t{}", unitig);
        callback(unitig, 0, 0, std::vector<int64_t>{}, std::vector<int64_t>{});
    });

    logger->info("Found {} superbubbles covering {} unitigs",
                  sdsl::util::cnt_one_bits(is_superbubble_start),
                  sdsl::util::cnt_one_bits(in_superbubble));

    // chain superbubbles
    sdsl::bit_vector chain_middle(num_unitigs, false);
    call_ones(is_superbubble_start, [&](size_t i) {
        const auto &[terminus, seens] = superbubbles[i];

        if (superbubbles[terminus].first != std::numeric_limits<size_t>::max())
            chain_middle[terminus] = true;
    });

    std::vector<std::vector<size_t>> chains;
    call_ones(is_superbubble_start, [&](size_t i) {
        if (chain_middle[i])
            return;

        // construct a chain
        VectorSet<size_t> chain;
        while (true) {
            const auto &[terminus, seens] = superbubbles[i];
            if (terminus == std::numeric_limits<size_t>::max() || !chain.emplace(i).second)
                break;

            i = terminus;
        }
        chains.emplace_back(const_cast<std::vector<size_t>&&>(chain.values_container()));
    });

    logger->info("Found {} chains", chains.size());
    for (size_t i = 0; i < chains.size(); ++i) {
        logger->info("Processing chain of size {}", chains[i].size());
        assert(chains[i].size());
        size_t term_id = cur_unitig_id + superbubbles[chains[i].front()].second.size();
        std::vector<tsl::hopscotch_set<int64_t>> distances_to_end(chains[i].size());
        distances_to_end.back().emplace(0);
        auto jt = distances_to_end.rbegin() + 1;
        for (auto it = chains[i].rbegin(); it + 1 != chains[i].rend(); ++it) {
            size_t jj = *it;
            const auto &[terminus, seens] = superbubbles[jj];
            term_id += seens.size() - 1;
            assert(jt != distances_to_end.rend());
            for (int64_t dd : seens.find(terminus)->second.first) {
                for (int64_t ww : *(jt - 1)) {
                    jt->emplace(dd + ww);
                }
            }
            ++jt;
        }
        assert(jt == distances_to_end.rend());

        tsl::hopscotch_set<int64_t> widths;
        widths.emplace(0);
        for (size_t jj : chains[i]) {
            --jt;
            const auto &[terminus, seens] = superbubbles[jj];
            size_t term_unitig_id = cur_unitig_id + seens.size();
            for (const auto &[j, d] : seens) {
                if (j == terminus) {
                    if (jt != distances_to_end.rbegin()) {
                        tsl::hopscotch_set<int64_t> new_widths;
                        for (int64_t ww : widths) {
                            for (int64_t dd : d.first) {
                                new_widths.emplace(ww + dd);
                            }
                        }
                        std::swap(widths, new_widths);
                    }
                    continue;
                }

                const auto &[unitig, path] = unitigs[j];
                ++cur_unitig_id;
                std::vector<int64_t> d_chain;
                d_chain.reserve(d.first.size() * widths.size());
                for (int64_t ww : widths) {
                    for (int64_t dd : d.first) {
                        d_chain.emplace_back(dd + ww);
                    }
                }
                std::sort(d_chain.begin(), d_chain.end());
                d_chain.erase(std::unique(d_chain.begin(), d_chain.end()), d_chain.end());

                std::vector<int64_t> de_chain;
                de_chain.reserve(d.second.size() * jt->size());
                for (int64_t ww : *jt) {
                    for (int64_t dd : d.second) {
                        de_chain.emplace_back(dd - ww);
                    }
                }
                std::sort(de_chain.begin(), de_chain.end());
                de_chain.erase(std::unique(de_chain.begin(), de_chain.end()), de_chain.end());
                // logger->info("U:\t{}\t{};{}\t{};{}", unitig, fmt::join(d.first, ","), fmt::join(d.second, ","), fmt::join(d_chain, ","), fmt::join(de_chain, ","));
                callback(unitig, term_unitig_id, term_id, d_chain, de_chain);
            }
        }

        const auto &[terminus, seens] = superbubbles[chains[i].back()];
        const auto &[unitig, path] = unitigs[terminus];
        const auto &d = seens.find(terminus)->second;
        ++cur_unitig_id;
        assert(cur_unitig_id == term_id);
        std::vector<int64_t> d_chain;
        d_chain.reserve(d.first.size() * widths.size());
        for (int64_t ww : widths) {
            for (int64_t dd : d.first) {
                d_chain.emplace_back(dd + ww);
            }
        }
        std::sort(d_chain.begin(), d_chain.end());
        d_chain.erase(std::unique(d_chain.begin(), d_chain.end()), d_chain.end());
        // logger->info("U:\t{}\t{};{}\t{};{}", unitig, fmt::join(d.first, ","), fmt::join(d.second, ","), fmt::join(d_chain, ","), fmt::join(d.second, ","));
        callback(unitig, cur_unitig_id, term_id, d_chain, d.second);
    }
}

std::string format_coords(const std::vector<int64_t> &dist_from_start,
                          const std::vector<int64_t> &dist_to_end) {
    return fmt::format("{};{}", fmt::join(dist_from_start, ","),
                                fmt::join(dist_to_end, ","));
}

std::string format_header(size_t superbubble_id,
                          size_t chain_id,
                          const std::vector<int64_t> &dist_from_start,
                          const std::vector<int64_t> &dist_to_end) {
    if (!superbubble_id)
        return "";

    if (!chain_id) {
        return fmt::format("{}\t{};{}",
                           superbubble_id,
                           format_coords(dist_from_start, dist_to_end));
    }

    return fmt::format("{};{}\t{}",
                       superbubble_id,
                       chain_id,
                       format_coords(dist_from_start, dist_to_end));
}

} // namespace graph
} // namespace mtg
