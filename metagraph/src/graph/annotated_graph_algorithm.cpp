#include "annotated_graph_algorithm.hpp"

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

SuperbubbleInfo assemble_with_coordinates(const SequenceGenerator &input_generator,
                                          const SequenceCoordCallback &callback) {
    DBGHashOrdered dbg;
    input_generator([&](const std::string &seq) { dbg.add_sequence(seq); });

    // assemble unitigs
    std::vector<size_t> unitig_ids(dbg.max_index() + 1);
    std::vector<std::string> unitigs;
    std::vector<std::tuple<size_t, node_index, node_index>> unitig_info;
    unitig_info.emplace_back();
    unitigs.emplace_back();
    std::mutex mu;
    dbg.call_unitigs([&](const std::string &seq, const auto &path) {
        size_t unitig_id = 0;
        {
            std::lock_guard<std::mutex> lock(mu);
            unitig_id = unitig_info.size();
            unitig_info.emplace_back(path.size(), path.front(), path.back());
            unitigs.emplace_back(seq);
        }

        for (node_index node : path) {
            unitig_ids[node] = unitig_id;
        }
    }, get_num_threads());

    size_t num_unitigs = unitig_info.size();

    // detect superbubbles
    SuperbubbleInfo ret_val(num_unitigs);

    std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_starts(num_unitigs);
    std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_termini(num_unitigs);
    std::vector<tsl::hopscotch_set<size_t>> superbubble_ends(num_unitigs);
    sdsl::bit_vector not_in_superbubble(num_unitigs, false);

    for (size_t i = 1; i < num_unitigs; ++i) {
        if (not_in_superbubble[i])
            continue;

        auto [i_length, i_front, i_back] = unitig_info[i];
        tsl::hopscotch_set<size_t> visited;
        VectorMap<size_t, tsl::hopscotch_set<size_t>> seen;
        tsl::hopscotch_map<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>> parents_children;
        std::vector<std::pair<size_t, size_t>> traversal_stack;
        traversal_stack.emplace_back(i, 0);
        seen[i].emplace(0);
        bool is_terminal_superbubble = false;
        size_t terminus = 0;
        size_t term_dist = 0;
        while (traversal_stack.size()) {
            auto [unitig_id, dist] = traversal_stack.back();
            traversal_stack.pop_back();
            assert(!visited.count(unitig_id));

            if (not_in_superbubble[unitig_id]) {
                is_terminal_superbubble = false;
                break;
            }

            bool has_cycle = false;
            visited.insert(unitig_id);
            size_t num_children = 0;
            auto [length, front, back] = unitig_info[unitig_id];
            auto &children = parents_children[unitig_id].second;
            dbg.call_outgoing_kmers(back, [&](node_index next, char) {
                if (next == front)
                    return;

                ++num_children;

                if (has_cycle)
                    return;

                if (next == i_front) {
                    has_cycle = true;
                    return;
                }

                size_t next_id = unitig_ids[next];
                children.emplace_back(next_id);

                bool add_parents = !seen.count(next_id);

                seen[next_id].emplace(dist + length);
                bool all_visited = true;
                dbg.call_incoming_kmers(next, [&](node_index sibling, char) {
                    size_t sibling_id = unitig_ids[sibling];
                    if (sibling_id == next_id)
                        return;

                    if (add_parents)
                        parents_children[next_id].first.emplace_back(sibling_id);

                    if (all_visited && !visited.count(sibling_id))
                        all_visited = false;
                });

                if (all_visited)
                    traversal_stack.emplace_back(next_id, dist + length);
            });

            if (has_cycle) {
                for (const auto &[u_id, d] : seen) {
                    not_in_superbubble[u_id] = true;
                }
                is_terminal_superbubble = false;
                break;
            }

            if (!num_children) {
                is_terminal_superbubble = true;
            }

            bool reached_end = (traversal_stack.size() == 1 && visited.size() + 1 == seen.size());
            if (reached_end) {
                auto [unitig_id, dist] = traversal_stack.back();
                traversal_stack.pop_back();

                terminus = unitig_id;
                term_dist = dist;
                for (const auto &[u_id, d] : seen) {
                    for (size_t dist : d) {
                        if (dist > term_dist) {
                            terminus = u_id;
                            term_dist = dist;
                        }
                    }
                }

                std::get<1>(ret_val[i]) = terminus;

                const auto &d = seen[terminus];
                superbubble_termini[i] = std::make_pair(terminus,
                    std::vector<size_t>(d.begin(), d.end()));
                std::sort(superbubble_termini[i].second.begin(),
                          superbubble_termini[i].second.end());
                std::get<2>(ret_val[i]) = superbubble_termini[i].second;

                for (const auto &[u_id, d] : seen) {
                    if (u_id == i)
                        continue;

                    std::vector<size_t> cur_d(d.begin(), d.end());
                    std::sort(cur_d.begin(), cur_d.end());

                    if (!superbubble_starts[u_id].first
                            || cur_d.back() < superbubble_starts[u_id].second.back()) {
                        superbubble_starts[u_id] = std::make_pair(i, std::move(cur_d));
                    }
                }
            }
        }

        if (reached_end) {
            std::vector<std::tuple<size_t, size_t, std::vector<size_t>>> back_traversal_stack;
            back_traversal_stack.reserve(seen.size());
            back_traversal_stack.emplace_back(terminus, 0,
                std::vector<size_t>(seen[terminus].begin(), seen[terminus].end()));
            while (back_traversal_stack.size()) {
                auto [cur_id, df, d] = back_traversal_stack.back();
                back_traversal_stack.pop_back();

                if (superbubble_starts[cur_id] == i) {
                    bool inserted = superbubble_ends[cur_id].emplace(df).second;
                    if (!inserted)
                        continue;
                }

                if (parents[cur_id].empty())
                    continue;

                for (auto &dd : d) {
                    --dd;
                }

                for (size_t parent : parents[cur_id]) {
                    auto &[p, p_df, p_d] = back_traversal_stack.emplace_back(parent, df + 1, std::vector<size_t>{});
                    for (auto dd : d) {
                        if (seen[parent].count(dd))
                            p_d.emplace_back(dd);
                    }

                    if (p_d.empty())
                        back_traversal_stack.pop_back();
                }
            }
        }
    }

    for (size_t i = 1; i < num_unitigs; ++i) {
        const auto &[superbubble_id, sb_d] = superbubble_starts[i];
        const auto &[ends] = superbubble_ends[i];
        std::vector<size_t> ends_v(ends.begin(), ends.end());
        std::sort(ends_v.begin(), ends_v.end());
        callback(superbubble_id, i, sb_d, ends_v);
    }

    // chain superbubbles
    tsl::hopscotch_map<size_t, std::vector<size_t>> superbubble_storage;

    sdsl::int_vector<> chain_parent(num_unitigs);
    for (size_t i = 1; i < num_unitigs; ++i) {
        size_t t = superbubble_termini[i].first;
        if (!t) {
            if (size_t sb = superbubble_starts[i].first)
                superbubble_storage[sb].emplace(i);
        } else {
            // this is a superbubble start
            superbubble_storage[i].emplace(i);
            if (size_t sb = superbubble_starts[i].first)
                superbubble_storage[sb].emplace(i);
        }

        if (!t || t == i)
            continue;

        // we are at a superbubble start
        size_t t2 = superbubble_termini[t].first;
        if (!t2)
            continue;

        // the superbubble starting at i is chained to the one starting at t
        assert(!chain_parent[t2]);
        chain_parent[t2] = t;
        chain_parent[t] = i;
    }

    // sdsl::int_vector<> superbubble_chain(num_unitigs * 2);
    // std::vector<bool> chain_bounds;
    // std::vector<uint64_t> chain_widths;
    // std::vector<size_t> chain_starts;
    size_t chain_i = 1;
    // size_t chain_widths_id = 0;
    for (size_t i = 1; i < num_unitigs; ++i) {
        if (chain_parent[i] && !superbubble_termini[i].first) {
            // assert(!superbubble_chain[i * 2]);

            // end of a superbubble chain
            size_t j = i;
            // size_t chain_term = i;
            // superbubble_chain[j * 2] = chain_term;
            // superbubble_chain[j * 2 + 1] = ++chain_widths_id;
            // chain_widths.emplace_back(0);
            // chain_bounds.emplace_back(true);

            tsl::hopscotch_set<size_t> widths;
            widths.emplace(0);
            while (chain_parent[j]) {
                std::get<0>(ret_val[j]) = chain_i;
                std::get<3>(ret_val[j]) = std::vector<int64_t>(widths.begin(), widths.end());
                std::sort(std::get<3>(ret_val[j]).begin(), std::get<3>(ret_val[j]).end());
                size_t old_j = j;
                std::ignore = old_j;
                j = chain_parent[j];
                size_t widths_start = chain_widths.size();
                // superbubble_chain[j - 1 * 2] = chain_term;
                // superbubble_chain[j - 1 * 2 + 1] = ++chain_widths_id;
                const auto &[t, d] = superbubble_termini[j];
                assert(t == old_j);
                tsl::hopscotch_set<size_t> next_widths;
                for (size_t width : widths) {
                    std::for_each(d.begin(), d.end(), [&](size_t dd) {
                        next_widths.emplace(dd + width);
                    });
                }
                std::swap(next_widths, widths);
                assert(widths.size());
                // for (size_t width : widths) {
                    // chain_widths.emplace_back(width);
                    // chain_bounds.emplace_back(false);
                // }
                std::sort(chain_widths.end() - widths.size(), chain_widths.end());
                // chain_bounds[widths_start] = true;
            }

            std::get<0>(ret_val[j]) = chain_i;
            std::get<3>(ret_val[j]) = std::vector<int64_t>(widths.begin(), widths.end());
            std::sort(std::get<3>(ret_val[j]).begin(), std::get<3>(ret_val[j]).end());
            // chain_starts.emplace_back(j);

            ++chain_i;
        }
    }

    return ret_val;
}

} // namespace graph
} // namespace mtg
