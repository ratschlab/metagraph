#include "annotated_graph_algorithm.hpp"

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/vector_map.hpp"


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

void assemble_min_path_cover(const StringGenerator &generate,
                             size_t k,
                             const std::function<void(const std::string&)> &callback,
                             DeBruijnGraph::Mode mode,
                             size_t num_threads) {
    // assemble graph
    if (mode == DeBruijnGraph::PRIMARY)
        mode = DeBruijnGraph::CANONICAL;

    DBGHashFast dbg(k, mode);
    generate([&dbg](std::string_view seq) { dbg.add_sequence(seq); });
    assemble_min_path_cover(dbg, callback, num_threads);
}

using unitig_index = node_index;
template <class T>
void traverse(node_index v,
              const T &adjacent_outgoing_unitigs,
              std::vector<size_t> &component,
              std::vector<size_t> &preorder,
              std::vector<std::vector<unitig_index>> &sccs,
              size_t &C,
              std::vector<unitig_index> &S,
              std::vector<unitig_index> &P) {
    preorder[v] = C++;
    S.push_back(v);
    P.push_back(v);
    adjacent_outgoing_unitigs(v, [&](unitig_index w) {
        if (preorder[w] == component.size()) {
            traverse(w, adjacent_outgoing_unitigs, component, preorder, sccs, C, S, P);
        } else if (component[w] == component.size()) {
            assert(P.size());
            while (preorder[P.back()] > preorder[w]) {
                P.pop_back();
                assert(P.size());
            }
        }
    });

    assert(P.size());
    if (P.back() == v) {
        assert(S.size());
        size_t cur_component = sccs.size();
        std::vector<unitig_index> &scc = sccs.emplace_back();
        while (S.back() != v) {
            component[S.back()] = cur_component;
            scc.emplace_back(S.back());
            S.pop_back();
            assert(S.size());
        }
        component[v] = cur_component;
        scc.emplace_back(v);
        S.pop_back();
        P.pop_back();
    }
}

void assemble_min_path_cover(const DeBruijnGraph &dbg,
                             const std::function<void(const std::string&)> &callback,
                             size_t num_threads) {
    // assemble unitigs
    std::vector<std::tuple<node_index, node_index, size_t, std::string>> unitigs;
    std::atomic<unitig_index> num_unitigs { 0 };
    std::vector<unitig_index> node_to_unitig(dbg.max_index() + 1);

    std::mutex mu;

    dbg.call_unitigs([&](const std::string &seq, const auto &path) {
        unitig_index unitig_id = num_unitigs.fetch_add(1, std::memory_order_relaxed);
        std::lock_guard<std::mutex> lock(mu);
        node_to_unitig[path.front()] = unitig_id;
        node_to_unitig[path.back()] = unitig_id;
        unitigs.emplace_back(path.front(), path.back(), path.size(), seq);
    }, num_threads);

    logger->trace("Topologically sorting {} unitigs", unitigs.size());

    auto adjacent_incoming_unitigs = [&](unitig_index i,
                                         const auto &callback) {
        const auto &[front, back, size, seq] = unitigs[i];
        dbg.call_incoming_kmers(front, [&](node_index prev, char c) {
            if (c != boss::BOSS::kSentinel)
                callback(node_to_unitig[prev]);
        });
    };

    auto adjacent_outgoing_unitigs = [&](unitig_index i,
                                         const auto &callback) {
        const auto &[front, back, size, seq] = unitigs[i];
        dbg.call_outgoing_kmers(back, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                callback(node_to_unitig[next]);
        });
    };

    // find strongly-connected components
    std::vector<size_t> component(unitigs.size(), unitigs.size());
    std::vector<std::vector<unitig_index>> sccs;
    {
        std::vector<size_t> preorder(unitigs.size(), unitigs.size());
        std::vector<unitig_index> S;
        std::vector<unitig_index> P;
        size_t C = 0;
        for (size_t i = 0; i < unitigs.size(); ++i) {
            if (component[i] == unitigs.size())
                traverse(i, adjacent_outgoing_unitigs, component, preorder, sccs, C, S, P);
        }
    }

    logger->trace("Found {} connected components", sccs.size());

    logger->trace("Discarding back edges");
    tsl::hopscotch_map<unitig_index, tsl::hopscotch_set<unitig_index>> discarded_edges;
    size_t num_discarded = 0;
    for (const auto &scc : sccs) {
        std::vector<unitig_index> stack;
        tsl::hopscotch_set<unitig_index> visited;
        stack.emplace_back(scc[0]);
        while (stack.size()) {
            unitig_index v = stack.back();
            stack.pop_back();
            visited.emplace(v);
            adjacent_outgoing_unitigs(v, [&](unitig_index w) {
                if (visited.count(w)) {
                    discarded_edges[v].emplace(w);
                    ++num_discarded;
                } else {
                    stack.emplace_back(w);
                }
            });
        }
    }

    logger->trace("Discarded {} edges", num_discarded);

    auto adjacent_outgoing_unitigs_dag = [&](unitig_index i, const auto &callback) {
        auto find = discarded_edges.find(i);
        if (find == discarded_edges.end()) {
            adjacent_outgoing_unitigs(i, callback);
        } else {
            adjacent_outgoing_unitigs(i, [&](unitig_index j) {
                if (!find->second.count(j))
                    callback(j);
            });
        }
    };

    auto adjacent_incoming_unitigs_dag = [&](unitig_index j, const auto &callback) {
        adjacent_incoming_unitigs(j, [&](unitig_index i) {
            auto find = discarded_edges.find(i);
            if (find == discarded_edges.end() || !find->second.count(j))
                callback(i);
        });
    };

    logger->trace("Topologically sorting unitigs");
    std::vector<unitig_index> s;

    for (unitig_index i = 0; i < unitigs.size(); ++i) {
        size_t num_incoming = 0;
        adjacent_incoming_unitigs_dag(i, [&](unitig_index) { ++num_incoming; });
        if (!num_incoming)
            s.emplace_back(i);
    }

    VectorMap<unitig_index, tsl::hopscotch_map<unitig_index, bool>> seen;
    while (s.size()) {
        unitig_index n = s.back();
        s.pop_back();

        auto &seen_bucket = seen[n];

        adjacent_outgoing_unitigs_dag(n, [&](unitig_index next) {
            seen_bucket.emplace(next, false);
            bool all_visited = true;

            adjacent_incoming_unitigs_dag(next, [&](unitig_index sibling) {
                if (!all_visited)
                    return;

                auto find = seen.find(sibling);
                // TODO: simplify this
                if (find != seen.end()) {
                    if (!find->second.count(next)) {
                        all_visited = false;
                    }
                } else {
                    all_visited = false;
                }
            });

            if (all_visited)
                s.emplace_back(next);
        });
    }

    if (seen.size() != unitigs.size()) {
        logger->error("Found a cycle, since {} != {}", seen.size(), unitigs.size());
        exit(1);
    }

    logger->trace("Finding min path cover");

    sdsl::bit_vector m(unitigs.size(), false);

    size_t npaths = 0;

#ifndef NDEBUG
    sdsl::bit_vector output(unitigs.size(), false);
#endif

    auto make_paths = [&](unitig_index start_node, auto &edges) {
        if (m[start_node])
            return false;

        unitig_index n = start_node;
        std::string spelling = std::get<3>(unitigs[n]);
#ifndef NDEBUG
        output[n] = true;
#endif
        while (n != m.size()) {
            size_t old_n = n;
            auto find = edges.find(n);
            n = m.size();
            bool found = false;
            if (find != edges.end()) {
                for (auto it = find.value().begin(); it != find.value().end(); ++it) {
                    if (!it.value()) {
                        found = true;
                        it.value() = true;
                        n = it->first;

                        spelling += std::get<3>(unitigs[n]).substr(dbg.get_k() - 1);
#ifndef NDEBUG
                        output[n] = true;
#endif

                        break;
                    }
                }
            }

            if (!found)
                m[old_n] = true;
        }

        callback(spelling);
        ++npaths;
        return true;
    };

    logger->trace("Outputting paths");
    bool added = true;
    while (added) {
        added = false;
        for (const auto &[start_node, edges] : seen) {
            added |= make_paths(start_node, seen);
        }
    }

    logger->trace("Found {} paths", npaths);

    logger->trace("Outputting discarded edges as extra paths");
    sdsl::util::set_to_value(m, false);
    VectorMap<unitig_index, tsl::hopscotch_map<unitig_index, bool>> rest_edges;
    for (const auto &[j, edges] : discarded_edges) {
        size_t num_incoming = 0;
        adjacent_incoming_unitigs(j, [&](unitig_index i) {
            auto find = discarded_edges.find(i);
            if (find != discarded_edges.end() && find->second.count(j))
                ++num_incoming;
        });

        if (!num_incoming) {
            for (unitig_index k : edges) {
                rest_edges[j].emplace(k, false);
            }
        }
    }

    for (const auto &[j, edges] : discarded_edges) {
        for (unitig_index k : edges) {
            rest_edges[j].try_emplace(k, false);
        }
    }

    added = true;
    while (added) {
        added = false;
        for (const auto &[start_node, edges] : rest_edges) {
            added |= make_paths(start_node, rest_edges);
        }
    }

    assert(std::find(output.begin(), output.end(), false) == output.end());

    for (const auto &[i, edges] : seen) {
        for (const auto &[j, used] : edges) {
            if (!used) {
                logger->error("Not all DAG edges covered");
                exit(1);
            }
        }
    }

    for (const auto &[i, edges] : rest_edges) {
        for (const auto &[j, used] : edges) {
            if (!used) {
                logger->error("Not all back edges covered");
                exit(1);
            }
        }
    }

    logger->trace("Found {} paths in total", npaths);
}

void assemble_superbubbles(const StringGenerator &generate,
                           size_t k,
                           const std::function<void(const std::string&, size_t)> &callback,
                           DeBruijnGraph::Mode mode,
                           size_t num_threads) {
    // assemble graph
    if (mode == DeBruijnGraph::PRIMARY)
        mode = DeBruijnGraph::CANONICAL;

    DBGHashFast dbg(k, mode);
    generate([&dbg](std::string_view seq) { dbg.add_sequence(seq); });
    assemble_superbubbles(dbg, callback, num_threads);
}

void assemble_superbubbles(const DeBruijnGraph &dbg,
                           const std::function<void(const std::string&, size_t)> &callback,
                           size_t num_threads,
                           bool ignore_nested_superbubbles) {
    using unitig_index = node_index;
    using superbubble_index = node_index;

    // assemble unitigs
    std::vector<std::tuple<node_index, node_index, size_t, std::string>> unitigs;
    std::atomic<unitig_index> num_unitigs { 0 };
    std::vector<unitig_index> node_to_unitig(dbg.max_index() + 1);
    std::mutex mu;

    dbg.call_unitigs([&](const std::string &seq, const auto &path) {
        unitig_index unitig_id = num_unitigs.fetch_add(1, std::memory_order_relaxed) + 1;
        for (node_index node : path) {
            node_to_unitig[node] = unitig_id;
        }

        std::lock_guard<std::mutex> lock(mu);
        unitigs.emplace_back(path.front(), path.back(), path.size(), seq);
    }, num_threads);

#ifndef NDEBUG
    sdsl::bit_vector output_unitig(unitigs.size(), false);
#endif

    logger->trace("Assembled {} unitigs", unitigs.size());

    // assemble superbubbles
    superbubble_index num_superbubbles = 0;
    std::vector<std::tuple<superbubble_index, superbubble_index, superbubble_index, uint64_t, uint64_t>> unitig_to_superbubble(unitigs.size() + 1);
    std::vector<std::pair<std::vector<unitig_index>, bool>> superbubbles;

    {
        ProgressBar progress_bar(unitigs.size(), "Finding superbubbles",
                                std::cerr, !common::get_verbose());

        #pragma omp parallel for num_threads(num_threads)
        for (size_t unitig_id = 1; unitig_id <= unitigs.size(); ++unitig_id) {
            ++progress_bar;

            // If this node is in the middle of a superbubble, then any superbubble
            // found here is a nested one -> ignore this unitig.
            if (ignore_nested_superbubbles && std::get<1>(unitig_to_superbubble[unitig_id]))
                continue;

            const auto &[front, back, size, seq] = unitigs[unitig_id - 1];

            if (dbg.outdegree(back) <= 1)
                continue;

            // start superbubble search
            bool has_tips = false;
            bool found_cycle = false;
            std::vector<node_index> S { unitig_id };

            std::vector<unitig_index> visited;
            tsl::hopscotch_map<unitig_index, bool> seen;

            while (S.size()) {
                unitig_index unitig = S.back();
                S.pop_back();

                visited.emplace_back(unitig);
                seen[unitig] = true;

                node_index node = std::get<1>(unitigs[unitig - 1]);

                bool has_children = false;
                dbg.call_outgoing_kmers(node, [&](node_index next_front, char c) {
                    if (c == boss::BOSS::kSentinel)
                        return;

                    has_children = true;
                    if (found_cycle)
                        return;

                    bool all_visited = true;
                    unitig_index next_unitig = node_to_unitig[next_front];

                    if (next_unitig == unitig_id) {
                        // cycle detected
                        found_cycle = true;
                        return;
                    }

                    auto [find, inserted] = seen.try_emplace(next_unitig, false);
                    if (find->second) {
                        visited.emplace_back(next_unitig);
                    } else {
                        dbg.call_incoming_kmers(next_front, [&](node_index sibling_back, char c) {
                            if (c == boss::BOSS::kSentinel)
                                return;

                            if (all_visited) {
                                unitig_index sibling_unitig = node_to_unitig[sibling_back];
                                auto find = seen.find(sibling_unitig);
                                if (find != seen.end() && !find->second)
                                    all_visited = false;
                            }
                        });

                        if (all_visited)
                            S.push_back(next_unitig);
                    }
                });

                if (found_cycle || !has_children)
                    break;

                if (S.size() == 1 && seen.size() == visited.size() + 1) {
                    dbg.adjacent_outgoing_nodes(std::get<1>(unitigs[S[0] - 1]), [&](node_index next_front) {
                        const auto &[front, back, size, seq] = unitigs[unitig_id - 1];
                        if (next_front == front)
                            found_cycle = true;
                    });

                    if (!found_cycle) {
                        // found a superbubble
                        #pragma omp critical
                        {
                        assert(visited[0] == unitig_id);
                        size_t superbubble_id = ++num_superbubbles;
                        uint64_t coord = 1;

                        assert(visited[0] < unitig_to_superbubble.size());
                        auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[visited[0]];
                        start_sb = superbubble_id;

                        for (size_t i = 1; i < visited.size(); ++i) {
                            coord += std::get<2>(unitigs[visited[i] - 1]);
                            assert(visited[i] < unitig_to_superbubble.size());
                            auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[visited[i]];
                            if (coord > mid_c) {
                                mid_sb = superbubble_id;
                                mid_c = coord;
                            }
                        }

                        {
                            assert(S[0] < unitig_to_superbubble.size());
                            auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[S[0]];
                            term_sb = superbubble_id;
                            term_c = coord;
                        }

                        superbubbles.emplace_back(std::move(visited), has_tips).first.emplace_back(S[0]);
                        assert(superbubbles.back().first.size() >= 4);
                        }
                    }

                    S.clear();
                }
            }
        }
    }

    logger->trace("Found {} (possibly nested) superbubbles", superbubbles.size());

    // find chains and print non-chained unitigs
    size_t cluster_id = 0;
    sdsl::bit_vector in_chain(superbubbles.size() + 1);
    std::vector<std::pair<std::vector<superbubble_index>, bool>> chains;

    {
        ProgressBar progress_bar(unitigs.size(), "Finding linear chains",
                                 std::cerr, !common::get_verbose());

        for (unitig_index unitig_id = 1; unitig_id < unitig_to_superbubble.size(); ++unitig_id) {
            ++progress_bar;

            const auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[unitig_id];
            // logger->info("Checking at {}\t{},{},{}", dbg.get_node_sequence(std::get<0>(unitigs[unitig_id - 1])),start_sb,mid_sb,term_sb);
            if (!mid_sb && start_sb && !term_sb && !in_chain[start_sb]) {
                // start of a chain
                // logger->info("\tStart!");
                auto &[chain, has_cycle] = chains.emplace_back();
                assert(!in_chain[start_sb]);
                in_chain[start_sb] = true;
                chain.emplace_back(start_sb);
                auto cur_start_sb = start_sb;
                while (cur_start_sb) {
                    const auto &[superbubble, has_tips] = superbubbles[cur_start_sb - 1];
                    cur_start_sb = 0;
                    assert(superbubble.back() < unitig_to_superbubble.size());
                    const auto &[next_start_sb, next_mid_sb, next_term_sb, next_mid_c, next_term_c] = unitig_to_superbubble[superbubble.back()];
                    if (next_start_sb && !next_mid_sb) {
                        chain.emplace_back(next_start_sb);

                        // break if we find a superbubble cycle
                        if (in_chain[next_start_sb]) {
                            // logger->info("\tCycle!");
                            has_cycle = true;
                            break;
                        }

                        assert(!in_chain[next_start_sb]);
                        in_chain[next_start_sb] = true;
                        cur_start_sb = next_start_sb;
                    }
                }
            }

            if (!start_sb && !mid_sb && !term_sb) {
                // logger->info("\tIsolated superbubble!");
                const auto &[front, back, size, seq] = unitigs[unitig_id - 1];
                assert(!output_unitig[unitig_id - 1]);
#ifndef NDEBUG
                output_unitig[unitig_id - 1] = true;
#endif
                callback(seq, cluster_id++);
            }
        }
    }

    {
        ProgressBar progress_bar(unitigs.size(), "Finding cyclic chains",
                                 std::cerr, !common::get_verbose());

        for (unitig_index unitig_id = 1; unitig_id < unitig_to_superbubble.size(); ++unitig_id) {
            ++progress_bar;

            const auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[unitig_id];
            // logger->info("Checking at {}\t{},{},{}", dbg.get_node_sequence(std::get<0>(unitigs[unitig_id - 1])),start_sb,mid_sb,term_sb);
            if (!mid_sb && start_sb && term_sb && !in_chain[start_sb]) {
                // start of a cycle chain
                // logger->info("\tStart!");
                auto &[chain, has_cycle] = chains.emplace_back();
                assert(!in_chain[start_sb]);
                in_chain[start_sb] = true;
                chain.emplace_back(start_sb);
                auto cur_start_sb = start_sb;
                while (cur_start_sb) {
                    const auto &[superbubble, has_tips] = superbubbles[cur_start_sb - 1];
                    cur_start_sb = 0;
                    assert(superbubble.back() < unitig_to_superbubble.size());
                    const auto &[next_start_sb, next_mid_sb, next_term_sb, next_mid_c, next_term_c] = unitig_to_superbubble[superbubble.back()];
                    if (next_start_sb && !next_mid_sb) {
                        chain.emplace_back(next_start_sb);

                        // break if we find a superbubble cycle
                        if (in_chain[next_start_sb]) {
                            has_cycle = true;
                            break;
                        }

                        assert(!in_chain[next_start_sb]);
                        in_chain[next_start_sb] = true;
                        cur_start_sb = next_start_sb;
                    }
                }

                assert(has_cycle);
            }
        }
    }

    // outputting superbubbles
    for (const auto &[chain, has_cycle] : chains) {
        std::for_each(chain.begin(), chain.end() - 1, [&](superbubble_index superbubble_id) {
            const auto &[unitig_ids, has_tips] = superbubbles[superbubble_id - 1];
            std::for_each(unitig_ids.begin(), unitig_ids.end() - 1,
                [&](unitig_index unitig_id) {
                    const auto &[front, back, size, seq] = unitigs[unitig_id - 1];
                    assert(!output_unitig[unitig_id - 1]);
#ifndef NDEBUG
                    output_unitig[unitig_id - 1] = true;
#endif
                    callback(seq, cluster_id);
                    if (has_cycle)
                        ++cluster_id;
                }
            );
        });

        if (has_cycle)
            continue;

        const auto &[unitig_ids, has_tips] = superbubbles[chain.back() - 1];
        std::for_each(unitig_ids.begin(), unitig_ids.end(),
            [&](unitig_index unitig_id) {
                const auto &[front, back, size, seq] = unitigs[unitig_id - 1];
                assert(!output_unitig[unitig_id - 1]);
#ifndef NDEBUG
                output_unitig[unitig_id - 1] = true;
#endif
                callback(seq, cluster_id);
            }
        );

        ++cluster_id;
    }

    ProgressBar progress_bar(unitigs.size(), "Outputting remaining superbubbles",
                             std::cerr, !common::get_verbose());

    for (unitig_index unitig_id = 1; unitig_id < unitig_to_superbubble.size(); ++unitig_id) {
        ++progress_bar;

        const auto &[start_sb, mid_sb, term_sb, mid_c, term_c] = unitig_to_superbubble[unitig_id];
        if (!mid_sb && start_sb && !in_chain[start_sb]) {
            in_chain[start_sb] = true;
            auto superbubble_id = start_sb;
            const auto &[unitig_ids, has_tips] = superbubbles[superbubble_id - 1];
            std::for_each(unitig_ids.begin(), unitig_ids.end(),
                [&](unitig_index unitig_id) {
                    const auto &[front, back, size, seq] = unitigs[unitig_id - 1];
                    assert(!output_unitig[unitig_id - 1]);
#ifndef NDEBUG
                    output_unitig[unitig_id - 1] = true;
#endif
                    callback(seq, cluster_id);
                }
            );

            ++cluster_id;
        }
    }

    logger->trace("Found {} unitig clusters", cluster_id);
}

} // namespace graph
} // namespace mtg
