#include "path_index.hpp"

#include <mutex>

#include <progress_bar.hpp>

#include "graph/annotated_dbg.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "annotation/annotation_converters.hpp"
#include "graph/alignment/alignment.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::graph {
using namespace annot;
using namespace annot::matrix;
using namespace annot::binmat;

using common::logger;

using Label = AnnotatedDBG::Label;
using Row = MultiIntMatrix::Row;
using Tuple = MultiIntMatrix::Tuple;

constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

static const std::vector<Label> DUMMY { Label(1, 1) };

void IPathIndex::call_dists(size_t path_id_1,
                            size_t path_id_2,
                            const std::function<void(size_t)> &callback,
                            size_t max_dist) const {
    if (path_id_1 == path_id_2) {
        callback(0);
        return;
    }

    if (!max_dist || !is_unitig(path_id_1) || !is_unitig(path_id_2) || path_length(path_id_1) > max_dist)
        return;

    size_t sb1 = path_id_1;
    size_t sb2 = path_id_2;
    std::vector<size_t> d1 { 0 };
    std::vector<size_t> d2 { 0 };

    bool is_source1 = get_superbubble_terminus(path_id_1).first;
    bool is_source2 = get_superbubble_terminus(path_id_2).first;

    if (!is_source1)
        std::tie(sb1, d1) = get_superbubble_and_dist(path_id_1);

    if (!is_source2)
        std::tie(sb2, d2) = get_superbubble_and_dist(path_id_2);

    // logger->info("checking dists for {}:{}:{}:{} -> {}:{}:{}\t{}",
    //     sb1, path_id_1, is_source1,path_length(path_id_1),
    //     sb2, path_id_2, is_source2,
    //     max_dist);

    if (!sb1) {
        // logger->info("\tfull traversal\tmax_dist: {}", max_dist);
        std::vector<std::pair<size_t, size_t>> path;
        path.emplace_back(path_id_1, 0);
        while (path.size()) {
            auto [uid, d] = path.back();
            path.pop_back();
            if (d > max_dist)
                continue;

            d += path_length(uid);
            adjacent_outgoing_unitigs(uid, [&](size_t next) {
                assert(d);
                if (max_dist > d && get_superbubble_terminus(next).first) {
                    call_dists(next, path_id_2, [&](size_t l) {
                        callback(l + d);
                    }, max_dist - d);
                } else {
                    path.emplace_back(next, d);
                }
            });
        }

        return;
    }

    if (!sb2) {
        if (!can_reach_superbubble_terminus(path_id_1))
            return;

        auto [t, d] = get_superbubble_terminus(sb1);
        if (!is_source1 && d.size() != 1) {
            //TODO
            return;
        }

        // logger->info("\tskip superbubble, then full traversal");

        if (is_source1) {
            if (d.size() == 1 && !d[0])
                return;

            size_t msize = d[1];
            call_dists(t, path_id_2, [&](size_t l) {
                for (size_t dd : d) {
                    callback(l + dd);
                }
            }, max_dist - msize);
            return;
        }

        size_t msize = std::numeric_limits<size_t>::max();
        for (size_t dd1 : d1) {
            for (size_t dd : d) {
                if (dd > dd1)
                    msize = std::min(msize, dd - dd1);
            }
        }
        assert(msize);
        if (max_dist > msize) {
            // logger->info("\t\tgo");
            call_dists(t, path_id_2, [&](size_t l) {
                for (size_t dd1 : d1) {
                    for (size_t dd : d) {
                        callback(l + dd - dd1);
                    }
                }
            }, max_dist - msize);
        }

        return;
    }

    if (sb1 == sb2) {
        if (is_source1) {
            std::for_each(d2.begin(), d2.end(), callback);
            return;
        }

        // logger->info("\tsame superbubble");

        if (!can_reach_superbubble_terminus(path_id_1)) {
            if (can_reach_superbubble_terminus(path_id_2))
                return;

            std::vector<std::pair<size_t, size_t>> path;
            path.emplace_back(path_id_1, 0);
            while (path.size()) {
                auto [uid, d] = path.back();
                path.pop_back();
                if (d > max_dist)
                    continue;

                size_t len = path_length(uid);
                adjacent_outgoing_unitigs(uid, [&](size_t next) {
                    if (next == path_id_2) {
                        callback(len + d);
                    } else {
                        path.emplace_back(next, len + d);
                    }
                });
            }

            return;
        }

        auto [t, d] = get_superbubble_terminus(sb1);
        if (path_id_2 == t) {
            if (d.size() != 1) {
                // TODO later
                return;
            }

            for (size_t dd1 : d1) {
                for (size_t dd : d) {
                    callback(dd - dd1);
                }
            }

            return;
        }

        // TODO
        return;
    }

    if (!can_reach_superbubble_terminus(path_id_1))
        return;

    auto [t, d] = get_superbubble_terminus(sb1);
    if (!is_source1 && d.size() != 1) {
        // TODO:
        return;
    }

    // logger->info("\tdifferent superbubble");

    if (is_source1) {
        if (d.size() == 1 && !d[0])
            return;

        size_t msize = d[1];
        call_dists(t, path_id_2, [&](size_t l) {
            for (size_t dd : d) {
                callback(l + dd);
            }
        }, max_dist - msize);
        return;
    }

    size_t msize = std::numeric_limits<size_t>::max();
    for (size_t dd : d) {
        for (size_t dd1 : d1) {
            if (dd > dd1)
                msize = std::min(msize, dd - dd1);
        }
    }
    assert(msize);
    if (max_dist > msize) {
        call_dists(t, path_id_2, [&](size_t l) {
            for (size_t dd : d) {
                for (size_t dd1 : d1) {
                    callback(l + dd - dd1);
                }
            }
        }, max_dist - msize);
    }
}

template <class Indicator, class Storage>
std::pair<size_t, std::vector<size_t>> get_range(const Indicator &indicator,
                                                 const Storage &storage,
                                                 size_t i) {
    assert(i);
    size_t num_bits = indicator.num_set_bits();
    assert(storage.size() >= num_bits);
    assert(i <= num_bits);
    size_t begin = indicator.select1(i) + 1;
    size_t end = i != num_bits ? indicator.next1(begin) : storage.size();
    assert(end <= storage.size());
    assert(begin <= end);
    std::vector<size_t> values;
    values.reserve(end - begin);
    std::copy(storage.begin() + begin, storage.begin() + end, std::back_inserter(values));
    return std::make_pair(storage[begin - 1], std::move(values));
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::adjacent_outgoing_unitigs(size_t path_id,
                            const std::function<void(size_t)> &callback) const {
    if (--path_id >= num_unitigs_)
        return;

    dbg_succ_->adjacent_outgoing_nodes(unitig_backs_[path_id], [&](node_index next) {
        if (size_t next_id = unitig_fronts_[next])
            callback(next_id);
    });
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
std::pair<size_t, std::vector<size_t>>
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_terminus(size_t path_id) const {
    if (path_id > num_unitigs_ || !path_id)
        return {};

    return get_range(superbubble_termini_b_, superbubble_termini_, path_id);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
std::pair<size_t, std::vector<size_t>>
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_and_dist(size_t path_id) const {
    if (path_id > num_unitigs_ || !path_id)
        return {};

    return get_range(superbubble_sources_b_, superbubble_sources_, path_id);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::can_reach_superbubble_terminus(size_t path_id) const {
    return --path_id < num_unitigs_ && can_reach_terminus_[path_id];
}


template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::load(const std::string &filename_base) {
    auto in = utils::open_ifstream(filename_base + kPathIndexExtension);
    if (!in->good())
        return false;

    num_unitigs_ = load_number(*in);
    num_superbubbles_ = load_number(*in);

    if (!paths_indices_.load(*in))
        return false;

    if (!path_boundaries_.load(*in))
        return false;

    logger->trace("Loaded {} paths", path_boundaries_.num_set_bits());

    try {
        unitig_backs_.load(*in);
        unitig_fronts_.load(*in);
        superbubble_sources_.load(*in);
        superbubble_termini_.load(*in);
    } catch (...) {
        return false;
    }

    if (!superbubble_sources_b_.load(*in))
        return false;

    if (!superbubble_termini_b_.load(*in))
        return false;

    if (!can_reach_terminus_.load(*in))
        return false;

    logger->trace("Loaded {} superbubbles", num_superbubbles_);

    return true;
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::serialize(const std::string &filename_base) const {
    std::ofstream fout(filename_base + kPathIndexExtension);
    serialize_number(fout, num_unitigs_);
    serialize_number(fout, num_superbubbles_);
    paths_indices_.serialize(fout);
    path_boundaries_.serialize(fout);
    unitig_backs_.serialize(fout);
    unitig_fronts_.serialize(fout);
    superbubble_sources_.serialize(fout);
    superbubble_termini_.serialize(fout);
    superbubble_sources_b_.serialize(fout);
    superbubble_termini_b_.serialize(fout);
    can_reach_terminus_.serialize(fout);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::set_graph(std::shared_ptr<const DBGSuccinct> graph) {
    dbg_succ_ = graph;
    if constexpr(std::is_base_of_v<IRowDiff, PathStorage>) {
        static_cast<IRowDiff&>(paths_indices_).set_graph(dbg_succ_.get());
    }
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::PathIndex(std::shared_ptr<const DBGSuccinct> graph,
            const std::string &graph_name,
            const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences) {
    const DBGSuccinct &dbg_succ = *graph;

    LabelEncoder<Label> label_encoder;
    label_encoder.insert_and_encode(DUMMY[0]);

    AnnotatedDBG anno_graph(std::const_pointer_cast<DBGSuccinct>(graph),
                            std::make_unique<ColumnCompressed<>>(dbg_succ.max_index()));

    auto &annotator = const_cast<ColumnCompressed<>&>(
        static_cast<const ColumnCompressed<>&>(anno_graph.get_annotator())
    );

    std::shared_ptr<const DeBruijnGraph> check_graph = graph;
    std::shared_ptr<const CanonicalDBG> canonical;
    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        canonical = std::make_shared<CanonicalDBG>(graph);
        check_graph = canonical;
    }

    std::vector<uint64_t> boundaries { 0 };
    std::vector<node_index> unitig_fronts;
    std::vector<node_index> unitig_backs;
    tsl::hopscotch_map<node_index, size_t> front_to_unitig_id;
    tsl::hopscotch_map<node_index, size_t> back_to_unitig_id;

    std::mutex mu;
    dbg_succ.call_unitigs([&](const auto &seq, const auto &path) {
        std::ignore = seq;
        auto rows = path;
        std::transform(rows.begin(), rows.end(), rows.begin(), AnnotatedDBG::graph_to_anno_index);

        std::lock_guard<std::mutex> lock(mu);
        front_to_unitig_id[path.front()] = unitig_fronts.size();
        back_to_unitig_id[path.back()] = unitig_fronts.size();
        unitig_fronts.emplace_back(path.front());
        unitig_backs.emplace_back(path.back());
        uint64_t coord = boundaries.back();
        annotator.add_labels(rows, DUMMY);
        for (auto row : rows) {
            annotator.add_label_coord(row, DUMMY, coord++);
        }

        boundaries.emplace_back(coord);
    }, get_num_threads());

    num_unitigs_ = boundaries.size() - 1;

    // enumerate superbubbles
    {
        sdsl::bit_vector can_reach_terminus(num_unitigs_, false);
        std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_starts(num_unitigs_);
        std::atomic<size_t> superbubble_start_size { num_unitigs_ };
        std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_termini(num_unitigs_);
        std::atomic<size_t> superbubble_termini_size { num_unitigs_ };

        std::atomic<size_t> num_terminal_superbubbles { 0 };
        std::atomic_thread_fence(std::memory_order_release);

        ProgressBar progress_bar(num_unitigs_, "Indexing superbubbles",
                                 std::cerr, !common::get_verbose());
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < num_unitigs_; ++i) {
            ++progress_bar;
            tsl::hopscotch_set<size_t> visited;
            VectorMap<size_t, tsl::hopscotch_set<size_t>> seen;
            tsl::hopscotch_map<size_t, std::vector<size_t>> parents;
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

                visited.insert(unitig_id);
                bool has_cycle = false;
                size_t num_children = 0;
                size_t length = boundaries[unitig_id + 1] - boundaries[unitig_id];
                dbg_succ.call_outgoing_kmers(unitig_backs[unitig_id], [&](node_index next, char c) {
                    if (c == boss::BOSS::kSentinel || next == unitig_fronts[unitig_id])
                        return;

                    ++num_children;

                    if (has_cycle)
                        return;

                    if (next == unitig_fronts[i]) {
                        has_cycle = true;
                        return;
                    }

                    assert(front_to_unitig_id.count(next));
                    size_t next_id = front_to_unitig_id[next];

                    bool add_parents = !seen.count(next_id);

                    seen[next_id].emplace(dist + length);
                    bool all_visited = true;
                    dbg_succ.call_incoming_kmers(next, [&](node_index sibling, char c) {
                        if (c != boss::BOSS::kSentinel) {
                            assert(back_to_unitig_id.count(sibling));
                            size_t sibling_id = back_to_unitig_id[sibling];
                            if (sibling_id == next_id)
                                return;

                            if (add_parents)
                                parents[next_id].emplace_back(sibling_id);

                            if (all_visited && !visited.count(sibling_id))
                                all_visited = false;
                        }
                    });

                    if (all_visited)
                        traversal_stack.emplace_back(next_id, dist + length);
                });

                bool reached_end = (traversal_stack.size() == 1 && visited.size() + 1 == seen.size());
                if (has_cycle) {
                    // logger->info("found cycle: num_children: {}\tv: {}\ts: {}", num_children, visited.size(), seen.size());
                    is_terminal_superbubble = false;
                    break;
                }

                if (!num_children) {
                    is_terminal_superbubble = true;
                }

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


                    set_bit(can_reach_terminus.data(), i, true, MO_RELAXED);
                    const auto &d = seen[terminus];

                    {
                        std::lock_guard<std::mutex> lock(mu);
                        superbubble_termini_size -= superbubble_termini[i].second.size();
                        superbubble_termini_size += d.size();
                        superbubble_termini[i] = std::make_pair(terminus + 1,
                            std::vector<size_t>(d.begin(), d.end()));
                        std::sort(superbubble_termini[i].second.begin(),
                                  superbubble_termini[i].second.end());
                    }

                    for (const auto &[u_id, d] : seen) {
                        if (u_id == i)
                            continue;

                        std::vector<size_t> cur_d(d.begin(), d.end());
                        std::sort(cur_d.begin(), cur_d.end());

                        if (!superbubble_starts[u_id].first
                                || cur_d.back() < superbubble_starts[u_id].second.back()) {
                            superbubble_start_size -= superbubble_starts[u_id].second.size();
                            superbubble_start_size += cur_d.size();
                            superbubble_starts[u_id] = std::make_pair(i + 1, std::move(cur_d));
                            can_reach_terminus[u_id] = !is_terminal_superbubble;
                        }
                    }
                }
            }

            if (is_terminal_superbubble && seen.size() > 1 && terminus) {
                sdsl::bit_vector found_map(seen.size(), false);
                std::vector<size_t> back_traversal_stack;
                back_traversal_stack.reserve(seen.size());
                back_traversal_stack.emplace_back(terminus);
                while (back_traversal_stack.size()) {
                    size_t cur_id = back_traversal_stack.back();
                    back_traversal_stack.pop_back();

                    found_map[seen.find(cur_id) - seen.begin()] = true;
                    for (size_t parent : parents[cur_id]) {
                        back_traversal_stack.emplace_back(parent);
                    }
                }

                auto it = found_map.begin();
                for (const auto &[cur_id, stuff] : seen) {
                    std::lock_guard<std::mutex> lock(mu);

                    // TODO: what if this is a smaller superbubble, but we can't reach the terminus?
                    if (superbubble_starts[cur_id].first == i + 1)
                        can_reach_terminus[cur_id] = *it;

                    ++it;
                }

                num_terminal_superbubbles.fetch_add(1, MO_RELAXED);
            }
        }

        std::atomic_thread_fence(std::memory_order_acquire);

        size_t num_in_superbubble = num_unitigs_;
        {
            sdsl::int_vector<> sources(superbubble_start_size);
            sdsl::bit_vector source_indicator(superbubble_start_size);
            auto it = sources.begin();
            auto jt = source_indicator.begin();
            size_t j = 0;
            for (const auto &[start, d] : superbubble_starts) {
                *jt = true;
                *it = start;
                ++it;
                ++jt;

                // TODO: this is not correct
                if (!start && !superbubble_termini[j].first)
                    --num_in_superbubble;

                std::copy(d.begin(), d.end(), it);
                it += d.size();
                jt += d.size();
                ++j;
            }
            superbubble_sources_ = SuperbubbleStorage(std::move(sources));
            superbubble_sources_b_ = SuperbubbleIndicator(std::move(source_indicator));
        }

        num_superbubbles_ = 0;
        size_t num_multiple_sizes = 0;
        {
            sdsl::int_vector<> termini(superbubble_termini_size);
            sdsl::bit_vector termini_indicator(superbubble_termini_size);
            auto it = termini.begin();
            auto jt = termini_indicator.begin();
            for (const auto &[start, d] : superbubble_termini) {
                *jt = true;
                *it = start;
                ++it;
                ++jt;
                num_superbubbles_ += !d.empty();
                num_multiple_sizes += (d.size() > 1);
                std::copy(d.begin(), d.end(), it);
                it += d.size();
                jt += d.size();
            }
            superbubble_termini_ = SuperbubbleStorage(std::move(termini));
            superbubble_termini_b_ = SuperbubbleIndicator(std::move(termini_indicator));
        }

        can_reach_terminus_ = SuperbubbleIndicator(std::move(can_reach_terminus));

        logger->info("{} / {} unitigs are in a superbubble. Indexed {} superbubbles, of which {} have dead ends and {} have multiple paths to the terminus.",
                     num_in_superbubble, num_unitigs_,
                     num_superbubbles_,
                     num_terminal_superbubbles,
                     num_multiple_sizes);

        sdsl::int_vector<> unitig_backs_sdsl(num_unitigs_);
        std::copy(unitig_backs.begin(), unitig_backs.end(), unitig_backs_sdsl.begin());
        unitig_backs_ = SuperbubbleStorage(std::move(unitig_backs_sdsl));

        sdsl::int_vector<> unitig_fronts_sdsl(dbg_succ.max_index() + 1);
        for (size_t i = 0; i < unitig_fronts.size(); ++i) {
            unitig_fronts_sdsl[unitig_fronts[i]] = i;
        }
        unitig_fronts_ = SuperbubbleStorage(std::move(unitig_fronts_sdsl));
    }

    size_t seq_count = 0;
    size_t total_seq_count = 0;

    generate_sequences([&](std::string_view seq) {
        total_seq_count += 1 + (dbg_succ.get_mode() != DeBruijnGraph::BASIC);
        auto nodes = map_to_nodes_sequentially(*check_graph, seq);

        if (std::any_of(nodes.begin(), nodes.end(), [&](node_index node) {
            return !node || check_graph->has_multiple_outgoing(node) || check_graph->indegree(node) > 1;
        })) {
            ++seq_count;
            uint64_t coord = boundaries.back();
            for (node_index &node : nodes) {
                if (node) {
                    if (canonical)
                        node = canonical->get_base_node(node);

                    annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);
                }
                ++coord;
            }
            boundaries.emplace_back(coord);

            if (canonical) {
                ++seq_count;
                uint64_t coord = boundaries.back();
                for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
                    if (*it)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(*it), DUMMY, coord);

                    ++coord;
                }
                boundaries.emplace_back(coord);
            } else if (dbg_succ.get_mode() == DeBruijnGraph::CANONICAL) {
                ++seq_count;
                uint64_t coord = boundaries.back();
                std::string seq_rc(seq);
                reverse_complement_seq_path(dbg_succ, seq_rc, nodes);
                for (node_index node : nodes) {
                    if (node)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);

                    ++coord;
                }
                boundaries.emplace_back(coord);
            }
        }
    });

    if (total_seq_count)
        logger->info("Indexed {} / {} sequences", seq_count, total_seq_count);

    assert(annotator.num_labels() <= 1);
    assert(std::adjacent_find(boundaries.begin(), boundaries.end()) == boundaries.end());

    path_boundaries_ = bit_vector_smart([&](const auto &callback) {
        std::for_each(boundaries.begin(), boundaries.end(), callback);
    }, boundaries.back() + 1, boundaries.size());

    logger->info("Indexed a total of {} paths", path_boundaries_.num_set_bits() - 1);

    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
    std::string out_path = tmp_dir/"test_col";
    annotator.serialize(out_path);

    std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
    if (!std::filesystem::exists(files[0])) {
        logger->error("Failed to serialize annotation to {}.", files[0]);
        std::exit(1);
    }

    if constexpr(std::is_same_v<PathStorage, ColumnCoordAnnotator::binary_matrix_type>) {
        paths_indices_ = const_cast<PathStorage&&>(load_coords(
            std::move(annotator),
            files
        ).get_matrix());

    } else if constexpr(std::is_same_v<PathStorage, RowDiffCoordAnnotator::binary_matrix_type>) {
        std::string graph_fname = graph_name;
        if (graph_fname.empty()) {
            graph->serialize(out_path);
            graph_fname = out_path + graph->file_extension();
        }

        if (!std::filesystem::exists(graph_fname)) {
            logger->error("Graph path incorrect: {}.", graph_fname);
            std::exit(1);
        }

        {
            std::filesystem::path swap_dir = utils::create_temp_dir("", "swap_col");
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(0), out_path + ".row_count", false, true);
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(1), out_path + ".row_reduction", false, true);
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(2), out_path, false, true);
        }

        const std::string anchors_file = graph_fname + kRowDiffAnchorExt;
        const std::string fork_succ_file = graph_fname + kRowDiffForkSuccExt;
        if (!std::filesystem::exists(anchors_file)) {
            logger->error("Anchor bitmap {} does not exist.", anchors_file);
            std::exit(1);
        }
        if (!std::filesystem::exists(fork_succ_file)) {
            logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
            std::exit(1);
        }

        std::unique_ptr<AnnotatedDBG::Annotator> annotator;
        size_t num_columns = anno_graph.get_annotator().get_label_encoder().size();
        auto diff_annotator = std::make_unique<ColumnCompressed<>>(0);
        if (!diff_annotator->merge_load(files)) {
            logger->error("Cannot load annotations from {}", files[0]);
            exit(1);
        }
        assert(diff_annotator->num_labels() == num_columns);
        std::vector<bit_vector_smart> delimiters;
        std::vector<sdsl::int_vector<>> column_values;

        typedef ColumnCoordAnnotator::binary_matrix_type CoordDiff;
        auto coords_fname = utils::remove_suffix(files[0], ColumnCompressed<>::kExtension)
                                                        + ColumnCompressed<>::kCoordExtension;
        std::ifstream in(coords_fname, std::ios::binary);
        try {
            CoordDiff::load_tuples(in, num_columns, [&](auto&& delims, auto&& values) {
                delimiters.emplace_back(std::move(delims));
                column_values.emplace_back(std::move(values));
            });
        } catch (const std::exception &e) {
            logger->error("Couldn't load coordinates from {}\nException: {}", coords_fname, e.what());
            exit(1);
        } catch (...) {
            logger->error("Couldn't load coordinates from {}", coords_fname);
            exit(1);
        }

        annotator.reset(new RowDiffCoordAnnotator(
            label_encoder,
            graph.get(),
            std::move(*diff_annotator->release_matrix()),
            std::move(delimiters), std::move(column_values)
        ));

        auto &row_diff = const_cast<PathStorage&>(dynamic_cast<const PathStorage&>(
            annotator->get_matrix()
        ));
        row_diff.load_anchor(anchors_file);
        row_diff.load_fork_succ(fork_succ_file);
        paths_indices_ = const_cast<PathStorage&&>(row_diff);

    } else {
        throw std::runtime_error("Only ColumnCoord and RowDiffCoord annotators supported");
    }

    set_graph(graph);
}

auto IPathIndex
::get_coords(const std::vector<node_index> &nodes) const -> std::vector<RowTuples> {
    sdsl::bit_vector picked(nodes.size(), true);

    std::vector<Row> rows;
    rows.reserve(nodes.size());
    tsl::hopscotch_map<size_t, std::vector<size_t>> path_id_to_nodes;

    for (size_t i = 0; i < nodes.size(); ++i) {
        if (!has_coord(nodes[i])) {
            picked[i] = false;
            continue;
        }

        rows.emplace_back(AnnotatedDBG::graph_to_anno_index(nodes[i]));
    }

    auto it = picked.begin();
    auto row_tuples = get_row_tuples(rows);
    std::vector<RowTuples> ret_val;
    ret_val.reserve(nodes.size());
    while (it != picked.end() && !*it) {
        ret_val.emplace_back();
        ++it;
    }

    for (auto &tuples : row_tuples) {
        VectorMap<size_t, Tuple> out_tuples;
        assert(tuples.size() <= 1);
        for (const auto &[c, tuple] : tuples) {
            assert(!c);
            assert(std::adjacent_find(tuple.begin(), tuple.end()) == tuple.end());
            for (auto coord : tuple) {
                size_t path_id = coord_to_path_id(coord);
                out_tuples[path_id].emplace_back(coord);
                path_id_to_nodes[path_id].emplace_back(it - picked.begin());
            }
        }
        ret_val.emplace_back(out_tuples.values_container().begin(),
                             out_tuples.values_container().end());
        ++it;
        while (it != picked.end() && !*it) {
            ret_val.emplace_back();
            ++it;
        }
    }

    return ret_val;
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::has_coord(node_index node) const {
    assert(dbg_succ_);
    return node != DeBruijnGraph::npos
        && dbg_succ_->get_node_sequence(node).find(boss::BOSS::kSentinel) == std::string::npos;
}

template class PathIndex<>;

} // namespace mtg::graph
