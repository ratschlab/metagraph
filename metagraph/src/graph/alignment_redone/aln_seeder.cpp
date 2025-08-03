#include "aln_seeder.hpp"

#include <tsl/hopscotch_set.h>

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::graph::align_redone {

std::vector<Anchor> ExactSeeder::get_anchors() const {
    const DeBruijnGraph &graph = query_.get_graph();
    std::vector<Anchor> anchors;

    if (query_.get_query().size() < config_.min_seed_length)
        return anchors;

    for (bool orientation : { false, true }) {
        std::string_view this_query = query_.get_query(orientation);
        auto nodes = map_to_nodes_sequentially(graph, this_query);
        for (size_t begin = 0; begin < nodes.size(); ++begin) {
            Match::node_index node = nodes[begin];
            if (node != DeBruijnGraph::npos) {
                anchors.emplace_back(this_query,
                                        begin, begin + graph.get_k(),
                                        orientation,
                                        std::vector<Match::node_index>{ node },
                                        config_);
            }
        }

        if (config_.min_seed_length < graph.get_k()) {
            const DBGSuccinct *dbg_succ;
            const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
            if (canonical) {
                dbg_succ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
            } else {
                dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
            }

            if (dbg_succ) {
                const auto &boss = dbg_succ->get_boss();
                auto encoded = boss.encode(this_query);

                size_t k = graph.get_k();
                for (size_t i = 0; i + config_.min_seed_length <= this_query.size(); ++i) {
                    if (std::find(encoded.begin() + i, encoded.begin() + i + config_.min_seed_length,
                                boss.alph_size) != encoded.begin() + i + config_.min_seed_length) {
                        continue;
                    }

                    auto [first, last, it] = boss.index_range(
                        encoded.begin() + i,
                        encoded.begin() + i + config_.min_seed_length
                    );

                    size_t match_size = it - (encoded.begin() + i);
                    if (match_size < config_.min_seed_length)
                        continue;

                    first = first == boss.select_last(1) ? 1 : boss.pred_last(first - 1) + 1;
                    // std::cerr << i << "\t" << std::string_view(this_query.begin() + i, match_size) << "\t" << first << "," << last << "\t" << std::flush << graph.get_node_sequence(first) << "\t" << graph.get_node_sequence(last) << "\n";
                    using edge_index = boss::BOSS::edge_index;

                    std::vector<std::tuple<edge_index, edge_index, std::string, bool>> traverse;
                    traverse.emplace_back(first, last, "", true);
                    while (traverse.size()) {
                        auto [first, last, suffix, is_exact_match] = traverse.back();
                        assert(boss.get_last(last));
                        traverse.pop_back();

                        size_t cur_match_size = match_size + suffix.size();

                        if (canonical) {
                            if (is_exact_match) {
                                if (dbg_succ->in_graph(last)) {
                                    dbg_succ->adjacent_incoming_nodes(last, [&](DeBruijnGraph::node_index node) {
                                        node = canonical->reverse_complement(node);
                                        size_t i_rc = this_query.size() - (i + cur_match_size);
                                        anchors.emplace_back(query_.get_query(!orientation),
                                                            i_rc, i_rc + cur_match_size,
                                                            !orientation,
                                                            std::vector<Match::node_index>{ node },
                                                            config_,
                                                            graph.get_node_sequence(node).substr(cur_match_size));
                                        assert(anchors.back().get_spelling().size() + anchors.back().get_end_trim() == graph.get_k());
                                        assert(anchors.back().is_spelling_valid(graph));
                                    });
                                }
                            } else if (nodes[i] != DeBruijnGraph::npos) {
                                break;
                            }
                        }

                        assert(cur_match_size < k);
                        if (cur_match_size == k - 1) {
                            assert(first == last || !boss.get_last(first));
                            assert(boss.succ_last(first) == last);
                            // std::cerr << "foo\t" << i << "\t" << std::string_view(this_query.begin() + i, match_size) << "\n";
                            for (edge_index node = first; node <= last; ++node) {
                                if (dbg_succ->in_graph(node)) {
                                    char c = boss.decode(boss.get_W(node) % boss.alph_size);
                                    if (c != boss::BOSS::kSentinel) {
                                        // std::cerr << "\t" << suffix << c << "\n";
                                        auto first_mm = std::mismatch(this_query.begin() + i + match_size, this_query.end(),
                                                                    suffix.begin(), suffix.end()).second;
                                        bool last_char_matches = i + k <= this_query.size() && first_mm == suffix.end() && c == this_query[i + k - 1];
                                        size_t full_match_size = match_size + (first_mm - suffix.begin()) + last_char_matches;

                                        is_exact_match &= last_char_matches;

                                        if (!is_exact_match) {
                                            anchors.emplace_back(this_query,
                                                                i, i + full_match_size,
                                                                orientation,
                                                                std::vector<Match::node_index>{ node },
                                                                config_,
                                                                (suffix + c).substr(full_match_size - match_size));
                                            assert(anchors.back().get_spelling().size() + anchors.back().get_end_trim() == graph.get_k());
                                            assert(anchors.back().is_spelling_valid(graph));
                                        }
                                    }
                                }
                            }

                            continue;
                        }

                        // for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
                        for (boss::BOSS::edge_index e = first; e <= last; ++e) {
                            boss::BOSS::TAlphabet s = boss.get_W(e) % boss.alph_size;
                            // if (s) {
                            boss::BOSS::edge_index next_first = first;
                            boss::BOSS::edge_index next_last = last;
                            // std::cerr << i << "\t" << std::string_view(this_query.begin() + i, match_size) << "\ttry " << suffix << boss.decode(s) << "\t" << graph.get_node_sequence(next_first) << "\t" << graph.get_node_sequence(next_last) << "\n";
                            if (s && boss.tighten_range(&next_first, &next_last, s)) {
                                // std::cerr << "\tworked!\n";
                                char c = boss.decode(s);
                                bool next_is_exact_match = is_exact_match && i + cur_match_size < this_query.size() && this_query[i + cur_match_size] == c;
                                traverse.emplace_back(next_first, next_last, suffix + c,
                                                    next_is_exact_match);
                            }
                        }
                    }
                }

                if (anchors.size() > 1 && canonical) {
                    std::sort(anchors.begin(), anchors.end(), [](const auto &a, const auto &b) {
                        return std::make_tuple(a.get_orientation(), a.get_clipping(), a.get_path()[0], b.get_end_clipping())
                            < std::make_tuple(b.get_orientation(), b.get_clipping(), b.get_path()[0], a.get_end_clipping());
                    });
                    for (auto it = anchors.begin(), jt = anchors.begin() + 1; jt != anchors.end(); ++it, ++jt) {
                        if (it->get_orientation() != jt->get_orientation() || it->get_clipping() != jt->get_clipping() || it->get_path()[0] != jt->get_path()[0])
                            continue;

                        assert(it->get_end_clipping() >= jt->get_end_clipping());
                        assert(it->get_end_trim() >= jt->get_end_trim());
                        // std::cerr << "redundant\t" << *it << "\t" << *jt << "\n";
                        *it = Anchor();
                    }
                    anchors.erase(std::remove_if(anchors.begin(), anchors.end(), [](const auto &a) { return a.get_path().empty(); }),
                                anchors.end());
                }
            } else if (const auto *sshash = dynamic_cast<const DBGSSHash*>(&graph)) {
                // TODO
            }
        }
    }

    assert(std::all_of(anchors.begin(), anchors.end(), [k=graph.get_k()](const auto &a) {
        return a.get_spelling().size() + a.get_end_trim() == k;
    }));

    return anchors;
}

std::vector<Anchor> LabeledSeeder::get_anchors() const {
    auto anchors = ExactSeeder::get_anchors();

    for (const auto &anchor : anchors) {
        anno_buffer_.queue_path(anchor.get_path());
    }

    anno_buffer_.fetch_queued_annotations();

    std::vector<Anchor> labeled_anchors;
    labeled_anchors.reserve(anchors.size());
    for (auto &anchor : anchors) {
        auto [labels_it, coord_it] = anno_buffer_.get_labels_and_coords(anchor.get_path()[0]);
        if (labels_it) {
            for (size_t i = 0; i < labels_it->size(); ++i) {
                anchor.set_label_class(anno_buffer_.cache_column_set(labels_it->begin() + i,
                                                                     labels_it->begin() + i + 1));
                if (coord_it) {
                    const auto &coords = (*coord_it)[i];
                    for (auto coord : coords) {
                        labeled_anchors.emplace_back(anchor).set_coord(coord);
                    }
                } else {
                    labeled_anchors.emplace_back(anchor);
                }
            }
        }
        anchor = Anchor();
    }

    return labeled_anchors;
}

template <typename T>
class OffsetVector {
  public:
    OffsetVector() : offset_(0) {}

    template <typename... Args>
    OffsetVector(size_t offset, Args&&... args) : offset_(offset), null_(T()), data_(std::forward<Args>(args)...) {}

    size_t size() const { return offset_ + data_.size(); }
    bool empty() const { return offset_ == 0 && data_.empty(); }

    void resize(size_t new_size, const T &val = T()) {
        if (new_size >= offset_) {
            data_.resize(new_size - offset_, val);
            return;
        }

        new_size -= data_.size();
        data_.clear();
        assert(new_size <= offset_);
        offset_ -= new_size;
    }

    size_t offset() const { return offset_; }

    template <typename... Args>
    T& set(size_t i, Args&&... args) {
        get(i) = T(std::forward<Args>(args)...);
        return data_[i - offset_];
    }

    T& get(size_t i) {
        if (i < offset_) {
            data_.insert(data_.begin(), offset_ - i, null_);
            offset_ = i;
        }

        if (i - offset_ >= data_.size())
            data_.resize(i - offset_ + 1);

        return data_[i - offset_];
    }

    T& operator[](size_t i) {
        assert(i >= offset_);
        return data_[i - offset_];
    }

    const T& operator[](size_t i) const {
        return i >= offset_ ? data_[i - offset_] : null_;
    }

    using const_iterator = typename std::vector<T>::const_iterator;
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend() const { return data_.cend(); }

  private:
    size_t offset_;
    const T null_;
    std::vector<T> data_;
};

template <typename Tuple>
using WaveFront = OffsetVector<VectorMap<DeBruijnGraph::node_index, Tuple>>; // S(cost)

template <typename Tuple>
using ScoreTable = std::vector<WaveFront<Tuple>>;

DBGAlignerConfig::score_t cost_to_score(size_t cost,
                                        size_t query_size,
                                        size_t match_size,
                                        DBGAlignerConfig::score_t match_score) {
    DBGAlignerConfig::score_t double_score = match_score * (query_size + match_size) - cost;
    assert(double_score % 2 == 0);
    return double_score / 2;
}

template <typename Tuple>
const Tuple* get_bucket(const ScoreTable<Tuple> &table,
                        size_t cost,
                        size_t query_dist,
                        DeBruijnGraph::node_index node) {
    if (cost >= table.size())
        return nullptr;

    if (query_dist >= table[cost].size() || query_dist < table[cost].offset())
        return nullptr;

    const auto &table_slot = table[cost][query_dist];
    auto it = table_slot.find(node);
    if (it == table_slot.end())
        return nullptr;

    return &it->second;
};

// dist, num_ops, last_node, last_op, last_char_of_this_node, num_matches
using SMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index, Cigar::Operator, char, size_t>;

// dist, num_ops
using EMap = std::tuple<size_t, size_t>;

// dist, num_ops, last_node
using FMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index, char>;

template <typename StrItr>
void align_impl(const std::function<size_t(DeBruijnGraph::node_index, size_t, size_t)> &has_single_incoming,
                const std::function<DeBruijnGraph::node_index(DeBruijnGraph::node_index, char, size_t, size_t)> &traverse_back,
                const std::function<void(DeBruijnGraph::node_index, const std::function<void(DeBruijnGraph::node_index, char)>&, size_t, size_t)> &call_incoming_kmers,
                const DBGAlignerConfig &config,
                const DeBruijnGraph::node_index start_node,
                size_t start_num_matches,
                const StrItr query_window_begin,
                const StrItr query_window_end,
                size_t max_dist,
                const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t cost)> &callback,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    using node_index = DeBruijnGraph::node_index;

    if (query_window_begin == query_window_end)
        return;

    size_t query_size = query_window_end - query_window_begin;

    const auto query_window_rbegin = std::make_reverse_iterator(query_window_end);
    const auto query_window_rend = std::make_reverse_iterator(query_window_begin);

    // derived from 2.4.1 in https://doi.org/10.1101/2022.01.12.476087 (v1)
    DBGAlignerConfig::score_t match_score = config.match_score("A");
    DBGAlignerConfig::score_t mismatch_score = config.score_sequences("A", "T");

    assert(config.gap_opening_penalty <= config.gap_extension_penalty);
    assert(config.gap_extension_penalty < 0);
    assert(mismatch_score < match_score);

    ssize_t mismatch_cost = (match_score - mismatch_score) * 2;
    ssize_t gap_ext = match_score - 2 * config.gap_extension_penalty;
    ssize_t gap_opn = match_score - 2 * config.gap_opening_penalty;

    assert(mismatch_cost > 0);
    assert(gap_ext > 0);
    assert(gap_opn >= gap_ext);

    assert(mismatch_cost % 2 == 0);
    assert(gap_ext % 2 == 0);
    assert(gap_opn % 2 == 0);

    // common::logger->info("x: {}\to: {}\te: {}", mismatch_cost, gap_opn, gap_ext);

    // S(cost, query_dist, node) = best_dist
    ScoreTable<SMap> S; // best cost
    ScoreTable<EMap> E; // best cost (last operation is insertion)
    ScoreTable<FMap> F; // best cost (last operation is deletion)

    auto set_value = [&](
            auto &table,
            size_t cost,
            size_t query_dist,
            node_index node,
            size_t dist,
            node_index last_node,
            size_t num_ops,
            Cigar::Operator last_op,
            char c,
            size_t num_matches) -> size_t {
        assert(last_op == Cigar::MATCH
            || last_op == Cigar::MISMATCH
            || last_op == Cigar::INSERTION
            || last_op == Cigar::DELETION);
        if (cost >= table.size())
            table.resize(cost + 1);

        if (query_dist >= table[cost].size())
            table[cost].resize(query_dist + 1);

        auto &bucket = table[cost].get(query_dist)[node];
        using table_t = std::decay_t<decltype(bucket)>;
        static_assert(std::is_same_v<table_t, FMap>
                        || std::is_same_v<table_t, EMap>
                        || std::is_same_v<table_t, SMap>);

        bool inserted = dist > std::get<0>(bucket);

        if constexpr(std::is_same_v<table_t, SMap>) {
            inserted |= std::get<3>(bucket) == Cigar::CLIPPED || std::make_pair(dist, num_matches) > std::make_pair(std::get<0>(bucket), std::get<5>(bucket));
        }

        if constexpr(std::is_same_v<table_t, FMap>)
            inserted |= std::get<2>(bucket) == DeBruijnGraph::npos;

        if constexpr(std::is_same_v<table_t, EMap>)
            inserted |= std::get<1>(bucket) == 0;

        if (inserted) {
            std::get<0>(bucket) = dist;
            std::get<1>(bucket) = num_ops;

            if constexpr(std::is_same_v<table_t, FMap>
                            || std::is_same_v<table_t, SMap>) {
                assert(dist >= num_ops);
                std::get<2>(bucket) = last_node;
            }

            if constexpr(std::is_same_v<table_t, FMap>) {
                assert(c != '\0');
                std::get<3>(bucket) = c;
            }

            if constexpr(std::is_same_v<table_t, SMap>) {
                std::get<3>(bucket) = last_op;
                std::get<4>(bucket) = c;
                std::get<5>(bucket) = num_matches;
            }

            inserted = true;
        }

        if constexpr(std::is_same_v<table_t, SMap>) {
            assert(std::get<3>(bucket) == Cigar::MATCH
                || std::get<3>(bucket) == Cigar::MISMATCH
                || std::get<3>(bucket) == Cigar::INSERTION
                || std::get<3>(bucket) == Cigar::DELETION);
        }

        return inserted;
    };

    size_t num_explored_nodes = 0;
    {
        size_t query_dist = 0;
        DeBruijnGraph::node_index node = start_node;
        SMap data(0, 0, start_node, Cigar::MATCH, *(query_window_rbegin - 1), start_num_matches);
        auto &[best_dist, last_ext, last_node, last_op, c, num_matches] = data;

        // std::cerr << "start ext" << std::endl;

        for (auto it = query_window_rbegin; it != query_window_rend; ++it) {
            // std::cerr << "\tcheck\t" << node << std::endl;
            if (terminate_branch(0, data, query_dist, node) || best_dist == max_dist || !has_single_incoming(node, best_dist, query_dist))
                break;

            if (auto prev = traverse_back(node, *it, best_dist, query_dist)) {
                // std::cerr << "\t" << node << " -> " << prev << std::endl;
                ++num_explored_nodes;
                ++best_dist;
                ++query_dist;
                ++last_ext;
                ++num_matches;
                node = prev;
            } else {
                break;
            }
        }

        // std::cerr << "end ext" << std::endl;

        set_value(
            S,
            0,
            query_dist,
            node,
            best_dist,
            start_node,
            best_dist,
            Cigar::MATCH,
            best_dist > 0 ? *query_window_rbegin : '\0',
            num_matches
        );

        // std::cerr << "init ext\n";
    }

    // size_t max_edit_cost = std::max({ gap_opn, gap_ext, mismatch_cost });

    // for (size_t cost = 0; cost < S.size() + max_edit_cost; ++cost) {
    for (size_t cost = 0; cost < S.size(); ++cost) {
        bool done = false;
        for (size_t query_dist = S[cost].offset(); query_dist < S[cost].size(); ++query_dist) {
            assert(query_dist <= query_size);

            size_t it_dist = 0;
            for (auto it = S[cost][query_dist].begin(); it != S[cost][query_dist].end(); ++it, ++it_dist) {
                node_index node = it->first;
                size_t best_dist = std::get<0>(it->second);
                Cigar::Operator last_op = std::get<3>(it->second);
                size_t num_matches = std::get<5>(it->second);
                assert(last_op != Cigar::CLIPPED);

                // std::cerr << "Exp: " << cost << "\t" << node << "\t" << best_dist << "," << query_dist << "\n";

                if (terminate(cost, it->second, query_dist, node)) {
                    done = true;
                    // continue;
                    // common::logger->info("Halted at cost {}", cost);
                    break;
                }

                if (terminate_branch(cost, it->second, query_dist, node))
                    continue;

                // forward creation of insertions
                if (query_dist < query_size) {
                    // extend a previous insertion
                    if (const auto *ins_ext_bucket = get_bucket(E, cost, query_dist, node)) {
                        auto [last_dist, last_num_ops] = *ins_ext_bucket;
                        assert(last_num_ops > 0);
                        assert(query_dist >= last_num_ops);
                        ssize_t next_ext_cost = cost + gap_ext;
                        if (!terminate_branch(next_ext_cost, SMap(last_dist, 0, node, Cigar::INSERTION, '\0', 0), query_dist + 1, node)
                                && set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION, '\0', 0)
                                && set_value(S, next_ext_cost, query_dist + 1, node, last_dist, node, 0, Cigar::INSERTION, '\0', 0)) {
                            it = S[cost][query_dist].begin() + it_dist;
                        }
                    }

                    // open an insertion
                    if (last_op != Cigar::DELETION) {
                        ssize_t next_opn_cost = cost + gap_opn;
                        if (!terminate_branch(next_opn_cost, SMap(best_dist, 0, node, Cigar::INSERTION, '\0', 0), query_dist + 1, node)
                                && set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION, '\0', 0)
                                && set_value(S, next_opn_cost, query_dist + 1, node, best_dist, node, 0, Cigar::INSERTION, '\0', 0)) {
                            it = S[cost][query_dist].begin() + it_dist;
                        }
                    }
                }

                if (best_dist < max_dist) {
                    std::vector<std::pair<node_index, char>> prevs;
                    call_incoming_kmers(node, [&](auto prev, char c) {
                        ++num_explored_nodes;
                        if (c != boss::BOSS::kSentinel)
                            prevs.emplace_back(prev, c);
                    }, best_dist, query_dist);

                    if (prevs.empty())
                        continue;

                    // deletion extension
                    if (const auto *del_ext_bucket = get_bucket(F, cost, query_dist, node)) {
                        auto [last_dist, last_num_ops, last_node, del_char] = *del_ext_bucket;
                        assert(last_node != DeBruijnGraph::npos);
                        assert(last_num_ops > 0);
                        assert(last_dist >= last_num_ops);
                        if (last_dist < max_dist) {
                            ssize_t next_ext_cost = cost + gap_ext;
                            std::vector<std::pair<DeBruijnGraph::node_index, char>> last_prevs;
                            if (best_dist == last_dist) {
                                last_prevs = prevs;
                            } else {
                                call_incoming_kmers(node, [&](auto prev, char c) {
                                    ++num_explored_nodes;
                                    if (c != boss::BOSS::kSentinel)
                                        last_prevs.emplace_back(prev, c);
                                }, last_dist, query_dist);
                            }

                            for (const auto &[prev, c] : last_prevs) {
                                if (!terminate_branch(next_ext_cost, SMap(last_dist + 1, 0, node, Cigar::DELETION, c, 0), query_dist, prev)
                                        && set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION, c, 0)
                                        && set_value(S, next_ext_cost, query_dist, prev, last_dist + 1, node, 0, Cigar::DELETION, c, 0)) {
                                    it = S[cost][query_dist].begin() + it_dist;
                                }
                            }
                        }
                    }

                    // deletion open
                    if (last_op != Cigar::INSERTION) {
                        ssize_t next_opn_cost = cost + gap_opn;
                        for (const auto &[prev, c] : prevs) {
                            if (!terminate_branch(next_opn_cost, SMap(best_dist + 1, 0, node, Cigar::DELETION, c, 0), query_dist, prev)
                                    && set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION, c, 0)
                                    && set_value(S, next_opn_cost, query_dist, prev, best_dist + 1, node, 0, Cigar::DELETION, c, 0)) {
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }

                    if (query_dist < query_size) {
                        // match
                        auto local_query_window_begin = query_window_begin;
                        auto local_query_window_end = query_window_end;
                        local_query_window_end -= query_dist;

                        auto local_query_window_rbegin = std::make_reverse_iterator(local_query_window_end);
                        auto local_query_window_rend = std::make_reverse_iterator(local_query_window_begin);

                        for (auto [prev, c] : prevs) {
                            SMap data(best_dist + 1, 1, node, Cigar::MATCH, c, num_matches + 1);
                            auto &[cur_best, num_ops, last_node, op, cur_c, cur_num_matches] = data;
                            size_t cur_query_dist = query_dist + 1;
                            size_t cur_cost = cost;
                            if (c != *local_query_window_rbegin) {
                                cur_cost += mismatch_cost;
                                op = Cigar::MISMATCH;
                                cur_num_matches = 0;
                            }

                            for (auto jt = local_query_window_rbegin + 1; jt != local_query_window_rend; ++jt) {
                                if (terminate_branch(cur_cost, data, cur_query_dist, prev) || cur_best == max_dist || !has_single_incoming(prev, cur_best, cur_query_dist))
                                    break;

                                if (auto pprev = traverse_back(prev, *jt, cur_best, cur_query_dist)) {
                                    ++num_explored_nodes;
                                    ++cur_best;
                                    ++cur_query_dist;
                                    ++cur_num_matches;
                                    prev = pprev;
                                } else {
                                    break;
                                }
                            }

                            if (set_value(
                                S,
                                cur_cost,
                                cur_query_dist,
                                prev,
                                cur_best,
                                node,
                                cur_query_dist - query_dist,
                                op,
                                c,
                                cur_num_matches
                            )) {
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }
                }
            }
            // if (done)
            //     break;
        }

        if (done)
            break;
    }

    // backtrack
    for (size_t start_cost = 0; start_cost < S.size(); ++start_cost) {
        std::vector<std::tuple<size_t, ssize_t, size_t, DeBruijnGraph::node_index>> starts;
        for (size_t start_query_dist = S[start_cost].offset(); start_query_dist < S[start_cost].size(); ++start_query_dist) {
            assert(start_query_dist <= query_size);
            const auto &bucket = S[start_cost][start_query_dist];
            for (const auto &[node, data] : bucket) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                assert(last_op != Cigar::CLIPPED);
                // std::cerr << "BTIter: " << start_cost << "\t" << node << "\t" << dist << "," << start_query_dist << std::endl;
                if (start_backtrack(start_cost, data, start_query_dist, node))
                    starts.emplace_back(dist, static_cast<ssize_t>(dist) - start_query_dist, start_query_dist, node);
            }
        }

        for (auto [dist, diag, query_dist, node] : starts) {
            size_t cost = start_cost;

            const auto *bt = get_bucket(S, cost, query_dist, node);
            assert(bt);

            std::vector<node_index> path;
            Cigar cigar;
            std::string spelling;

            do {
                assert(std::get<0>(*bt) == dist);

                if (std::get<3>(*bt) == Cigar::INSERTION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *et = get_bucket(E, cost, query_dist, node);
                    assert(et);

                    auto [ins_dist, num_ins] = *et;
                    assert(query_dist >= num_ins);
                    assert(num_ins > 0);
                    assert(ins_dist == dist);

                    size_t prev_cost = cost - gap_opn - (num_ins - 1) * gap_ext;
                    assert(prev_cost < cost);
                    size_t prev_query_dist = query_dist - num_ins;

                    const auto *check_bt = get_bucket(S, prev_cost, prev_query_dist, node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);

                    query_dist = prev_query_dist;
                    cost = prev_cost;
                    bt = check_bt;
                    cigar.append(Cigar::INSERTION, num_ins);
                } else if (std::get<3>(*bt) == Cigar::DELETION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *ft = get_bucket(F, cost, query_dist, node);
                    assert(ft);

                    auto [del_dist, num_del, prev_node, del_char] = *ft;
                    assert(prev_node != DeBruijnGraph::npos);
                    assert(del_char != '\0');
                    assert(del_dist == dist);
                    assert(dist >= num_del);

                    cigar.append(Cigar::DELETION, num_del);

                    size_t prev_cost = cost;
                    size_t prev_dist = dist;
                    while (num_del > 1) {
                        path.emplace_back(node);
                        spelling += del_char;

                        prev_cost -= gap_ext;
                        --prev_dist;
                        node = prev_node;
                        --num_del;

                        const auto *prev_ft = get_bucket(F, prev_cost, query_dist, node);
                        assert(prev_ft);
                        assert(std::get<0>(*prev_ft) == prev_dist);
                        assert(std::get<1>(*prev_ft) == num_del);
                        ft = prev_ft;
                        prev_node = std::get<2>(*ft);
                        del_char = std::get<3>(*ft);
                    }

                    path.emplace_back(node);
                    spelling += del_char;

                    prev_cost -= gap_opn;
                    --prev_dist;

                    assert(prev_cost < cost);

                    const auto *check_bt = get_bucket(S, prev_cost, query_dist, prev_node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == prev_dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);

                    dist = prev_dist;
                    cost = prev_cost;
                    node = prev_node;
                    bt = check_bt;
                } else if (std::get<3>(*bt) == Cigar::MATCH || std::get<3>(*bt) == Cigar::MISMATCH) {
                    auto [cur_dist, num_match, last_node, last_op, mismatch_char, num_matches] = *bt;

                    if (!num_match) {
                        assert(last_op == Cigar::MATCH);
                        bt = nullptr;
                        break;
                    }

                    assert(query_dist >= num_match);
                    assert(dist >= num_match);
                    assert(cur_dist == dist);

                    node_index traverse_node = last_node;

                    cigar.append(Cigar::MATCH, num_match - 1);
                    cigar.append(last_op);
                    size_t prev_dist = dist - num_match;
                    size_t prev_query = query_dist - num_match;

                    // reconstruct the backwards traversal from last_node to node
                    auto it = query_window_rbegin + prev_query;
                    assert(it != query_window_rend);
                    size_t cur_path_size = path.size();
                    size_t num_traversed = 0;
                    if (std::get<3>(*bt) == Cigar::MISMATCH) {
                        assert(mismatch_char != *it);
                        traverse_node = traverse_back(traverse_node, mismatch_char, prev_dist, prev_query);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        spelling += mismatch_char;
                        ++it;
                        --num_match;
                        ++num_traversed;
                    } else {
                        assert(mismatch_char == *it);
                    }

                    while (num_match) {
                        assert(it != query_window_rend);
                        traverse_node = traverse_back(traverse_node, *it, prev_dist + num_traversed, prev_query + num_traversed);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        spelling += *it;
                        ++it;
                        --num_match;
                        ++num_traversed;
                    }
                    assert(traverse_node == node);

                    std::reverse(path.begin() + cur_path_size, path.end());
                    std::reverse(spelling.begin() + cur_path_size, spelling.end());

                    size_t prev_cost = cost - (last_op == Cigar::MATCH ? 0 : mismatch_cost);
                    assert(prev_cost <= cost);

                    const auto *check_bt = bt;
                    check_bt = nullptr;
                    if (prev_cost || prev_query || prev_dist) {
                        check_bt = get_bucket(S, prev_cost, prev_query, last_node);
                        assert(check_bt);
                        assert(std::get<0>(*check_bt) == prev_dist);
                        assert(std::get<3>(*check_bt) != Cigar::CLIPPED);
                    } else {
                        assert(last_node == start_node);
                    }

                    dist = prev_dist;
                    cost = prev_cost;
                    query_dist = prev_query;
                    node = last_node;

                    if (!check_bt)
                        break;

                    bt = check_bt;
                } else {
                    assert(false && "Different operation found");
                    throw std::runtime_error("Inf. loop");
                }

            } while (bt);

            assert(dist == 0);
            assert(cost == 0);
            assert(query_dist == 0);
            // std::cerr << "Ext: " << start_cost << "\t" << cigar.to_string() << "\n";
            callback(std::move(spelling), std::move(path), std::move(cigar), start_cost);
        }
    }
    common::logger->trace("Explored {} nodes", num_explored_nodes);
}

template <class Graph>
void align_bwd(const Graph &graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    align_impl(
        [&graph](DeBruijnGraph::node_index node, size_t, size_t) { return graph.has_single_incoming(node); },
        [&graph](DeBruijnGraph::node_index node, char c, size_t, size_t) { return graph.traverse_back(node, c); },
        [&graph](DeBruijnGraph::node_index node, const auto &callback, size_t, size_t) {
            graph.call_incoming_kmers(node, callback);
        },
        config,
        start_node,
        num_matches,
        query_window.begin(), query_window.end(),
        max_dist,
        callback,
        start_backtrack,
        terminate_branch,
        terminate
    );
}

template <class Graph>
void align_fwd(const Graph &graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate,
               std::string_view suffix = "") {
    align_impl(
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, size_t dist, size_t) {
            std::ignore = start_node;
            assert(dist >= suffix.size() || node == start_node);
            return dist < suffix.size() || graph.has_single_outgoing(node);
        },
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, char c, size_t dist, size_t) {
            if (dist < suffix.size()) {
                std::ignore = start_node;
                assert(node == start_node);
                return c == suffix[dist] ? node : DeBruijnGraph::npos;
            } else {
                return graph.traverse(node, c);
            }
        },
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, const auto &callback, size_t dist, size_t) {
            if (dist < suffix.size()) {
                std::ignore = start_node;
                assert(node == start_node);
                callback(node, suffix[dist]);
            } else {
                graph.call_outgoing_kmers(node, callback);
            }
        },
        config,
        start_node,
        num_matches,
        query_window.rbegin(), query_window.rend(),
        max_dist,
        [&](auto&& spelling, auto&& path, auto&& cigar, size_t cost) {
            if (path.size() >= suffix.size()) {
                path.resize(path.size() - suffix.size());
                std::reverse(spelling.begin(), spelling.end());
                std::reverse(path.begin(), path.end());
                std::reverse(cigar.data().begin(), cigar.data().end());
                callback(std::move(spelling), std::move(path), std::move(cigar), cost);
            }
        },
        start_backtrack,
        terminate_branch,
        terminate
    );
}

class AlignmentGraph {
    using node_index = DeBruijnGraph::node_index;
    using label_class_t = Anchor::label_class_t;

  public:
    AlignmentGraph(const DeBruijnGraph &graph)
          : graph_(graph), anno_buffer_(nullptr), target_(Anchor::nlabel) {}

    AlignmentGraph(const DeBruijnGraph &graph, AnnotationBuffer &anno_buffer, label_class_t target)
          : graph_(graph), anno_buffer_(&anno_buffer), target_(target) {}

    node_index traverse(node_index node, char c) const {
        node_index next = graph_.traverse(node, c);
        return has_labels(next) ? next : DeBruijnGraph::npos;
    }

    void traverse(node_index node, const char *begin, const char *end, const std::function<void(node_index)> &callback) const {
        if (target_ == Anchor::nlabel) {
            graph_.traverse(node, begin, end, callback);
            return;
        }

        for (auto it = begin; it != end && (node = traverse(node, *it)); ++it) {
            callback(node);
        }
    }

    bool has_single_outgoing(node_index node) const {
        if (target_ == Anchor::nlabel)
            return graph_.has_single_outgoing(node);

        size_t outdegree = 0;
        adjacent_outgoing_nodes(node, [&](node_index) { ++outdegree; });
        return outdegree == 1;
    }

    void call_outgoing_kmers(node_index node, const std::function<void(node_index, char)> &callback) const {
        if (target_ == Anchor::nlabel) {
            graph_.call_outgoing_kmers(node, callback);
        } else {
            graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
                if (has_labels(next))
                    callback(next, c);
            });
        }
    }

    void adjacent_outgoing_nodes(node_index node, const std::function<void(node_index)> &callback) const {
        if (target_ == Anchor::nlabel) {
            graph_.adjacent_outgoing_nodes(node, callback);
        } else {
            graph_.adjacent_outgoing_nodes(node, [&](node_index next) {
                if (has_labels(next))
                    callback(next);
            });
        }
    }

    node_index traverse_back(node_index node, char c) const {
        node_index prev = graph_.traverse_back(node, c);
        return has_labels(prev) ? prev : DeBruijnGraph::npos;
    }

    bool has_single_incoming(node_index node) const {
        // std::cerr << "hsi\t" << node << "\t" << target_ << std::endl;
        if (target_ == Anchor::nlabel) {
            return graph_.has_single_incoming(node);
        } else {
            size_t indegree = 0;
            adjacent_incoming_nodes(node, [&](node_index) { ++indegree; });
            return indegree == 1;
        }
    }

    void call_incoming_kmers(node_index node, const std::function<void(node_index, char)> &callback) const {
        // std::cerr << "cik\t" << node << "\t" << target_ << std::endl;
        if (target_ == Anchor::nlabel) {
            graph_.call_incoming_kmers(node, callback);
        } else {
            graph_.call_incoming_kmers(node, [&](node_index prev, char c) {
                if (has_labels(prev))
                    callback(prev, c);
            });
        }
    }

    void adjacent_incoming_nodes(node_index node, const std::function<void(node_index)> &callback) const {
        if (target_ == Anchor::nlabel) {
            graph_.adjacent_incoming_nodes(node, callback);
        } else {
            graph_.adjacent_incoming_nodes(node, [&](node_index prev) {
                if (has_labels(prev))
                    callback(prev);
            });
        }
    }

  private:
    const DeBruijnGraph &graph_;
    AnnotationBuffer* const anno_buffer_;
    label_class_t target_;

    bool has_labels(node_index node) const {
        if (!anno_buffer_ || target_ == Anchor::nlabel)
            return true;

        anno_buffer_->queue_path(std::vector<node_index>{ node });
        anno_buffer_->fetch_queued_annotations();

        const auto *node_labels = anno_buffer_->get_labels(node);
        assert(node_labels);

        const auto &target_labels = anno_buffer_->get_cached_column_set(target_);

        return utils::count_intersection(node_labels->begin(), node_labels->end(),
                                         target_labels.begin(), target_labels.end()) == target_labels.size();
    }
};

void Extender::extend(const Alignment &aln, const std::function<void(Alignment&&)> &callback, bool no_bwd, bool no_fwd) const {
    // std::cerr << "Base " << aln << "\n";
    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    std::vector<Alignment> fwd_exts;
    if (no_fwd || !aln.get_end_clipping()) {
        // std::cerr << "skip fwd\t" << aln << "\t" << aln.get_label_classes().back() << "\n";
        fwd_exts.emplace_back(aln);
    } else {
        // std::cerr << "fwd extending\t" << aln << "\t" << aln.get_label_classes().back() << "\n";
        std::string_view query_window = aln.get_query();
        query_window.remove_prefix(aln.get_clipping() + aln.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return aln.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.right_end_bonus : 0);
        };

        auto it = std::make_reverse_iterator(aln.get_cigar().data().end());
        auto end = std::make_reverse_iterator(aln.get_cigar().data().begin());
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;
        size_t max_dist = 0;
        DBGAlignerConfig::score_t best_score = aln.get_score();
        size_t max_query_dist = 0;
        align_fwd(
            make_aln_graph(aln.get_label_classes().back()), config_, aln.get_path().back(), num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& spelling, auto&& path, auto&& cigar, size_t) {
                auto ext_path = aln.get_path();
                ext_path.insert(ext_path.end(), path.begin(), path.end());

                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                cigar.append(Cigar::CLIPPED, query_window.size() - query_dist);

                Cigar ext_cigar = aln.get_cigar();
                ext_cigar.trim_end_clipping();
                ext_cigar.append(std::move(cigar));

                fwd_exts.emplace_back(
                    query_.get_graph(),
                    aln.get_query(),
                    aln.get_orientation(),
                    std::move(ext_path),
                    config_,
                    std::move(ext_cigar),
                    std::string(aln.get_spelling()) + spelling,
                    aln.get_end_trim() - std::min(aln.get_end_trim(), spelling.size()),
                    aln.get_label_classes().back()
                );

                // std::cerr << "\tend" << fwd_exts.back() << "\t" << fwd_exts.back().get_label_classes().back() << std::endl;

                // assert(fwd_exts.back().get_score() == best_score);
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if (!num_matches || dist < aln.get_end_trim() || (dist == 0 && query_dist == 0))
                // if (!num_matches || (dist == 0 && query_dist == 0))
                if (dist == 0 && query_dist == 0)
                    return false;

                // if (best_score == get_score(cost, dist, query_dist))
                    // std::cerr << "BT?: " << get_score(cost, dist, query_dist) << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                if (query_dist == query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                    // std::cerr << "TB?: " << score << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\n";
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            },
            aln.get_trim_spelling()
        );

        if (fwd_exts.empty())
            fwd_exts.emplace_back(aln);
    }

    for (auto &fwd_ext : fwd_exts) {
        fwd_ext.trim_end();
        if (no_bwd || !fwd_ext.get_clipping()) {
            // std::cerr << "\tskipbwd\t" << fwd_ext << "\n";
            callback(std::move(fwd_ext));
            continue;
        }

        std::string_view query_window = fwd_ext.get_query();
        query_window.remove_suffix(fwd_ext.get_end_clipping() + fwd_ext.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return fwd_ext.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.left_end_bonus : 0);
        };

        DBGAlignerConfig::score_t best_score = fwd_ext.get_score();

        auto it = aln.get_cigar().data().begin();
        auto end = aln.get_cigar().data().end();
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;

        bool found = false;
        size_t max_dist = 0;
        size_t max_query_dist = 0;
        // std::cerr << "bwd extending\t" << fwd_ext << "\t" << fwd_ext.get_label_classes()[0] << "\n";
        align_bwd(
            make_aln_graph(fwd_ext.get_label_classes()[0]), config_, fwd_ext.get_path()[0], num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& spelling, auto&& path, auto&& cigar, size_t) {
                path.insert(path.end(), fwd_ext.get_path().begin(), fwd_ext.get_path().end());
                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                Cigar ext_cigar(Cigar::CLIPPED, query_window.size() - query_dist);
                ext_cigar.append(std::move(cigar));

                Cigar aln_cigar = fwd_ext.get_cigar();
                aln_cigar.trim_clipping();
                ext_cigar.append(std::move(aln_cigar));

                Alignment next_aln(query_.get_graph(),
                    fwd_ext.get_query(),
                    fwd_ext.get_orientation(),
                    std::move(path),
                    config_,
                    std::move(ext_cigar),
                    spelling + fwd_ext.get_path_spelling(),
                    fwd_ext.get_end_trim(),
                    fwd_ext.get_label_classes()[0]
                );
                next_aln.trim_end();

                // assert(next_aln.get_score() == best_score);
                // std::cerr << "found\t" << next_aln << "\n";
                callback(std::move(next_aln));
                found = true;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if ((dist == 0 && query_dist == 0) || !num_matches)
                //     return false;

                if (dist == 0 && query_dist == 0)
                    return false;

                // if (best_score == get_score(cost, dist, query_dist))
                //     std::cerr << "BT?: " << get_score(cost, dist, query_dist) << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                if (query_dist == query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\t" << score << " vs. " << best_score << "\t" << Cigar::opt_to_char(std::get<3>(data)) << "\n";
                // return false;
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            }
        );

        if (!found)
            callback(std::move(fwd_ext));
    }
}

AlignmentGraph ExactSeeder::make_aln_graph(Anchor::label_class_t) const {
    return AlignmentGraph(query_.get_graph());
}

AlignmentGraph LabeledSeeder::make_aln_graph(Anchor::label_class_t target) const {
    return AlignmentGraph(query_.get_graph(), anno_buffer_, target);
}

AlignmentGraph Extender::make_aln_graph(Anchor::label_class_t) const {
    return AlignmentGraph(query_.get_graph());
}

AlignmentGraph LabeledExtender::make_aln_graph(Anchor::label_class_t target) const {
    return AlignmentGraph(query_.get_graph(), anno_buffer_, target);
}

std::vector<Alignment> ExactSeeder::get_inexact_anchors(bool align) const {
    common::logger->trace("Seeding");
    std::vector<Anchor> anchors = get_anchors();

    const auto &graph = query_.get_graph();

    auto full_it = std::partition(anchors.begin(), anchors.end(), [](const auto &a) { return a.get_end_trim(); });
    std::sort(full_it, anchors.end(), AnchorLess<Anchor>());
    if (anchors.end() - full_it > 1) {
        assert(full_it->get_spelling().size() + full_it->get_end_trim() == graph.get_k());
        for (auto it = full_it; it + 1 != anchors.end(); ++it) {
            auto &a_i = *it;
            auto &a_j = *(it + 1);
            assert(!a_i.get_end_trim());
            assert(!a_j.get_end_trim());
            assert(a_j.get_spelling().size() + a_j.get_end_trim() == graph.get_k());
            if (a_i.get_orientation() != a_j.get_orientation() || a_i.get_label_class() != a_j.get_label_class())
                continue;

            if (a_i.get_seed().size() < config_.max_seed_length
                    && a_i.get_end_clipping() == a_j.get_end_clipping() + 1
                    && graph.has_single_outgoing(a_i.get_path().back())
                    && graph.has_single_incoming(a_j.get_path()[0])) {
                // merge
                a_i.append(a_j, config_, &graph);
                std::swap(a_i, a_j);
                a_i = Anchor();
            }
        }

        full_it = std::remove_if(full_it, anchors.end(), [](const auto &a) { return a.get_path().empty(); });
        anchors.erase(full_it, anchors.end());
    }

    using AnchorIt = std::vector<Anchor>::iterator;
    std::sort(anchors.begin(), anchors.end(), AnchorLess<Anchor>());

    // std::cerr << "Anchors\n";
    // for (const auto &a : anchors) {
    //     std::cerr << "\t" << a << "\t" << a.get_path_spelling() << "\n";
    // }

    std::vector<Alignment> alignments;

    if (anchors.empty())
        return alignments;

    ssize_t k = graph.get_k();

    assert(!dynamic_cast<const LabeledSeeder*>(this) || std::all_of(anchors.begin(), anchors.end(),
                                          [](const auto &a) { return a.get_label_class() != Anchor::nlabel; }));

    // sdsl::bit_vector selected(anchors.size(), false);

    bool first_chain = true;
    DBGAlignerConfig::score_t first_chain_score = DBGAlignerConfig::ninf;
    size_t num_chains = 0;

    std::string dummy(query_.get_query().size(), '$');

    tsl::hopscotch_set<Anchor::label_class_t> found_labels;

    bool ext_success = false;
    Anchor::score_t last_chain_score = 0;
    common::logger->trace("Chaining");
    chain_anchors<AnchorIt>(query_, config_, anchors.begin(), anchors.end(),
        [this,k,&dummy](
                const Anchor &a_j,
                ssize_t max_dist,
                AnchorIt rbegin,
                AnchorIt rend,
                typename ChainScores<AnchorIt>::iterator rchain_scores,
                const ScoreUpdater<AnchorIt> &update_score) {
            auto rchain_scores_end = rchain_scores + (rend - rbegin);

            auto begin = std::make_reverse_iterator(rend);
            auto end = std::make_reverse_iterator(rbegin);
            auto chain_scores = std::make_reverse_iterator(rchain_scores_end);

            std::string_view query_j = a_j.get_seed();
            // const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));
            const DBGAlignerConfig::score_t &score_j = std::get<0>(*(rchain_scores_end));
            // const DBGAlignerConfig::score_t dist_j = std::get<2>(*(rchain_scores_end));

            auto try_chain = [&](auto it, auto chain_scores) {
                const Anchor &a_i = *it;
                assert(a_i.get_label_class() == a_j.get_label_class());
                assert(a_i.get_orientation() == a_j.get_orientation());

                auto graph = make_aln_graph(a_i.get_label_class());

                std::string_view query_i = a_i.get_seed();

                DBGAlignerConfig::score_t dist = query_j.end() - query_i.end();

                if (dist < 0 || query_i.begin() >= query_j.begin())
                    return;

                DBGAlignerConfig::score_t base_score = std::get<0>(*chain_scores);

                // std::cerr << "chh\t" << a_i << "\t" << a_i.get_path_spelling() << " -> " << a_j << "\t" << a_j.get_path_spelling() << "\t" << dist << "\n";
                if (dist == 0 && a_j.get_end_trim()) {
                    assert(a_j.get_path().size() == 1);
                    std::string_view a_j_trim_spelling = a_j.get_trim_spelling();
                    // std::cerr << "extleft" << a_i << "\t" << a_i.get_path_spelling() << " -> " << a_j << "\t" << a_j.get_path_spelling() << "\t" << a_j_trim_spelling.size() << "\n";
                    if (a_i.get_end_trim()) {
                        std::string_view a_i_trim_spelling = a_i.get_trim_spelling();
                        if (!std::equal(a_i_trim_spelling.begin(), a_i_trim_spelling.end(),
                                        a_j_trim_spelling.begin(), a_j_trim_spelling.begin() + a_i_trim_spelling.size()))
                            return;
                        a_j_trim_spelling = a_j_trim_spelling.substr(a_i_trim_spelling.size());
                    }
                    size_t traverse = 0;
                    DeBruijnGraph::node_index node = a_i.get_path().back();
                    graph.traverse(a_i.get_path().back(), a_j_trim_spelling.data(), a_j_trim_spelling.data() + a_j_trim_spelling.size(),
                                   [&](DeBruijnGraph::node_index next) {
                                       ++traverse;
                                       node = next;
                                   });
                    if (traverse == a_j_trim_spelling.size()) {
                        assert(node == a_j.get_path()[0]);
                        if (base_score > score_j) {
                            // std::cerr << "\tworked!\t" << base_score << "\n";
                            update_score(base_score, it, 0);
                        }
                    }
                }

                auto added_seq_start = std::max(query_j.begin(), query_i.end());
                std::string_view added_seq(added_seq_start, query_j.end() - added_seq_start);

                DBGAlignerConfig::score_t score = base_score + config_.match_score(added_seq)
                                                + (config_.right_end_bonus * !a_j.get_end_clipping());

                // if (score < score_j || (score == score_j && dist >= dist_j))
                if (score <= score_j)
                    return;
                // std::cerr << "\tchhs\n";

                size_t i_seed_size = k - a_i.get_end_trim();
                auto last_node_i_query = query_i.end() - i_seed_size;

                if (last_node_i_query >= query_j.begin()) {
                    // there should be overlapping nodes
                    // std::cerr << "olN" << a_i << " -> " << a_j << "\n";
                    size_t j_start = last_node_i_query - query_j.begin();
                    if (a_i.get_path().size() >= j_start
                            && a_j.get_path().size() >= j_start
                            && std::equal(a_j.get_path().begin(), a_j.get_path().begin() + j_start,
                                        a_i.get_path().end() - j_start, a_i.get_path().end())) {
                        // std::cerr << "\t\tworked! " << score << "\n";
                        update_score(score, it, dist);
                    }
                    return;
                } else if (query_i.end() > query_j.begin()) {
                    // overlapping nucleotides
                    // std::cerr << "oln" << a_i << " -> " << a_j << "\n";
                    assert(query_j.data() > query_i.data());
                    if (!a_i.get_end_trim() && !a_j.get_end_trim() && query_i.end() - query_j.begin() == k - 1) {
                        // no need to check
                        update_score(score, it, dist);
                    }
                    return;
                    // TODO: this check is expensive, skip it
                    // std::string_view left(query_i.data(), query_j.data() - query_i.data());
                    // DeBruijnGraph::node_index node = a_j.get_path()[0];
                    // ssize_t num_matches = 0;
                    // for (auto ait = left.rbegin(); ait != left.rend() && node != DeBruijnGraph::npos; ++ait) {
                    //     node = graph.traverse_back(node, *ait);
                    //     ++num_matches;
                    // }
                    // if (num_matches == static_cast<ssize_t>(left.size()) && node == a_i.get_path()[0]) {
                    //     assert(dist == num_matches - static_cast<ssize_t>(query_i.size()) + static_cast<ssize_t>(query_j.size()));
                    //     update_score(score, it, dist);
                    // }
                    // return;
                } else if (a_i.get_end_trim() && a_i.get_clipping() + a_i.get_seed().size() + a_i.get_end_trim() >= a_j.get_clipping()) {
                    // overlap in end parts
                    size_t olap = std::min(a_i.get_query().size(),
                                           a_i.get_clipping() + a_i.get_seed().size() + a_i.get_end_trim()) - a_j.get_clipping();
                    std::string_view a_i_trim = a_i.get_trim_spelling();
                    auto a_j_spelling = a_j.get_path_spelling();
                    assert(a_i_trim.size() >= olap);
                    for (size_t del = 0; del < olap; ++del) {
                        // std::cerr << "olap o:" << olap - del << "\td: " << del << "\t" << a_i << " " << a_i.get_path_spelling() << " -> " << a_j << " " << a_j.get_path_spelling() << std::endl;
                        if (std::equal(a_i_trim.end() - olap + del, a_i_trim.end(), a_j_spelling.begin(), a_j_spelling.begin() + olap - del)) {
                            // no need for indels
                            // std::cerr << "\tfound\n";
                            DeBruijnGraph::node_index node = a_j.get_path()[0];
                            size_t traversed = 0;
                            for (auto bit = a_i_trim.rbegin() + olap - del; bit != a_i_trim.rend() && (node = graph.traverse_back(node, *bit)); ++bit) {
                                ++traversed;
                            }
                            if (node != DeBruijnGraph::npos) {
                                // std::cerr << "\t\tstep2\n";
                                for (auto bit = a_i.get_seed().rbegin(); bit != a_i.get_seed().rend() && (node = graph.traverse_back(node, *bit)); ++bit) {
                                    ++traversed;
                                }
                                if (node == a_i.get_path().back()) {
                                    traversed += a_i.get_path().size() - 1;
                                    size_t nmismatch = query_j.begin() - query_i.end();
                                    assert(traversed >= del);
                                    assert(static_cast<ssize_t>(traversed - del) == query_j.begin() - query_i.begin());
                                    DBGAlignerConfig::score_t cur_score = score;
                                    DBGAlignerConfig::score_t cur_mismatch = config_.score_sequences(std::string_view(&*query_i.end(), nmismatch),
                                                                            std::string_view(dummy.c_str(), nmismatch));
                                    cur_score += cur_mismatch;
                                    if (del)
                                        cur_score += config_.gap_opening_penalty + (del - 1) * config_.gap_extension_penalty;
                                    if (cur_score > score_j) {
                                        // std::cerr << "\t\t\tworked! " << score << "\n";
                                        update_score(cur_score, it, traversed - query_i.size() + query_j.size());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    return;
                } else if (!a_i.get_end_trim()) {
                    // detect insertions
                    // std::cerr << "ins\tq: " << query_dist << "\td: " << traversed << "\tg: " << gap << "\t" << a_i << " " << a_i.get_path_spelling() << "\t" << a_j << " " << a_j.get_path_spelling() << "\n";
                    size_t traversed = 0;
                    DeBruijnGraph::node_index node = a_j.get_path()[0];
                    auto end = query_i.rend() - a_i.get_path().size() + 1;
                    for (auto ait = query_i.rbegin(); ait != end && (node = graph.traverse_back(node, *ait)); ++ait) {
                        ++traversed;
                    }
                    assert((node == DeBruijnGraph::npos && traversed < query_i.size()) || node == a_i.get_path().back());
                    if (traversed == query_i.size() - a_i.get_path().size() + 1) {
                        assert(node == a_i.get_path().back());
                        size_t query_dist = a_j.get_clipping() - (a_i.get_clipping() + a_i.get_path().size() - 1);
                        assert(query_dist >= traversed);
                        size_t gap = query_dist - traversed;
                        if (gap > 0)
                            score += config_.gap_opening_penalty + (gap - 1) * config_.gap_extension_penalty;

                        if (score > score_j) {
                            // std::cerr << "\t\tworked! " << score << "\n";
                            update_score(score, it, traversed - (query_i.size() - a_i.get_path().size() + 1) + query_j.size());
                        }
                    }
                }

                // TODO: if a_i.get_end_trim() != 0, then we need to deal with
                //       a combination of D from the trim, X to the trim, then I of the query
            };

            // for (auto it = rbegin; it != rend; ++it, ++rchain_scores) {
            //     try_chain(it, rchain_scores);
            // }

            // std::cerr << "round 1\n";
            for (auto rit = begin; rit != end; ++rit, ++chain_scores) {
                try_chain(rit.base() - 1, chain_scores);
            }

            // std::cerr << "round 2\n";

        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            assert(chain.size());

            ext_success = false;

            bool ret_val = (!chain[0].first->get_clipping() && !chain[0].first->get_end_clipping())
                            || chain.size() > 1
                            || std::any_of(chain.begin(), chain.end(),
                                           [](const auto &a) { return a.first->get_path().size() > 1; });

            // ret_val &= std::all_of(chain.begin(), chain.end(),
            //                        [&](const auto &a) { return !selected[a.first - anchors.begin()]
            //                                                     && (a.first->get_label_class() == Anchor::nlabel
            //                                                             || !found_labels.count(a.first->get_label_class())); });

            ret_val &= std::all_of(chain.begin(), chain.end(),
                                   [&](const auto &a) { return a.first->get_label_class() == Anchor::nlabel
                                                                        || !found_labels.count(a.first->get_label_class()); });

            if (!ret_val)
                return false;

            // common::logger->info("Chain\t{}\t{}", score_traceback.back(), ret_val);
            // for (const auto &[it, dist] : chain) {
            //     std::cerr << "\t" << *it << "\t" << dist << "\t" << it->get_label_class() << "\t" << it->get_path_spelling() << "\n";
            // }

            last_chain_score = score_traceback.back();
            ret_val = score_traceback.back() >= first_chain_score * config_.rel_score_cutoff;

            // if (ret_val) {
            //     ++num_chains;
            //     for (const auto &[it, dist] : chain) {
            //         // selected[it - anchors.begin()] = true;
            //         if (it->get_label_class() != Anchor::nlabel)
            //             found_labels.emplace(it->get_label_class());
            //     }
            // }

            return ret_val;
        },
        [this,align,match_score=config_.match_score("A")](
                AnchorIt last,
               AnchorIt next,
               Alignment&& aln,
               size_t next_to_last_dist,
               DBGAlignerConfig::score_t score_up_to_now,
               const AlignmentCallback &callback) {
            size_t traversal_dist = next_to_last_dist - last->get_seed().size() + next->get_seed().size();
            const auto &graph = query_.get_graph();
            if (next->get_seed().begin() + 1 == last->get_seed().begin() && traversal_dist == 1) {
                // simple extension
                // std::cerr << "Simple extension: " << *next << "\t->\t" << aln << std::endl;
#ifndef NDEBUG
                size_t k = graph.get_k();
                if (graph.traverse(next->get_path()[0], aln.get_path_spelling()[k-1]) != aln.get_path()[0]) {
                    std::cerr << "\t" << next->get_path()[0] << " > " << aln.get_path_spelling()[k-1] << "\t" << graph.traverse(next->get_path()[0], aln.get_path_spelling()[k-1]) << " vs. " << aln.get_path()[0] << std::endl;
                    assert(false);
                }
#endif
                // common::logger->info("Simple extension {} -> {}", *next, aln);
                std::vector<DeBruijnGraph::node_index> path;
                path.resize(aln.get_path().size() + 1);
                assert(next->get_path().size() == 1);
                path[0] = next->get_path()[0];
                std::copy(aln.get_path().begin(), aln.get_path().end(), path.begin() + 1);

                Cigar cigar(Cigar::CLIPPED, next->get_clipping());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();

                cigar.append(Cigar::MATCH, 1);
                cigar.append(std::move(aln_cigar));

                Alignment next_aln(graph,
                                aln.get_query(),
                                aln.get_orientation(),
                                std::move(path),
                                config_,
                                std::move(cigar),
                                std::string(1, next->get_seed()[0]) + aln.get_path_spelling(),
                                aln.get_end_trim(),
                                next->get_label_class());
                // std::cerr << "\tnew aln: " << next_aln << std::endl;
                callback(std::move(next_aln));
                return;
            }

            if (!align) {
                // connect anchors without actually aligning or traversing
                size_t query_dist = last->get_clipping() - next->get_clipping();
                std::vector<DeBruijnGraph::node_index> path;
                path.resize(aln.get_path().size() + traversal_dist);
                std::copy(aln.get_path().rbegin(), aln.get_path().rend(), path.rbegin());

                Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                std::string spelling;

                if (next->get_seed().end() > last->get_seed().begin()) {
                    // overlap
                    // std::cerr << "foo\t" << next->get_seed().end() - last->get_seed().begin() << "\t" << *next << " " << next->get_path_spelling() << " -> " << *last << " " << last->get_path_spelling() << std::endl;
                    assert(traversal_dist == query_dist);
                    // TODO: don't bother traversing, trust that the connection worked.
                    //       keep intermediate nodes as DeBruijnGraph::npos
#ifndef NDEBUG
                    DeBruijnGraph::node_index node = *(path.rbegin() + aln.get_path().size() - 1);
                    // auto jt = path.rbegin() + aln.get_path().size();
                    // for (auto it = next->get_seed().rend() - query_dist; it != next->get_seed().rend(); ++it, ++jt) {
                    for (auto it = next->get_seed().rend() - query_dist; it != next->get_seed().rend(); ++it) {
                        // *jt = graph.traverse_back(*(jt - 1), *it);
                        // assert(graph.in_graph(*jt));
                        node = graph.traverse_back(node, *it);
                        assert(graph.in_graph(node));
                    }
                    assert(node == next->get_path()[0]);
#endif
                    // assert(path[0] == next->get_path()[0]);
                    cigar.append(Cigar::MATCH, query_dist);
                    assert(next->get_query().begin() + next->get_clipping() + query_dist <= next->get_query().end());
                    spelling += std::string(next->get_seed().data(), next->get_seed().data() + query_dist);
                } else {
                    cigar.append(Cigar::MATCH, next->get_seed().size());
                    spelling += next->get_seed();
                    if (traversal_dist == query_dist) {
                        // no gap
                        assert(query_dist > next->get_seed().size());
                        cigar.append(Cigar::MISMATCH, query_dist - next->get_seed().size());
                        // assert(query_dist - next->get_seed().size() <= next->get_end_trim());
                        // spelling += next->get_trim_spelling().substr(0, query_dist - next->get_seed().size());
                        spelling += std::string(query_dist - next->get_seed().size(), '$');
                    } else if (traversal_dist > query_dist) {
                        // deletion
                        // TODO: which first, mismatch or deletion?
                        cigar.append(Cigar::DELETION, traversal_dist - query_dist);
                        spelling += std::string(traversal_dist - query_dist, '$');
                        if (query_dist > next->get_seed().size()) {
                            cigar.append(Cigar::MISMATCH, query_dist - next->get_seed().size());
                            spelling += std::string(query_dist - next->get_seed().size(), '$');
                        }
                    } else {
                        // insertion
                        // TODO: which first, mismatch or insert?
                        cigar.append(Cigar::INSERTION, query_dist - traversal_dist);
                        if (traversal_dist > next->get_seed().size()) {
                            cigar.append(Cigar::MISMATCH, traversal_dist - next->get_seed().size());
                            spelling += std::string(traversal_dist - next->get_seed().size(), '$');
                        }
                    }

                    if (next->get_end_trim()) {
                        assert(next->get_path().size() == 1);
                        std::string_view next_trim_spelling = next->get_trim_spelling();
                        auto end = next_trim_spelling.begin() + std::min(next_trim_spelling.size(), spelling.size() - next->get_seed().size());
                        assert(next->get_seed().size() + (end - next_trim_spelling.begin()) <= spelling.size());
                        std::copy(next_trim_spelling.begin(), end,
                                  spelling.begin() + next->get_seed().size());
                    }
                }

                std::copy(next->get_path().begin(), next->get_path().end(), path.begin());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();
                cigar.append(std::move(aln_cigar));

                Alignment next_aln(graph,
                                aln.get_query(),
                                aln.get_orientation(),
                                std::move(path),
                                config_,
                                std::move(cigar),
                                spelling + aln.get_path_spelling(),
                                aln.get_end_trim(),
                                next->get_label_class());
                callback(std::move(next_aln));
                return;
            }

            assert(traversal_dist > 0);

            std::string_view query_window(
                next->get_seed().data(),
                aln.get_seed().begin() - next->get_seed().begin()
            );

            bool aln_found = false;

            auto start_backtracking = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                // common::logger->info("Checkbt: n: {}, d: {}, qd: {} / {}, s: {}, nm: {}",
                //                      query_.get_graph().get_node_sequence(node), dist, query_dist, query_window.size(),
                //                      get_score(cost, query_dist, dist), num_matches);

                if (aln_found || dist == 0 || query_dist == 0 || num_matches < next->get_seed().size())
                    return false;

                return query_dist == query_window.size() && node == next->get_path()[0];
            };

            auto terminate_branch = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                // common::logger->info("Check: n: {}, d: {}, qd: {} / {}, s: {}",
                //                      query_.get_graph().get_node_sequence(node), dist, query_dist, query_window.size(),
                //                      get_score(cost, query_dist, dist));

                if (dist == 0 && query_dist == 0)
                    return false;

                if (query_dist == query_window.size()) {
                    aln_found |= start_backtracking(cost, data, query_dist, node);
                    return true;
                }

                if (query_window.size() - query_dist <= next->get_seed().size() && num_matches < next->get_seed().size()) {
                    // sanity checks once we've hit the target seed

                    // cut off the branch if we did not match some of the target seed
                    // characters along the way
                    if (query_window.size() - query_dist < next->get_seed().size() - num_matches) {
                        // std::cerr << "ttt: " << query_window.size() << "," << query_dist << "\t" << next->get_seed().size() << "," << num_matches << "\n";
                        return true;
                    }

                    // cut off the branch if we've hit the target seed on an incorrect
                    // diagonal (i.e., not enough indels)
                    assert(traversal_dist >= dist);
                    if (query_window.size() - query_dist != traversal_dist - dist)
                        return true;
                }

                return false;
            };


            align_bwd(
                make_aln_graph(last->get_label_class()),
                config_,
                aln.get_path()[0],
                last->get_seed().size(),
                query_window,
                traversal_dist,
                [&](auto&& spelling, auto&& bt_path, auto&& bt_cigar, size_t cost) { // callback
                    assert(bt_path[0] == next->get_path()[0]);

                    Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                    cigar.append(std::move(bt_cigar));
                    // std::cerr << "\textension: " << *next << "\t" << cigar.to_string() << "\t" << aln << std::endl;

                    Cigar aln_cigar = aln.get_cigar();
                    aln_cigar.trim_clipping();
                    cigar.append(std::move(aln_cigar));
                    bt_path.insert(bt_path.end(), aln.get_path().begin(), aln.get_path().end());

                    Alignment next_aln(query_.get_graph(),
                                    aln.get_query(),
                                    aln.get_orientation(),
                                    std::move(bt_path),
                                    config_,
                                    std::move(cigar),
                                    spelling + aln.get_path_spelling(),
                                    aln.get_end_trim(),
                                    last->get_label_class());
                    // std::cerr << "\t\t" << next_aln << std::endl;
                    aln_found = true;
                    callback(std::move(next_aln));
                },
                start_backtracking,
                terminate_branch,
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) { // terminate
                    std::ignore = cost;
                    std::ignore = data;
                    std::ignore = query_dist;
                    // const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                    // common::logger->info("Checkt: n: {}, d: {}, qd: {} / {}, s: {}",
                    //                  graph.get_node_sequence(node), dist, query_dist, query_window.size(),
                    //                  get_score(cost, query_dist, dist));
                    return aln_found || (query_dist == query_window.size() && node == next->get_path()[0]);
                }
            );
            assert(aln_found);
        },
        [&num_chains,&alignments,&ext_success,&last_chain_score,&first_chain,&first_chain_score,&found_labels](Alignment&& aln) {
            aln.trim_end();
            if (!ext_success) {
                ++num_chains;
                ext_success = true;
                if (first_chain) {
                    first_chain_score = last_chain_score;
                    first_chain = false;
                }
                for (auto label : aln.get_label_classes()) {
                    if (label != Anchor::nlabel)
                        found_labels.emplace(label);
                }
            }
            // std::cerr << "\tInit aln\t" << aln << "\t" << aln.get_label_classes()[0] << std::endl;
            alignments.emplace_back(std::move(aln));
        }
    );

    // size_t num_selected = sdsl::util::cnt_one_bits(selected);
    // common::logger->info("Selected {} / {} anchors from {} chains. Produced {} alignments",
    //                         num_selected,
    //                         selected.size(),
    //                         num_chains, alignments.size());

    return alignments;
}


} // namespace mtg::graph::align
