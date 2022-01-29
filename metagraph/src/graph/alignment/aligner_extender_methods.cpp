#include "aligner_extender_methods.hpp"

#include <map>

#include <tsl/hopscotch_set.h>

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/rc_dbg.hpp"


namespace mtg {
namespace graph {
namespace align {

using score_t = DBGAlignerConfig::score_t;
constexpr score_t ninf = DBGAlignerConfig::ninf;

GreedyExtender::GreedyExtender(const DeBruijnGraph &graph,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               const std::vector<Alignment> &seeds)
      : graph_(&graph), config_(config), query_(query) {
    std::ignore = seeds;
    if (config_.alignment_match_score % 2) {
        throw std::runtime_error("Only supported with even match score");
        exit(1);
    }
    config_.alignment_mm_transition_score = config_.alignment_mm_transversion_score;
    config_.gap_extension_penalty = -config_.alignment_mm_transition_score
        - config_.alignment_match_score / 2;
    config_.gap_opening_penalty = config_.gap_extension_penalty;

    for (const auto &seed : seeds) {
        if (!seed.empty())
            seeds_[std::make_pair(seed.get_nodes().back(), seed.get_query().size() - 1 + seed.get_cigar().get_clipping())] = seed.get_query().size();
    }
}

bool GreedyExtender::set_seed(const Alignment &seed) {
    seed_ = seed;
    return true;
}

bool GreedyExtender::terminate(uint64_t query_i, uint64_t) const {
    return query_i == query_.size();
}

bool GreedyExtender::check_seed(const Alignment &seed) const {
    return !std::all_of(seed.get_nodes().begin(), seed.get_nodes().end(),
                        [this](const auto &node) { return nodes_.count(node); });
}

std::vector<Alignment> GreedyExtender::extend(score_t min_path_score, bool force_fixed_seed) {
    std::ignore = min_path_score;
    std::ignore = force_fixed_seed;
    ++num_extensions_;
    num_explored_nodes_ += seed_.get_nodes().size();

    using dist_t = uint32_t;

    std::string_view window = query_.substr(seed_.get_cigar().get_clipping());
    score_t hmat = config_.alignment_match_score / 2;
    score_t dmatmis = config_.alignment_match_score + config_.alignment_mm_transition_score;
    dist_t d_offset = (config_.xdrop + hmat) / dmatmis + 1;

    auto get_score = [&](score_t ij, score_t d) -> score_t {
        return ij * hmat - d * dmatmis;
    };


    using diag_t = int32_t;
    using aln_tree_node_t = uint32_t;
    using limits_t = std::tuple<diag_t, diag_t, diag_t>;
    std::vector<std::tuple<node_index, aln_tree_node_t, dist_t, char>> aln_tree;
    aln_tree.emplace_back(seed_.get_nodes().back(),
                          std::numeric_limits<aln_tree_node_t>::max(),
                          seed_.get_query().size() - 1,
                          seed_.get_query().back());

    diag_t max_i = window.size();

    using diag_v_t = std::tuple<limits_t, diag_t, std::deque<diag_t>, bool>;

    using score_v = std::pair<score_t, diag_v_t>;
    std::vector<std::map<aln_tree_node_t, score_v>> RT;
    score_t best_score = config_.ninf;

    using queue_t = std::vector<std::tuple<bool, score_t, aln_tree_node_t>>;
    queue_t queue;

    RT.emplace_back();
    std::vector<aln_tree_node_t> traversal_stack;
    traversal_stack.emplace_back(0);
    while (traversal_stack.size()) {
        aln_tree_node_t cur = traversal_stack.back();
        traversal_stack.pop_back();
        auto [node, parent, dist, c_prev] = aln_tree[cur];
        ++dist;
        if (dist == window.size()) {
            score_t score = get_score(dist + dist, 0);
            best_score = std::max(best_score, score);
            if (terminate(dist, dist)) {
                queue.clear();
                break;
            }

            auto [it, inserted] = RT[0].emplace(cur, score_v{});
            it->second.first = score;
            auto &[limits, k_min, dists, has_seed] = it->second.second;
            has_seed = true;
            k_min = 0;
            dists.emplace_back(dist);
            limits = limits_t{ 0, 0, 0 };
            continue;
        }

        graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
            ++num_explored_nodes_;
            aln_tree_node_t aln_tree_next = aln_tree.size();
            aln_tree.emplace_back(next, cur, dist, c);
            if (c == window[dist]) {
                size_t global_query_i = dist + seed_.get_cigar().get_clipping();
                auto find_seed = seeds_.find(std::make_pair(next, global_query_i));
                if (find_seed != seeds_.end()) {
                    // common::logger->trace(
                    //     "found next seed: {}, {}\t{}-{}",
                    //     next, find_seed->second,
                    //     global_query_i + 1 - find_seed->second,
                    //     global_query_i + 1
                    // );
                    RT[0].clear();
                    queue.clear();
                }
                traversal_stack.emplace_back(aln_tree_next);
            } else {
                // std::cerr << "foo\t" << dist << "\n";
                auto [it, inserted] = RT[0].emplace(aln_tree_next, score_v{});
                it->second.first = get_score(dist + dist, 0);
                best_score = std::max(best_score, it->second.first);
                auto &[limits, k_min, dists, has_seed] = it->second.second;
                has_seed = true;
                k_min = 0;
                dists.emplace_back(dist);
                limits = limits_t{ 0, 0, 0 };
                queue.emplace_back(true, it->second.first, aln_tree_next);
            }
        });
    }

    for (dist_t d = 1; queue.size(); ++d) {
        assert(RT.size());
        std::make_heap(queue.begin(), queue.end());
        queue_t next_queue;
        RT.emplace_back();
        while (queue.size()) {
            std::pop_heap(queue.begin(), queue.end());
            auto [last_found_seed, score_prev, cur] = queue.back();
            queue.pop_back();

            auto it = RT[d - 1].find(cur);
            assert(it != RT[d - 1].end());
            assert(score_prev == it->second.first);
            // const auto &score_prev = it->second.first;
            // if (score_p < score_prev)
            //     continue;

            // std::cerr << "foo\t" << d << "\t" << score_prev << "\t\t" << num_explored_nodes_ << "\n";

            const auto &[limits_prev, k_prev, vals_prev, has_seed_prev] = it->second.second;
            const auto &[L, U, Uv2] = limits_prev;
            if (L > std::min(U, Uv2))
                continue;

            score_t xdrop_cutoff = config_.ninf;
            for (const auto &[n, sv] : RT[d - std::min(d, d_offset)]) {
                xdrop_cutoff = std::max(xdrop_cutoff, sv.first - config_.xdrop);
            }

            std::vector<std::pair<diag_t, diag_t>> to_update;

            for (diag_t k = L - 1; k <= U + 1; ++k) {
                diag_t query_i = config_.ninf;
                if (L < k && vals_prev[k - k_prev - 1] < max_i && vals_prev[k - k_prev - 1] + 1 > query_i) {
                    query_i = vals_prev[k - k_prev - 1] + 1;
                }

                if (L <= k && k <= U && vals_prev[k - k_prev] < max_i && vals_prev[k - k_prev] + 1 > query_i) {
                    query_i = vals_prev[k - k_prev] + 1;
                }

                if (k < U && vals_prev[k - k_prev + 1] <= max_i && vals_prev[k - k_prev + 1] > query_i) {
                    query_i = vals_prev[k - k_prev + 1];
                }

                assert(query_i <= config_.ninf + 1 || query_i >= 0);

                if (query_i <= config_.ninf + 1 || get_score(query_i + query_i - k, d) < xdrop_cutoff)
                    continue;

                to_update.emplace_back(query_i - k - std::get<2>(aln_tree[cur]), k);
            }

            if (to_update.empty())
                continue;

            std::sort(to_update.begin(), to_update.end());
            std::vector<std::vector<aln_tree_node_t>> local_tree;
            local_tree.emplace_back();
            local_tree.back().emplace_back(it->first);
            size_t start = 0;
            if (to_update.front().first < 0) {
                for (diag_t i = 0; i < -to_update.front().first; ++i) {
                    if (local_tree.back().empty()) {
                        local_tree.emplace_back();
                        continue;
                    }

                    local_tree.emplace_back();
                    local_tree.back().emplace_back(std::get<1>(aln_tree[local_tree[local_tree.size() - 2].back()]));
                    if (local_tree.back().back() == std::numeric_limits<aln_tree_node_t>::max())
                        local_tree.back().clear();
                }
                std::reverse(local_tree.begin(), local_tree.end());
                start = local_tree.size() - 1;
            }
            assert(local_tree[start].size() == 1);
            assert(local_tree[start][0] == it->first);
            for (diag_t i = 0; i < to_update.back().first; ++i) {
                local_tree.emplace_back();
                for (aln_tree_node_t ref_j : local_tree[i + start]) {
                    auto [node, parent, dist, c_prev] = aln_tree[ref_j];
                    ++dist;
                    graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
                        ++num_explored_nodes_;
                        local_tree.back().emplace_back(aln_tree.size());
                        aln_tree.emplace_back(next, ref_j, dist, c);
                    });
                }
            }
            diag_t local_tree_offset = std::min(std::get<0>(to_update[0]), 0);

            tsl::hopscotch_set<aln_tree_node_t> pushed;

            for (auto [ref_offset, k] : to_update) {
                diag_t query_i = ref_offset + k + std::get<2>(aln_tree[cur]);
                assert(query_i >= k);
                assert(static_cast<dist_t>(ref_offset - local_tree_offset) < local_tree.size());

                auto push_score = [&](dist_t query_i, aln_tree_node_t ref_j, bool found_seed) {
                    auto [it, inserted] = RT[d].emplace(
                        ref_j,
                        score_v{ score_prev, diag_v_t{} }
                    );
                    score_t score = get_score(query_i + std::get<2>(aln_tree[ref_j]), d);
                    pushed.emplace(ref_j);

                    it->second.first = std::max(it->second.first, score);
                    best_score = std::max(best_score, score);

                    auto &[newlimits, k_min, dists, has_seed] = it->second.second;
                    has_seed |= found_seed;

                    if (inserted || has_seed) {
                        k_min = k;
                        dists.clear();
                        dists.emplace_back(query_i);
                    } else {
                        if (k >= k_min) {
                            // append to the end of the array
                            dists.resize(std::max(
                                dists.size(), static_cast<size_t>(k - k_min + 1)
                            ), config_.ninf);
                            dists[k - k_min] = query_i;
                        } else {
                            dists.insert(dists.begin(), k_min - k, config_.ninf);
                            k_min = k;
                            dists[0] = query_i;
                        }
                    }
                    auto &[newL, newU, newUv2] = newlimits;
                    if (inserted || has_seed) {
                        newL = k;
                        newU = k;
                    } else {
                        newL = std::min(newL, k);
                        newU = std::max(newU, k);
                    }

                    if (query_i == window.size()) {
                        if (inserted || has_seed) {
                            newUv2 = k - 2;
                        } else {
                            newUv2 = std::min(newUv2, k - 2);
                        }
                    }
                };

                tsl::hopscotch_map<aln_tree_node_t, bool> prevs;
                for (auto ref_j : local_tree[ref_offset - local_tree_offset]) {
                    prevs.emplace(std::get<1>(aln_tree[ref_j]), false);
                }
                bool reached_end = false;
                dist_t num_exact_match = 1;
                for (dist_t j = ref_offset - local_tree_offset; !reached_end; ++j, ++query_i, ++num_exact_match) {
                    assert(query_i <= max_i);
                    assert(j <= local_tree.size());
                    if (aln_tree.size() > config_.max_nodes_per_seq_char * window.size()) {
                        pushed.clear();
                        next_queue.clear();
                        reached_end = true;
                        break;
                    }

                    if (j == local_tree.size()) {
                        local_tree.emplace_back();
                        for (auto ref_j : local_tree[j - 1]) {
                            auto [node, parent, dist, c_prev] = aln_tree[ref_j];

                            ++dist;
                            graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
                                ++num_explored_nodes_;
                                aln_tree_node_t aln_tree_next = aln_tree.size();
                                aln_tree.emplace_back(next, ref_j, dist, c);
                                local_tree.back().emplace_back(aln_tree_next);
                            });
                        }
                    }

                    tsl::hopscotch_map<aln_tree_node_t, bool> prevs_next;
                    reached_end = true;
                    for (auto ref_j : local_tree[j]) {
                        assert(std::get<2>(aln_tree[ref_j]) == static_cast<dist_t>(query_i - k));
                        auto [node, parent, dist, c] = aln_tree[ref_j];
                        auto pt = prevs.find(parent);
                        if (pt == prevs.end())
                            continue;

                        bool found_seed = pt->second;

                        if (query_i == max_i) {
                            push_score(query_i, ref_j, found_seed);
                            if (terminate(query_i, ref_j)) {
                                pushed.clear();
                                queue.clear();
                                break;
                            }
                        } else if (window[query_i] != c) {
                            push_score(query_i, ref_j, found_seed);
                            if (found_seed) {
                                pushed.clear();
                                queue.clear();
                                pushed.emplace(ref_j);
                                break;
                            }
                        } else {
                            if (!found_seed) {
                                size_t global_query_i = query_i + seed_.get_cigar().get_clipping();
                                auto find_seed = seeds_.find(std::make_pair(node, global_query_i));
                                found_seed |= find_seed != seeds_.end() && find_seed->second <= num_exact_match;
                                // if (found_seed) {
                                //     common::logger->trace(
                                //         "found next seed, terminating: {}, {}\t{}-{}",
                                //         node, find_seed->second,
                                //         global_query_i + 1 - find_seed->second,
                                //         global_query_i + 1
                                //     );
                                // }
                            }
                            reached_end = false;
                            prevs_next.emplace(ref_j, found_seed);
                        }
                    }

                    std::swap(prevs, prevs_next);
                }

                if (aln_tree.size() > config_.max_nodes_per_seq_char * window.size()) {
                    pushed.clear();
                    next_queue.clear();
                    reached_end = true;
                    break;
                }
            }

            for (auto ref_j : pushed) {
                const auto &bucket = RT[d][ref_j];
                next_queue.emplace_back(std::get<3>(bucket.second), bucket.first, ref_j);
            }
        }

        std::swap(queue, next_queue);

        if (aln_tree.size() > config_.max_nodes_per_seq_char * window.size()) {
            common::logger->trace("Explored too much, giving up");
            break;
        }
    }

    // common::logger->trace("Best score: {}, num_explored: {}, num_extensions: {}",
    //                       best_score, num_explored_nodes_, num_extensions_);

    for (auto&& [node, parent, dist, c] : aln_tree) {
        nodes_.emplace(node);
    }

    return {};
}

// to ensure that SIMD operations on arrays don't read out of bounds
constexpr size_t kPadding = 5;


DefaultColumnExtender::DefaultColumnExtender(const DeBruijnGraph &graph,
                                             const DBGAlignerConfig &config,
                                             std::string_view query)
      : SeedFilteringExtender(config, query), graph_(&graph), query_(query) {
    // compute exact-match scores for all suffixes of the query
    partial_sums_.reserve(query_.size() + 1);
    partial_sums_.resize(query_.size(), 0);
    std::transform(query_.begin(), query_.end(), partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query_) == partial_sums_.front());
    assert(config_.get_row(query_.back())[query_.back()] == partial_sums_.back());
    partial_sums_.push_back(0);

    // precompute profiles to store match/mismatch scores and Cigar::Operators
    // in contiguous arrays
    for (char c : graph_->alphabet()) {
        auto &p_score_row = profile_score_.emplace(c, query_.size() + kPadding).first.value();
        auto &p_op_row = profile_op_.emplace(c, query_.size() + kPadding).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = kCharToOp[c];

        // the first cell in a DP table row is one position before the first character,
        // so we need to shift the indices of profile_score_ and profile_op_
        std::transform(query_.begin(), query_.end(), p_score_row.begin() + 1,
                       [&row](char q) { return row[q]; });

        std::transform(query_.begin(), query_.end(), p_op_row.begin() + 1,
                       [&op_row](char q) { return op_row[q]; });
    }
}

bool SeedFilteringExtender::check_seed(const Alignment &seed) const {
    if (seed.empty())
        return false;

    assert(seed.get_nodes().back());
    assert(seed.get_cigar().size());

    node_index node = seed.get_nodes().back();
    if (dynamic_cast<const RCDBG*>(&get_graph()))
        node += get_graph().max_index();

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end())
        return true;

    size_t pos = seed.get_query().size() + seed.get_clipping() - 1;
    const auto &[start, vec] = it->second;
    if (pos < start || pos - start >= vec.size()
            || vec[pos - start] < seed.get_score()) {
        return true;
    }

    return false;
}

bool SeedFilteringExtender::set_seed(const Alignment &seed) {
    seed_ = &seed;
    clear_conv_checker();
    return seed_;
}

score_t SeedFilteringExtender::update_seed_filter(node_index node,
                                                  size_t query_start,
                                                  const score_t *s_begin,
                                                  const score_t *s_end) {
    if (dynamic_cast<const RCDBG*>(&get_graph()))
        node += get_graph().max_index();

    assert(s_end >= s_begin);
    assert(query_start + (s_end - s_begin) <= query_size_);

    size_t size = s_end - s_begin;

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end()) {
        conv_checker_.emplace(node, ScoreVec(query_start, { s_begin, s_end }));
        return *std::max_element(s_begin, s_end);
    }

    auto &[start, vec] = it.value();
    if (query_start + size <= start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        std::copy(s_begin, s_end, vec.begin());
        start = query_start;
        return *std::max_element(s_begin, s_end);
    }

    if (query_start >= start + vec.size()) {
        vec.reserve(query_start + size - start);
        vec.insert(vec.end(), query_start - start - vec.size(), ninf);
        vec.insert(vec.end(), s_begin, s_end);
        return *std::max_element(s_begin, s_end);
    }

    // overlap
    if (query_start < start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        start = query_start;
    }

    if (query_start + size > start + vec.size())
        vec.resize(query_start + size - start, ninf);

    score_t max_changed_value = ninf;
    score_t *v = vec.data() + query_start - start;
    for (size_t j = 0; j < size; ++j) {
        if (s_begin[j] > v[j] * config_.rel_score_cutoff) {
            v[j] = std::max(v[j], s_begin[j]);
            max_changed_value = std::max(max_changed_value, v[j]);
        }
    }

    return max_changed_value;
}

bool SeedFilteringExtender
::filter_nodes(node_index node, size_t query_start, size_t query_end) {
    assert(query_end >= query_start);
    assert(query_end <= query_size_);
    constexpr score_t mscore = -ninf;
    size_t size = query_end - query_start;

    auto it = conv_checker_.find(node);
    if (it == conv_checker_.end()) {
        conv_checker_.emplace(
            node, ScoreVec(query_start, AlignedVector<score_t>(size, mscore))
        );
        return true;
    }

    auto &[start, vec] = it.value();
    if (query_start + size <= start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        std::fill(vec.begin(), vec.begin() + size, mscore);
        start = query_start;
        return true;
    }

    if (query_start >= start + vec.size()) {
        vec.reserve(query_start + size - start);
        vec.insert(vec.end(), query_start - start - vec.size(), ninf);
        vec.insert(vec.end(), size, mscore);
        return true;
    }

    // overlap
    if (query_start < start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        start = query_start;
    }

    if (query_start + size > start + vec.size())
        vec.resize(query_start + size - start, ninf);

    bool converged = true;
    score_t *v = vec.data() + query_start - start;
    for (size_t j = 0; j < size; ++j) {
        if (mscore > v[j]) {
            converged = false;
            v[j] = mscore;
        }
    }

    return !converged;
}

void update_column(size_t prev_end,
                   const score_t *S_prev_v,
                   const score_t *F_prev_v,
                   AlignedVector<score_t> &S_v,
                   AlignedVector<score_t> &E_v,
                   AlignedVector<score_t> &F_v,
                   const score_t *profile_scores,
                   score_t xdrop_cutoff,
                   const DBGAlignerConfig &config_,
                   score_t init_score,
                   size_t offset) {
#ifndef __SSE4_1__
    for (size_t j = 0; j < prev_end; ++j) {
        score_t match = j ? (S_prev_v[j - 1] + profile_scores[j] + init_score) : ninf;
        if (offset > 1) {
            F_v[j] = std::max(S_prev_v[j] + init_score + config_.gap_opening_penalty,
                              F_prev_v[j] + init_score + config_.gap_extension_penalty);
        }

        match = std::max(F_v[j], match);

        if (j + 1 < prev_end)
            E_v[j + 1] = match + config_.gap_opening_penalty;

        if (match >= xdrop_cutoff)
            S_v[j] = std::max(match, E_v[j]);
    }
#else
    const __m128i gap_open = _mm_set1_epi32(config_.gap_opening_penalty);
    const __m128i gap_extend = _mm_set1_epi32(config_.gap_extension_penalty);
    const __m128i xdrop_v = _mm_set1_epi32(xdrop_cutoff - 1);
    const __m128i ninf_v = _mm_set1_epi32(ninf);
    const __m128i prev_end_v = _mm_set1_epi32(prev_end);
    const __m128i score_v = _mm_set1_epi32(init_score);
    __m128i j_v = _mm_set_epi32(3, 2, 1, 0);
    for (size_t j = 0; j < prev_end; j += 4) {
        // match = j ? S_prev_v[j - 1] + profile_scores[j] : ninf;
        __m128i match;
        if (j) {
            match = _mm_add_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j - 1]),
                                  _mm_loadu_si128((__m128i*)&profile_scores[j]));
            match = _mm_add_epi32(match, score_v);
        } else {
            // rotate elements to the right, then insert ninf in first cell
            match = _mm_shuffle_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), 0b10010000);
            match = _mm_add_epi32(match, _mm_loadu_si128((__m128i*)&profile_scores[j]));
            match = _mm_add_epi32(match, score_v);
            match = _mm_insert_epi32(match, ninf, 0);
        }

        // del_score = std::max(del_open, del_extend);
        __m128i del_score;
        if (offset > 1) {
            del_score = _mm_max_epi32(
                _mm_add_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), gap_open),
                _mm_add_epi32(_mm_loadu_si128((__m128i*)&F_prev_v[j]), gap_extend)
            );
            del_score = _mm_add_epi32(del_score, score_v);
        } else {
            del_score = ninf_v;
        }

        // F_v[j] = del_score
        _mm_store_si128((__m128i*)&F_v[j], del_score);

        // match = max(match, del_score)
        match = _mm_max_epi32(match, del_score);

        // match >= xdrop_cutoff
        __m128i mask = _mm_cmpgt_epi32(match, xdrop_v);

        // j < prev_end
        __m128i bound = _mm_cmpgt_epi32(prev_end_v, j_v);
        j_v = _mm_add_epi32(j_v, _mm_set1_epi32(4));
        mask = _mm_and_si128(mask, bound);
        match = _mm_blendv_epi8(ninf_v, match, mask);

        // ins_open_next = S[j] + gap_open
        __m128i ins_open_next = _mm_add_epi32(match, gap_open);

        // E_v[j + 1] = ins_open_next
        _mm_storeu_si128((__m128i*)&E_v[j + 1], ins_open_next);

        // load E_v[j] vector by rotating elements of ins_open_next right, then inserting E_v[j]
        __m128i ins_open = _mm_shuffle_epi32(ins_open_next, 0b10010000);
        ins_open = _mm_insert_epi32(ins_open, E_v[j], 0);

        // E_v[j] >= xdrop_cutoff
        ins_open = _mm_blendv_epi8(ninf_v, ins_open, mask);

        // S_v[j] = max(match, E_v[j])
        match = _mm_max_epi32(match, ins_open);
        _mm_store_si128((__m128i*)&S_v[j], match);
    }

#endif

    if (S_v.size() > prev_end) {
        size_t j = S_v.size() - 1;
        score_t match = std::max(S_prev_v[j - 1] + init_score + profile_scores[j], E_v[j]);
        if (match >= xdrop_cutoff)
            S_v[j] = match;
    }
}

// update insertion extension scores
void update_ins_extension(AlignedVector<score_t> &S,
                          AlignedVector<score_t> &E,
                          score_t xdrop_cutoff,
                          const DBGAlignerConfig &config_) {
    // elements are dependent on the previous one, so this can't be vectorized easily
    for (size_t j = 1; j < S.size(); ++j) {
        score_t ins_extend = E[j - 1] + config_.gap_extension_penalty;
        if (ins_extend > std::max(E[j], xdrop_cutoff - 1)) {
            E[j] = ins_extend;
            S[j] = std::max(S[j], ins_extend);
        }
    }
}

// add insertions to the end of the array until the score drops too low
void extend_ins_end(AlignedVector<score_t> &S,
                    AlignedVector<score_t> &E,
                    AlignedVector<score_t> &F,
                    size_t max_size,
                    score_t xdrop_cutoff,
                    const DBGAlignerConfig &config_) {
    while (S.back() >= xdrop_cutoff && S.size() < max_size) {
        score_t ins_score = std::max(
            S.back() + config_.gap_opening_penalty,
            E.back() + config_.gap_extension_penalty
        );

        if (ins_score < xdrop_cutoff)
            break;

        S.push_back(ins_score);
        E.push_back(ins_score);
        F.push_back(ninf);
    }

    // allocate and initialize enough space to allow the SIMD code to access these
    // vectors in 16 byte blocks without reading out of bounds
    S.reserve(S.size() + kPadding);
    E.reserve(E.size() + kPadding);
    F.reserve(F.size() + kPadding);

    std::fill(S.data() + S.size(), S.data() + S.capacity(), ninf);
    std::fill(E.data() + E.size(), E.data() + E.capacity(), ninf);
    std::fill(F.data() + F.size(), F.data() + F.capacity(), ninf);
}

void DefaultColumnExtender
::call_outgoing(node_index node,
                size_t /* max_prefetch_distance */,
                const std::function<void(node_index, char, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    assert(node == std::get<3>(table[table_i]));
    size_t next_offset = std::get<6>(table[table_i]) + 1;
    size_t seed_pos = next_offset - this->seed_->get_offset();
    bool in_seed = seed_pos < this->seed_->get_sequence().size();

    // Get the next node(s) from the graph. If the current node is
    // part of the seed and the arguments tell us to force using
    // the full seed, then pick the next node from the seed.
    if (in_seed && next_offset < graph_->get_k()) {
        assert(this->seed_->get_nodes().front());
        callback(this->seed_->get_nodes().front(),
                 this->seed_->get_sequence()[seed_pos],
                 0);
    } else if (in_seed && (force_fixed_seed || this->fixed_seed())) {
        size_t node_i = next_offset - graph_->get_k() + 1;
        assert(node == this->seed_->get_nodes()[node_i - 1]);
        node_index next_node = this->seed_->get_nodes()[node_i];
        char next_c = this->seed_->get_sequence()[seed_pos];
        callback(next_node, next_c, next_node
            ? 0
            : (!node ? config_.gap_extension_penalty : config_.gap_opening_penalty));
        assert(!node || graph_->traverse(node, next_c) == next_node);
    } else {
        assert(node);
        graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                callback(next, c, 0);
        });
    }
}

// allocate and initialize with padding to ensure that SIMD operations don't
// read/write out of bounds
template <class Column, typename... RestArgs>
Column alloc_column(size_t size, RestArgs... args) {
    Column column { {}, {}, {}, args... };
    auto &[S, E, F, node, i_prev, c, offset, max_pos, trim,
           xdrop_cutoff_i, last_fork_i, score] = column;

    // allocate and initialize enough space to allow the SIMD code to access these
    // vectors in 16 byte blocks without reading out of bounds
    S.reserve(size + kPadding);
    E.reserve(size + kPadding);
    F.reserve(size + kPadding);

    // the size is set properly to allow for AlignedVector methods (size(), push_back())
    // to function properly
    S.resize(size, ninf);
    E.resize(size, ninf);
    F.resize(size, ninf);

    std::fill(S.data() + S.size(), S.data() + S.capacity(), ninf);
    std::fill(E.data() + E.size(), E.data() + E.capacity(), ninf);
    std::fill(F.data() + F.size(), F.data() + F.capacity(), ninf);

    return column;
}

std::vector<Alignment> DefaultColumnExtender
::extend(score_t min_path_score, bool force_fixed_seed) {
    assert(this->seed_);

    ++num_extensions_;

    min_path_score = std::max(0, min_path_score);

    table.clear();
    prev_starts.clear();

    assert(config_.xdrop > 0);

    xdrop_cutoffs_.assign(1, std::make_pair(0u, std::max(-config_.xdrop, ninf + 1)));
    assert(xdrop_cutoffs_[0].second < 0);

    if (!config_.global_xdrop)
        scores_reached_.assign(1, 0);

    typedef typename decltype(xdrop_cutoffs_)::value_type xdrop_cutoff_v;
    table_size_bytes_ = sizeof(table) + sizeof(xdrop_cutoffs_) + sizeof(xdrop_cutoff_v)
        + sizeof(scores_reached_) + sizeof(score_t);

    size_t start = this->seed_->get_clipping();

    // the sequence to align (the suffix of the query starting from the seed)
    std::string_view window(this->seed_->get_query().data(),
                            query_.data() + query_.size() - this->seed_->get_query().data());
    assert(partial_sums_.at(start) == config_.match_score(window));

    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset()) - 1;

    // initialize the root of the tree
    table.emplace_back(
        alloc_column<Column>(1, this->seed_->get_nodes().front(), static_cast<size_t>(-1),
                             '\0', seed_offset, 0, 0, 0u, 0u, 0)
    );

    {
        auto &[S, E, F, node, i_prev, c, offset, max_pos, trim,
               xdrop_cutoff_i, last_fork_i, score] = table[0];
        S[0] = config_.left_end_bonus && !seed_->get_clipping() ? config_.left_end_bonus : 0;
        extend_ins_end(S, E, F, window.size() + 1 - trim,
                       xdrop_cutoffs_[xdrop_cutoff_i].second, config_);

        static_assert(std::is_same_v<decltype(table)::value_type, Column>);
        static_assert(std::is_same_v<decltype(S)::value_type, score_t>);
        static_assert(std::is_same_v<decltype(E)::value_type, score_t>);
        static_assert(std::is_same_v<decltype(F)::value_type, score_t>);

        table_size_bytes_ = sizeof(Column) * table.capacity()
                                + S.capacity() * sizeof(score_t) * 3;
    }

    // The nodes in the traversal (with corresponding score columns) are sorted by
    // 1) their score (higher is better), then by
    // 2) the absolute distance of their highest scoring index from the score
    //    matrix diagonal (lower is better), finally by
    // 3) Their index in the table vector (higher is better, for better cache locality)
    using TableIt = std::tuple<score_t, /* max changed value */
                               ssize_t, /* negative off_diag */
                               size_t, /* table idx */
                               score_t /* max score */>;
    TableIt best_score { 0, 0, 0, 0 };

    min_cell_score_ = 0;

    std::priority_queue<TableIt> queue;

    // Initialize the node traversal heap with the root.
    queue.emplace(best_score);

    while (queue.size()) {
        std::vector<TableIt> next_nodes{ queue.top() };
        queue.pop();

        // try all paths which have the same best partial alignment score
        // (this performs a BFS-like search)
        while (queue.size() && std::get<0>(queue.top()) == std::get<0>(next_nodes.back())) {
            next_nodes.push_back(queue.top());
            queue.pop();
        }

        while (next_nodes.size()) {
            size_t i = std::get<2>(next_nodes.back());
            next_nodes.pop_back();

            std::vector<std::tuple<node_index, char, score_t>> outgoing;
            size_t next_offset = -1;

            size_t prev_begin = 0;
            size_t prev_end = window.size() + 1;
            size_t prev_xdrop_cutoff_i = 0;
            size_t prev_last_fork_i = 0;
            score_t prev_xdrop_cutoff = 0;

            bool in_seed;

            {
                const auto &[S, E, F, node, i_prev, c, offset, max_pos, trim,
                             xdrop_cutoff_i, last_fork_i, score] = table[i];
                next_offset = offset + 1;
                prev_xdrop_cutoff_i = xdrop_cutoff_i;
                prev_last_fork_i = last_fork_i;
                prev_xdrop_cutoff = xdrop_cutoffs_[xdrop_cutoff_i].second;

                size_t node_counter = config_.global_xdrop
                    ? table.size()
                    : offset - seed_offset + 1;

                // if we are currently not on the optimal path, check early cutoff criteria
                if (S[max_pos - trim] < std::get<3>(best_score)) {
                    // if too many nodes have been explored, give up
                    if (static_cast<double>(node_counter) / window.size()
                            >= config_.max_nodes_per_seq_char) {
                        DEBUG_LOG("Alignment node limit reached, stopping this branch");
                        if (config_.global_xdrop) {
                            queue = std::priority_queue<TableIt>();
                            next_nodes.clear();
                        }

                        continue;
                    }

                    if (static_cast<double>(table_size_bytes_) / 1'000'000
                            > config_.max_ram_per_alignment) {
                        DEBUG_LOG("Alignment RAM limit reached, stopping extension");
                        queue = std::priority_queue<TableIt>();
                        next_nodes.clear();
                        continue;
                    }
                }

                // determine maximal range within the xdrop score cutoff
                auto in_range = [prev_xdrop_cutoff](score_t s) {
                    return s >= prev_xdrop_cutoff;
                };

                prev_begin = std::find_if(S.begin(), S.end(), in_range) - S.begin() + trim;
                prev_end = std::find_if(S.rbegin(), S.rend(), in_range).base() - S.begin() + trim;

                if (prev_end <= prev_begin)
                    continue;

                in_seed = next_offset - this->seed_->get_offset()
                                < this->seed_->get_sequence().size();

                call_outgoing(node, window.size() + 1 - offset - S.size(),
                              [&](node_index next, char c, score_t s) {
                                  outgoing.emplace_back(next, c, s);
                              }, i, force_fixed_seed);
            }

            ssize_t begin = prev_begin;
            size_t end = std::min(prev_end, window.size()) + 1;

            std::vector<TableIt> to_push;
            for (const auto &[next, c, score] : outgoing) {
                bool forked = outgoing.size() > 1;
                bool forked_xdrop = !config_.global_xdrop && forked;
                size_t xdrop_cutoffs_sizediff = xdrop_cutoffs_.capacity();
                if (forked_xdrop) {
                    xdrop_cutoffs_.emplace_back(table.size(), prev_xdrop_cutoff);
                    xdrop_cutoffs_sizediff = xdrop_cutoffs_.capacity() - xdrop_cutoffs_sizediff;
                } else {
                    xdrop_cutoffs_sizediff = 0;
                }

                size_t table_sizediff = table.capacity();
                table.emplace_back(alloc_column<Column>(
                    end - begin, next, i, c,
                    static_cast<ssize_t>(next_offset),
                    begin, begin,
                    forked_xdrop ? xdrop_cutoffs_.size() - 1 : prev_xdrop_cutoff_i,
                    forked ? table.size() : prev_last_fork_i,
                    score
                ));

                const auto &[S_prev, E_prev, F_prev, node_prev, i_prev, c_prev,
                             offset_prev, max_pos_prev, trim_prev,
                             xdrop_cutoff_i_prev, last_fork_i_prev, score_prev] = table[i];

                auto &[S, E, F, node_cur, i_cur, c_stored, offset, max_pos, trim,
                       xdrop_cutoff_i, last_fork_i, score_cur] = table.back();

                score_t &xdrop_cutoff = xdrop_cutoffs_[xdrop_cutoff_i].second;

                assert(i_cur == i);
                assert(node_cur == next);
                assert(c_stored == c);
                assert(offset == offset_prev + 1);
                assert(!node_cur || c == graph_->get_node_sequence(node_cur)[std::min(
                    static_cast<ssize_t>(graph_->get_k()) - 1, offset
                )]);

                // compute column scores
                update_column(prev_end - trim,
                              S_prev.data() + trim - trim_prev,
                              F_prev.data() + trim - trim_prev,
                              S, E, F,
                              profile_score_[c].data() + start + trim,
                              xdrop_cutoff, config_, score, offset);

                update_ins_extension(S, E, xdrop_cutoff, config_);
                extend_ins_end(S, E, F, window.size() + 1 - trim, xdrop_cutoff, config_);

                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                // find the maximal scoring position which is closest to the diagonal
                // TODO: this can be done with SIMD, but it's not a bottleneck
                ssize_t cur_offset = begin;
                ssize_t diag_i = offset - seed_offset;
                bool has_extension = in_seed;
                const score_t *partial_sums = &partial_sums_[start + trim];
                score_t extension_cutoff = std::get<3>(best_score) * config_.rel_score_cutoff;
                score_t max_diff = ninf;

                size_t scores_reached_sizediff = 0;
                bool scores_reached_cutoff = true;
                if (!config_.global_xdrop) {
                    scores_reached_sizediff = scores_reached_.capacity();
                    scores_reached_.resize(S.size() + trim + 1, ninf);
                    scores_reached_sizediff = scores_reached_.capacity()
                                                - scores_reached_sizediff;
                }

                for (size_t j = 0; j < S.size(); ++j, ++cur_offset) {
                    if (S[j] != ninf)
                        min_cell_score_ = std::min(min_cell_score_, S[j]);

                    if (std::make_pair(S[j], std::abs(max_pos - diag_i))
                            > std::make_pair(S[max_pos - begin], std::abs(cur_offset - diag_i))) {
                        max_pos = j + begin;
                    }

                    if (!config_.global_xdrop) {
                        scores_reached_[trim + j] = std::max(scores_reached_[trim + j], S[j]);
                        scores_reached_cutoff = (S[j] >= scores_reached_[trim + j] * config_.rel_score_cutoff);
                    }

                    // check if this node can be extended to get a better alignment
                    assert(partial_sums[j] == config_.match_score(window.substr(j + trim)));
                    if (!has_extension && S[j] + partial_sums[j] >= extension_cutoff
                            && scores_reached_cutoff)
                        has_extension = true;

                    if (static_cast<size_t>(trim - trim_prev) < S_prev.size()
                            && S[j] - S_prev[j + trim - trim_prev] > max_diff) {
                        max_diff = S[j] - S_prev[j + trim - trim_prev];
                    }
                }

                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                score_t max_val = S[max_pos - trim];

                if (static_cast<size_t>(offset - seed_offset) < config_.target_distance + 1) {
                    has_extension = true;
                } else if (config_.target_distance) {
                    has_extension = false;
                }

                if (!in_seed && (max_val < xdrop_cutoff || !has_extension)) {
                    pop(table.size() - 1);
                    if (forked_xdrop)
                        xdrop_cutoffs_.pop_back();

                    continue;
                }

                table_sizediff = table.capacity() - table_sizediff;

                static_assert(std::is_same_v<decltype(table)::value_type, Column>);
                static_assert(std::is_same_v<decltype(S)::value_type, score_t>);
                static_assert(std::is_same_v<decltype(E)::value_type, score_t>);
                static_assert(std::is_same_v<decltype(F)::value_type, score_t>);
                static_assert(std::is_same_v<decltype(scores_reached_)::value_type, score_t>);

                table_size_bytes_ += sizeof(Column) * table_sizediff
                    + S.capacity() * sizeof(score_t) * 3
                    + sizeof(score_t) * scores_reached_sizediff
                    + sizeof(xdrop_cutoff_v) * xdrop_cutoffs_sizediff;

                if (max_val - xdrop_cutoff > config_.xdrop)
                    xdrop_cutoff = max_val - config_.xdrop;

                // if the best score in this column is above the xdrop score
                // then check if the extension can continue
                TableIt next_score { 0, -std::abs(max_pos - diag_i),
                                     table.size() - 1, max_val };

                if (max_val > std::get<3>(best_score))
                    best_score = next_score;

                // skip the first index since it corresponds to the position
                // before the query start
                size_t vec_offset = start + begin - static_cast<bool>(begin);
                score_t *s_begin = S.data() + !begin;
                score_t *s_end = S.data() + S.size();

                assert(s_begin <= s_end);
                assert(vec_offset + (s_end - s_begin) <= query_.size());

                // if this node has not been reached by a different
                // alignment with a better score, continue
                std::get<0>(next_score) = this->update_seed_filter(
                    next, vec_offset, s_begin, s_end
                );
                if (std::get<0>(next_score) != ninf) {
                    // if this next node is the only next option, or if it's
                    // better than all other options, take it without pushing
                    // to the queue
                    if (next_nodes.size()
                            && std::get<0>(next_score) == std::get<0>(next_nodes[0])) {
                        next_nodes.emplace_back(std::move(next_score));
                    } else {
                        queue.emplace(std::move(next_score));
                    }
                } else {
                    DEBUG_LOG("Dropped due to convergence");
                }
            }
        }
    }

    if (config_.no_backtrack)
        return { *seed_ };

    return backtrack(min_path_score, window);
}

Alignment DefaultColumnExtender::construct_alignment(Cigar cigar,
                                                     size_t clipping,
                                                     std::string_view window,
                                                     std::vector<node_index> final_path,
                                                     std::string match,
                                                     score_t score,
                                                     size_t offset,
                                                     score_t extra_penalty) const {
    assert(final_path.size());
    cigar.append(Cigar::CLIPPED, clipping);

    std::reverse(cigar.data().begin(), cigar.data().end());
    std::reverse(final_path.begin(), final_path.end());
    std::reverse(match.begin(), match.end());

    Alignment extension(window, std::move(final_path), std::move(match), score,
                        std::move(cigar), 0, this->seed_->get_orientation(), offset);
    extension.extend_query_begin(query_.data());
    extension.extend_query_end(query_.data() + query_.size());
    extension.extra_penalty = extra_penalty;
    assert(extension.is_valid(*this->graph_, &config_));

    if (config_.trim_offset_after_extend) {
        extension.trim_offset();
        assert(extension.is_valid(*this->graph_, &config_));
    }

    return extension;
}

std::vector<Alignment> DefaultColumnExtender
::backtrack(score_t min_path_score, std::string_view window) {
    std::vector<Alignment> extensions;
    size_t seed_clipping = this->seed_->get_clipping();
    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset() - 1);
    ssize_t k_minus_1 = graph_->get_k() - 1;
    ssize_t last_pos = window.size();
    ssize_t seed_dist = std::max(graph_->get_k(), this->seed_->get_sequence().size()) - 1;
    score_t min_start_score = config_.semiglobal ? ninf : min_path_score;
    size_t min_trace_length = this->graph_->get_k() - this->seed_->get_offset();

    std::vector<std::tuple<score_t, ssize_t, ssize_t, ssize_t>> indices;
    indices.reserve(table.size());
    for (size_t i = 1; i < table.size(); ++i) {
        auto check_and_add_pos = [&](ssize_t start_pos) {
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                         xdrop_cutoff_i, last_fork_i, score] = table[i];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                         xdrop_cutoff_i_p, last_fork_i_p, score_p] = table[j_prev];
            if (start_pos < trim_p + 1)
                return;

            size_t pos = start_pos - trim;
            size_t pos_p = start_pos - trim_p - 1;

            if (S[pos] == ninf || S_p[pos_p] == ninf)
                return;

            score_t end_bonus = start_pos == last_pos ? config_.right_end_bonus : 0;

            if (config_.semiglobal) {
                if (node == config_.terminal_node
                        && S[pos] == S_p[pos_p] + score + profile_score_.find(c)->second[seed_clipping + start_pos]) {
                    indices.emplace_back(S[pos] + end_bonus, -std::abs(start_pos - offset + seed_offset),
                                         -static_cast<ssize_t>(i), start_pos);
                }
            } else if (S[pos] + end_bonus >= min_start_score) {
                bool is_match = S[pos] == S_p[pos_p] + score + profile_score_.find(c)->second[seed_clipping + start_pos]
                    && profile_op_.find(c)->second[seed_clipping + start_pos] == Cigar::MATCH;
                if (is_match || start_pos == last_pos) {
                    indices.emplace_back(
                        S[pos] + end_bonus, -std::abs(start_pos - offset + seed_offset),
                        -static_cast<ssize_t>(i), start_pos
                    );
                }
            }
        };

        const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                     xdrop_cutoff_i, last_fork_i, score] = table[i];

        if (offset < seed_dist)
            continue;

        if (!config_.semiglobal)
            check_and_add_pos(max_pos);

        if (S.size() + trim == window.size() + 1
                && (config_.semiglobal || max_pos != last_pos)) {
            check_and_add_pos(last_pos);
        }
    }

    // find highest scoring which is closest to the diagonal
    // use heap sort to make this run in O(n + (num_alternative_paths) * log(n)) time
    std::make_heap(indices.begin(), indices.end());

    score_t best_score = std::numeric_limits<score_t>::min();

    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        std::pop_heap(indices.begin(), it.base());
        const auto &[start_score, neg_off_diag, neg_j_start, start_pos] = *it;

        if (terminate_backtrack_start(extensions))
            break;

        size_t j = -neg_j_start;

        if (skip_backtrack_start(j))
            continue;

        std::vector<DeBruijnGraph::node_index> path;
        std::vector<size_t> trace;
        Cigar ops;
        std::string seq;
        score_t score = start_score;

        if (score - min_cell_score_ < best_score)
            break;

        size_t dummy_counter = 0;

        ssize_t pos = start_pos;
        ssize_t end_pos = pos;
        size_t align_offset = this->seed_->get_offset();

        auto append_node = [&](node_index node, ssize_t offset, Cigar::Operator op) {
            ops.append(op);
            if (offset >= k_minus_1) {
                path.emplace_back(node);
                if (!node) {
                    ++dummy_counter;
                } else if (dummy_counter) {
                    ops.append(Cigar::NODE_INSERTION, dummy_counter);
                    dummy_counter = 0;
                }
            }
        };

        score_t extra_penalty = 0;

        while (j && !terminate_backtrack()) {
            assert(j != static_cast<size_t>(-1));
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                         xdrop_cutoff_i, last_fork_i, score_cur] = table[j];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                         xdrop_cutoff_i_p, last_fork_i_p, score_cur_p] = table[j_prev];

            assert(pos >= trim);
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);
            assert(!node || c == graph_->get_node_sequence(node)[std::min(k_minus_1, offset)]);

            align_offset = std::min(offset, k_minus_1);

            if (pos == max_pos)
                prev_starts.emplace(j);

            if (S[pos - trim] == ninf) {
                j = 0;

            } else if (pos && S[pos - trim] == E[pos - trim] && (ops.empty()
                    || ops.data().back().first != Cigar::DELETION)) {
                // prefer insertions over match/mismatch to prevent premature
                // clipping of the alignment beginning
                Cigar::Operator last_op = Cigar::INSERTION;
                while (last_op == Cigar::INSERTION) {
                    ops.append(last_op);

                    assert(E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        || E[pos - trim] == S[pos - trim - 1] + config_.gap_opening_penalty);

                    last_op = E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        ? Cigar::INSERTION
                        : Cigar::MATCH;

                    --pos;
                }
            } else if (pos && pos >= trim_p + 1
                    && S[pos - trim] == S_p[pos - trim_p - 1] + score_cur
                        + profile_score_.find(c)->second[seed_clipping + pos]) {
                // match/mismatch
                trace.emplace_back(j);

                seq += c;
                append_node(node, offset, profile_op_.find(c)->second[seed_clipping + pos]);
                --pos;
                assert(j_prev != static_cast<size_t>(-1));
                j = j_prev;
                extra_penalty += score_cur;

            } else if (S[pos - trim] == F[pos - trim] && (ops.empty()
                    || ops.data().back().first != Cigar::INSERTION)) {
                // deletion
                Cigar::Operator last_op = Cigar::DELETION;
                while (last_op == Cigar::DELETION && j) {
                    const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                                 xdrop_cutoff_i, last_fork_i, score_cur] = table[j];
                    const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                                 xdrop_cutoff_i_p, last_fork_i_p, score_cur_p] = table[j_prev];

                    assert(pos >= trim_p);

                    assert(F[pos - trim] == F_p[pos - trim_p] + score_cur + config_.gap_extension_penalty
                        || F[pos - trim] == S_p[pos - trim_p] + score_cur + config_.gap_opening_penalty);

                    last_op = F[pos - trim] == F_p[pos - trim_p] + score_cur + config_.gap_extension_penalty
                        ? Cigar::DELETION
                        : Cigar::MATCH;

                    trace.emplace_back(j);
                    append_node(node, offset, Cigar::DELETION);
                    extra_penalty += score_cur;

                    seq += c;
                    assert(j_prev != static_cast<size_t>(-1));
                    j = j_prev;
                }

            } else {
                DEBUG_LOG("Backtracking failed, trying next start point");
                assert(false);
                break;
            }

            if (trace.size() >= min_trace_length && path.size() && path.back()) {
                const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                             xdrop_cutoff_i, last_fork_i, score_cur] = table[j];

                best_score = std::max(best_score, score - S[pos - trim]);
                if (score - min_cell_score_ < best_score)
                    break;

                call_alignments(S[pos - trim], score, min_start_score, path, trace, j,
                                ops, pos, align_offset, window.substr(pos, end_pos - pos),
                                seq, extra_penalty,
                                [&](Alignment&& alignment) {
                    extensions.emplace_back(std::move(alignment));
                });
            }
        }
    }

    if (extensions.empty() && this->seed_->get_score() >= min_path_score) {
        extensions.emplace_back(*this->seed_);
        if (config_.trim_offset_after_extend) {
            extensions.back().trim_offset();
            extensions.back().extend_query_end(query_.data() + query_.size());
        } else {
            extensions.back().get_cigar().trim_clipping();
        }
    }

    return extensions;
}

} // namespace align
} // namespace graph
} // namespace mtg
