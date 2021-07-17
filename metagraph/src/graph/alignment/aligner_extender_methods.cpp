#include "aligner_extender_methods.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

using score_t = DBGAlignerConfig::score_t;
constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

// to ensure that SIMD operations on arrays don't read out of bounds
constexpr size_t kPadding = 5;


DefaultColumnExtender::DefaultColumnExtender(const DeBruijnGraph &graph,
                                             const DBGAlignerConfig &config,
                                             std::string_view query)
      : SeedFilteringExtender(query),
        graph_(&graph), config_(config), query_(query) {
    assert(config_.check_config_scores());

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

bool SeedFilteringExtender::set_seed(const Alignment &seed) {
    assert(seed.size());
    assert(seed.get_cigar().size());
    assert(seed.get_cigar().back().first == Cigar::MATCH
        || seed.get_cigar().back().first == Cigar::MISMATCH);

    seed_ = nullptr;

    auto it = conv_checker_.find(seed.back());

    if (it != conv_checker_.end()) {
        size_t pos = seed.get_query().size() + seed.get_clipping() - 1;
        const auto &[start, vec] = it->second;
        if (pos < start || pos - start >= vec.size() || vec[pos - start] < seed.get_score())
            it = conv_checker_.end();
    }

    if (it == conv_checker_.end()) {
        seed_ = &seed;
    } else {
        DEBUG_LOG("Skipping seed: {}", seed);
    }

    return seed_;
}

bool SeedFilteringExtender::update_seed_filter(node_index node,
                                               size_t query_start,
                                               const score_t *s_begin,
                                               const score_t *s_end) {
    assert(s_end >= s_begin);
    assert(query_start + (s_end - s_begin) <= query_size_);

    size_t size = s_end - s_begin;

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end()) {
        conv_checker_.emplace(node, ScoreVec(query_start, { s_begin, s_end }));
        return true;

    } else {
        auto &[start, vec] = it.value();
        if (query_start + size <= start) {
            vec.insert(vec.begin(), start - query_start, ninf);
            std::copy(s_begin, s_end, vec.begin());
            start = query_start;
            return true;

        } else if (query_start >= start + vec.size()) {
            vec.reserve(query_start + size - start);
            vec.insert(vec.end(), query_start - start - vec.size(), ninf);
            vec.insert(vec.end(), s_begin, s_end);
            return true;

        } else {
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
                if (s_begin[j] > v[j]) {
                    converged = false;
                    v[j] = s_begin[j];
                }
            }

            return !converged;
        }
    }
}

bool SeedFilteringExtender::filter_nodes(node_index node,
                                         size_t query_start,
                                         size_t query_end) {
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

    } else {
        auto &[start, vec] = it.value();
        if (query_start + size <= start) {
            vec.insert(vec.begin(), start - query_start, ninf);
            std::fill(vec.begin(), vec.begin() + size, mscore);
            start = query_start;
            return true;

        } else if (query_start >= start + vec.size()) {
            vec.reserve(query_start + size - start);
            vec.insert(vec.end(), query_start - start - vec.size(), ninf);
            vec.insert(vec.end(), size, mscore);
            return true;

        } else {
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
    }
}

void update_column(size_t prev_end,
                   const score_t *S_prev_v,
                   const score_t *F_prev_v,
                   AlignedVector<score_t> &S_v,
                   AlignedVector<score_t> &E_v,
                   AlignedVector<score_t> &F_v,
                   const score_t *profile_scores,
                   score_t xdrop_cutoff,
                   const DBGAlignerConfig &config_) {
#ifndef __SSE4_1__
    for (size_t j = 0; j < prev_end; ++j) {
        score_t match = j ? (S_prev_v[j - 1] + profile_scores[j]) : ninf;
        F_v[j] = std::max(S_prev_v[j] + config_.gap_opening_penalty,
                          F_prev_v[j] + config_.gap_extension_penalty);

        match = std::max(del_score, match);

        if (match >= xdrop_cutoff) {
            S_v[j] = match;
            score_t ins_score = match + config_.gap_opening_penalty;
            if (j + 1 < prev_end)
                E_v[j + 1] = ins_score;
        }
    }
#else
    const __m128i gap_open = _mm_set1_epi32(config_.gap_opening_penalty);
    const __m128i gap_extend = _mm_set1_epi32(config_.gap_extension_penalty);
    const __m128i xdrop_v = _mm_set1_epi32(xdrop_cutoff - 1);
    const __m128i ninf_v = _mm_set1_epi32(ninf);
    const __m128i prev_end_v = _mm_set1_epi32(prev_end);
    __m128i j_v = _mm_set_epi32(3, 2, 1, 0);
    for (size_t j = 0; j < prev_end; j += 4) {
        // match = j ? S_prev_v[j - 1] + profile_scores[j] : ninf;
        __m128i match;
        if (j) {
            match = _mm_add_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j - 1]),
                                  _mm_loadu_si128((__m128i*)&profile_scores[j]));
        } else {
            // rotate elements to the right, then insert ninf in first cell
            match = _mm_shuffle_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), 0b10010000);
            match = _mm_add_epi32(match, _mm_loadu_si128((__m128i*)&profile_scores[j]));
            match = _mm_insert_epi32(match, ninf, 0);
        }

        // del_score = std::max(del_open, del_extend);
        __m128i del_score = _mm_max_epi32(
            _mm_add_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), gap_open),
            _mm_add_epi32(_mm_loadu_si128((__m128i*)&F_prev_v[j]), gap_extend)
        );

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

        // S_v[j] = match
        _mm_store_si128((__m128i*)&S_v[j], match);

        // ins_open = S[j] + gap_open
        __m128i ins_open = _mm_add_epi32(match, gap_open);

        // E_v[j + 1] = ins_open
        _mm_storeu_si128((__m128i*)&E_v[j + 1], ins_open);
    }

#endif

    if (S_v.size() > prev_end) {
        size_t j = S_v.size() - 1;
        score_t match = S_prev_v[j - 1] + profile_scores[j];
        if (match >= xdrop_cutoff)
            S_v[j] = match;
    }
}

void update_ins(AlignedVector<score_t> &S,
                AlignedVector<score_t> &E,
                score_t xdrop_cutoff,
                const DBGAlignerConfig &config_) {
    // update insertion scores
    // elements are dependent on the previous one, so this can't be vectorized easily
    // this takes 15% of the run time when aligning long reads...
#ifndef __SSE4_1__
    for (size_t j = 1; j < S.size(); ++j) {
        E[j] = std::max(E[j - 1] + config_.gap_extension_penalty, E[j]);
        if (E[j] >= xdrop_cutoff)
            S[j] = std::max(S[j], E[j]);
    }
#else
    // first update the best score vector with insert open penalties
    const __m128i xdrop_v = _mm_set1_epi32(xdrop_cutoff - 1);
    const __m128i ninf_v = _mm_set1_epi32(ninf);
    for (size_t j = 0; j < S.size(); j += 4) {
        __m128i ins_open = _mm_load_si128((__m128i*)&E[j]);
        __m128i mask = _mm_cmpgt_epi32(ins_open, xdrop_v);
        ins_open = _mm_blendv_epi8(ninf_v, ins_open, mask);
        ins_open = _mm_max_epi32(_mm_load_si128((__m128i*)&S[j]), ins_open);
        _mm_store_si128((__m128i*)&S[j], ins_open);
    }

    // then compute insert extension penalties
    for (size_t j = 1; j < S.size(); ++j) {
        score_t ins_extend = E[j - 1] + config_.gap_extension_penalty;
        if (ins_extend > std::max(E[j], xdrop_cutoff - 1)) {
            E[j] = ins_extend;
            S[j] = std::max(S[j], ins_extend);
        }
    }
#endif
}

// add insertions to the end of the array until the score drops too low
void extend_ins(AlignedVector<score_t> &S,
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

        if (ins_score >= xdrop_cutoff) {
            S.push_back(ins_score);
            E.push_back(ins_score);
            F.push_back(ninf);
        } else {
            break;
        }
    }
}

void DefaultColumnExtender
::call_outgoing(node_index node,
                size_t /* max_prefetch_distance */,
                const std::function<void(node_index, char)> &callback) {
    graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
        if (c != boss::BOSS::kSentinel)
            callback(next, c);
    });
}

// allocate and initialize with padding to ensure that SIMD operations don't
// read/write out of bounds
template <class Column, typename... RestArgs>
Column alloc_column(size_t size, RestArgs... args) {
    Column column { {}, {}, {}, args... };
    auto &[S, E, F, node, i_prev, c, offset, max_pos, trim] = column;

    S.reserve(size + kPadding);
    E.reserve(size + kPadding);
    F.reserve(size + kPadding);

    S.resize(size, ninf);
    E.resize(size, ninf);
    F.resize(size, ninf);

    std::fill(S.data() + S.size(), S.data() + S.capacity(), ninf);
    std::fill(E.data() + E.size(), E.data() + E.capacity(), ninf);
    std::fill(F.data() + F.size(), F.data() + F.capacity(), ninf);

    return column;
}

auto DefaultColumnExtender::extend(score_t min_path_score) -> std::vector<Alignment> {
    assert(this->seed_);

    table.clear();
    prev_starts.clear();

    size_t start = this->seed_->get_clipping();

    // the sequence to align (the suffix of the query starting from the seed)
    std::string_view window(this->seed_->get_query().data(),
                            query_.data() + query_.size() - this->seed_->get_query().data());

    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset()) - 1;

    // initialize the root of the tree
    table.emplace_back(alloc_column<Column>(window.size() + 1, this->seed_->front(),
                                            static_cast<size_t>(-1), '\0', seed_offset, 0, 0));
    std::fill(std::get<0>(table[0]).begin(), std::get<0>(table[0]).end(), 0);

    score_t xdrop_cutoff = std::max(-config_.xdrop, ninf + 1);
    assert(config_.xdrop > 0);
    assert(xdrop_cutoff < 0);

    using TableIt = std::tuple<score_t,
                               ssize_t, /* negative off_diag */
                               size_t /* table idx */>;
    TableIt best_score { 0, 0, 0 };

    // node traversal heap
    std::priority_queue<TableIt> queue;
    queue.emplace(best_score);

    assert(partial_sums_.at(start) == config_.match_score(window));

    while (queue.size()) {
        size_t i = std::get<2>(queue.top());
        queue.pop();

        // If a node has a single outgoing edge, which has the best score, then
        // that node is taken without updating the heap. This loop ensures that
        // that happens
        bool nofork = true;
        while (nofork) {
            nofork = false;

            std::vector<std::pair<node_index, char>> outgoing;
            size_t next_offset = -1;

            size_t prev_begin = 0;
            size_t prev_end = window.size() + 1;

            {
                const auto &[S, E, F, node, i_prev, c, offset, max_pos, trim] = table[i];
                next_offset = offset + 1;

                // if too many nodes have been explored, give up
                if (static_cast<double>(table.size()) / window.size() >= config_.max_nodes_per_seq_char)
                    continue;

                // determine maximal range within the xdrop score cutoff
                auto in_range = [xdrop_cutoff](score_t s) { return s >= xdrop_cutoff; };

                prev_begin = std::find_if(S.begin(), S.end(), in_range) - S.begin() + trim;
                prev_end = std::find_if(S.rbegin(), S.rend(), in_range).base() - S.begin() + trim;

                if (prev_end <= prev_begin)
                    continue;

                // check if this node can be extended to get a better alignment
                bool has_extension = false;
                for (size_t j = prev_begin; j < prev_end; ++j) {
                    assert(partial_sums_.at(start + j) == config_.match_score(window.substr(j)));
                    score_t ext_score = S[j - trim] + partial_sums_.at(start + j);
                    if ((config_.num_alternative_paths == 1 && ext_score > std::get<0>(best_score))
                            || ext_score >= min_path_score) {
                        has_extension = true;
                        break;
                    }
                }

                if (!has_extension)
                    continue;

                // Get the next node(s) from the graph. If the current node is
                // part of the seed, then pick the next node from the seed.
                if (next_offset - this->seed_->get_offset() < this->seed_->get_sequence().size()) {
                    if (next_offset < graph_->get_k()) {
                        outgoing.emplace_back(
                            this->seed_->front(),
                            this->seed_->get_sequence()[next_offset - this->seed_->get_offset()]
                        );
                    } else {
                        outgoing.emplace_back(
                            (*this->seed_)[next_offset - graph_->get_k() + 1],
                            this->seed_->get_sequence()[next_offset - this->seed_->get_offset()]
                        );
                        assert(graph_->traverse(node, outgoing.back().second) == outgoing.back().first);
                    }
                } else {
                    call_outgoing(node, window.size() + 1 - offset - S.size(),
                                  [&](node_index next, char c) { outgoing.emplace_back(next, c); });
                }
            }

            ssize_t begin = prev_begin;
            size_t end = std::min(prev_end, window.size()) + 1;

            for (const auto &[next, c] : outgoing) {
                assert(std::get<0>(best_score) > xdrop_cutoff);

                table.emplace_back(alloc_column<Column>(
                    end - begin, next, i, c,
                    static_cast<ssize_t>(next_offset),
                    begin, begin
                ));

                const auto &[S_prev, E_prev, F_prev, node_prev, i_prev, c_prev,
                             offset_prev, max_pos_prev, trim_prev] = table[i];

                auto &[S, E, F, node_cur, i_cur, c_stored, offset, max_pos, trim]
                    = table.back();

                if (next != node_prev && !trim && !trim_prev)
                    S[0] = S_prev[0];

                assert(i_cur == i);
                assert(node_cur == next);
                assert(c_stored == c);
                assert(offset == offset_prev + 1);
                assert(c == graph_->get_node_sequence(node_cur)[std::min(static_cast<ssize_t>(graph_->get_k()) - 1, offset)]);

                // compute column scores
                update_column(prev_end - trim,
                              S_prev.data() + trim - trim_prev,
                              F_prev.data() + trim - trim_prev,
                              S, E, F,
                              profile_score_[c].data() + start + trim,
                              xdrop_cutoff, config_);

                update_ins(S, E, xdrop_cutoff, config_);
                extend_ins(S, E, F, window.size() + 1 - trim, xdrop_cutoff, config_);

                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                // find the maximal scoring position which is closest to the diagonal
                // TODO: this can be done with SIMD, but it's not a bottleneck
                ssize_t cur_offset = begin;
                ssize_t diag_i = offset - seed_offset;
                for (size_t j = 0; j < S.size(); ++j, ++cur_offset) {
                    if (std::make_pair(S[j], std::abs(max_pos - diag_i))
                            > std::make_pair(S[max_pos - begin], std::abs(cur_offset - diag_i))) {
                        max_pos = j + begin;
                    }
                }
                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                // if the best score in this column is above the xdrop score
                // then check if the extension can continue
                score_t max_val = S[max_pos - trim];
                if (max_val >= xdrop_cutoff) {
                    TableIt next_score { max_val, -std::abs(max_pos - diag_i),
                                         table.size() - 1 };

                    if (max_val - xdrop_cutoff > config_.xdrop)
                        xdrop_cutoff = max_val - config_.xdrop;

                    if (max_val > std::get<0>(best_score))
                        best_score = next_score;

                    size_t vec_offset = start + begin;
                    score_t *s_begin = S.data();
                    score_t *s_end = S.data() + S.size();

                    // skip the first index since it corresponds to the position
                    // before the query start
                    if (!begin) {
                        ++s_begin;
                    } else {
                        --vec_offset;
                    }

                    assert(s_begin <= s_end);
                    assert(vec_offset + (s_end - s_begin) <= query_.size());

                    // if this node has not been reached by a different
                    // alignment with a better score, continue
                    if (this->update_seed_filter(next, vec_offset, s_begin, s_end)) {
                        // if there is only one outgoing edge, which is the best,
                        // don't update the node traversal heap
                        if (outgoing.size() == 1 && max_val >= std::get<0>(queue.top())) {
                            nofork = true;
                            i = table.size() - 1;
                        } else {
                            queue.emplace(next_score);
                        }
                    }

                } else {
                    table.pop_back();
                }
            }
        }
    }

    return backtrack(min_path_score, window);
}

auto DefaultColumnExtender::construct_alignment(Cigar cigar,
                                                size_t clipping,
                                                std::string_view window,
                                                std::vector<node_index> final_path,
                                                std::string match,
                                                score_t score,
                                                size_t offset) const -> Alignment {
    assert(final_path.size());
    cigar.append(Cigar::CLIPPED, clipping);

    std::reverse(cigar.begin(), cigar.end());
    std::reverse(final_path.begin(), final_path.end());
    std::reverse(match.begin(), match.end());

    Alignment extension(window, std::move(final_path), std::move(match), score,
                           std::move(cigar), 0, this->seed_->get_orientation(), offset);
    assert(extension.is_valid(*this->graph_, &config_));

    extension.trim_offset();
    extension.extend_query_begin(query_.data());
    extension.extend_query_end(query_.data() + query_.size());

    assert(extension.is_valid(*this->graph_, &config_));

    return extension;
}

auto DefaultColumnExtender
::backtrack(score_t min_path_score, std::string_view window) -> std::vector<Alignment> {
    if (table.empty())
        return {};

    init_backtrack();

    std::vector<Alignment> extensions;

    size_t seed_clipping = this->seed_->get_clipping();
    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset() - 1);
    ssize_t k_minus_1 = graph_->get_k() - 1;

    std::vector<std::tuple<score_t, ssize_t, size_t>> indices;
    indices.reserve(table.size());
    for (size_t i = 1; i < table.size(); ++i) {
        const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[i];
        const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

        if (max_pos < trim_p + 1)
            continue;

        size_t pos = max_pos - trim;
        size_t pos_p = max_pos - trim_p - 1;
        if (S[pos] >= min_path_score
                && offset >= k_minus_1
                && S[pos] == S_p[pos_p] + profile_score_.find(c)->second[seed_clipping + max_pos]
                && profile_op_.find(c)->second[seed_clipping + max_pos] == Cigar::MATCH) {
            indices.emplace_back(-S[pos], std::abs(max_pos - offset + seed_offset), i);
        }
    }

    // find highest scoring which is closest to the diagonal
    std::sort(indices.begin(), indices.end());

    for (const auto &[neg_score, off_diag, j_start] : indices) {
        if (terminate_backtrack_start(extensions))
            break;

        if (skip_backtrack_start(j_start))
            continue;

        std::vector<DeBruijnGraph::node_index> path;
        std::vector<size_t> trace;
        Cigar ops;
        std::string seq;
        score_t score = -neg_score;

        size_t j = j_start;
        ssize_t pos = std::get<7>(table[j]);
        ssize_t end_pos = pos;
        size_t align_offset = this->seed_->get_offset();

        while (j && !terminate_backtrack()) {
            assert(j != static_cast<size_t>(-1));
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

            assert(pos >= trim);
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);
            assert(c == graph_->get_node_sequence(node)[std::min(k_minus_1, offset)]);

            align_offset = std::min(offset, k_minus_1);

            if (pos == max_pos)
                prev_starts.emplace(j);

            if (S[pos - trim] == ninf) {
                j = 0;

            } else if (pos && pos >= trim_p + 1
                    && S[pos - trim] == S_p[pos - trim_p - 1]
                        + profile_score_.find(c)->second[seed_clipping + pos]) {
                // match/mismatch
                trace.emplace_back(j);
                if (offset >= k_minus_1)
                    path.emplace_back(node);

                seq += c;
                ops.append(profile_op_.find(c)->second[seed_clipping + pos]);
                --pos;
                assert(j_prev != static_cast<size_t>(-1));
                j = j_prev;

            } else if (S[pos - trim] == F[pos - trim] && ops.size() && ops.back().first != Cigar::INSERTION) {
                // deletion
                Cigar::Operator last_op = Cigar::DELETION;
                while (last_op == Cigar::DELETION && j) {
                    const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];
                    const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

                    assert(pos >= trim_p);

                    assert(F[pos - trim] == F_p[pos - trim_p] + config_.gap_extension_penalty
                        || F[pos - trim] == S_p[pos - trim_p] + config_.gap_opening_penalty);

                    last_op = F[pos - trim] == F_p[pos - trim_p] + config_.gap_extension_penalty
                        ? Cigar::DELETION
                        : Cigar::MATCH;

                    trace.emplace_back(j);
                    if (offset >= k_minus_1)
                        path.emplace_back(node);

                    seq += c;
                    ops.append(Cigar::DELETION);
                    assert(j_prev != static_cast<size_t>(-1));
                    j = j_prev;
                }
            } else if (pos && S[pos - trim] == E[pos - trim] && ops.size() && ops.back().first != Cigar::DELETION) {
                // insertion
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
#ifndef NDEBUG
            } else {
                assert(false && "Failure to backtrack. One of the above should apply");
#endif
            }

            if (trace.size() >= this->graph_->get_k()) {
                const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];

                call_alignments(S[pos - trim], score, min_path_score, path, trace, ops,
                                pos, align_offset, window.substr(pos, end_pos - pos), seq,
                                [&](Alignment&& alignment) {
                    extensions.emplace_back(std::move(alignment));
                });
            }
        }
    }

    return extensions;
}

} // namespace align
} // namespace graph
} // namespace mtg
