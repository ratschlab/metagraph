#ifndef __ALIGN_CHAIN__
#define __ALIGN_CHAIN__

#include "graph/alignment/alignment.hpp"

namespace mtg::graph::align {

template <typename Anchor>
using ChainScores = std::vector<std::tuple<score_t, const Anchor*, size_t>>;

using AlignmentCallback = std::function<void(Alignment&&)>;

template <typename Anchor>
using AnchorConnector = std::function<void(const Anchor&, // start anchor
                                           ssize_t b,     // max dist
                                           const Anchor*, // target anchors begin
                                           const Anchor*, // target anchors end
                                           typename ChainScores<Anchor>::pointer,
                                           const std::function<bool(score_t,       // connect score
                                                                    const Anchor*, // last
                                                                    size_t         // distance
                                                                   )>&
                                          )>;

template <typename Anchor>
using AnchorExtender = std::function<void(const Anchor*, // first ptr,
                                          Alignment&&,   // target,
                                          size_t,        // distance,
                                          score_t,       // chain score up to this point
                                          const AlignmentCallback&)>;

template <typename Anchor>
struct AnchorChain {
    using storage_t = std::vector<std::pair<const Anchor*, size_t>>;
    using value_type = typename storage_t::value_type;
    using reference = typename storage_t::reference;
    using const_reference = typename storage_t::const_reference;
    using iterator = typename storage_t::iterator;
    using const_iterator = typename storage_t::const_iterator;
    using const_reverse_iterator = typename storage_t::const_reverse_iterator;
    using size_type = typename storage_t::size_type;

    template <typename... Args>
    AnchorChain(score_t score, Args&&... args)
          : data_(std::forward<Args>(args)...), score_(score) {}

    template <typename... Args>
    constexpr reference& emplace_back(Args&&... args) {
        return data_.emplace_back(std::forward<Args>(args)...);
    }

    constexpr iterator begin() { return data_.begin(); }
    constexpr iterator end() { return data_.end(); }

    constexpr const_iterator begin() const { return data_.begin(); }
    constexpr const_iterator end() const { return data_.end(); }

    constexpr const_reverse_iterator rbegin() const { return data_.rbegin(); }
    constexpr const_reverse_iterator rend() const { return data_.rend(); }

    constexpr const_reference front() const { return data_.front(); }
    constexpr const_reference back() const { return data_.back(); }

    constexpr const_reference operator[](size_t i) const { return data_[i]; }

    constexpr size_type size() const { return data_.size(); }

    constexpr score_t get_score() const { return score_; }

    constexpr size_t get_clipping() const { return front().first->get_clipping(); }
    constexpr size_t get_end_clipping() const { return back().first->get_end_clipping(); }

    constexpr std::string_view get_query_view() const {
        if (data_.empty())
            return {};

        return std::string_view(
            front().first->get_query_view().data(),
            back().first->get_query_view().end() - front().first->get_query_view().begin()
        );
    }

    constexpr bool get_orientation() const {
        return data_.size() ? front().first->get_orientation() : false;
    }

    std::vector<std::pair<const Anchor*, size_t>> data_;
    score_t score_;
};

template <typename Anchor>
using BacktrackStarter = std::function<bool(const AnchorChain<Anchor>&, // chain
                                            const std::vector<score_t>& // score traceback
                                           )>;

template <typename Anchor>
void extend_chain(const AnchorChain<Anchor> &chain,
                  const std::vector<score_t> &score_traceback,
                  const AnchorExtender<Anchor> &anchor_extender,
                  const AlignmentCallback &callback,
                  const std::function<bool()> &terminate = []() { return false; }) {
    auto jt = score_traceback.rbegin();
    std::vector<Alignment> alns;
    anchor_extender(chain.back().first, Alignment(), 0, *jt,
                    [&](Alignment&& aln) { alns.emplace_back(aln); });
    ++jt;

    for (auto it = chain.rbegin(); it + 1 != chain.rend(); ++it, ++jt) {
        assert(jt != score_traceback.rend());
        std::vector<Alignment> next_alns;
        for (auto&& aln : alns) {
            anchor_extender((it + 1)->first, std::move(aln), it->second, *jt,
                [&](Alignment&& next_aln) {
                    next_alns.emplace_back(std::move(next_aln));
                }
            );
        }
        std::swap(next_alns, alns);
    }

    assert(jt == score_traceback.rend());

    for (auto&& aln : alns) {
        if (terminate())
            return;

        callback(std::move(aln));
    }
}

template <typename Anchor>
void chain_anchors(const DBGAlignerConfig &config,
                   const Anchor *anchors_begin,
                   const Anchor *anchors_end,
                   const AnchorConnector<Anchor> &anchor_connector,
                   const BacktrackStarter<Anchor> &start_backtrack
                       = [](const AnchorChain<Anchor>&, const std::vector<score_t>&) { return true; },
                   bool extend_anchors = true,
                   const AnchorExtender<Anchor> &anchor_extender
                       = [](const Anchor*,
                            Alignment&&,
                            size_t,
                            score_t,
                            const AlignmentCallback&) {},
                   const AlignmentCallback &callback = [](Alignment&&) {},
                   const std::function<bool()> &terminate = []() { return false; }) {
    if (terminate() || anchors_begin == anchors_end)
        return;

    ssize_t query_size = anchors_begin->get_clipping() + anchors_begin->get_end_clipping()
                            + anchors_begin->get_query_view().size();

    assert(std::is_sorted(anchors_begin, anchors_end, [&](const Anchor &a,
                                                          const Anchor &b) {
        return std::make_pair(b.get_orientation(), a.get_query_view().end())
             > std::make_pair(a.get_orientation(), b.get_query_view().end());
    }));

    const Anchor *orientation_change = anchors_end;
    ChainScores<Anchor> chain_scores;
    chain_scores.reserve(anchors_end - anchors_begin);
    for (const Anchor *it = anchors_begin; it != anchors_end; ++it) {
        chain_scores.emplace_back(it->get_score(config),
                                  anchors_end,
                                  std::numeric_limits<size_t>::max());
        if (it != anchors_begin && (it - 1)->get_orientation() != it->get_orientation()) {
            assert(it->get_orientation());
            orientation_change = it;
        }
    }

    // forward pass
    ssize_t max_gap_between_anchors = std::min(
        static_cast<ssize_t>(config.max_dist_between_seeds),
        query_size
    );

    auto forward_pass = [&](const Anchor *anchors_begin,
                            const Anchor *anchors_end,
                            auto *chain_scores) {
        if (anchors_begin == anchors_end)
            return;

        auto make_anchor_connector = [&](const Anchor *i) {
            return [&,i](score_t score, const Anchor* last, size_t dist) {
                assert(last != i);
                auto &[max_score, best_last, best_dist] = chain_scores[i - anchors_begin];
                if (std::tie(score, best_dist) > std::tie(max_score, dist)) {
                    max_score = score;
                    best_last = last;
                    best_dist = dist;
                    return true;
                } else {
                    return false;
                }
            };
        };

        if (static_cast<double>(anchors_end - anchors_begin) / query_size
                <= config.chaining_algorithm_switch_cutoff) {
            // if there are fewer seeds, this algorithm is faster
            ssize_t b = max_gap_between_anchors;
            ssize_t b_last;
            score_t best_score = std::get<0>(chain_scores[0]);
            do {
                const Anchor *j = anchors_begin;
                for (const Anchor *i = anchors_begin + 1; i != anchors_end; ++i) {
                    auto end = i->get_query_view().end();
                    j = std::find_if(j, anchors_end, [&](const Anchor &s_j) {
                        return s_j.get_query_view().end() - end <= b;
                    });

                    // align anchor i forwards
                    anchor_connector(*i, b, j, i, chain_scores + (j - anchors_begin),
                                     make_anchor_connector(i));
                    best_score = std::max(best_score, std::get<0>(chain_scores[i - anchors_begin]));
                }
                b_last = b;
                b *= config.max_gap_shrinking_factor;
            } while (best_score < query_size - b_last / 2);
        } else {
            // otherwise, use this algorithm
            for (const Anchor *i = anchors_begin + 1; i != anchors_end; ++i) {
                const Anchor *j = std::max(i - config.chaining_bandwidth, anchors_begin);
                anchor_connector(*i, std::numeric_limits<ssize_t>::max(), j, i,
                                 chain_scores + (j - anchors_begin),
                                 make_anchor_connector(i));
            }
        }
    };

    size_t num_forward = orientation_change - anchors_begin;

    forward_pass(anchors_begin, orientation_change, chain_scores.data());
    forward_pass(orientation_change, anchors_end, chain_scores.data() + num_forward);

    // backtracking
    std::vector<std::tuple<score_t, bool, size_t>> best_chains;
    best_chains.reserve(chain_scores.size());
    for (size_t i = 0; i < chain_scores.size(); ++i) {
        const auto &[score, last, dist] = chain_scores[i];

        if (score > 0)
            best_chains.emplace_back(-score, i >= num_forward, i);
    }

    std::sort(best_chains.begin(), best_chains.end());

    sdsl::bit_vector used(chain_scores.size());
    for (auto [nscore, orientation, i] : best_chains) {
        if (terminate())
            return;

        if (used[i])
            continue;

        AnchorChain<Anchor> chain(-nscore);
        std::vector<score_t> scores;
        const Anchor *last_anchor = anchors_begin + i;
        chain.emplace_back(last_anchor, 0);
        auto [score, last, dist] = chain_scores[i];
        assert(score == -nscore);
        scores.emplace_back(score);
        while (last != anchors_end) {
            last_anchor = last;
            size_t to_traverse = dist;
            std::tie(score, last, dist) = chain_scores[last - anchors_begin];
            chain.emplace_back(last_anchor, to_traverse);
            scores.emplace_back(score);
        }

        if (start_backtrack(chain, scores)) {
            for (const auto &[a_ptr, dist] : chain) {
                used[a_ptr - anchors_begin] = true;
            }

            if (extend_anchors) {
                extend_chain<Anchor>(chain, scores, anchor_extender,
                                     callback, terminate);
            }
        }
    }
}

} // namespace mtg::graph::align

#endif // __ALIGN_CHAIN__
