#pragma once

#include "aligner_config.hpp"
#include "aln_query.hpp"
#include "aln_match.hpp"

namespace mtg::graph::align {

template <typename AnchorIt>
using ChainScores = std::vector<std::tuple<DBGAlignerConfig::score_t, AnchorIt, size_t>>;

template <typename AnchorIt>
using ScoreUpdater = std::function<bool(DBGAlignerConfig::score_t,       // connect score
                                        AnchorIt, // last
                                        size_t         // distance
                                       )>;

template <typename AnchorIt>
using AnchorChain = std::vector<std::pair<AnchorIt, size_t>>;

template <typename AnchorIt>
using ExtensionStarter = std::function<bool(const AnchorChain<AnchorIt>&, // chain
                                            const std::vector<DBGAlignerConfig::score_t>& // score traceback
                                           )>;


template <typename AnchorIt>
using AnchorConnector = std::function<void(const Anchor&, // start anchor
                                           ssize_t b,     // max dist
                                           AnchorIt, // target anchors begin
                                           AnchorIt, // target anchors end
                                           typename ChainScores<AnchorIt>::iterator,
                                           const ScoreUpdater<AnchorIt>&
                                          )>;

using AlignmentCallback = std::function<void(Alignment&&)>;

template <typename AnchorIt>
using AnchorExtender = std::function<void(AnchorIt, // first ptr,
                                          AnchorIt, // prev ptr
                                          Alignment&&,   // target,
                                          size_t,        // distance,
                                          DBGAlignerConfig::score_t,       // chain score up to this point
                                          const AlignmentCallback&)>;

template <typename AnchorIt>
void extend_chain(const AnchorChain<AnchorIt> &chain,
                  const std::vector<DBGAlignerConfig::score_t> &score_traceback,
                  const AnchorExtender<AnchorIt> &anchor_extender,
                  const AlignmentCallback &callback,
                  const std::function<bool()> &terminate = []() { return false; }) {
    std::vector<Alignment> alns;
    alns.emplace_back(*chain.back().first);

    for (auto it = chain.rbegin() + 1; it != chain.rend(); ++it) {
        std::vector<Alignment> next_alns;
        for (auto&& aln : alns) {
            size_t dist = (it - 1)->second;
            anchor_extender((it - 1)->first, it->first, std::move(aln), dist, 0,
                [&](Alignment&& next_aln) { next_alns.emplace_back(std::move(next_aln)); }
            );
        }
        std::swap(next_alns, alns);
    }
    // auto jt = score_traceback.begin();
    // std::vector<Alignment> alns;
    // alns.emplace_back(*chain[0].first);
    // ++jt;

    // for (auto it = chain.begin() + 1; it != chain.end(); ++it, ++jt) {
    //     assert(jt != score_traceback.end());
    //     std::vector<Alignment> next_alns;
    //     for (auto&& aln : alns) {
    //         anchor_extender(it->first, (it - 1)->first, std::move(aln), it->second, *jt,
    //             [&](Alignment&& next_aln) {
    //                 next_alns.emplace_back(std::move(next_aln));
    //             }
    //         );
    //     }
    //     std::swap(next_alns, alns);
    // }

    assert(jt == score_traceback.end());

    for (auto&& aln : alns) {
        if (terminate())
            return;

        callback(std::move(aln));
    }
}

template <typename AnchorIt>
void chain_anchors(const Query &query,
                   const DBGAlignerConfig &config,
                   AnchorIt begin,
                   AnchorIt end,
                   const AnchorConnector<AnchorIt> &anchor_connector,
                   const ExtensionStarter<AnchorIt> &extension_starter,
                   const AnchorExtender<AnchorIt> &anchor_extender,
                   const AlignmentCallback &callback = [](Alignment&&) {},
                   const std::function<bool()> &terminate = []() { return false; }) {
    if (terminate() || begin == end)
        return;

    ssize_t query_size = query.get_query().size();

    std::sort(begin, end, [&](const Anchor &a, const Anchor &b) {
        return std::make_tuple(a.get_orientation(), a.get_label_class(), b.get_end_clipping())
             < std::make_tuple(b.get_orientation(), b.get_label_class(), a.get_end_clipping());
    });

    ChainScores<AnchorIt> chain_scores;
    chain_scores.reserve(end - begin);
    std::transform(begin, end, std::back_inserter(chain_scores), [&](const Anchor &a) {
        return std::make_tuple(a.get_score(), end, 0);
    });

    // forward pass
    ssize_t max_gap_between_anchors = std::min(
        static_cast<ssize_t>(config.max_dist_between_seeds),
        query_size
    );

    auto forward_pass = [&](AnchorIt begin,
                            AnchorIt end,
                            typename ChainScores<AnchorIt>::iterator chain_scores) {
        if (begin == end)
            return;

        auto make_anchor_connector = [&](const AnchorIt i) {
            return [&,i](DBGAlignerConfig::score_t score, const AnchorIt last, size_t dist) {
                assert(last != i);
                auto &[max_score, best_last, best_dist] = *(chain_scores + (i - begin));
                if (score > max_score) {
                    max_score = score;
                    best_last = last;
                    best_dist = dist;
                    return true;
                } else {
                    return false;
                }
            };
        };

        if (static_cast<double>(end - begin) / query_size
                <= config.chaining_algorithm_switch_cutoff) {
            // if there are fewer seeds, this algorithm is faster
            ssize_t b = max_gap_between_anchors;
            ssize_t b_last;
            DBGAlignerConfig::score_t best_score = std::get<0>(chain_scores[0]);
            ssize_t best_cost = std::numeric_limits<ssize_t>::max();
            DBGAlignerConfig::score_t match_score = config.match_score("A");
            do {
                AnchorIt j = begin;
                for (AnchorIt i = begin; i != end; ++i) {
                    j = std::find_if(j, end, [&](const Anchor &a) {
                        return a.get_end_clipping() <= b + i->get_end_clipping();
                    });

                    // align anchors [j, i) to i
                    anchor_connector(*i, b, j, i, chain_scores + (j - begin),
                                     make_anchor_connector(i));
                    best_score = std::max(best_score, std::get<0>(chain_scores[i - begin]));
                }
                b_last = b;
                b *= config.gap_shrinking_factor;

                // best_score = 1/2(match_score(query_size + match_spelling_size - best_cost))
                // => 2*best_score/match_score - query_size - match_spelling_size = -best_cost
                // => best_cost = match_spelling_size + query_size - 2*best_score/match_score

                // if we assume that match_spelling_size == query_size
                best_cost = 2*(query_size - best_score/match_score);
            } while (best_cost > b_last);
        } else {
            // otherwise, use this algorithm
            for (AnchorIt i = begin + 1; i != end; ++i) {
                AnchorIt j = i - std::min(static_cast<size_t>(i - begin), config.chaining_bandwidth);
                anchor_connector(*i, std::numeric_limits<ssize_t>::max(), j, i,
                                 chain_scores + (j - begin),
                                 make_anchor_connector(i));
            }
        }
    };

    auto cur_begin = begin;
    auto jt = chain_scores.begin();
    for (auto it = begin + 1; it != end; ++it) {
        if (it->get_orientation() != (it - 1)->get_orientation()
                || it->get_label_class() != (it - 1)->get_label_class()) {
            forward_pass(cur_begin, it, jt);
            jt += it - cur_begin;
            cur_begin = it;
        }
    }

    forward_pass(cur_begin, end, jt);

    // backtracking
    std::vector<std::tuple<DBGAlignerConfig::score_t, bool, size_t>> best_chains;
    best_chains.reserve(chain_scores.size());
    for (size_t i = 0; i < chain_scores.size(); ++i) {
        const auto &[score, last, dist] = chain_scores[i];

        if (score > 0)
            best_chains.emplace_back(-score, (begin + i)->get_orientation(), i);
    }

    std::sort(best_chains.begin(), best_chains.end());

    sdsl::bit_vector used(chain_scores.size());
    for (auto [nscore, orientation, i] : best_chains) {
        if (terminate())
            return;

        if (used[i])
            continue;

        AnchorChain<AnchorIt> chain;
        std::vector<DBGAlignerConfig::score_t> scores;
        AnchorIt last_anchor = begin + i;
        auto [score, last, dist] = chain_scores[i];
        chain.emplace_back(last_anchor, dist);
        assert(score == -nscore);
        scores.emplace_back(score);
        while (last != end) {
            last_anchor = last;
            std::tie(score, last, dist) = chain_scores[last - begin];
            chain.emplace_back(last_anchor, dist);
            scores.emplace_back(score);
            if (used[last_anchor - begin]) {
                chain.clear();
                scores.clear();
                break;
            }
        }

        std::reverse(chain.begin(), chain.end());
        std::reverse(scores.begin(), scores.end());

        if (chain.size() && extension_starter(chain, scores)) {
            for (const auto &[a_ptr, dist] : chain) {
                used[a_ptr - begin] = true;
            }

            extend_chain<AnchorIt>(chain, scores, anchor_extender, callback, terminate);
        }
    }
}

} // namespace mtg::graph::align