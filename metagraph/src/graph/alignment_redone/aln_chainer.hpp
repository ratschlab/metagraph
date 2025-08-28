#pragma once

#include <iostream>

#include "aligner_config.hpp"
#include "aln_query.hpp"
#include "aln_match.hpp"
#include "annotation_buffer.hpp"

#include "common/logger.hpp"

namespace mtg::graph::align_redone {

template <typename AnchorIt>
using ChainScores = std::vector<std::tuple<DBGAlignerConfig::score_t, AnchorIt, size_t, size_t, size_t>>;

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
struct ChainHash {
    inline std::size_t operator()(const AnchorChain<AnchorIt> &chain) const {
        uint64_t hash = 0;

        auto update_hash = [&](auto val) {
            hash ^= val + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        };

        for (const auto &[anchor_it, dist] : chain) {
            std::for_each(anchor_it->get_path().begin(), anchor_it->get_path().end(),
                          update_hash);
            update_hash(dist);
            update_hash(anchor_it->get_coord());
        }
        return hash;
    }
};

template <typename AnchorIt>
struct ChainEqual {
    inline std::size_t operator()(const AnchorChain<AnchorIt> &a,
                                  const AnchorChain<AnchorIt> &b) const {
        if (a.size() != b.size())
            return false;

        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].second != b[i].second)
                return false;
        }

        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].first->get_coord() != b[i].first->get_coord())
                return false;

            const Match &a_i = static_cast<const Match&>(*a[i].first);
            const Match &b_i = static_cast<const Match&>(*b[i].first);
            if (!(a_i == b_i))
                return false;
        }

        return true;
    }
};

template <typename AnchorIt>
bool extend_chain(const AnchorChain<AnchorIt> &chain,
                  const AnchorExtender<AnchorIt> &anchor_extender,
                  const AlignmentCallback &callback,
                  const std::function<bool()> &terminate = []() { return false; }) {
    assert(chain.size());

    std::vector<Alignment> alns;
    alns.emplace_back(*chain.back().first);

    for (auto it = chain.rbegin(); it + 1 != chain.rend(); ++it) {
        std::vector<Alignment> next_alns;
        for (auto&& aln : alns) {
            assert(it->first->get_clipping() == aln.get_clipping());

            size_t dist = it->second;
            anchor_extender(it->first, (it + 1)->first, std::move(aln), dist, 0,
                [&](Alignment&& next_aln) { next_alns.emplace_back(std::move(next_aln)); }
            );
        }

        std::swap(next_alns, alns);
    }

    size_t num_alns = 0;

    for (auto&& aln : alns) {
        if (terminate())
            return num_alns;

        callback(std::move(aln));
        ++num_alns;
    }

    return num_alns;
}

template <typename Anchor>
struct AnchorLess {
    bool operator()(const Anchor &a, const Anchor &b) const {
        return std::make_tuple(a.get_orientation(), a.get_label_class(), b.get_end_clipping(), a.get_clipping())
             < std::make_tuple(b.get_orientation(), b.get_label_class(), a.get_end_clipping(), b.get_clipping());
    }
};

template <typename AnchorIt>
AnchorIt chain_anchors(const Query &query,
                       const DBGAlignerConfig &config,
                       AnchorIt begin,
                       AnchorIt end,
                       const AnchorConnector<AnchorIt> &anchor_connector,
                       const ExtensionStarter<AnchorIt> &extension_starter,
                       const AnchorExtender<AnchorIt> &anchor_extender,
                       const AlignmentCallback &callback = [](Alignment&&) {},
                       const std::function<bool()> &terminate = []() { return false; },
                       AnnotationBuffer *anno_buffer = nullptr) {
    if (terminate() || begin == end)
        return end;

    ssize_t query_size = query.get_query().size();

    assert(std::is_sorted(begin, end, AnchorLess<Anchor>()));

    ChainScores<AnchorIt> chain_scores;
    chain_scores.reserve(end - begin);
    std::transform(begin, end, std::back_inserter(chain_scores), [&](const Anchor &a) {
        return std::make_tuple(score_match(a, config),
                               end, 0, a.get_spelling().size(), a.get_seed().size());
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
                auto &[max_score, best_last, best_dist, dist_traversed, queried] = *(chain_scores + (i - begin));
                if (std::tie(score, best_dist) > std::tie(max_score, dist)) {
                    if (last != end) {
                        dist_traversed = std::get<3>(*(chain_scores + (last - begin))) + dist;
                        assert(i->get_seed().end() >= last->get_seed().end());
                        queried = std::get<4>(*(chain_scores + (last - begin))) + (i->get_seed().end() - last->get_seed().end());
                    }

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
            size_t traversed = std::get<3>(chain_scores[0]);
            size_t queried = begin->get_seed().size();
            ssize_t best_cost = std::numeric_limits<ssize_t>::max();
            DBGAlignerConfig::score_t match_score = config.match_score("A");
            do {
                AnchorIt j = begin;
                for (AnchorIt i = begin; i != end; ++i) {
                    j = std::find_if(j, end, [&](const Anchor &a) {
                        return i->get_seed().end() - a.get_seed().end() <= b;
                    });

                    // align anchors [j, i) to i
                    anchor_connector(*i, b, j, i, chain_scores + (j - begin),
                                     make_anchor_connector(i));

                    if (std::get<0>(chain_scores[i - begin]) > best_score) {
                        AnchorIt last;
                        size_t added_dist;
                        std::tie(best_score, last, added_dist, traversed, queried)
                            = chain_scores[i - begin];
                    }
                }
                b_last = b;
                b *= config.gap_shrinking_factor;

                // best_score = 1/2(match_score(query_size + match_spelling_size - best_cost))
                // => 2*best_score/match_score - query_size - match_spelling_size = -best_cost
                // => best_cost = match_spelling_size + query_size - 2*best_score/match_score

                best_cost = traversed + queried - 2 * best_score / match_score;
            } while (best_cost > b_last);
        } else {
            common::logger->trace("Using fallback chaining with bandwidth {}", config.chaining_bandwidth);
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
    std::vector<std::tuple<DBGAlignerConfig::score_t, size_t, bool, size_t>> best_chains;
    best_chains.reserve(chain_scores.size());
    for (size_t i = 0; i < chain_scores.size(); ++i) {
        const auto &[score, last, dist, traversed, queried] = chain_scores[i];

        if (score > 0) {
            assert(last <= end);
            assert(last >= begin);
            best_chains.emplace_back(
                -score,
                (begin + i)->get_end_clipping(),
                (begin + i)->get_orientation(),
                i
            );
        }
    }

    std::sort(best_chains.begin(), best_chains.end());

    sdsl::bit_vector used(chain_scores.size());

    tsl::hopscotch_set<Anchor::label_class_t> label_classes;
    AnchorChain<AnchorIt> last_chain;

    auto flush_chains = [&]() {
        if (last_chain.size()) {
            if (label_classes.size() == 1) {
                extend_chain<AnchorIt>(last_chain, anchor_extender, callback, terminate);
            } else {
                AnnotationBuffer::Columns merged_columns;
                for (auto label_class : label_classes) {
                    if (merged_columns.empty()) {
                        merged_columns = anno_buffer->get_cached_column_set(label_class);
                    } else {
                        const auto &next_columns = anno_buffer->get_cached_column_set(label_class);
                        AnnotationBuffer::Columns next_merged;
                        std::set_union(merged_columns.begin(), merged_columns.end(),
                                        next_columns.begin(), next_columns.end(),
                                        std::back_inserter(next_merged));
                        std::swap(next_merged, merged_columns);
                    }
                }
                auto merged_label_class = anno_buffer->cache_column_set(std::move(merged_columns));
                std::vector<Anchor> relabeled_anchors;
                relabeled_anchors.reserve(last_chain.size());
                for (const auto &[it, dist] : last_chain) {
                    relabeled_anchors.emplace_back(*it).set_label_class(merged_label_class);
                }
                for (size_t i = 0; i < last_chain.size(); ++i) {
                    static_assert(std::is_same_v<AnchorIt, std::vector<Anchor>::iterator>);
                    last_chain[i].first = relabeled_anchors.begin() + i;
                }
                extend_chain<AnchorIt>(last_chain, anchor_extender, callback, terminate);
            }
        }

        last_chain.clear();
        label_classes.clear();
    };

    for (auto [nscore, end_clipping, orientation, i] : best_chains) {
        if (terminate())
            return end;

        if (used[i])
            continue;

        AnchorChain<AnchorIt> chain;
        std::vector<DBGAlignerConfig::score_t> scores;
        AnchorIt last_anchor = begin + i;
        auto [score, last, dist, traversed, queried] = chain_scores[i];
        chain.emplace_back(last_anchor, dist);
        assert(score == -nscore);
        scores.emplace_back(score);
        while (last != end) {
            last_anchor = last;
            std::tie(score, last, dist, traversed, queried) = chain_scores[last - begin];
            assert(last_anchor != end);
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
            if (!anno_buffer) {
                extend_chain<AnchorIt>(chain, anchor_extender, callback, terminate);
            } else {
                if (label_classes.empty() || !ChainEqual<AnchorIt>()(last_chain, chain)) {
                    flush_chains();
                    last_chain = chain;
                }
                label_classes.emplace(chain[0].first->get_label_class());
            }
            for (const auto &[a_ptr, dist] : chain) {
                used[a_ptr - begin] = true;
            }
        }
    }

    flush_chains();

    return end;
}

} // namespace mtg::graph::align