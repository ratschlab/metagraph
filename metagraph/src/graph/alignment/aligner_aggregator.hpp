#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__


#include <priority_deque.hpp>

#include "aligner_helper.hpp"

namespace mtg {
namespace graph {
namespace align {

template <typename NodeType, class AlignmentCompare>
class AlignmentAggregator {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;

    AlignmentAggregator(const std::string_view query, const std::string_view rc_query,
                        const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), graph_(graph), config_(config),
            best_score_(config_.chain_alignments ? query.size() : 1,
                        config_.min_path_score) {
        assert(config_.num_alternative_paths);
    }

    inline void add_alignment(DBGAlignment&& alignment);

    inline score_t get_min_path_score(const DBGAlignment &seed) const;

    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback);

    size_t size() const { return path_queue_.size(); }
    bool empty() const { return path_queue_.empty(); }

  private:
    const std::string_view query_;
    const std::string_view rc_query_;
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    boost::container::priority_deque<DBGAlignment,
                                     std::vector<DBGAlignment>,
                                     AlignmentCompare> path_queue_;
    std::vector<score_t> best_score_;

    void call_alignment_chains(const std::vector<DBGAlignment> &paths,
                               const std::function<void(std::vector<DBGAlignment>&&,
                                                        score_t)> &callback,
                               size_t i = 0,
                               std::vector<DBGAlignment>&& chain = {},
                               score_t cur_score = 0,
                               std::shared_ptr<std::vector<score_t>> best_score = {}) const;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    if (!config_.chain_alignments) {
        // no alignment chaining
        if (path_queue_.size() < config_.num_alternative_paths) {
            path_queue_.emplace(std::move(alignment));
        } else if (alignment.get_score() > path_queue_.minimum().get_score()) {
            path_queue_.update(path_queue_.begin(), std::move(alignment));
        }
    } else {
        // alignment chaining
        bool updated = false;
        if (!alignment.get_orientation()) {
            AlignmentSuffix suffix(alignment, config_, graph_.get_k());
            const std::string_view aln_query = alignment.get_query();
            auto start = best_score_.begin() + (aln_query.data() - query_.data());
            auto end = start + aln_query.size();
            for (auto it = start; it != end; ++it) {
                if (suffix.get_score() > *it) {
                    *it = suffix.get_score();
                    updated = true;
                }
                ++suffix;
                while (!suffix.eof() && suffix.get_front_op() == Cigar::Operator::INSERTION) {
                    ++suffix;
                }
            }
        } else {
            AlignmentSuffix suffix(alignment, config_, graph_.get_k());
            const std::string_view aln_query = alignment.get_query();
            auto start = best_score_.rbegin() + (aln_query.data() - rc_query_.data());
            auto end = start + aln_query.size();
            for (auto it = start; it != end; ++it) {
                if (suffix.get_score() > *it) {
                    *it = suffix.get_score();
                    updated = true;
                }
                ++suffix;
                while (!suffix.eof() && suffix.get_front_op() == Cigar::Operator::INSERTION) {
                    ++suffix;
                }
            }
        }

        if (updated)
            path_queue_.emplace(std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(const DBGAlignment &seed) const -> score_t {
    if (!config_.chain_alignments) {
        return path_queue_.size() ? path_queue_.minimum().get_score()
                                  : config_.min_path_score;
    } else if (!seed.get_orientation()) {
        const std::string_view seed_query = seed.get_query();
        auto start = best_score_.begin() + (seed_query.data() - query_.data());
        return *std::min_element(start, start + seed_query.size());
    } else {
        const std::string_view seed_query = seed.get_query();
        auto start = best_score_.rbegin() + (seed_query.data() - rc_query_.data());
        return *std::min_element(start, start + seed_query.size());
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback) {
    if (path_queue_.empty())
        return;

    if (!config_.chain_alignments || path_queue_.size() == 1) {
        // no chaining
        while (path_queue_.size()) {
            callback(DBGAlignment(path_queue_.maximum()));
            path_queue_.pop_maximum();
        }

        return;
    }

    mtg::common::logger->trace("Chaining {} local alignment(s)", path_queue_.size());

    // heap sort
    std::vector<DBGAlignment> cur_paths;
    DBGAlignment dummy;
    const auto *last_path = &dummy;

    while (path_queue_.size()) {
        const auto &top_path = path_queue_.maximum();
        mtg::common::logger->trace("Alignment: {}\tOrientation: {}\tOffset: {}",
                                   top_path.get_cigar().to_string(),
                                   top_path.get_orientation(),
                                   top_path.get_offset());

        if (!top_path.get_clipping() && !top_path.get_end_clipping()) {
            callback(DBGAlignment(top_path));
            return;
        }

        if (top_path == *last_path) {
            mtg::common::logger->trace("\tSkipping duplicate alignment");
        } else if (top_path.get_orientation()) {
            // TODO
            mtg::common::logger->trace("\tSkipping reverse alignment");
        } else if (top_path.get_query().size() < graph_.get_k()) {
            mtg::common::logger->trace("\tSkipping short alignment");
        } else {
            cur_paths.emplace_back(top_path);
            last_path = &cur_paths.back();
        }

        path_queue_.pop_maximum();
    }

    if (cur_paths.size() <= 1) {
        mtg::common::logger->trace("Nothing to chain: {} path(s) selected",
                                   cur_paths.size());
        for (auto&& path : cur_paths) {
            callback(std::move(path));
        }

        return;
    }

    std::sort(cur_paths.begin(), cur_paths.end(),
              [&](const auto &a, const auto &b) {
                  return a.get_query_end() < b.get_query_end()
                      || (a.get_query_end() == b.get_query_end() && a.get_query().data() < b.get_query().data());
              });

    std::pair<score_t, std::vector<DBGAlignment>> cur_chain(config_.min_path_score, {});
    call_alignment_chains(cur_paths, [&](std::vector<DBGAlignment>&& chain,
                                         score_t score) {
        if (score > cur_chain.first) {
            for (auto &path : chain) {
                path.extend_query_end(query_.data() + query_.size());
            }
            cur_chain = std::make_pair(score, std::move(chain));
        }
    });

    if (cur_chain.first <= config_.min_path_score) {
        mtg::common::logger->trace("No good chain found");
    } else {
        mtg::common::logger->trace("Best chain score: {}", cur_chain.first);
        assert(cur_chain.second.size());

        for (auto&& path : cur_chain.second) {
            callback(std::move(path));
        }
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignment_chains(const std::vector<DBGAlignment> &paths,
                        const std::function<void(std::vector<DBGAlignment>&&,
                                                 score_t)> &callback,
                        size_t i,
                        std::vector<DBGAlignment>&& chain,
                        score_t cur_score,
                        std::shared_ptr<std::vector<score_t>> best_score) const {
    if (i >= paths.size()) {
        callback(std::move(chain), cur_score);
        return;
    }

    if (!best_score)
        best_score = std::make_shared<std::vector<score_t>>(paths.size(), cur_score);

    if (chain.empty()) {
        assert(!i);
        for (auto it = paths.begin(); it != paths.end(); ++it) {
            if (cur_score + it->get_score() > best_score->at(it - paths.begin())) {
                best_score->at(it - paths.begin()) = cur_score + it->get_score();
                call_alignment_chains(paths, callback, it - paths.begin() + 1,
                                      { *it }, cur_score + it->get_score(), best_score);
            }
        }

        return;
    }

    auto it = std::lower_bound(paths.begin() + i, paths.end(), chain.back(),
                               [&](const auto &a, const auto &b) {
                                   return a.get_query_end() < b.get_query_end();
                               });

    for ( ; it != paths.end(); ++it) {
        assert(it->get_query_end() >= chain.back().get_query_end());

        // skip this alignment if it covers the same region as the last one
        if (it->get_query().data() == chain.back().get_query().data()
                && it->get_query_end() == chain.back().get_query_end())
            continue;

        if (it->get_query().data() < chain.back().get_query_end()) {
            // the alignments overlap
            auto [first, second] = DBGAlignment::get_best_overlap(
                chain.back(), *it, graph_, config_
            );

            if (first.empty() || second.empty())
                continue;

            mtg::common::logger->trace("Chain: {}\n{}", first, second);

            assert(first.is_valid(graph_, &config_));
            assert(second.is_valid(graph_, &config_));

            score_t score = cur_score - chain.back().get_score()
                + first.get_score() + second.get_score();

            if (score <= best_score->at(it - paths.begin()))
                continue;

            best_score->at(it - paths.begin()) = score;

            auto next_chain = const_cast<const std::vector<DBGAlignment>&>(chain);
            next_chain.back() = std::move(first);
            next_chain.push_back(std::move(second));

            call_alignment_chains(paths, callback, it - paths.begin() + 1,
                                  std::move(next_chain), score, best_score);

        } else {
            // select this alignment if it leads to a better score
            score_t score = cur_score + it->get_score();
            if (score > best_score->at(it - paths.begin())) {
                auto next_chain = const_cast<const std::vector<DBGAlignment>&>(chain);
                next_chain.push_back(*it);
                best_score->at(it - paths.begin()) = score;
                call_alignment_chains(paths, callback, it - paths.begin() + 1,
                                      std::move(next_chain), score, best_score);
            }
        }
    }

    if (cur_score > best_score->back()) {
        best_score->back() = cur_score;
        call_alignment_chains(paths, callback, paths.size(), std::move(chain),
                              cur_score, best_score);
    }
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
