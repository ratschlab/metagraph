#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__


#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "common/vector_map.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType, class AlignmentCompare>
class AlignmentAggregator {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;
    typedef boost::container::priority_deque<DBGAlignment,
                                             std::vector<DBGAlignment>,
                                             AlignmentCompare> PathQueue;

    AlignmentAggregator(std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(DBGAlignment&& alignment);

    score_t get_min_path_score() const;
    score_t get_max_path_score() const;

    score_t get_min_path_score(const DBGAlignment &) const {
        return get_min_path_score();
    }

    score_t get_max_path_score(const DBGAlignment &) const {
        return get_max_path_score();
    }


    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                         const std::function<bool()> &terminate = []() { return false; });

    size_t size() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    void clear() { path_queue_.clear(); }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    PathQueue path_queue_;
    AlignmentCompare cmp_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    if (std::find(path_queue_.begin(), path_queue_.end(), alignment) != path_queue_.end())
        return;

    if (path_queue_.size() < config_.num_alternative_paths) {
        path_queue_.emplace(std::move(alignment));
    } else if (!cmp_(alignment, path_queue_.minimum())) {
        path_queue_.update(path_queue_.begin(), std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score() const -> score_t {
    return path_queue_.size() < config_.num_alternative_paths
        ? config_.min_path_score
        : std::max(static_cast<score_t>(path_queue_.maximum().get_score() * config_.fraction_of_top),
                   path_queue_.minimum().get_score());
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_max_path_score() const -> score_t {
    return path_queue_.size() ? path_queue_.maximum().get_score() : config_.min_path_score;
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                  const std::function<bool()> &terminate) {
    while (!terminate() && path_queue_.size()) {
        callback(DBGAlignment(path_queue_.maximum()));
        path_queue_.pop_maximum();
    }
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
