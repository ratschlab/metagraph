#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__


#include <priority_deque.hpp>

#include "aligner_alignment.hpp"

namespace mtg {
namespace graph {
namespace align {

template <typename NodeType, class AlignmentCompare>
class AlignmentAggregator {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;

    AlignmentAggregator(std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config) {
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
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    boost::container::priority_deque<DBGAlignment,
                                     std::vector<DBGAlignment>,
                                     AlignmentCompare> path_queue_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    if (path_queue_.size() < config_.num_alternative_paths) {
        path_queue_.emplace(std::move(alignment));
    } else if (!AlignmentCompare()(alignment, path_queue_.minimum())) {
        path_queue_.update(path_queue_.begin(), std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(const DBGAlignment &) const -> score_t {
    return path_queue_.size() ? path_queue_.minimum().get_score()
                              : config_.min_path_score;
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback) {
    if (path_queue_.empty())
        return;

    while (path_queue_.size()) {
        callback(DBGAlignment(path_queue_.maximum()));
        path_queue_.pop_maximum();
    }
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
