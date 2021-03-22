#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__


#include <priority_deque.hpp>
#include <tsl/hopscotch_map.h>

#include "aligner_alignment.hpp"

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

    inline void add_alignment(DBGAlignment&& alignment);

    inline score_t get_min_path_score(const DBGAlignment &seed) const;

    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback);

    size_t size() const {
        size_t size = 0;
        for (const auto &[target, queue] : path_queue_)
            size += queue.size();

        return size;
    }

    size_t num_targets() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    const tsl::hopscotch_map<uint64_t, PathQueue>& data() const { return path_queue_; }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;

    tsl::hopscotch_map<uint64_t, PathQueue> path_queue_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    auto &cur_queue = path_queue_[alignment.target_column];

    if (cur_queue.size() < config_.num_alternative_paths) {
        cur_queue.emplace(std::move(alignment));
    } else if (!AlignmentCompare()(alignment, cur_queue.minimum())) {
        cur_queue.update(cur_queue.begin(), std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(const DBGAlignment &alignment) const -> score_t {
    auto find = path_queue_.find(alignment.target_column);
    return find == path_queue_.end() || find->second.empty()
        ? config_.min_path_score
        : find->second.minimum().get_score();
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback) {
    for (auto it = path_queue_.begin(); it != path_queue_.end(); ++it) {
        while (it.value().size()) {
            callback(DBGAlignment(it.value().maximum()));
            it.value().pop_maximum();
        }
    }
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
