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

    inline void add_alignment(DBGAlignment&& alignment);

    inline score_t get_min_path_score(const DBGAlignment &seed) const;

    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                         const std::function<bool()> &terminate = []() { return false; });

    size_t size() const {
        size_t size = 0;
        for (const auto &[target, queue] : path_queue_)
            size += queue.size();

        return size;
    }

    size_t num_targets() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    VectorMap<uint64_t, PathQueue>& data() { return path_queue_; }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;

    VectorMap<uint64_t, PathQueue> path_queue_;
    AlignmentCompare cmp_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    auto add_to_target = [&](uint64_t target) {
        DBGAlignment this_alignment(alignment);
        this_alignment.target_columns.assign(1, target);

        auto &cur_queue = path_queue_[target];

        if (cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(std::move(this_alignment));
        } else if (!cmp_.operator()(alignment, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), std::move(this_alignment));
        }
    };

    if (alignment.target_columns.empty()) {
        add_to_target(std::numeric_limits<uint64_t>::max());
    } else {
        std::for_each(alignment.target_columns.begin(),
                      alignment.target_columns.end(),
                      add_to_target);
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(const DBGAlignment &alignment) const -> score_t {
    if (alignment.target_columns.empty()) {
        auto find = path_queue_.find(std::numeric_limits<uint64_t>::max());
        return find == path_queue_.end() || find->second.empty()
            ? config_.min_path_score
            : find->second.minimum().get_score();
    }

    score_t min_score = std::numeric_limits<score_t>::max();
    for (uint64_t target : alignment.target_columns) {
        auto find = path_queue_.find(target);
        if (find == path_queue_.end() || find->second.empty())
            return config_.min_path_score;

        min_score = std::min(min_score, find->second.minimum().get_score());
    }

    return min_score;
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                  const std::function<bool()> &terminate) {
    typedef std::pair<uint64_t, PathQueue> queue_value;
    auto queues = const_cast<std::vector<queue_value>&&>(path_queue_.values_container());

    if (queues.empty())
        return;

    auto cmp = [acmp=AlignmentCompare()](const queue_value &a, const queue_value &b) {
        return b.second.size()
            && (a.second.empty()
                || acmp.operator()(a.second.maximum(), b.second.maximum()));
    };

    std::make_heap(queues.begin(), queues.end(), cmp);

    auto begin = queues.begin();
    auto end = queues.end();
    while (!terminate() && queues.size() && queues[0].second.size()) {
        callback(DBGAlignment(queues[0].second.maximum()));
        queues[0].second.pop_maximum();
        std::pop_heap(queues.begin(), queues.end(), cmp);
        if (queues.back().second.empty()) {
            queues.pop_back();
            if (queues.empty())
                break;

            begin = queues.begin();
            end = queues.end();
        }

        if (--end == begin && begin->second.size()) {
            end = queues.end();
            std::make_heap(begin, end, cmp);
        }
    }

    path_queue_.clear();
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
