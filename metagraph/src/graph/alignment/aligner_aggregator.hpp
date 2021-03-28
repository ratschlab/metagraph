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
    struct SharedPtrCmp {
        bool operator()(const std::shared_ptr<Alignment<NodeType>> &a,
                        const std::shared_ptr<Alignment<NodeType>> &b) const {
            return base_cmp_(*a, *b);
        }

        AlignmentCompare base_cmp_;
    };

  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;
    typedef boost::container::priority_deque<std::shared_ptr<DBGAlignment>,
                                             std::vector<std::shared_ptr<DBGAlignment>>,
                                             SharedPtrCmp> PathQueue;

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
    SharedPtrCmp cmp_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    auto packaged_alignment = std::make_shared<DBGAlignment>(std::move(alignment));

    auto add_to_target = [&](uint64_t target) {
        auto &cur_queue = path_queue_[target];

        if (cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(packaged_alignment);
        } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), packaged_alignment);
        }
    };

    if (packaged_alignment->target_columns.empty()) {
        add_to_target(std::numeric_limits<uint64_t>::max());
    } else {
        std::for_each(packaged_alignment->target_columns.begin(),
                      packaged_alignment->target_columns.end(),
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
            : find->second.minimum()->get_score();
    }

    score_t min_score = std::numeric_limits<score_t>::max();
    for (uint64_t target : alignment.target_columns) {
        auto find = path_queue_.find(target);
        if (find == path_queue_.end() || find->second.empty())
            return config_.min_path_score;

        min_score = std::min(min_score, find->second.minimum()->get_score());
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

    auto cmp = [this](const queue_value &a, const queue_value &b) {
        return b.second.size()
            && (a.second.empty() || cmp_(a.second.maximum(), b.second.maximum()));
    };

    std::make_heap(queues.begin(), queues.end(), cmp);

    auto begin = queues.begin();
    auto end = queues.end();
    std::shared_ptr<DBGAlignment> last_alignment;
    while (!terminate() && queues.size() && queues[0].second.size()) {
        if (queues[0].second.maximum().get() != last_alignment.get()) {
            if (last_alignment)
                callback(std::move(*last_alignment));

            last_alignment = queues[0].second.maximum();
        }

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

    if (last_alignment)
        callback(std::move(*last_alignment));

    path_queue_.clear();
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
