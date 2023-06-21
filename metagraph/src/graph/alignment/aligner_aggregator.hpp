#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>
#include <tsl/hopscotch_map.h>

#include "alignment.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename T, class Container = std::vector<T>, class Compare = std::less<T>>
class PriorityDeque : public boost::container::priority_deque<T, Container, Compare> {
  public:
    Container& data() { return this->sequence(); }
    const Container& data() const { return this->sequence(); }
};


template <class AlignmentCompare>
class AlignmentAggregator {
    typedef std::shared_ptr<Alignment> value_type;

    struct ValCmp {
        bool operator()(const value_type &a, const value_type &b) const {
            return base_cmp_(*a, *b);
        }

        AlignmentCompare base_cmp_;
    };

  public:
    typedef Alignment::score_t score_t;
    typedef Alignment::Column Column;
    typedef PriorityDeque<value_type, std::vector<value_type>, ValCmp> PathQueue;

    explicit AlignmentAggregator(const DBGAlignerConfig &config) : config_(config) {
        assert(config_.num_alternative_paths);
    }

    bool add_alignment(Alignment&& alignment);

    score_t get_global_cutoff() const;

    std::vector<Alignment> get_alignments();

    size_t num_aligned_labels() const { return path_queue_.size(); }

    void clear() {
        path_queue_.clear();
        best_alignment_.reset();
    }

  private:
    const DBGAlignerConfig &config_;
    tsl::hopscotch_map<Column, PathQueue> path_queue_;
    value_type best_alignment_;
    ValCmp cmp_;
};

// return true if the alignment was added
template <class AlignmentCompare>
inline bool AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    // first, wrap the alignment so that duplicates are not stored in each per-label queue
    auto a = std::make_shared<Alignment>(std::move(alignment));
    if (!best_alignment_ || cmp_(best_alignment_, a))
        best_alignment_ = a;

    if (a->label_columns.empty()) {
        path_queue_[std::numeric_limits<Column>::max()].emplace(a);

    } else {
        for (Column column : a->label_columns) {
            path_queue_[column].emplace(a);
        }
    }

    return true;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_global_cutoff() const -> score_t {
    if (!best_alignment_)
        return config_.ninf;

    score_t cur_max = best_alignment_->get_score();
    return cur_max > 0 ? cur_max * config_.rel_score_cutoff : cur_max;
}

template <class AlignmentCompare>
inline std::vector<Alignment> AlignmentAggregator<AlignmentCompare>::get_alignments() {
    if (!best_alignment_) {
        assert(path_queue_.empty());
        return {};
    }

    std::vector<value_type> alignment_ptrs;
    size_t max_num_alignments = config_.post_chain_alignments
        ? std::numeric_limits<size_t>::max()
        : config_.num_alternative_paths;

    for (auto it = path_queue_.begin(); it != path_queue_.end(); ++it) {
        auto &queue = it.value();
        size_t added = 0;
        while (queue.size() && added < max_num_alignments) {
            alignment_ptrs.emplace_back(queue.maximum());
            queue.pop_maximum();
            ++added;
        }
    }

    std::sort(alignment_ptrs.begin(), alignment_ptrs.end(), cmp_);

    std::vector<Alignment> alignments;
    std::for_each(alignment_ptrs.rbegin(), alignment_ptrs.rend(), [&](value_type &aln_ptr) {
        assert(aln_ptr);
        if (!aln_ptr->empty()) {
            alignments.emplace_back(std::move(*aln_ptr));
            *aln_ptr = Alignment();
        }
    });

    return alignments;
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
