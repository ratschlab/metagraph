#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "alignment.hpp"
#include "common/algorithms.hpp"
#include "common/vector_map.hpp"
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
    struct ValCmp {
        bool operator()(const std::shared_ptr<Alignment> &a,
                        const std::shared_ptr<Alignment> &b) const {
            return base_cmp_(*a, *b);
        }

        AlignmentCompare base_cmp_;
    };

  public:
    typedef Alignment::score_t score_t;
    typedef Alignment::Column Column;
    typedef PriorityDeque<std::shared_ptr<Alignment>,
                          std::vector<std::shared_ptr<Alignment>>, ValCmp> PathQueue;

    explicit AlignmentAggregator(const DBGAlignerConfig &config) : config_(config) {
        assert(config_.num_alternative_paths);
    }

    bool add_alignment(Alignment&& alignment);

    score_t get_global_cutoff() const;
    score_t get_score_cutoff(const Vector<Column> &labels) const;

    std::vector<Alignment> get_alignments();

    size_t num_aligned_labels() const { return path_queue_.size(); }

    void clear() { path_queue_.clear(); unlabeled_.clear(); }

  private:
    const DBGAlignerConfig &config_;
    VectorMap<Column, PathQueue> path_queue_;
    PathQueue unlabeled_;
    ValCmp cmp_;

    score_t get_label_cutoff(Column label) const;
};

// return true if the alignment was added
template <class AlignmentCompare>
inline bool AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    // first, wrap the alignment so that duplicates are not stored in each per-label queue
    auto a = std::make_shared<Alignment>(std::move(alignment));

    // if nothing has been added to the queue so far, add the alignment
    if (unlabeled_.empty()) {
        unlabeled_.emplace(a);
        for (Column c : a->label_columns) {
            path_queue_[c].emplace(a);
        }
        return true;
    }

    // if the score is less than the cutoff, don't add it
    if (a->get_score() < get_global_cutoff())
        return false;

    // helper for adding alignments to the queue
    auto push_to_queue = [&](auto &queue) {
        // check for duplicates
        for (const auto &aln : queue) {
            if (*a == *aln)
                return false;
        }
        // If post-alignment chaining is requested, never skip any alignments
        if (config_.post_chain_alignments || queue.size() < config_.num_alternative_paths) {
            queue.emplace(a);
            return true;
        }
        // the queue is full
        assert(queue.size() == config_.num_alternative_paths);
        if (cmp_(a, queue.minimum()))
            return false;

        queue.update(queue.begin(), a);
        return true;
    };

    // if we are in the unlabeled case, only consider the global queue
    if (a->label_columns.empty())
        return push_to_queue(unlabeled_);

    // if an incoming alignment has labels, and we haven't encountered a labeled
    // alignment yet, we only need the ncol queue for fetching the global minimum,
    // so shrink it to only one element
    if (path_queue_.empty()) {
        if (unlabeled_.size() > 1) {
            // maximum is stored at begin+1
            auto max = std::move(*(unlabeled_.begin() + 1));
            unlabeled_.clear();
            unlabeled_.push(std::move(max));
        }
    }
    assert(unlabeled_.size() == 1);

    // add the alignment to its labeled queues
    bool added = false;
    for (Column c : a->label_columns) {
        added |= push_to_queue(path_queue_[c]);
    }

    if (!added)
        return false;

    // TODO: maintain a pointer to the best alignment
    // if this is the best alignment so far, update the global queue
    if (!cmp_(a, unlabeled_.maximum()))
        unlabeled_.update(unlabeled_.begin(), a);

    return true;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_global_cutoff() const -> score_t {
    if (unlabeled_.empty())
        return config_.ninf;

    score_t cur_max = unlabeled_.maximum()->get_score();

    return cur_max > 0 ? cur_max * config_.rel_score_cutoff : cur_max;
}

// TODO: define it the same way as in get_global_cutoff()?
template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_score_cutoff(const Vector<Column> &labels) const -> score_t {
    assert(labels.size());

    score_t global_min = get_global_cutoff();

    score_t min_score = std::numeric_limits<score_t>::max();
    for (Column label : labels) {
        min_score = std::min(min_score, get_label_cutoff(label));
        if (min_score < global_min)
            return global_min;
    }
    return min_score;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_label_cutoff(Column label) const -> score_t {
    auto find = path_queue_.find(label);
    return find == path_queue_.end()
            || find->second.size() < config_.num_alternative_paths
            || config_.post_chain_alignments
        ? config_.ninf
        : find->second.minimum()->get_score();
}

template <class AlignmentCompare>
inline std::vector<Alignment> AlignmentAggregator<AlignmentCompare>::get_alignments() {
    // move all alignments to one vector
    std::vector<std::shared_ptr<Alignment>> ptrs;
    for (const auto &[_, alns] : path_queue_) {
        std::copy(alns.begin(), alns.end(), std::back_inserter(ptrs));
    }
    std::copy(unlabeled_.begin(), unlabeled_.end(), std::back_inserter(ptrs));
    clear();
    // sort by value (not by pointer value)
    std::sort(ptrs.begin(), ptrs.end(), cmp_);
    // transform pointers to objects
    std::vector<Alignment> alignments;
    alignments.reserve(ptrs.size());
    for (auto it = ptrs.rbegin(); it != ptrs.rend(); ++it) {
        // make sure this alignment hasn't been moved yet
        if ((*it)->size()) {
            alignments.emplace_back(std::move(**it));
            **it = Alignment();
        }
    }

    return alignments;
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
