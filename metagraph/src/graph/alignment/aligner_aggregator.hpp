#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
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

    static constexpr Column ncol = std::numeric_limits<Column>::max();

    explicit AlignmentAggregator(const DBGAlignerConfig &config) : config_(config) {
        assert(config_.num_alternative_paths);
    }

    bool add_alignment(Alignment&& alignment);

    score_t get_min_path_score(const Vector<Column> &labels = {}) const;
    score_t get_max_path_score(const Vector<Column> &labels = {}) const;

    std::vector<Alignment> get_alignments();

    const VectorMap<Column, PathQueue>& data() const { return path_queue_; }
    VectorMap<Column, PathQueue>& data() { return path_queue_; }

    void clear() { path_queue_.clear(); }

    size_t num_labels() const {
        // don't count the ncol queue
        return path_queue_.size() ? path_queue_.size() - 1 : 0;
    }

  private:
    const DBGAlignerConfig &config_;
    VectorMap<Column, PathQueue> path_queue_;
    ValCmp cmp_;

    score_t get_label_min_path_score(Column label = ncol) const;
    score_t get_label_max_path_score(Column label = ncol) const;

    score_t get_global_min() const {
        auto find = path_queue_.find(ncol);
        if (find == path_queue_.end()) {
            assert(path_queue_.empty());
            return config_.ninf;
        }

        score_t global_min = find->second.maximum()->get_score();
        if (global_min > 0)
            global_min *= config_.rel_score_cutoff;

        return global_min;
    }
};


template <class AlignmentCompare>
inline bool AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    // first, wrap the alignment so that duplicates are not stored in each per-label queue
    auto a = std::make_shared<Alignment>(std::move(alignment));

    // if nothing has been added to the queue so far, add the alignment
    if (path_queue_.empty()) {
        path_queue_[ncol].emplace(a);
        for (Column c : a->label_columns) {
            path_queue_[c].emplace(a);
        }
        return true;
    }

    // if the score is less than the global relative score cutoff, don't add it
    if (a->get_score() < get_global_min())
        return false;

    // check for duplicates
    auto has_aln = [&](auto find) {
        if (find == path_queue_.end())
            return false;

        for (const auto &aln : find->second) {
            if (*a == *aln)
                return true;
        }

        return false;
    };

    auto find = path_queue_.find(ncol);
    if (has_aln(find))
        return false;

    for (Column c : a->label_columns) {
        if (has_aln(path_queue_.find(c)))
            return false;
    }

    // helper for adding alignments to the queue
    auto add_alignment_to_label = [&](auto &queue) {
        // If the queue is not at capacity, or if post-alignment chaining is
        // requested, add the alignment. Otherwise, replace the minimum element
        // in the queue if the current alignment is better.
        if (queue.size() < config_.num_alternative_paths
                || config_.post_chain_alignments) {
            queue.emplace(a);
            return true;
        } else if (cmp_(a, queue.minimum())) {
            return false;
        }

        queue.update(queue.begin(), a);
        return true;
    };

    // if we are in the unlabeled case, only consider the ncol queue
    if (a->label_columns.empty() && path_queue_.size() == 1)
        return add_alignment_to_label(find.value());

    // if an incoming alignment has labels, and we haven't encountered a labeled
    // alignment yet, we only need the ncol queue for fetching the global minimum,
    // so shrink it to only one element
    if (path_queue_.size() == 1) {
        auto &data = find.value().data();
        if (data.size() > 1)
            data.erase(data.begin(), data.begin() + 1);

        data.erase(data.begin() + 1, data.end());
    }
    assert(find.value().size() == 1);

    // add the alignment to its labeled queues
    bool added = false;
    for (Column c : a->label_columns) {
        added |= add_alignment_to_label(path_queue_[c]);
    }

    if (!added)
        return false;

    // if this is the best alignment so far, update the ncol queue
    auto &nqueue = path_queue_[ncol];
    if (!cmp_(a, nqueue.maximum()))
        nqueue.update(nqueue.begin(), a);

    return true;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_min_path_score(const Vector<Column> &labels) const -> score_t {
    score_t global_min = get_global_min();

    if (labels.empty())
        return std::max(global_min, get_label_min_path_score());

    score_t min_score = std::numeric_limits<score_t>::max();
    for (Column label : labels) {
        if (min_score < global_min)
            break;

        min_score = std::min(min_score, get_label_min_path_score(label));
    }

    return std::max(global_min, min_score);
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_max_path_score(const Vector<Column> &labels) const -> score_t {
    if (labels.empty())
        return get_label_max_path_score();

    score_t max_score = std::numeric_limits<score_t>::min();
    for (Column label : labels) {
        max_score = std::max(max_score, get_label_max_path_score(label));
    }

    return max_score;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_label_min_path_score(Column label) const -> score_t {
    auto find = path_queue_.find(label);
    return find == path_queue_.end()
            || find->second.size() < config_.num_alternative_paths
            || config_.post_chain_alignments
        ? config_.ninf
        : find->second.minimum()->get_score();
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_label_max_path_score(Column label) const -> score_t {
    auto find = path_queue_.find(label);
    return find == path_queue_.end() ? config_.ninf : find->second.maximum()->get_score();
}

template <class AlignmentCompare>
inline std::vector<Alignment> AlignmentAggregator<AlignmentCompare>::get_alignments() {
    // move all alignments to one vector
    std::vector<std::shared_ptr<Alignment>> ptrs;
    for (const auto &[_, alns] : path_queue_) {
        std::copy(alns.begin(), alns.end(), std::back_inserter(ptrs));
    }
    path_queue_.clear();
    // sort and remove duplicates
    std::sort(ptrs.begin(), ptrs.end(), cmp_);
    ptrs.erase(std::unique(ptrs.begin(), ptrs.end()), ptrs.end());
    // transform pointers to objects
    std::vector<Alignment> alignments;
    alignments.reserve(ptrs.size());
    std::transform(ptrs.rbegin(), ptrs.rend(), std::back_inserter(alignments),
                   [](auto &a) { return std::move(*a); });
    return alignments;
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
