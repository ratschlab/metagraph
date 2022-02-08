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

    score_t get_min_path_score(const Vector<Column> &labels) const;
    score_t get_max_path_score(const Vector<Column> &labels) const;

    score_t get_min_path_score(Column label = ncol) const;
    score_t get_max_path_score(Column label = ncol) const;

    std::vector<Alignment> get_alignments();

    const VectorMap<Column, PathQueue>& data() const { return path_queue_; }
    VectorMap<Column, PathQueue>& data() { return path_queue_; }

    void clear() { path_queue_.clear(); }

  private:
    const DBGAlignerConfig &config_;
    VectorMap<Column, PathQueue> path_queue_;
    ValCmp cmp_;
};


template <class AlignmentCompare>
inline bool AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    auto a = std::make_shared<Alignment>(std::move(alignment));

    if (path_queue_.empty()) {
        path_queue_[ncol].emplace(a);
        for (Column c : a->label_columns) {
            path_queue_[c].emplace(a);
        }
        return true;
    }

    if (a->label_columns.empty()) {
        // If other alignments have been added to the aggregator which have labels,
        // then we're not interested in storing an unlabeled alignment
        if (path_queue_.size() != 1)
            return false;

        auto &nqueue = path_queue_[ncol];
        if (nqueue.size() < config_.num_alternative_paths
                || config_.post_chain_alignments) {
            for (const auto &aln : nqueue) {
                if (*a == *aln)
                    return false;
            }

            nqueue.emplace(a);
            return true;
        }

        if (cmp_(a, nqueue.minimum()))
            return false;

        nqueue.update(nqueue.begin(), a);
        return true;
    }

    // note: Each reference made to path_queue_ has the potential to increase
    //       the size of the underlying storage, and hence, reallocate it.
    //       So, storing a reference to path_queue_[ncol] may lead to segfaults.

    for (const auto &aln : path_queue_[ncol]) {
        if (*a == *aln)
            return false;
    }

    if (path_queue_[ncol].size() < config_.num_alternative_paths
            || config_.post_chain_alignments) {
        for (Column c : a->label_columns) {
            auto &cur_queue = path_queue_[c];
            for (const auto &aln : cur_queue) {
                if (*a == *aln)
                    return false;
            }

            cur_queue.emplace(a);
        }

        path_queue_[ncol].emplace(a);
        return true;
    }

    if (cmp_(a, path_queue_[ncol].minimum())
            && a->get_score()
                < path_queue_[ncol].maximum()->get_score() * config_.rel_score_cutoff) {
        return false;
    }

    bool added = false;
    for (Column c : a->label_columns) {
        auto &cur_queue = path_queue_[c];
        if (cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(a);
            added = true;
        } else if (!cmp_(a, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), a);
            added = true;
        }
    }
    if (!added)
        return false;

    if (!cmp_(a, path_queue_[ncol].minimum())) {
        path_queue_[ncol].update(path_queue_[ncol].begin(), a);
    } else {
        path_queue_[ncol].emplace(a);
    }
    return true;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_min_path_score(const Vector<Column> &labels) const -> score_t {
    score_t global_min = get_max_path_score() * config_.rel_score_cutoff;

    if (labels.empty())
        return std::max(global_min, get_min_path_score());

    score_t min_score = std::numeric_limits<score_t>::max();
    for (Column label : labels) {
        if (min_score < global_min)
            break;

        min_score = std::min(min_score, get_min_path_score(label));
    }

    return std::max(global_min, min_score);
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_max_path_score(const Vector<Column> &labels) const -> score_t {
    if (labels.empty())
        return get_max_path_score();

    score_t max_score = std::numeric_limits<score_t>::min();
    for (Column label : labels) {
        max_score = std::max(max_score, get_max_path_score(label));
    }

    return max_score;
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_min_path_score(Column label) const -> score_t {
    auto find = path_queue_.find(label);
    return find == path_queue_.end()
            || find->second.size() < config_.num_alternative_paths
            || config_.post_chain_alignments
        ? config_.min_path_score
        : std::max(static_cast<score_t>(find->second.maximum()->get_score()
                                            * config_.rel_score_cutoff),
                   find->second.minimum()->get_score());
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_max_path_score(Column label) const -> score_t {
    auto find = path_queue_.find(label);
    return find == path_queue_.end() ? config_.min_path_score
                                     : find->second.maximum()->get_score();
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
