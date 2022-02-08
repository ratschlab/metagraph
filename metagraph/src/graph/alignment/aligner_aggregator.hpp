#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/algorithms.hpp"
#include "common/vector_map.hpp"
#include "common/utils/template_utils.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

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
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef PriorityDeque<std::shared_ptr<Alignment>,
                          std::vector<std::shared_ptr<Alignment>>, ValCmp> PathQueue;

    static constexpr Column ncol = std::numeric_limits<Column>::max();

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
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

    std::string_view get_query(bool is_reverse_complement) const {
        return is_reverse_complement ? rc_query_ : query_;
    }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    const DeBruijnGraph &graph_;
    VectorMap<Column, PathQueue> path_queue_;
    ValCmp cmp_;
};


template <class AlignmentCompare>
inline bool AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    auto packaged_alignment = std::make_shared<Alignment>(std::move(alignment));

    if (path_queue_.empty()) {
        path_queue_[ncol].emplace(packaged_alignment);
        for (Column c : packaged_alignment->label_columns) {
            path_queue_[c].emplace(packaged_alignment);
        }
        return true;
    }

    if (packaged_alignment->label_columns.empty()) {
        if (path_queue_.size() != 1)
            return false;

        auto &queue = path_queue_[ncol];
        if (queue.size() < config_.num_alternative_paths
                || config_.post_chain_alignments) {
            for (const auto &aln : queue) {
                if (*packaged_alignment == *aln)
                    return false;
            }

            queue.emplace(packaged_alignment);
            return true;
        }

        if (cmp_(packaged_alignment, queue.minimum()))
            return false;

        queue.update(queue.begin(), packaged_alignment);
        return true;
    }

    for (const auto &aln : path_queue_[ncol]) {
        if (*packaged_alignment == *aln)
            return false;
    }

    if (path_queue_[ncol].size() < config_.num_alternative_paths
            || config_.post_chain_alignments) {
        for (Column c : packaged_alignment->label_columns) {
            auto &cur_queue = path_queue_[c];
            for (const auto &aln : cur_queue) {
                if (*packaged_alignment == *aln)
                    return false;
            }

            cur_queue.emplace(packaged_alignment);
        }

        path_queue_[ncol].emplace(packaged_alignment);
        return true;
    }

    if (cmp_(packaged_alignment, path_queue_[ncol].minimum())
            && (packaged_alignment->get_score()
                < path_queue_[ncol].maximum()->get_score() * config_.rel_score_cutoff)) {
        return false;
    }

    bool added = false;
    for (Column c : packaged_alignment->label_columns) {
        auto &cur_queue = path_queue_[c];
        if (cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(packaged_alignment);
            added = true;
        } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), packaged_alignment);
            added = true;
        }
    }
    if (!added)
        return false;

    auto &cur_queue = path_queue_[ncol];
    if (!cmp_(packaged_alignment, cur_queue.minimum())) {
        cur_queue.update(cur_queue.begin(), packaged_alignment);
    } else {
        cur_queue.emplace(packaged_alignment);
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
    return find == path_queue_.end() || find->second.size() < config_.num_alternative_paths
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
    typedef std::pair<Column, PathQueue> queue_value;
    auto queues = const_cast<std::vector<queue_value>&&>(path_queue_.values_container());
    path_queue_.clear();

    std::vector<std::shared_ptr<Alignment>> alignment_ptrs;
    for (auto &[label, queue] : queues) {
        std::vector<std::shared_ptr<Alignment>> merged;
        boost::heap::sort_interval_heap(queue.data().begin(), queue.data().end(), cmp_);
        std::set_union(alignment_ptrs.begin(), alignment_ptrs.end(),
                       queue.data().begin(), queue.data().end(),
                       std::back_inserter(merged), cmp_);
        std::swap(merged, alignment_ptrs);
    }

    std::vector<Alignment> alignments;
    alignments.reserve(alignment_ptrs.size());

    std::transform(alignment_ptrs.rbegin(), alignment_ptrs.rend(),
                   std::back_inserter(alignments),
                   [](auto a) { return std::move(*a); });

    return alignments;
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
