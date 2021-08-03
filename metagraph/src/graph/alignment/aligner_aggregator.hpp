#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <algorithm>

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/algorithms.hpp"
#include "common/vector_map.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename Type, typename Sequence, typename Compare>
class PriorityDeque : public boost::container::priority_deque<Type, Sequence, Compare> {
  public:
    Sequence& data() { return this->sequence(); }
    Compare& cmp() { return this->compare(); }
};


template <class AlignmentCompare>
class AlignmentAggregator {
  public:
    struct SharedPtrCmp {
        bool operator()(const std::shared_ptr<Alignment> &a,
                        const std::shared_ptr<Alignment> &b) const {
            return base_cmp_(*a, *b);
        }

        AlignmentCompare base_cmp_;
    };

    typedef Alignment::node_index node_index;
    typedef Alignment::score_t score_t;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef PriorityDeque<std::shared_ptr<Alignment>,
                          std::vector<std::shared_ptr<Alignment>>,
                          SharedPtrCmp> PathQueue;

    static constexpr Column ncol = std::numeric_limits<Column>::max();

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(Alignment&& alignment);

    score_t get_min_path_score(const Vector<Column> &labels) const;
    score_t get_max_path_score(const Vector<Column> &labels) const;

    score_t get_min_path_score(Column label = ncol) const;
    score_t get_max_path_score(Column label = ncol) const;

    score_t get_min_path_score(const Alignment &seed) const {
        return get_min_path_score(seed.label_columns);
    }

    score_t get_max_path_score(const Alignment &seed) const {
        return get_max_path_score(seed.label_columns);
    }

    const Alignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    std::vector<Alignment> get_alignments();

    size_t size() const {
        size_t size = 0;
        for (const auto &[label, queue] : path_queue_)
            size += queue.size();

        return size;
    }

    size_t num_labels() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

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
    SharedPtrCmp cmp_;

    void construct_alignment_chains();
    void construct_alignment_chain(std::string_view query,
                                   Alignment&& chain,
                                   typename std::vector<Alignment>::iterator begin,
                                   typename std::vector<Alignment>::iterator end,
                                   std::vector<score_t> &best_score,
                                   const std::function<void(Alignment&&)> &callback);
};


template <class AlignmentCompare>
inline void AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    auto packaged_alignment = std::make_shared<Alignment>(std::move(alignment));

    auto add_to_label = [&](Column label) {
        auto &cur_queue = path_queue_[label];

        for (const auto &aln : cur_queue) {
            if (*packaged_alignment == *aln)
                return;
        }

        if (config_.chain_alignments || cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(packaged_alignment);
        } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), packaged_alignment);
        }
    };

    add_to_label(ncol);
    std::for_each(packaged_alignment->label_columns.begin(),
                  packaged_alignment->label_columns.end(),
                  add_to_label);
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>
::get_min_path_score(const Vector<Column> &labels) const -> score_t {
    score_t global_min = !config_.chain_alignments
        ? get_max_path_score() * config_.rel_score_cutoff
        : std::numeric_limits<score_t>::min();

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
    return config_.chain_alignments || find == path_queue_.end()
            || find->second.size() < config_.num_alternative_paths
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
    if (config_.chain_alignments)
        construct_alignment_chains();

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

template <class AlignmentCompare>
inline void AlignmentAggregator<AlignmentCompare>::construct_alignment_chains() {
    if (path_queue_.empty())
        return;

    std::vector<Alignment> alignments[2];
    for (const auto &[label, queue] : path_queue_) {
        for (std::shared_ptr<Alignment> aln : queue) {
            if (aln->size()) {
                auto &bucket = alignments[aln->get_orientation()];
                bucket.emplace_back(std::move(*aln));

                // once an alignment is used from one label queue, clear it so
                // it can't be fetched from another queue
                *aln = Alignment();
            }
        }
    }

    if (alignments[0].empty() && alignments[1].empty())
        return;

    path_queue_.clear();

    auto push_to_queue = [&](Alignment&& chain) {
        auto packaged_alignment = std::make_shared<Alignment>(std::move(chain));

        auto add_to_label = [&](Column label) {
            auto &cur_queue = path_queue_[label];

            for (const auto &aln : cur_queue) {
                if (*packaged_alignment == *aln)
                    return;
            }

            if (cur_queue.size() < config_.num_alternative_paths) {
                cur_queue.emplace(packaged_alignment);
            } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
                cur_queue.update(cur_queue.begin(), packaged_alignment);
            }
        };

        add_to_label(ncol);
        std::for_each(packaged_alignment->label_columns.begin(),
                      packaged_alignment->label_columns.end(),
                      add_to_label);
    };

    for (bool orientation : { false, true }) {
        auto &aln = alignments[orientation];

        // sort by endpoint (using beginning point and scores as tie-breakers)
        std::sort(aln.begin(), aln.end(), [](const auto &a, const auto &b) {
            return std::make_tuple(a.get_clipping() + a.get_query().size(),
                                   a.get_clipping(),
                                   b.get_score(),
                                   a.get_sequence().size())
                < std::make_tuple(b.get_clipping() + b.get_query().size(),
                                  b.get_clipping(),
                                  a.get_score(),
                                  b.get_sequence().size());
        });

        // recursively construct chains
        std::string_view this_query = get_query(orientation);
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = aln.begin(); it != aln.end(); ++it) {
            size_t end_pos = it->get_query().data() + it->get_query().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain(this_query, Alignment(*it),
                                          it + 1, aln.end(),
                                          best_score, push_to_queue);
            }
        }
    }
}

struct ConsecutiveCoords {
    template <class InputIt1, class InputIt2, class OutputIt>
    void operator()(InputIt1 a_begin,
                    InputIt1 a_end,
                    InputIt2 b_begin,
                    InputIt2 b_end,
                    OutputIt out) const {
        while (a_begin + 1 < a_end) {
            *out = *a_begin;
            ++out;
            ++a_begin;
        }

        if (a_begin != a_end) {
            auto a_last = *a_begin;
            ++a_begin;
            *out = std::move(a_last);
            ++out;

            // find the first set of coordinates which is disjoint from the first range
            b_begin = std::lower_bound(b_begin, b_end, a_last,
                [](const auto &a, const auto &b) {
                    return a.first < b.first || a.second < b.second;
                }
            );
        }

        while (b_begin != b_end) {
            *out = *b_begin;
            ++out;
            ++b_begin;
        }
    }
};

template <class AlignmentCompare>
inline void AlignmentAggregator<AlignmentCompare>
::construct_alignment_chain(std::string_view query,
                            Alignment&& chain,
                            typename std::vector<Alignment>::iterator begin,
                            typename std::vector<Alignment>::iterator end,
                            std::vector<score_t> &best_score,
                            const std::function<void(Alignment&&)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query().data();
    const char *chain_end = chain_begin + chain.get_query().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    size_t k = graph_.get_k();
    score_t score = chain.get_score();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        if (it->get_offset())
            continue;

        const char *next_begin = it->get_query().data();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        const char *next_end = next_begin + it->get_query().size();

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        Alignment::LabelSet label_columns;
        Alignment::CoordinateSet label_coordinates;
        if (chain.label_coordinates.size()) {
            // only put together a chain if the labels intersect and if there exist
            // disjoint coordinates in the incoming alignment
            assert(chain.label_columns.size() == chain.label_coordinates.size());
            utils::indexed_set_op<std::vector<std::pair<uint64_t, uint64_t>>,
                                  ConsecutiveCoords>(
                chain.label_columns.begin(), chain.label_columns.end(),
                chain.label_coordinates.begin(),
                it->label_columns.begin(), it->label_columns.end(),
                it->label_coordinates.begin(),
                std::back_inserter(label_columns),
                std::back_inserter(label_coordinates)
            );

            if (label_columns.empty())
                continue;

        } else if (chain.label_columns.size()) {
            std::set_intersection(it->label_columns.begin(), it->label_columns.end(),
                                  chain.label_columns.begin(),
                                  chain.label_columns.end(),
                                  std::back_inserter(label_columns));

            if (label_columns.empty())
                continue;
        }

        Alignment aln(*it);

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, graph_, config_);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().data().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query().data(), graph_, config_)
            );

            if (aln.empty() || aln.get_sequence().size() < graph_.get_k()
                    || aln.get_cigar().data().begin()->first != Cigar::MATCH)
                continue;

            assert(aln.get_query().data() == chain.get_query().data() + chain.get_query().size());

            if (overlap < k - 1)
                aln.insert_gap_prefix(-overlap, graph_, config_);
        }

        score_t next_score = score + aln.get_score();
        if (next_score > best_score[next_end - query.data()]) {
            best_score[next_end - query.data()] = next_score;

            Alignment next_chain(chain);
            next_chain.trim_end_clipping();
            next_chain.append(std::move(aln));
            assert(next_chain.get_score() == next_score);
            assert(next_chain.is_valid(graph_, &config_));
            if (next_chain.size()) {
                if (label_columns.size() == chain.label_columns.size())
                    called = true;

                next_chain.label_columns = std::move(label_columns);
                next_chain.label_coordinates = std::move(label_coordinates);
                construct_alignment_chain(query, std::move(next_chain), it + 1,
                                          end, best_score, callback);
            }
        }
    }

    if (!called)
        callback(std::move(chain));
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
