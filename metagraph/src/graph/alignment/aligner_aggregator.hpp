#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/vector_map.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType, class AlignmentCompare>
class AlignmentAggregator {
    template <typename Type, typename Sequence, typename Compare>
    class PriorityDeque : public boost::container::priority_deque<Type, Sequence, Compare> {
      public:
        void sort_interval_heap() {
            boost::heap::sort_interval_heap(this->sequence().begin(), this->sequence().end(),
                                            this->compare());
        }

        Sequence& data() { return this->sequence(); }
    };


  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;
    typedef PriorityDeque<DBGAlignment, std::vector<DBGAlignment>, AlignmentCompare> PathQueue;

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(DBGAlignment&& alignment);

    score_t get_min_path_score() const;
    score_t get_max_path_score() const;

    score_t get_min_path_score(const DBGAlignment &) const {
        return get_min_path_score();
    }

    score_t get_max_path_score(const DBGAlignment &) const {
        return get_max_path_score();
    }


    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                         const std::function<bool()> &terminate);

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback);

    size_t size() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    void clear() { path_queue_.clear(); }

    std::string_view get_query(bool is_reverse_complement) {
        return is_reverse_complement ? rc_query_ : query_;
    }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    const DeBruijnGraph &graph_;
    PathQueue path_queue_;
    AlignmentCompare cmp_;

    void construct_alignment_chains();
    void construct_alignment_chain(std::string_view query,
                                   DBGAlignment&& chain,
                                   typename std::vector<DBGAlignment>::iterator begin,
                                   typename std::vector<DBGAlignment>::iterator end,
                                   std::vector<score_t> &best_score,
                                   const std::function<void(DBGAlignment&&)> &callback);
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    if (std::find(path_queue_.begin(), path_queue_.end(), alignment) != path_queue_.end())
        return;

    if (config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths) {
        path_queue_.emplace(std::move(alignment));
    } else if (!cmp_(alignment, path_queue_.minimum())) {
        path_queue_.update(path_queue_.begin(), std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score() const -> score_t {
    return config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths
        ? config_.min_path_score
        : std::max(static_cast<score_t>(path_queue_.maximum().get_score() * config_.fraction_of_top),
                   path_queue_.minimum().get_score());
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_max_path_score() const -> score_t {
    return path_queue_.size() ? path_queue_.maximum().get_score() : config_.min_path_score;
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                  const std::function<bool()> &terminate) {
    if (config_.chain_alignments)
        construct_alignment_chains();

    while (!terminate() && path_queue_.size()) {
        callback(DBGAlignment(path_queue_.maximum()));
        path_queue_.pop_maximum();
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback) {
    if (config_.chain_alignments)
        construct_alignment_chains();

    path_queue_.sort_interval_heap();
    std::vector<DBGAlignment> paths(std::move(path_queue_.data()));
    std::for_each(std::make_move_iterator(paths.rbegin()),
                  std::make_move_iterator(paths.rend()),
                  callback);

    path_queue_.clear();
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::construct_alignment_chains() {
    if (path_queue_.empty())
        return;

    std::vector<DBGAlignment> alignments[2];
    for (const auto &alignment : path_queue_) {
        alignments[alignment.get_orientation()].push_back(alignment);
    }

    if (alignments[0].empty() && alignments[1].empty())
        return;

    path_queue_.clear();

    auto push_to_queue = [&](DBGAlignment&& chain) {
        if (std::find(path_queue_.begin(), path_queue_.end(), chain) != path_queue_.end())
            return;

        if (path_queue_.size() < config_.num_alternative_paths) {
            path_queue_.emplace(std::move(chain));
        } else if (!cmp_(chain, path_queue_.minimum())) {
            path_queue_.update(path_queue_.begin(), std::move(chain));
        }
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
                construct_alignment_chain(this_query, DBGAlignment(*it),
                                          it + 1, aln.end(),
                                          best_score, push_to_queue);
            }
        }
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::construct_alignment_chain(std::string_view query,
                            DBGAlignment&& chain,
                            typename std::vector<DBGAlignment>::iterator begin,
                            typename std::vector<DBGAlignment>::iterator end,
                            std::vector<score_t> &best_score,
                            const std::function<void(DBGAlignment&&)> &callback) {
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

        DBGAlignment aln(*it);

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, graph_, config_);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query().data(), graph_, config_)
            );

            if (aln.empty() || aln.get_sequence().size() < graph_.get_k()
                    || aln.get_cigar().front().first != Cigar::MATCH)
                continue;

            assert(aln.get_query().data() == chain.get_query().data() + chain.get_query().size());

            if (overlap < k - 1)
                aln.insert_gap_prefix(-overlap, graph_, config_);
        }

        score_t next_score = score + aln.get_score();
        if (next_score > best_score[next_end - query.data()]) {
            best_score[next_end - query.data()] = next_score;

            DBGAlignment next_chain(chain);
            next_chain.trim_end_clipping();
            next_chain.append(std::move(aln));
            assert(next_chain.get_score() == next_score);
            assert(next_chain.is_valid(graph_, &config_));
            if (next_chain.size()) {
                called = true;
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
