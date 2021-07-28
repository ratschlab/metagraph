#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"

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
    typedef Alignment::node_index node_index;
    typedef Alignment::score_t score_t;
    typedef PriorityDeque<Alignment, std::vector<Alignment>, AlignmentCompare> PathQueue;

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(Alignment&& alignment);

    score_t get_min_path_score() const;
    score_t get_max_path_score() const;

    score_t get_min_path_score(const Alignment &) const { return get_min_path_score(); }
    score_t get_max_path_score(const Alignment &) const { return get_max_path_score(); }

    const Alignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    std::vector<Alignment> get_alignments();

    size_t size() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    void clear() { path_queue_.clear(); }

    std::string_view get_query(bool is_reverse_complement) const {
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
                                   Alignment&& chain,
                                   typename std::vector<Alignment>::iterator begin,
                                   typename std::vector<Alignment>::iterator end,
                                   std::vector<score_t> &best_score,
                                   const std::function<void(Alignment&&)> &callback);
};


template <class AlignmentCompare>
inline void AlignmentAggregator<AlignmentCompare>::add_alignment(Alignment&& alignment) {
    if (std::find(path_queue_.begin(), path_queue_.end(), alignment) != path_queue_.end())
        return;

    if (config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths) {
        path_queue_.emplace(std::move(alignment));
    } else if (!cmp_(alignment, path_queue_.minimum())) {
        path_queue_.update(path_queue_.begin(), std::move(alignment));
    }
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>::get_min_path_score() const -> score_t {
    return config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths
        ? config_.min_path_score
        : std::max(static_cast<score_t>(path_queue_.maximum().get_score() * config_.rel_score_cutoff),
                   path_queue_.minimum().get_score());
}

template <class AlignmentCompare>
inline auto AlignmentAggregator<AlignmentCompare>::get_max_path_score() const -> score_t {
    return path_queue_.size() ? path_queue_.maximum().get_score() : config_.min_path_score;
}

template <class AlignmentCompare>
inline std::vector<Alignment> AlignmentAggregator<AlignmentCompare>::get_alignments() {
    if (config_.chain_alignments)
        construct_alignment_chains();

    std::vector<Alignment> data(std::move(path_queue_.data()));
    path_queue_.clear();

    // Pop off the min element to the back of the range each time. This results
    // in the vector being in non-increasing order
    for (auto it = data.rbegin(); it != data.rend(); ++it) {
        boost::heap::pop_interval_heap_min(data.begin(), it.base(), path_queue_.cmp());
    }

    return data;
}

template <class AlignmentCompare>
inline void AlignmentAggregator<AlignmentCompare>::construct_alignment_chains() {
    if (path_queue_.empty())
        return;

    std::vector<Alignment> alignments[2];
    for (auto&& alignment : path_queue_.data()) {
        std::vector<Alignment> &bucket = alignments[alignment.get_orientation()];
        bucket.push_back(std::forward<decltype(alignment)>(alignment));
    }

    if (alignments[0].empty() && alignments[1].empty())
        return;

    path_queue_.clear();

    auto push_to_queue = [&](Alignment&& chain) {
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
                construct_alignment_chain(this_query, Alignment(*it),
                                          it + 1, aln.end(),
                                          best_score, push_to_queue);
            }
        }
    }
}

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
