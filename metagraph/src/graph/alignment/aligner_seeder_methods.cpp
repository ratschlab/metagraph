#include "aligner_methods.hpp"

#include "common/algorithms.hpp"
#include "dbg_succinct.hpp"


template <typename NodeType>
std::vector<Alignment<NodeType>> ExactSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    assert(seed_end >= seed_begin);

    const size_t length = seed_end - seed_begin;

    if (length < graph_.get_k()
            || config_.min_seed_length > graph_.get_k()
            || config_.max_seed_length < graph_.get_k())
        return {};

    auto exact = graph_.kmer_to_node(std::string(seed_begin, seed_begin + graph_.get_k()));

    if (exact == DeBruijnGraph::npos)
        return {};

    auto match_score = config_.match_score(seed_begin, seed_begin + graph_.get_k());

    if (match_score <= config_.min_cell_score)
        return {};

    return { Alignment<NodeType>(seed_begin,
                                 seed_begin + graph_.get_k(),
                                 { exact },
                                 match_score,
                                 clipping,
                                 orientation) };
}

template <typename NodeType>
std::vector<Alignment<NodeType>> SuffixSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    if (!dynamic_cast<const DBGSuccinct*>(&get_graph()))
        throw std::runtime_error("Only implemented for DBGSuccinct");

    const auto &graph = dynamic_cast<const DBGSuccinct&>(get_graph());
    const auto &config = get_config();

    assert(seed_end >= seed_begin);
    const size_t length = seed_end - seed_begin;

    if (config.min_seed_length > std::min(graph.get_k(), length))
        return {};

    auto seeds = exact_seeder_(seed_begin, seed_end, clipping, orientation);

    if (seeds.size())
        return seeds;

    auto max_seed_length = std::min(
        length, std::min(config.max_seed_length, graph.get_k())
    );

    graph.call_nodes_with_suffix(
        seed_begin,
        seed_begin + max_seed_length,
        [&](auto node, uint64_t seed_length) {
            assert(node != DeBruijnGraph::npos);

            auto match_score = config.match_score(seed_begin, seed_begin + seed_length);

            if (match_score > get_config().min_cell_score) {
                seeds.emplace_back(seed_begin,
                                   seed_begin + seed_length,
                                   std::vector<NodeType>{ node },
                                   match_score,
                                   clipping,
                                   orientation,
                                   graph.get_k() - seed_length);
            }
        },
        config.min_seed_length,
        config.max_num_seeds_per_locus
    );

    assert(seeds.size() <= config.max_num_seeds_per_locus);

    return seeds;
}

// TODO: make this work with unmasked DBGSuccinct
template <typename NodeType>
std::vector<Alignment<NodeType>> MEMSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    if (query_nodes_.empty())
        throw std::runtime_error("MEMSeeder uninitialized");

    if (orientation != orientation_)
        throw std::runtime_error("wrong orientation passed");

    assert(seed_end >= seed_begin);
    const size_t length = seed_end - seed_begin;

    if (!length)
        return {};

    assert(length + clipping == query_nodes_.size() + graph_.get_k() - 1);

    auto start = std::find_if(
        query_nodes_.begin() + clipping,
        query_nodes_.end(),
        [](auto node) { return node != DeBruijnGraph::npos; }
    );

    if (start != query_nodes_.begin() + clipping)
        return {};

    assert(*start != DeBruijnGraph::npos);
    auto next = std::find_if(
        start,
        query_nodes_.end(),
        [&](auto i) { return i == DeBruijnGraph::npos || (*is_mem_terminus_)[i]; }
    );

    if (UNLIKELY(next != query_nodes_.end() && *next != DeBruijnGraph::npos))
        next++;

    assert(next != start);

    auto end_it = next == query_nodes_.end()
        ? seed_end
        : seed_begin + (next - query_nodes_.begin() - clipping) + graph_.get_k() - 1;

    assert(end_it >= seed_begin);
    assert(end_it <= seed_end);
    assert(static_cast<size_t>(end_it - seed_begin) >= graph_.get_k());

    auto match_score = config_.match_score(seed_begin, end_it);

    if (match_score <= config_.min_cell_score)
        return {};

    return { Alignment<NodeType>(seed_begin,
                                 end_it,
                                 std::vector<NodeType>(start, next),
                                 config_.match_score(seed_begin, end_it),
                                 clipping,
                                 orientation) };
}


template class ExactSeeder<>;
template class SuffixSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;
