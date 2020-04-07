#include "aligner_methods.hpp"

#include "common/algorithms.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


template <typename NodeType>
void ExactSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    query_ = query;
    orientation_ = orientation;
    query_nodes_.clear();
    query_nodes_.reserve(query_.size() - graph_.get_k() + 1);
    graph_.map_to_nodes_sequentially(query_, [&](auto i) { query_nodes_.push_back(i); });
    assert(query_nodes_.size() == query_.size() - graph_.get_k() + 1);

    partial_sum.resize(query_.size() + 1);
    std::transform(query_.begin(), query_.end(),
                   partial_sum.begin() + 1,
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sum.begin(), partial_sum.end(), partial_sum.begin());
    assert(config_.match_score(query_) == partial_sum.back());
    assert(config_.get_row(query_.front())[query_.front()] == partial_sum[1]);
    assert(!partial_sum.front());
}

template <typename NodeType>
void ExactSeeder<NodeType>
::call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const {
    size_t k = graph_.get_k();
    for (size_t i = 0; i < query_nodes_.size(); ++i) {
        assert(query_.begin() + i + k <= query_.end());
        if (query_nodes_[i] == DeBruijnGraph::npos)
            continue;

        score_t match_score = partial_sum[i + k] - partial_sum[i];

        if (match_score <= config_.min_cell_score)
            continue;

        callback(Alignment<NodeType>(std::string_view(query_.data() + i, k),
                                     { query_nodes_[i] },
                                     match_score,
                                     i,
                                     orientation_));

        // if this k-mer was matched, then the next one will be accessible
        // by an extend call, so there's no need to seed from it
        ++i;
    }
}

template <typename NodeType>
SuffixSeeder<NodeType>::SuffixSeeder(const DeBruijnGraph &graph,
                                     const DBGAlignerConfig &config)
      : exact_seeder_(graph, config) {
    if (!dynamic_cast<const DBGSuccinct*>(&graph))
        throw std::runtime_error("Only implemented for DBGSuccinct");
}

template <typename NodeType>
void SuffixSeeder<NodeType>
::call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const {
    // TODO: handle cases where
    // 1) the query is less than k
    const auto &graph_ = dynamic_cast<const DBGSuccinct&>(get_graph());
    const auto &config_ = exact_seeder_.get_config();
    auto query_ = exact_seeder_.get_query();
    const auto &query_nodes_ = exact_seeder_.get_query_nodes();
    const auto &partial_sum = exact_seeder_.get_partial_sums();
    bool orientation_ = exact_seeder_.get_orientation();

    if (config_.min_seed_length > std::min(graph_.get_k(), query_.size()))
        return;

    size_t max_seed_length = std::min(
        query_.size(),
        std::min(config_.max_seed_length, graph_.get_k())
    );

    for (size_t i = 0; i < query_nodes_.size(); ++i) {
        assert(query_.begin() + i + graph_.get_k() <= query_.end());
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            score_t match_score = partial_sum[i + graph_.get_k()] - partial_sum[i];

            if (match_score <= config_.min_cell_score)
                continue;

            callback(Alignment<NodeType>(
                std::string_view(query_.data() + i, graph_.get_k()),
                { query_nodes_[i] },
                match_score,
                i,
                orientation_)
            );

            // if this k-mer was matched, then the next one will be accessible
            // by an extend call, so there's no need to seed from it
            ++i;
        } else {
            graph_.call_nodes_with_suffix(
                query_.substr(i, max_seed_length),
                [&](NodeType node, size_t seed_length) {
                    assert(node != DeBruijnGraph::npos);

                    score_t match_score = partial_sum[i + seed_length] - partial_sum[i];

                    if (match_score > config_.min_cell_score) {
                        callback(Alignment<NodeType>(
                            std::string_view(query_.data() + i, seed_length),
                            std::vector<NodeType>{ node },
                            match_score,
                            i,
                            orientation_,
                            graph_.get_k() - seed_length)
                        );
                    }
                },
                config_.min_seed_length,
                config_.max_num_seeds_per_locus
            );
        }
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    query_ = query;
    orientation_ = orientation;
    query_nodes_.clear();
    query_nodes_.reserve(query.size() - graph_.get_k() + 1);
    graph_.map_to_nodes_sequentially(query, [&](auto i) { query_nodes_.push_back(i); });

    partial_sum.resize(query_.size() + 1);
    std::transform(query_.begin(), query_.end(),
                   partial_sum.begin() + 1,
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sum.begin(), partial_sum.end(), partial_sum.begin());
    assert(config_.match_score(query_) == partial_sum.back());
    assert(config_.get_row(query_.front())[query_.front()] == partial_sum[1]);
    assert(!partial_sum.front());
}

// TODO: make this work with unmasked DBGSuccinct
template <typename NodeType>
void MEMSeeder<NodeType>
::call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const {
    size_t k = graph_.get_k();

    // find start of MEM
    auto it = query_nodes_.begin();

    while ((it = std::find_if(it, query_nodes_.end(),
                              [](NodeType node) { return node != DeBruijnGraph::npos; }))
            < query_nodes_.end()) {
        // find end of MEM
        auto next = std::find_if(
            it, query_nodes_.end(),
            [&](NodeType node) { return node == DeBruijnGraph::npos || (*is_mem_terminus_)[node]; }
        );

        if (next != query_nodes_.end() && *next != DeBruijnGraph::npos)
            ++next;

        assert(next > it);
        assert(next <= query_nodes_.end());

        // compute the correct string offsets
        const char *begin_it = query_.data() + (it - query_nodes_.begin());
        const char *end_it = (next == query_nodes_.end())
            ? query_.data() + query_.size()
            : query_.data() + (next - query_nodes_.begin()) + k - 1;

        score_t match_score = partial_sum[end_it - query_.data()]
            - partial_sum[begin_it - query_.data()];

        if (match_score > config_.min_cell_score) {
            callback(Alignment<NodeType>(std::string_view(begin_it, end_it - begin_it),
                                         std::vector<NodeType>(it, next),
                                         match_score,
                                         begin_it - query_.data(),
                                         orientation_));
        }

        it = next;
    }
}


template class ExactSeeder<>;
template class SuffixSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;
