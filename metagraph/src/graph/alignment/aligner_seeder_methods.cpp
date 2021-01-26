#include "aligner_methods.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef DBGAlignerConfig::score_t score_t;

template <typename NodeType>
ExactMapSeeder<NodeType>::ExactMapSeeder(const DeBruijnGraph &graph,
                                         std::string_view query,
                                         bool orientation,
                                         std::vector<NodeType>&& nodes,
                                         const DBGAlignerConfig &config)
      : graph_(graph),
        query_(query),
        orientation_(orientation),
        query_nodes_(std::move(nodes)),
        config_(config),
        num_matching_(0) {
    assert(config_.check_config_scores());

    // count the number of matching nucleotides
    size_t last_match_count = 0;
    for (auto it = query_nodes_.begin(); it != query_nodes_.end(); ++it) {
        if (*it) {
            auto jt = std::find(it + 1, query_nodes_.end(), NodeType());
            num_matching_ += graph_.get_k() + std::distance(it, jt) - 1 - last_match_count;
            last_match_count = graph_.get_k();
            it = jt - 1;
        } else if (last_match_count) {
            --last_match_count;
        }
    }
    assert(num_matching_ <= query_.size());

    partial_sum_.resize(query_.size() + 1);
    std::transform(query_.begin(), query_.end(),
                   partial_sum_.begin() + 1,
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sum_.begin(), partial_sum_.end(), partial_sum_.begin());
    assert(config_.match_score(query_) == partial_sum_.back());
    assert(config_.get_row(query_.front())[query_.front()] == partial_sum_[1]);
    assert(!partial_sum_.front());
}

template <typename NodeType>
void ExactSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const DeBruijnGraph &graph = this->graph_;
    size_t k = graph.get_k();

    const DBGAlignerConfig &config = this->config_;
    assert(k >= config.min_seed_length);

    const std::vector<NodeType> &query_nodes = this->query_nodes_;
    const std::vector<score_t> &partial_sum = this->partial_sum_;
    std::string_view query = this->query_;
    bool orientation = this->orientation_;

    if (this->num_matching_ < config.min_exact_match * query.size())
        return;

    for (size_t i = 0; i < query_nodes.size(); ++i) {
        if (query_nodes[i] != DeBruijnGraph::npos) {
            assert(i + k <= query.size());

            score_t match_score = partial_sum[i + k] - partial_sum[i];

            if (match_score > config.min_cell_score) {
                Seed seed(std::string_view(query.data() + i, k), { query_nodes[i] },
                          match_score, i, orientation);

                assert(seed.is_valid(graph, &config));
                callback(std::move(seed));
            }
        }
    }
}

template <class BaseSeeder>
void SuffixSeeder<BaseSeeder>::call_seeds(std::function<void(Seed&&)> callback) const {
    this->BaseSeeder::call_seeds(callback);

    std::string_view query = this->query_;
    const DBGAlignerConfig &config = this->config_;
    const std::vector<node_index> &query_nodes = this->query_nodes_;

    if (this->num_matching_ < config.min_exact_match * query.size())
        return;

    size_t k = dbg_succ_.get_k();
    bool orientation = this->orientation_;

    auto call_suffix_seed = [&](size_t i, node_index alt_node, size_t seed_length) {
        std::string_view seed_seq(query.data() + i, seed_length);
        DBGAlignerConfig::score_t match_score = config.match_score(seed_seq);

        if (match_score > config.min_cell_score) {
            Seed seed(seed_seq, { alt_node }, match_score, i,
                      orientation, k - seed_length);

            assert(seed.is_valid(dbg_succ_, &config));
            callback(std::move(seed));
        }
    };

    for (size_t i = 0; i + config.min_seed_length <= query.size(); ++i) {
        if (i >= query_nodes.size() || !query_nodes[i]) {
            size_t max_seed_length = std::min({
                config.max_seed_length, k, query.size() - i
            });
            dbg_succ_.call_nodes_with_suffix_matching_longest_prefix(
                std::string_view(query.data() + i, max_seed_length),
                [&](node_index alt_node, size_t seed_length) {
                    call_suffix_seed(i, alt_node, seed_length);
                },
                config.min_seed_length,
                config.max_num_seeds_per_locus
            );
        }
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const DeBruijnGraph &graph = this->graph_;
    size_t k = graph.get_k();

    auto query = this->query_;
    const std::vector<NodeType> &query_nodes = this->query_nodes_;
    const DBGAlignerConfig &config = this->config_;
    bool orientation = this->orientation_;
    const std::vector<score_t> &partial_sum = this->partial_sum_;

    if (this->num_matching_ < config.min_exact_match * query.size())
        return;

    std::vector<uint8_t> query_node_flags(query_nodes.size(), 0);
    for (size_t i = 0; i < query_node_flags.size(); ++i) {
        if (query_nodes[i] != DeBruijnGraph::npos) {
            // the second bit indicates that a node has been found, while the
            // first bit indicates if the node is a maximal exact match terminus
            query_node_flags[i] = 2 | get_mem_terminator()[query_nodes[i]];
        }
    }

    // find start of MEM
    auto it = query_node_flags.begin();
    while ((it = std::find_if(it, query_node_flags.end(),
                              [](uint8_t flags) { return flags & 2; }))
            != query_node_flags.end()) {
        // find end of MEM
        auto next = std::find_if(
            it, query_node_flags.end(),
            [](uint8_t flags) { return (flags & 1) == 1 || (flags & 2) == 0; }
        );

        if (next != query_node_flags.end() && ((*next) & 2))
            ++next;

        assert(next > it);
        assert(next <= query_node_flags.end());

        size_t i = it - query_node_flags.begin();
        assert(it == query_node_flags.end() || query_nodes[i] != DeBruijnGraph::npos);

        size_t mem_length = (next - it) + k - 1;
        assert(i + mem_length <= query.size());

        const char *begin_it = query.data() + i;
        const char *end_it = begin_it + mem_length;

        assert(end_it >= begin_it + config.min_seed_length);

        score_t match_score = partial_sum[end_it - query.data()] - partial_sum[i];

        auto node_begin_it = query_nodes.begin() + i;
        auto node_end_it = node_begin_it + (next - it);
        assert(std::find(node_begin_it, node_end_it, DeBruijnGraph::npos) == node_end_it);

        if (match_score > config.min_cell_score) {
            Seed seed(std::string_view(begin_it, mem_length),
                      { node_begin_it, node_end_it },
                      match_score, i, orientation);
            assert(seed.is_valid(graph, &config));
            callback(std::move(seed));
        }

        it = next;
    }
}


template class ExactMapSeeder<>;
template class ExactSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;
template class SuffixSeeder<ExactSeeder<>>;
template class SuffixSeeder<UniMEMSeeder<>>;

} // namespace align
} // namespace graph
} // namespace mtg
