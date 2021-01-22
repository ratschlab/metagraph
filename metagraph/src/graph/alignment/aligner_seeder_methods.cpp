#include "aligner_methods.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename NodeType>
void ExactMapSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    query_ = query;
    orientation_ = orientation;
    query_nodes_.clear();

    if (query_.size() < this->get_graph().get_k())
        return;

    query_nodes_ = map_sequence_to_nodes(this->get_graph(), query_);

    partial_sum_.resize(query_.size() + 1);
    std::transform(query_.begin(), query_.end(),
                   partial_sum_.begin() + 1,
                   [&](char c) { return this->get_config().get_row(c)[c]; });

    std::partial_sum(partial_sum_.begin(), partial_sum_.end(), partial_sum_.begin());
    assert(this->get_config().match_score(query_) == partial_sum_.back());
    assert(this->get_config().get_row(query_.front())[query_.front()] == partial_sum_[1]);
    assert(!partial_sum_.front());
}

template <typename NodeType>
void ExactSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const auto &graph = this->get_graph();
    size_t k = graph.get_k();

    const auto &config = this->get_config();
    assert(k >= config.min_seed_length);

    const auto &query_nodes = this->get_query_nodes();
    const auto &partial_sum = this->get_partial_sums();
    auto query = this->get_query();
    bool orientation = this->get_orientation();

    // if node suffixes of length < k were matched, then query_nodes.size() may
    // be greater than query.size() - k + 1
    if (this->get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

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

template <typename NodeType>
SuffixSeeder<NodeType>::SuffixSeeder(const DeBruijnGraph &graph,
                                     const DBGAlignerConfig &config) {
    if (!dynamic_cast<const DBGSuccinct*>(&graph))
        throw std::runtime_error("Only implemented for DBGSuccinct");

    if (config.max_seed_length > graph.get_k()) {
        base_seeder_ = std::make_unique<UniMEMSeeder<NodeType>>(graph, config);
    } else {
        base_seeder_ = std::make_unique<ExactSeeder<NodeType>>(graph, config);
    }
}

template <typename NodeType>
void SuffixSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    base_seeder_->call_seeds(callback);

    const auto query = get_query();
    const auto &config = get_config();
    const auto &query_nodes = get_query_nodes();
    const auto &graph = dynamic_cast<const DBGSuccinct&>(get_graph());

    size_t k = graph.get_k();
    if (this->get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

    bool orientation = this->get_orientation();

    for (size_t i = 0; i + config.min_seed_length <= query.size(); ++i) {
        if (i >= query_nodes.size() || !query_nodes[i]) {
            size_t max_seed_length = std::min(config.max_seed_length,
                std::min(k, query.size() - i)
            );
            graph.call_nodes_with_suffix_matching_longest_prefix(
                std::string_view(query.data() + i, max_seed_length),
                [&](NodeType alt_node, size_t seed_length) {
                    std::string_view seed_seq(query.data() + i, seed_length);
                    DBGAlignerConfig::score_t match_score = config.match_score(seed_seq);

                    if (match_score > config.min_cell_score) {
                        Seed seed(seed_seq, { alt_node }, match_score, i,
                                  orientation, k - seed_length);

                        assert(seed.is_valid(graph, &config));
                        callback(std::move(seed));
                    }
                },
                config.min_seed_length,
                config.max_num_seeds_per_locus
            );
        }
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const auto &graph = this->get_graph();
    size_t k = graph.get_k();

    auto query = this->get_query();
    const auto &query_nodes = this->get_query_nodes();
    const auto &config = this->get_config();
    bool orientation = this->get_orientation();
    const auto &partial_sum = this->get_partial_sums();

    if (this->get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

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
            < query_node_flags.end()) {
        // find end of MEM
        auto next = std::find_if(it, query_node_flags.end(),
                                 [](uint8_t flags) { return (flags & 1) == 1 || (flags & 2) == 0; });

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

        score_t match_score = partial_sum[end_it - query.data()]
            - partial_sum[begin_it - query.data()];

        auto node_begin_it = query_nodes.begin() + (it - query_node_flags.begin());
        auto node_end_it = node_begin_it + (next - it);
        assert(std::find(node_begin_it, node_end_it, DeBruijnGraph::npos) == node_end_it);

        if (match_score > config.min_cell_score) {
            Seed seed(std::string_view(begin_it, end_it - begin_it),
                      std::vector<NodeType>(node_begin_it, node_end_it),
                      match_score,
                      begin_it - query.data(),
                      orientation);
            assert(seed.is_valid(graph, &config));
            callback(std::move(seed));
        }

        it = next;
    }
}


template class ExactSeeder<>;
template class SuffixSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;

} // namespace align
} // namespace graph
} // namespace mtg
