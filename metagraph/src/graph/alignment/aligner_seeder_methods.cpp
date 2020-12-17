#include "aligner_methods.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename NodeType>
void ExactMapSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    query_ = query;
    orientation_ = orientation;
    query_nodes_.clear();
    num_matching_kmers_ = 0;

    if (query_.size() < config_.min_seed_length)
        return;

    query_nodes_ = map_sequence_to_nodes(graph_, query_);
    num_matching_kmers_ = query_nodes_.size()
        - std::count(query_nodes_.begin(), query_nodes_.end(), NodeType());

    offsets_.resize(query_nodes_.size());

    for (size_t i = 0; i < query_nodes_.size(); ++i) {
        if (!query_nodes_[i]) {
            offsets_[i] = 0;
        } else {
            size_t node_length = graph_.get_node_length(query_nodes_[i]);
            assert(i + node_length <= query_.size());

            size_t offset = graph_.get_k() - node_length;
            assert(graph_.get_node_sequence(query_nodes_[i]).substr(offset)
                == std::string(query_.data() + i, query_.data() + i + node_length));

            if (offset)
                --num_matching_kmers_;

            if (node_length >= config_.min_seed_length) {
                offsets_[i] = offset;
            } else {
                query_nodes_[i] = 0;
                offsets_[i] = 0;
            }
        }
    }

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
void ExactMapSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const auto &graph = this->get_graph();
    size_t k = graph.get_k();
    const auto &config = this->get_config();
    const auto &query_nodes = this->get_query_nodes();
    const auto &partial_sum = this->get_partial_sums();
    const auto &offsets = this->get_offsets();
    auto query = this->get_query();
    bool orientation = this->get_orientation();

    // if node suffixes of length < k were matched, then query_nodes.size() may
    // be greater than query.size() - k + 1
    if (this->get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

    for (size_t i = 0; i < query_nodes.size(); ++i) {
        if (query_nodes[i] == DeBruijnGraph::npos)
            continue;

        size_t seed_length = k - offsets[i];

        assert(i + seed_length <= query.size());
        assert(seed_length >= config.min_seed_length);

        score_t match_score = partial_sum[i + seed_length] - partial_sum[i];

        if (match_score <= config.min_cell_score)
            continue;

        Seed seed(std::string_view(query.data() + i, seed_length),
                  { query_nodes[i] },
                  match_score,
                  i,
                  orientation,
                  offsets[i]);

        assert(seed.is_valid(graph, &config));
        callback(std::move(seed));
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    ExactMapSeeder<NodeType>::initialize(query, orientation);

    const auto &query_nodes = this->get_query_nodes();
    auto &offsets = this->get_offsets();

    query_node_flags_.assign(query_nodes.size(), 0);

    for (size_t i = 0; i < query_nodes.size(); ++i) {
        if (query_nodes[i] != DeBruijnGraph::npos) {
            query_node_flags_[i] = 2
                | ((offsets[i] > 0)
                || (i + 1 < query_nodes.size() && offsets[i + 1] > 0)
                || (i && offsets[i - 1] > 0)
                || (*is_mem_terminus_)[query_nodes[i]]);
        }
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    const auto &graph = this->get_graph();
    size_t k = graph.get_k();

    auto query = this->get_query();
    const auto &query_nodes = this->get_query_nodes();
    const auto &offsets = this->get_offsets();
    const auto &config = this->get_config();
    bool orientation = this->get_orientation();
    const auto &partial_sum = this->get_partial_sums();

    if (this->get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

    // find start of MEM
    auto it = query_node_flags_.begin();
    while ((it = std::find_if(it, query_node_flags_.end(),
                              [](uint8_t flags) { return flags & 2; }))
            < query_node_flags_.end()) {
        // find end of MEM
        auto next = std::find_if(it, query_node_flags_.end(),
                                 [](uint8_t flags) { return (flags & 1) == 1 || (flags & 2) == 0; });

        if (next != query_node_flags_.end() && ((*next) & 2))
            ++next;

        assert(next > it);
        assert(next <= query_node_flags_.end());

        size_t i = it - query_node_flags_.begin();
        assert(it == query_node_flags_.end() || query_nodes[i] != DeBruijnGraph::npos);

        size_t mem_length = (next - it) + k - 1 - offsets.at(i);

        assert(i + mem_length <= query.size());

        const char *begin_it = query.data() + i;
        const char *end_it = begin_it + mem_length;

        assert(end_it >= begin_it + config.min_seed_length);

        score_t match_score = partial_sum[end_it - query.data()]
            - partial_sum[begin_it - query.data()];

        auto node_begin_it = query_nodes.begin() + (it - query_node_flags_.begin());
        auto node_end_it = node_begin_it + (next - it);
        assert(std::find(node_begin_it, node_end_it, DeBruijnGraph::npos) == node_end_it);

        if (match_score > config.min_cell_score) {
            Seed seed(std::string_view(begin_it, end_it - begin_it),
                      std::vector<NodeType>(node_begin_it, node_end_it),
                      match_score,
                      begin_it - query.data(),
                      orientation,
                      offsets.at(i));
            assert(seed.is_valid(graph, &config));
            callback(std::move(seed));
        }

        it = next;
    }
}


template class ExactMapSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;

} // namespace align
} // namespace graph
} // namespace mtg
