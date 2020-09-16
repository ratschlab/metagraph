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
    num_matching_kmers_ = 0;

    if (query_.size() < graph_.get_k())
        return;

    query_nodes_ = map_sequence_to_nodes(graph_, query_);
    num_matching_kmers_ = query_nodes_.size()
        - std::count(query_nodes_.begin(), query_nodes_.end(), NodeType());

    offsets_.assign(query_nodes_.size(), 0);

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
void SuffixSeeder<NodeType>
::initialize(std::string_view query, bool orientation) {
    const auto &graph = dynamic_cast<const DBGSuccinct&>(get_graph());
    const auto &config = get_config();
    auto &query_nodes = get_query_nodes();
    auto &offsets = get_offsets();
    alt_query_nodes.clear();

    base_seeder_->initialize(query, orientation);
    size_t k = graph.get_k();

    // if node suffixes of length < k were matched, then query_nodes.size() may
    // be greater than query.size() - k + 1
    if (get_num_matching_kmers()
            < config.exact_kmer_match_fraction * (query.size() - k + 1)) {
        return;
    }

    if (query_nodes.empty())
        query_nodes.assign(1, DeBruijnGraph::npos);

    size_t max_seed_length = std::min(query.size(), std::min(config.max_seed_length, k));

    if (config.max_num_seeds_per_locus > 1)
        alt_query_nodes.resize(query_nodes.size());

    auto *unimem_seeder = dynamic_cast<MEMSeeder<NodeType>*>(base_seeder_.get());
    auto *query_node_flags = unimem_seeder ? &unimem_seeder->get_query_node_flags() : nullptr;

    for (size_t i = 0; i < query_nodes.size(); ++i) {
        if (query_nodes[i] != DeBruijnGraph::npos)
            continue;

        bool is_alt = false;
        graph.call_nodes_with_suffix(
            std::string_view(query.data() + i, max_seed_length),
            [&](NodeType alt_node, size_t seed_length) {
                assert(alt_node != DeBruijnGraph::npos);
                assert(seed_length < k);
                assert(i + seed_length + 1 == query.size()
                    || graph.traverse(alt_node, query[i + seed_length])
                        == DeBruijnGraph::npos);

                if (config.max_num_seeds_per_locus != 1) {
                    for (size_t j = 1; seed_length + j < k && i >= j; ++j) {
                        if (alt_node == query_nodes[i - j]
                                && k - offsets[i - j] - j == seed_length) {
                            return;
                        }
                    }
                }

                if (is_alt) {
                    alt_query_nodes.at(i).push_back(alt_node);
                    assert(offsets.at(i) == k - seed_length);
                } else {
                    query_nodes.at(i) = alt_node;
                    offsets.at(i) = k - seed_length;

                    if (query_node_flags) {
                        query_node_flags->at(i) = 3;

                        if (i)
                            query_node_flags->at(i - 1) |= 1;

                        if (i + 1 < query_node_flags->size())
                            query_node_flags->at(i + 1) |= 1;
                    }
                }

                is_alt = true;
            },
            config.min_seed_length,
            config.max_num_seeds_per_locus
        );
    }
}

template <typename NodeType>
void SuffixSeeder<NodeType>::call_seeds(std::function<void(Seed&&)> callback) const {
    auto query = get_query();
    base_seeder_->call_seeds([&](auto&& alignment) {
        size_t i = alignment.get_query().data() - query.data();

        assert(get_config().max_num_seeds_per_locus > alt_query_nodes.size());
        if (alt_query_nodes.size()) {
            for (auto alt_node : alt_query_nodes.at(i)) {
                auto alt_alignment = alignment;
                alt_alignment.nodes_.front() = alt_node;
                assert(alt_alignment.is_valid(get_graph(), &get_config()));
                callback(std::move(alt_alignment));
            }
        }

        callback(std::move(alignment));
    });

    const auto &query_nodes = base_seeder_->get_query_nodes();
    if (query_nodes.back())
        return;

    const auto &config = get_config();
    const auto &graph = dynamic_cast<const DBGSuccinct&>(get_graph());
    size_t k = graph.get_k();
    bool orientation = this->get_orientation();

    for (size_t i = query.size() - k + 1; i + config.min_seed_length <= query.size(); ++i) {
        graph.call_nodes_with_suffix(
            std::string_view(query.data() + i, query.size() - i),
            [&](NodeType alt_node, size_t seed_length) {
                assert(alt_node != DeBruijnGraph::npos);
                assert(seed_length < k);
                assert(i + seed_length + 1 == query.size()
                    || graph.traverse(alt_node, query[i + seed_length])
                        == DeBruijnGraph::npos);

                std::string_view seed_seq(query.data() + i, seed_length);
                DBGAlignerConfig::score_t match_score = config.match_score(seed_seq);

                if (match_score <= config.min_cell_score)
                    return;

                Seed seed(seed_seq, { alt_node }, match_score, i,
                          orientation, k - seed_length);

                assert(seed.is_valid(graph, &config));
                callback(std::move(seed));
            },
            config.min_seed_length,
            config.max_num_seeds_per_locus
        );
    }
}

template <typename NodeType>
void MEMSeeder<NodeType>::initialize(std::string_view query, bool orientation) {
    ExactMapSeeder<NodeType>::initialize(query, orientation);

    const auto &query_nodes = this->get_query_nodes();

    query_node_flags_.assign(query_nodes.size(), 0);

    for (size_t i = 0; i < query_nodes.size(); ++i) {
        if (query_nodes[i] != DeBruijnGraph::npos) {
            query_node_flags_[i] = 2 | (*is_mem_terminus_)[query_nodes[i]];
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


template class ExactSeeder<>;
template class SuffixSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;

} // namespace align
} // namespace graph
} // namespace mtg
