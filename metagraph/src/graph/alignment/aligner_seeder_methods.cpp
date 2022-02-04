#include "aligner_seeder_methods.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/utils/template_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace align {

ExactSeeder::ExactSeeder(const DeBruijnGraph &graph,
                         std::string_view query,
                         bool orientation,
                         const std::vector<node_index> &nodes,
                         const DBGAlignerConfig &config)
      : graph_(graph),
        query_(query),
        orientation_(orientation),
        query_nodes_(nodes),
        config_(config),
        num_matching_(num_exact_matching()) {
    assert(config_.check_config_scores());

    partial_sum_.resize(query_.size() + 1);
    std::transform(query_.begin(), query_.end(),
                   partial_sum_.begin() + 1,
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sum_.begin(), partial_sum_.end(), partial_sum_.begin());
    assert(config_.match_score(query_) == partial_sum_.back());
    assert(config_.get_row(query_.front())[query_.front()] == partial_sum_[1]);
    assert(!partial_sum_.front());
}

size_t ExactSeeder::num_exact_matching() const {
    size_t num_matching = 0;
    size_t last_match_count = 0;
    for (auto it = query_nodes_.begin(); it != query_nodes_.end(); ++it) {
        if (*it) {
            auto jt = std::find(it + 1, query_nodes_.end(), node_index());
            num_matching += graph_.get_k() + std::distance(it, jt) - 1 - last_match_count;
            last_match_count = graph_.get_k();
            it = jt - 1;
        } else if (last_match_count) {
            --last_match_count;
        }
    }
    assert(num_matching <= query_.size());

    return num_matching;
}

auto ExactSeeder::get_seeds() const -> std::vector<Seed> {
    size_t k = graph_.get_k();
    assert(k >= config_.min_seed_length);

    if (num_matching_ < config_.min_exact_match * query_.size())
        return {};

    std::vector<Seed> seeds;

    if (config_.max_seed_length < k)
        return seeds;

    const char *query_end = query_.data() + query_.size();

    for (size_t i = 0; i < query_nodes_.size(); ++i) {
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            assert(i + k <= query_.size());

            score_t match_score = partial_sum_[i + k] - partial_sum_[i]
                + !i * config_.left_end_bonus
                + (i + 1 == query_nodes_.size()) * config_.right_end_bonus;

            if (match_score > config_.min_cell_score) {
                seeds.emplace_back(query_.substr(i, k),
                                   std::vector<node_index>{ query_nodes_[i] },
                                   std::string(query_.substr(i, k)), match_score,
                                   i, orientation_);
                seeds.back().extend_query_end(query_end);
                assert(seeds.back().is_valid(graph_, &config_));
            }
        }
    }

    return seeds;
}

template <class BOSSEdgeRange>
void suffix_to_prefix(const DBGSuccinct &dbg_succ,
                      const BOSSEdgeRange &index_range,
                      const std::function<void(DBGSuccinct::node_index)> &callback) {
    const auto &boss = dbg_succ.get_boss();
    assert(std::get<2>(index_range));
    assert(std::get<2>(index_range) < dbg_succ.get_k());

    auto call_nodes_in_range = [&](const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        assert(seed_length == boss.get_k());
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            DBGSuccinct::node_index node = dbg_succ.boss_to_kmer_index(i);
            if (node)
                callback(node);
        }
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }

    std::vector<BOSSEdgeRange> range_stack{ index_range };

    while (range_stack.size()) {
        BOSSEdgeRange cur_range = std::move(range_stack.back());
        range_stack.pop_back();
        assert(std::get<2>(cur_range) < boss.get_k());
        ++std::get<2>(cur_range);

        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
            auto next_range = cur_range;
            auto &[first, last, seed_length] = next_range;

            if (boss.tighten_range(&first, &last, s)) {
                if (seed_length == boss.get_k()) {
                    call_nodes_in_range(next_range);
                } else {
                    range_stack.emplace_back(std::move(next_range));
                }
            }
        }
    }
}

template <class BaseSeeder>
void SuffixSeeder<BaseSeeder>::generate_seeds() {
    // this method assumes that seeds from the BaseSeeder are exact match only
    static_assert(std::is_base_of_v<ExactSeeder, BaseSeeder>);

    if (this->config_.min_seed_length >= this->graph_.get_k()) {
        seeds_ = this->BaseSeeder::get_seeds();
        return;
    }

    const DBGSuccinct &dbg_succ = get_base_dbg_succ(this->graph_);

    if (this->query_.size() < this->config_.min_seed_length)
        return;

    std::vector<std::vector<Seed>> suffix_seeds(
        this->query_.size() - this->config_.min_seed_length + 1
    );

    std::vector<size_t> min_seed_length(
        this->query_.size() - this->config_.min_seed_length + 1,
        this->config_.min_seed_length
    );

    for (auto&& seed : this->BaseSeeder::get_seeds()) {
        assert(seed.get_query().size() >= this->config_.min_seed_length);

        size_t i = seed.get_clipping();
        assert(i + seed.size() <= min_seed_length.size());

        for (size_t j = 0; j < seed.size(); ++j)
            min_seed_length[i + j] = this->graph_.get_k();

        if (i + seed.size() < min_seed_length.size())
            min_seed_length[i + seed.size()] = this->graph_.get_k();

        suffix_seeds[i].emplace_back(std::move(seed));
    }

    const char *query_end = this->query_.data() + this->query_.size();

    // when a seed is found, append it to the seed vector
    auto append_suffix_seed = [&](size_t i, node_index alt_node, size_t seed_length) {
        assert(i < suffix_seeds.size());

        std::string_view seed_seq = this->query_.substr(i, seed_length);
        score_t match_score = this->config_.match_score(seed_seq)
            + !i * this->config_.left_end_bonus
            + (i + seed_seq.size() == this->query_.size()) * this->config_.right_end_bonus;

        if (match_score <= this->config_.min_cell_score)
            return;

        if (seed_length > min_seed_length[i])
            suffix_seeds[i].clear();

        min_seed_length[i] = seed_length;

        assert(seed_length == min_seed_length[i]);
        suffix_seeds[i].emplace_back(seed_seq, std::vector<node_index>{ alt_node },
                                     std::string(seed_seq), match_score, i,
                                     this->orientation_,
                                     this->graph_.get_k() - seed_length);
        suffix_seeds[i].back().extend_query_end(query_end);
        assert(suffix_seeds[i].back().is_valid(this->graph_, &this->config_));

        for (++i; i < min_seed_length.size() && seed_length > min_seed_length[i]; ++i) {
            min_seed_length[i] = seed_length--;
            suffix_seeds[i].clear();
        }
    };

    // find sub-k matches in the forward orientation
    size_t last_full_id = this->query_.size() >= this->graph_.get_k()
        ? this->query_.size() - this->graph_.get_k() + 1
        : min_seed_length.size();
    for (size_t i = 0; i < min_seed_length.size(); ++i) {
        size_t max_seed_length = std::min({ this->config_.max_seed_length,
                                            this->graph_.get_k() - 1,
                                            this->query_.size() - i });
        size_t seed_length = 0;
        std::vector<node_index> alt_nodes;

        dbg_succ.call_nodes_with_suffix_matching_longest_prefix(
            this->query_.substr(i, max_seed_length),
            [&](node_index alt_node, size_t len) {
                seed_length = len;
                alt_nodes.push_back(alt_node);
            },
            min_seed_length[i]
        );

        if (i >= last_full_id && alt_nodes.size() == 1
                && min_seed_length[last_full_id - 1] == this->graph_.get_k()
                && suffix_seeds[last_full_id - 1].size() == 1
                && alt_nodes[0] == suffix_seeds[last_full_id - 1][0].get_nodes()[0])
            continue;

        for (node_index alt_node : alt_nodes) {
            append_suffix_seed(i, alt_node, seed_length);
        }
    }

    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_)) {
        // find sub-k matches in the reverse complement
        // TODO: find sub-k seeds which are sink tips in the underlying graph
        std::string query_rc(this->query_);
        reverse_complement(query_rc.begin(), query_rc.end());

        // matching is done query prefix -> node suffix, so the query index of
        // a match to the reverse complement is not known beforehand
        // e.g.,
        // k = 6;
        // rev: rev_end_pos = 8
        //     j
        //     ****--      <-- start position in forward depends on match length
        // GCTAGCATCTGAGAGGGGA fwd
        // TCCCCTCTCAGATGCTAGC rc
        //          --****
        //            i    <-- match returned from call
        for (size_t i = 0; i + this->config_.min_seed_length <= query_rc.size(); ++i) {
            // initial estimate of the max seed length
            size_t max_seed_length = std::min({ this->config_.max_seed_length,
                                                this->graph_.get_k() - 1,
                                                this->query_.size() - i });

            // the reverse complement of the sub-k match will fall somewhere in this range
            size_t j_min = query_rc.size() - i - max_seed_length;
            size_t j_max = query_rc.size() - i - this->config_.min_seed_length;

            // skip over positions which have better matches
            while (j_min <= j_max && min_seed_length[j_min] > max_seed_length) {
                ++j_min;
                --max_seed_length;
            }

            if (j_min > j_max)
                continue;

            const auto &boss = dbg_succ.get_boss();

            auto encoded = boss.encode({ query_rc.data() + i, max_seed_length });
            auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());

            size_t seed_length = end - encoded.begin();
            size_t j = query_rc.size() - i - seed_length;

            assert(seed_length < this->config_.min_seed_length
                || j < min_seed_length.size());

            if (seed_length < this->config_.min_seed_length
                    || seed_length < min_seed_length[j]) {
                continue;
            }

            // e.g., matched: ***ATG, want ATG***
            suffix_to_prefix(
                dbg_succ,
                std::make_tuple(boss.pred_last(first - 1) + 1, last, seed_length),
                [&](node_index match) {
                    append_suffix_seed(j, canonical->reverse_complement(match), seed_length);
                }
            );
        }
    }

    // aggregate all seeds
    seeds_.clear();
    this->num_matching_ = 0;
    size_t last_end = 0;
    for (size_t i = 0; i < suffix_seeds.size(); ++i) {
        std::vector<Seed> &pos_seeds = suffix_seeds[i];
        if (pos_seeds.empty())
            continue;

        // all seeds should have the same properties, but they will be at different
        // graph nodes
        assert(std::equal(pos_seeds.begin() + 1, pos_seeds.end(), pos_seeds.begin(),
                          [](const Seed &a, const Seed &b) {
            return a.get_orientation() == b.get_orientation()
                && a.get_offset() == b.get_offset()
                && a.get_score() == b.get_score()
                && a.get_query() == b.get_query()
                && a.get_sequence() == b.get_sequence()
                && a.get_cigar() == b.get_cigar();
        }));

        if (!pos_seeds[0].get_offset()) {
            assert(min_seed_length[i] == this->graph_.get_k());
            assert(pos_seeds.size() == 1);
            seeds_.emplace_back(std::move(pos_seeds[0]));
        } else {
            assert(min_seed_length[i] == this->graph_.get_k() - pos_seeds[0].get_offset());
            if (pos_seeds.size() <= this->config_.max_num_seeds_per_locus) {
                for (auto&& seed : pos_seeds) {
                    seeds_.emplace_back(std::move(seed));
                }
            }
        }
        if (!pos_seeds[0].get_offset()
                || pos_seeds.size() <= this->config_.max_num_seeds_per_locus) {
            size_t begin = seeds_.back().get_clipping();
            size_t end = begin + seeds_.back().get_query().size();
            if (begin < last_end) {
                this->num_matching_ += end - begin - (last_end - begin);
            } else {
                this->num_matching_ += end - begin;
            }
            last_end = end;
        }
    }
}

auto MEMSeeder::get_seeds() const -> std::vector<Seed> {
    size_t k = graph_.get_k();

    if (k >= config_.max_seed_length)
        return ExactSeeder::get_seeds();

    if (num_matching_ < config_.min_exact_match * query_.size())
        return {};

    std::vector<uint8_t> query_node_flags(query_nodes_.size(), 0);
    for (size_t i = 0; i < query_node_flags.size(); ++i) {
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            // the second bit indicates that a node has been found, while the
            // first bit indicates if the node is a maximal exact match terminus
            query_node_flags[i] = 2 | get_mem_terminator()[query_nodes_[i]];
        }
    }

    std::vector<Seed> seeds;

    const char *query_end = query_.data() + query_.size();

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
        assert(it == query_node_flags.end()
                || query_nodes_[i] != DeBruijnGraph::npos);

        size_t mem_length = (next - it) + k - 1;
        assert(i + mem_length <= query_.size());

        if (mem_length >= config_.min_seed_length) {
            const char *begin_it = query_.data() + i;
            const char *end_it = begin_it + mem_length;

            score_t match_score = partial_sum_[end_it - query_.data()] - partial_sum_[i]
                + (!i) * config_.left_end_bonus
                + (i + mem_length == this->query_.size()) * config_.right_end_bonus;

            auto node_begin_it = query_nodes_.begin() + i;
            auto node_end_it = node_begin_it + (next - it);
            assert(std::find(node_begin_it, node_end_it, DeBruijnGraph::npos)
                    == node_end_it);

            if (match_score > config_.min_cell_score) {
                seeds.emplace_back(std::string_view(begin_it, mem_length),
                                   std::vector<node_index>{ node_begin_it, node_end_it },
                                   std::string(begin_it, begin_it + mem_length),
                                   match_score, i, orientation_);
                seeds.back().extend_query_end(query_end);
                assert(seeds.back().is_valid(graph_, &config_));
            }
        }

        it = next;
    }

    return seeds;
}

template <class BaseSeeder>
const DBGSuccinct& SuffixSeeder<BaseSeeder>
::get_base_dbg_succ(const DeBruijnGraph &graph) {
    try {
        if (const auto *wrapper = dynamic_cast<const DBGWrapper<>*>(&graph))
            return dynamic_cast<const DBGSuccinct&>(wrapper->get_graph());

        return dynamic_cast<const DBGSuccinct&>(graph);

    } catch (const std::bad_cast &e) {
        common::logger->error(
            "SuffixSeeder can be used only with succinct graph representation"
        );
        throw e;
    }
}

template class SuffixSeeder<ExactSeeder>;
template class SuffixSeeder<UniMEMSeeder>;

} // namespace align
} // namespace graph
} // namespace mtg
