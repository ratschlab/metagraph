#ifndef __ALIGNER_SEEDER_METHODS_HPP__
#define __ALIGNER_SEEDER_METHODS_HPP__

#include "alignment.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

// The symbols suffix_to_prefix tries at every depth unless its caller names others: every
// symbol of the graph's alphabet except the sentinel, s = 1 .. alph_size - 1 (on a DNA5
// build this includes N). This is the loop suffix_to_prefix always had; the SuffixSeeder
// relies on it unchanged (DESIGN-pattern-search.md §11).
struct NonSentinelSymbols {
    template <class BOSSEdgeRange, class TrySymbol>
    void operator()(const boss::BOSS &boss, const BOSSEdgeRange &, const TrySymbol &try_symbol) const {
        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
            try_symbol(s);
        }
    }
};

/**
 * Starting from a range of nodes sharing a suffix of length std::get<2>(index_range) (in
 * [1, k - 1]), extend the suffix symbol by symbol (depth first) up to whole (k-1)-mer nodes
 * and call every valid edge leaving those nodes: the k-mers having the suffix as a prefix.
 *
 * |symbols| chooses the symbols tried at each step: it is called once per range taken from
 * the stack, as symbols(boss, range, try_symbol), where |range| has its length already
 * incremented to the length of the ranges it is about to produce, and calls
 * try_symbol(s) for each symbol s to append, in the order the ranges are to be pushed.
 * The default (NonSentinelSymbols) tries every non-sentinel symbol. The pattern search
 * (pattern_search.cpp) passes the symbols a pattern allows at each depth.
 */
template <class BOSSEdgeRange, class SymbolSet = NonSentinelSymbols>
void suffix_to_prefix(const DBGSuccinct &dbg_succ,
                      const BOSSEdgeRange &index_range,
                      const std::function<void(DBGSuccinct::node_index)> &callback,
                      const SymbolSet &symbols = SymbolSet()) {
    const auto &boss = dbg_succ.get_boss();
    assert(std::get<2>(index_range));
    assert(std::get<2>(index_range) < dbg_succ.get_k());

    auto call_nodes_in_range = [&](const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        assert(seed_length == boss.get_k());
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            DBGSuccinct::node_index node = dbg_succ.validate_edge(i);
            if (node)
                callback(node);
        }
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }

    std::vector<BOSSEdgeRange> range_stack { index_range };

    while (range_stack.size()) {
        BOSSEdgeRange cur_range = std::move(range_stack.back());
        range_stack.pop_back();
        assert(std::get<2>(cur_range) < boss.get_k());
        ++std::get<2>(cur_range);

        symbols(boss, cur_range, [&](boss::BOSS::TAlphabet s) {
            auto next_range = cur_range;
            auto &[first, last, seed_length] = next_range;

            if (boss.tighten_range(&first, &last, s)) {
                if (seed_length == boss.get_k()) {
                    call_nodes_in_range(next_range);
                } else {
                    range_stack.emplace_back(std::move(next_range));
                }
            }
        });
    }
}

class ISeeder {
  public:
    virtual ~ISeeder() {}

    virtual const DBGAlignerConfig& get_config() const = 0;
    virtual std::vector<Seed> get_seeds() const = 0;
    virtual size_t get_num_matches() const = 0;

    virtual std::vector<Alignment> get_alignments() const {
        std::vector<Alignment> alignments;
        std::vector<Seed> seeds = get_seeds();
        alignments.reserve(seeds.size());
        for (const Seed &seed : seeds) {
            alignments.emplace_back(seed, get_config());
            alignments.back().trim_offset();
        }
        return alignments;
    }
};

class ManualMatchingSeeder : public ISeeder {
  public:
    ManualMatchingSeeder(std::vector<Seed>&& seeds,
                         size_t num_matching,
                         const DBGAlignerConfig &config)
          : config_(config), seeds_(std::move(seeds)), num_matching_(num_matching) {}

    virtual ~ManualMatchingSeeder() {}

    std::vector<Seed> get_seeds() const override { return seeds_; }
    const DBGAlignerConfig& get_config() const override { return config_; }
    size_t get_num_matches() const override final { return num_matching_; }
    std::vector<Seed>& data() { return seeds_; }

  private:
    const DBGAlignerConfig &config_;
    std::vector<Seed> seeds_;
    size_t num_matching_;
};

class ManualSeeder : public ISeeder {
  public:
    ManualSeeder(std::vector<Alignment>&& seeds = {}, size_t num_matching = 0)
        : seeds_(std::move(seeds)), num_matching_(num_matching) {}

    virtual ~ManualSeeder() {}

    std::vector<Seed> get_seeds() const override {
        throw std::runtime_error("Not implemented");
    }

    const DBGAlignerConfig& get_config() const override {
        throw std::runtime_error("Not implemented");
    }

    std::vector<Alignment> get_alignments() const override { return seeds_; }
    size_t get_num_matches() const override final { return num_matching_; }

    std::vector<Alignment>& data() { return seeds_; }

  private:
    std::vector<Alignment> seeds_;
    size_t num_matching_;
};

class ExactSeeder : public ISeeder {
  public:
    typedef DeBruijnGraph::node_index node_index;

    ExactSeeder(const DeBruijnGraph &graph,
                std::string_view query,
                bool orientation,
                std::vector<node_index>&& nodes,
                const DBGAlignerConfig &config);

    virtual ~ExactSeeder() {}

    const DBGAlignerConfig& get_config() const override { return config_; }
    std::vector<Seed> get_seeds() const override;
    size_t get_num_matches() const override final { return num_matching_; }

  protected:
    const DeBruijnGraph &graph_;
    std::string_view query_;
    bool orientation_;
    std::vector<node_index> query_nodes_;
    const DBGAlignerConfig &config_;
    size_t num_matching_;

    size_t num_exact_matching() const;
};

class MEMSeeder : public ExactSeeder {
  public:
    template <typename... Args>
    MEMSeeder(Args&&... args) : ExactSeeder(std::forward<Args>(args)...) {}

    virtual ~MEMSeeder() {}

    std::vector<Seed> get_seeds() const override;

    virtual const bitmap& get_mem_terminator() const = 0;
};

class UniMEMSeeder : public MEMSeeder {
  public:
    template <typename... Args>
    UniMEMSeeder(Args&&... args)
          : MEMSeeder(std::forward<Args>(args)...),
            is_mem_terminus_([&](auto i) {
                                 return graph_.has_multiple_outgoing(i)
                                     || !graph_.has_single_incoming(i);
                             },
                             graph_.max_index() + 1) {
        assert(is_mem_terminus_.size() == graph_.max_index() + 1);
    }

    virtual ~UniMEMSeeder() {}

    const bitmap& get_mem_terminator() const override { return is_mem_terminus_; }

  private:
    bitmap_lazy is_mem_terminus_;
};

template <class BaseSeeder>
class SuffixSeeder : public BaseSeeder {
  public:
    template <typename... Args>
    SuffixSeeder(Args&&... args) : BaseSeeder(std::forward<Args>(args)...) {
        generate_seeds();
    }

    virtual ~SuffixSeeder() {}

    std::vector<Seed> get_seeds() const override { return seeds_; }

  protected:
    void generate_seeds();

  private:
    std::vector<Seed> seeds_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_SEEDER_METHODS_HPP__
