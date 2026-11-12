#ifndef __BOSS_PATTERN_HPP__
#define __BOSS_PATTERN_HPP__

#include <chrono>
#include <functional>
#include <string_view>
#include <vector>

#include "boss.hpp"


namespace mtg {
namespace graph {
namespace boss {

// A set of symbols: bit c is set when the symbol with code c is in the set
typedef uint32_t SymbolSet;
// One set per position. The sentinel is in no set.
typedef std::vector<SymbolSet> Pattern;

// Read a pattern. A DNA graph takes the 15 IUPAC codes in any case. On a case-sensitive
// DNA graph a code stands for bases of its own case. A protein graph takes its residues,
// and X for any. Throws std::invalid_argument on another character.
Pattern parse_pattern(const BOSS &boss, std::string_view sequence);

#if ! _PROTEIN_GRAPH
Pattern reverse_complement_pattern(Pattern pattern);
#endif

// The bounds of a search. A step is one tighten_range or one edge visited.
class SearchBudget {
  public:
    SearchBudget(uint64_t max_steps, std::chrono::milliseconds max_time)
          : max_steps_(max_steps),
            deadline_(std::chrono::steady_clock::now() + max_time) {}

    // count a step, and return false once a bound is reached
    bool step();

    // the bound that was reached ("max_steps" or "time_budget"), or NULL
    const char* stop() const { return stop_; }

  private:
    uint64_t max_steps_;
    std::chrono::steady_clock::time_point deadline_;
    uint64_t steps_ = 0;
    const char *stop_ = nullptr;
};

// Every symbol but the sentinel
SymbolSet all_symbols(const BOSS &boss);

// The edges of [rl, ru] whose symbol, plain or marked, is in |symbols|
uint64_t count_edges(const BOSS &boss, BOSS::edge_index rl, BOSS::edge_index ru, SymbolSet symbols);
// Call them symbol by symbol, one step each. Returns false if |budget| stopped the calls.
bool call_edges(const BOSS &boss, BOSS::edge_index rl, BOSS::edge_index ru, SymbolSet symbols,
                const std::function<void(BOSS::edge_index)> &callback,
                SearchBudget *budget);

/**
 * Find the k-mers that hold |pattern|. A k-mer is an edge: the get_k() symbols of its
 * source node and the edge symbol. The pattern has 1 to get_k() + 1 symbols.
 * Every k-mer that holds the pattern at |offset| is an edge of exactly one range called
 * with that offset, and its symbol is in |symbols|. All such edges of a range hold the
 * pattern, except that the symbols left of the pattern were not read when |offset| > 0:
 * they may hold the sentinel of a dummy source k-mer.
 * Returns false if |budget| stopped the search.
 */
bool call_pattern_ranges(const BOSS &boss,
                         const Pattern &pattern,
                         const std::function<void(BOSS::edge_index /* rl */,
                                                  BOSS::edge_index /* ru */,
                                                  size_t /* offset */,
                                                  SymbolSet /* symbols */)> &callback,
                         SearchBudget *budget);

} // namespace boss
} // namespace graph
} // namespace mtg

#endif // __BOSS_PATTERN_HPP__
