#ifndef __METAGRAPH_TESTS_PATTERN_TEST_SUPPORT_HPP__
#define __METAGRAPH_TESTS_PATTERN_TEST_SUPPORT_HPP__

// Helpers of the pattern-search tests that the engine and the route do not use

#include <cassert>
#include <cstdint>
#include <limits>
#include <string>

#include <sdsl/int_vector.hpp>

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace test {

// A deadline that never expires: an infinite time budget and no reserve (tests that do not
// test time)
inline graph::pattern::Deadline unbounded_deadline() {
    return graph::pattern::Deadline(graph::pattern::Deadline::Clock::now(),
                                    std::numeric_limits<double>::infinity(), 0);
}

// The three bases of a codon index (graph::pattern::codon_index's order), upper case
inline std::string codon_string(int index) {
    static constexpr char kBases[] = "ACGT";
    return { kBases[(index >> 4) & 3], kBases[(index >> 2) & 3], kBases[index & 3] };
}

/**
 * f over every entry with W != $ (graph::pattern::RealFraction): lower == upper == value, the
 * source dummies found by BOSS's own traversal of the dummy tree (BOSS::mark_source_dummy_edges,
 * as `stats --count-dummy`), not by the test sample_real_fraction draws with. O(edges) bits.
 */
inline graph::pattern::RealFraction exact_real_fraction(const graph::DBGSuccinct &graph) {
    const graph::boss::BOSS &boss = graph.get_boss();
    // the edges, the W = $ entries and the seed as sample_real_fraction counts them
    graph::pattern::RealFraction f = graph::pattern::sample_real_fraction(graph, 0);
    f.samples = f.edges - f.sentinel_edges;
    sdsl::bit_vector source_dummies(boss.get_W().size(), false);
    boss.mark_source_dummy_edges(&source_dummies, 1);
    uint64_t dummies = 0;
    for (uint64_t e = 1; e < source_dummies.size(); ++e) {
        // (W modulo alph_size is 0 for $, plain or marked: BOSS::kSentinelCode)
        if (source_dummies[e] && boss.get_W(e) % boss.alph_size)
            ++dummies;
    }
    assert(dummies <= f.samples);
    f.real = f.samples - dummies;
    if (f.samples) {
        f.value = static_cast<double>(f.real) / static_cast<double>(f.samples);
        f.lower = f.upper = f.value;
    }
    return f;
}

} // namespace test
} // namespace mtg

#endif // __METAGRAPH_TESTS_PATTERN_TEST_SUPPORT_HPP__
