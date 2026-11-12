#include "gtest/gtest.h"

#include <random>
#include <set>

#include "graph/representation/succinct/boss_pattern.hpp"


#if ! _PROTEIN_GRAPH

namespace {

using namespace mtg::graph::boss;

typedef BOSS::edge_index edge_index;

// (edge, offset)
typedef std::set<std::pair<edge_index, size_t>> Contexts;

const std::chrono::milliseconds kNoTimeLimit(1'000'000);

std::string random_dna(std::mt19937 &gen, size_t length) {
    std::string s(length, 'A');
    for (char &c : s) {
        c = "ACGT"[gen() % 4];
    }
    return s;
}

// a graph of random sequences, some of them shorter than a k-mer
std::unique_ptr<BOSS> random_graph(std::mt19937 &gen, size_t boss_k) {
    auto boss = std::make_unique<BOSS>(boss_k);
    for (size_t i = 0; i < 6; ++i) {
        boss->add_sequence(random_dna(gen, 1 + gen() % (4 * boss_k + 4)));
    }
    return boss;
}

bool matches(const BOSS &boss, const Pattern &pattern, std::string_view kmer, size_t offset) {
    for (size_t i = 0; i < pattern.size(); ++i) {
        // the sentinel has code 0, which is in no set
        if (!((pattern[i] >> boss.encode(kmer[offset + i])) & 1))
            return false;
    }
    return true;
}

// scan the spelling of every edge
Contexts brute_force(const BOSS &boss, const Pattern &pattern) {
    Contexts contexts;
    for (edge_index i = 1; i <= boss.num_edges(); ++i) {
        if (!(boss.get_W(i) % boss.alph_size))
            continue;

        std::string kmer = boss.get_node_str(i) + boss.decode(boss.get_W(i) % boss.alph_size);
        for (size_t offset = 0; offset + pattern.size() <= kmer.size(); ++offset) {
            if (matches(boss, pattern, kmer, offset))
                contexts.emplace(i, offset);
        }
    }
    return contexts;
}

Contexts search(const BOSS &boss, const Pattern &pattern) {
    Contexts contexts;
    SearchBudget budget(1'000'000, kNoTimeLimit);
    EXPECT_TRUE(call_pattern_ranges(boss, pattern,
        [&](edge_index rl, edge_index ru, size_t offset, SymbolSet symbols) {
            uint64_t called = 0;
            EXPECT_TRUE(call_edges(boss, rl, ru, symbols, [&](edge_index i) {
                EXPECT_TRUE(rl <= i && i <= ru);
                EXPECT_TRUE(contexts.emplace(i, offset).second) << "an edge called twice";
                ++called;
            }, &budget));
            EXPECT_EQ(called, count_edges(boss, rl, ru, symbols));
        },
        &budget
    ));
    EXPECT_EQ(nullptr, budget.stop());
    return contexts;
}

void test_patterns(const std::string &alphabet) {
    std::mt19937 gen(1);
    for (size_t boss_k : { 1, 2, 3, 5, 8 }) {
        for (size_t i = 0; i < 5; ++i) {
            auto boss = random_graph(gen, boss_k);
            for (BOSS::State state : { BOSS::State::DYN, BOSS::State::STAT,
                                       BOSS::State::SMALL, BOSS::State::FAST }) {
                boss->switch_state(state);
                // every length from 1 to k, where k = boss_k + 1
                for (size_t length = 1; length <= boss_k + 1; ++length) {
                    for (size_t j = 0; j < 5; ++j) {
                        std::string sequence(length, 'A');
                        for (char &c : sequence) {
                            c = alphabet[gen() % alphabet.size()];
                        }
                        Pattern pattern = parse_pattern(*boss, sequence);
                        EXPECT_EQ(brute_force(*boss, pattern), search(*boss, pattern))
                            << "k: " << boss_k + 1 << ", pattern: " << sequence
                            << ", state: " << state;
                    }
                }
            }
        }
    }
}

TEST(BOSSPattern, Exact) {
    test_patterns("ACGT");
}

TEST(BOSSPattern, Iupac) {
    test_patterns("ACGTRYSWKMBDHVNacgtn");
}

TEST(BOSSPattern, AllN) {
    test_patterns("N");
}

TEST(BOSSPattern, ExactFindsEveryKmerOfASequence) {
    BOSS boss(4);
    const std::string sequence = "ACGTTGCAAGGCT";
    boss.add_sequence(sequence);
    for (size_t i = 0; i + 5 <= sequence.size(); ++i) {
        Contexts contexts = search(boss, parse_pattern(boss, sequence.substr(i, 5)));
        ASSERT_EQ(1u, contexts.size());
        EXPECT_EQ(0u, contexts.begin()->second);
        EXPECT_EQ(sequence.substr(i, 4), boss.get_node_str(contexts.begin()->first));
    }
}

TEST(BOSSPattern, Parse) {
    BOSS boss(3);
    EXPECT_EQ(parse_pattern(boss, "ACGTN"), parse_pattern(boss, "acgtn"));
    EXPECT_EQ(parse_pattern(boss, "AG"), reverse_complement_pattern(parse_pattern(boss, "CT")));
    // R = A|G and Y = C|T, B = not A and V = not T
    EXPECT_EQ(parse_pattern(boss, "RYNBV"),
              reverse_complement_pattern(parse_pattern(boss, "BVNRY")));
    EXPECT_NE(parse_pattern(boss, "R"), parse_pattern(boss, "Y"));
    EXPECT_TRUE(parse_pattern(boss, "").empty());
    EXPECT_THROW(parse_pattern(boss, "AC$"), std::invalid_argument);
    EXPECT_THROW(parse_pattern(boss, "ACX"), std::invalid_argument);
    EXPECT_THROW(parse_pattern(boss, "A-C"), std::invalid_argument);
}

TEST(BOSSPattern, StepBound) {
    std::mt19937 gen(2);
    auto boss = random_graph(gen, 8);
    Pattern pattern = parse_pattern(*boss, "NNN");
    auto ignore = [](edge_index, edge_index, size_t, SymbolSet) {};

    SearchBudget unbounded(1'000'000, kNoTimeLimit);
    EXPECT_TRUE(call_pattern_ranges(*boss, pattern, ignore, &unbounded));
    EXPECT_EQ(nullptr, unbounded.stop());

    for (uint64_t max_steps : { 0, 1, 10 }) {
        SearchBudget budget(max_steps, kNoTimeLimit);
        EXPECT_FALSE(call_pattern_ranges(*boss, pattern, ignore, &budget));
        EXPECT_STREQ("max_steps", budget.stop());
        EXPECT_FALSE(budget.step());
    }

    // the edges of a range are steps too
    SearchBudget budget(3, kNoTimeLimit);
    uint64_t called = 0;
    EXPECT_FALSE(call_edges(*boss, 1, boss->num_edges(), all_symbols(*boss),
                            [&](edge_index) { ++called; }, &budget));
    EXPECT_EQ(3u, called);
}

TEST(BOSSPattern, TimeBound) {
    SearchBudget budget(1'000'000, std::chrono::milliseconds(-1));
    // the clock is read every 1024 steps
    for (size_t i = 0; i < 1023; ++i) {
        EXPECT_TRUE(budget.step());
    }
    EXPECT_FALSE(budget.step());
    EXPECT_STREQ("time_budget", budget.stop());
}

} // namespace

#endif // ! _PROTEIN_GRAPH
