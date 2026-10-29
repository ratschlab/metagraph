#include "graph/representation/succinct/suffix_matching.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace {

using namespace mtg::graph;
using Node = DBGSuccinct::node_index;

auto prefix_range(const DBGSuccinct &graph, const std::string &prefix) {
    const auto &boss = graph.get_boss();
    auto encoded = boss.encode(prefix);
    auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
    EXPECT_EQ(end, encoded.end());
    first = first == boss.select_last(1) ? 1 : boss.pred_last(first - 1) + 1;
    return std::make_tuple(first, last, prefix.size());
}

std::vector<Node> matching_nodes(const DBGSuccinct &graph, const std::string &prefix) {
    std::vector<Node> nodes;
    // The range helper honours the validity mask. Without one, it can return
    // sink edges, which call_nodes() omits when it enumerates real sequences.
    for (Node node = 1; node <= graph.max_index(); ++node) {
        if (graph.in_graph(node) && graph.get_node_sequence(node).starts_with(prefix))
            nodes.push_back(node);
    }
    std::sort(nodes.begin(), nodes.end());
    return nodes;
}

TEST(SuffixMatching, MatchesSpelledPrefixes) {
    for (size_t k = 2; k <= 6; ++k) {
        for (bool mask : {false, true}) {
            DBGSuccinct graph(k);
            graph.add_sequence("ACGACCCGTTACACGGATACGCGTTT");
            graph.add_sequence("ACCATGACGCATATACCGGTC");
            if (graph.alphabet().find('N') != std::string::npos)
                graph.add_sequence("ACNACNNGACNGT");
            if (mask)
                graph.mask_dummy_kmers(1, false);

            std::set<std::string> prefixes;
            graph.call_nodes([&](Node node) {
                const auto spelling = graph.get_node_sequence(node);
                for (size_t length = 1; length < k; ++length) {
                    auto prefix = spelling.substr(0, length);
                    if (prefix.find('$') == std::string::npos)
                        prefixes.insert(std::move(prefix));
                }
            });

            for (auto state : {boss::BOSS::DYN, boss::BOSS::STAT,
                               boss::BOSS::FAST, boss::BOSS::SMALL}) {
                graph.switch_state(state);
                for (const auto &prefix : prefixes) {
                    SCOPED_TRACE(::testing::Message() << k << '/' << mask << '/' << state
                                                     << '/' << prefix);
                    std::vector<Node> found;
                    suffix_to_prefix(graph, prefix_range(graph, prefix), [&](Node node) {
                        found.push_back(node);
                        return true;
                    });
                    std::sort(found.begin(), found.end());
                    EXPECT_EQ(found, matching_nodes(graph, prefix));
                }
            }
        }
    }
}

TEST(SuffixMatching, StopsAcrossRangesAndWithinFinalRange) {
    DBGSuccinct graph(4);
    for (const auto &sequence : {"AAAAC", "AACGT", "ACGTA", "AGGTC", "ATGAC"})
        graph.add_sequence(sequence);
    graph.mask_dummy_kmers(1, false);

    for (const std::string prefix : {"A", "AA", "AAA"}) {
        const auto range = prefix_range(graph, prefix);
        std::vector<Node> all;
        suffix_to_prefix(graph, range, [&](Node node) {
            all.push_back(node);
            return true;
        });
        ASSERT_FALSE(all.empty());
        for (size_t limit = 1; limit <= all.size(); ++limit) {
            std::vector<Node> found;
            suffix_to_prefix(graph, range, [&](Node node) {
                found.push_back(node);
                return found.size() < limit;
            });
            EXPECT_EQ(found, (std::vector<Node>(all.begin(), all.begin() + limit)));
        }
    }
}

TEST(SuffixMatching, SymbolPolicyComposesWithEarlyStop) {
    DBGSuccinct graph(5);
    for (const auto &sequence : {"ACGCA", "ACGCT", "AGGCA", "AGGCT", "ATGCA", "ACACA"})
        graph.add_sequence(sequence);
    graph.mask_dummy_kmers(1, false);
    const auto range = prefix_range(graph, "A");

    // Permit C/G, then G, then C after the initial A. The last edge symbol
    // remains unrestricted, as in the existing suffix-to-prefix enumeration.
    auto symbols = [&](const boss::BOSS &boss, const auto &next, const auto &emit) {
        switch (std::get<2>(next)) {
            case 2: emit(boss.encode('C')); emit(boss.encode('G')); break;
            case 3: emit(boss.encode('G')); break;
            case 4: emit(boss.encode('C')); break;
            default: ADD_FAILURE() << "Unexpected range depth";
        }
    };
    std::vector<Node> all;
    suffix_to_prefix(graph, range, [&](Node node) {
        all.push_back(node);
        return true;
    }, symbols);
    auto expected = matching_nodes(graph, "ACGC");
    auto more = matching_nodes(graph, "AGGC");
    expected.insert(expected.end(), more.begin(), more.end());
    std::sort(expected.begin(), expected.end());
    auto sorted = all;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_EQ(sorted, expected);
    ASSERT_GT(all.size(), 1u);

    std::vector<Node> first;
    suffix_to_prefix(graph, range, [&](Node node) {
        first.push_back(node);
        return false;
    }, symbols);
    EXPECT_EQ(first, (std::vector<Node>{all.front()}));

    suffix_to_prefix(graph, range, [&](Node) {
        ADD_FAILURE() << "An empty symbol policy must prune the range";
        return true;
    }, [](const auto &, const auto &, const auto &) {});
}

} // namespace
