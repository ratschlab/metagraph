#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/resource.h>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;

const std::string test_input = test_dump_dir() + "/sshash_input.fa";

// Inputs SSHash cannot parse used to send its build into an endless loop writing
// temp files. Cap the size of written files for the duration of a test so that a
// regression kills the test binary with SIGXFSZ instead of filling the disk.
class DBGSSHashInput : public ::testing::Test {
  protected:
    void SetUp() override {
        ASSERT_EQ(0, getrlimit(RLIMIT_FSIZE, &old_limit_));
        rlimit limit = old_limit_;
        limit.rlim_cur = 64 << 20;
        if (old_limit_.rlim_cur != RLIM_INFINITY && old_limit_.rlim_cur < limit.rlim_cur)
            limit.rlim_cur = old_limit_.rlim_cur;
        ASSERT_EQ(0, setrlimit(RLIMIT_FSIZE, &limit));
    }

    void TearDown() override {
        ASSERT_EQ(0, setrlimit(RLIMIT_FSIZE, &old_limit_));
    }

  private:
    rlimit old_limit_;
};

#if _DNA_CASE_SENSITIVE_GRAPH
TEST_F(DBGSSHashInput, RejectsCaseSensitiveAlphabet) {
    EXPECT_THROW(DBGSSHash(3), std::domain_error);
}
#else
void write_input(const std::string &contents) {
    std::ofstream out(test_input);
    out << contents;
}

void expect_rejected(const std::string &contents,
                     const std::string &reason,
                     size_t k = 21,
                     DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    write_input(contents);
    try {
        DBGSSHash graph(test_input, k, mode);
        ADD_FAILURE() << "input accepted";
    } catch (const std::exception &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find(reason)) << e.what();
    }
}

const std::string seq_30 = "ACGTTGCAAGGCTTACCATGGTACCAGTTA";

#if _DNA5_GRAPH
TEST_F(DBGSSHashInput, ChoosesSmallestPackedWord) {
    constexpr auto k64 = DBGSSHash::kmer_t<uint64_t>::max_k;
    constexpr auto k128 = DBGSSHash::kmer_t<__uint128_t>::max_k;
    EXPECT_EQ(0u, DBGSSHash(k64).data().index());
    EXPECT_EQ(1u, DBGSSHash(k64 + 1).data().index());
    EXPECT_EQ(1u, DBGSSHash(k128).data().index());
    EXPECT_EQ(2u, DBGSSHash(k128 + 1).data().index());
}

#endif

TEST_F(DBGSSHashInput, RejectsSequenceShorterThanK) {
    expect_rejected(">a\n" + seq_30 + "\n>b\nACGTA\n",
                    "sequence 2 has length 5, shorter than k = 21");
}

TEST_F(DBGSSHashInput, AcceptsWrappedFasta) {
    // a wrapped record whose tail is shorter than k, like the
    // records in tests/data/transcripts_100.fa
    write_input(">a\n" + seq_30 + "\nTTGACCATGCAGTCAGGTACCATGACTTGA\nAC\n");
    DBGSSHash graph(test_input, 21, DeBruijnGraph::BASIC);
    EXPECT_EQ(62u - 21u + 1u, graph.num_nodes());
}

TEST_F(DBGSSHashInput, RejectsFastq) {
    expect_rejected("@a\n" + seq_30 + "\n+\n" + std::string(30, 'I') + "\n",
                    "expected a FASTA header starting with '>' before sequence 1");
}

TEST_F(DBGSSHashInput, RejectsMissingTrailingNewline) {
    expect_rejected(">a\n" + seq_30, "sequence 1 is not terminated by a newline");
}

TEST_F(DBGSSHashInput, RejectsHeaderWithoutSequence) {
    expect_rejected(">a\n" + seq_30 + "\n>b\n", "record 2 has a header but no sequence");
}

TEST_F(DBGSSHashInput, RejectsInvalidCharacter) {
#if _DNA_GRAPH
    const char invalid = 'N';
#else
    const char invalid = '#';
#endif
    std::string seq = seq_30;
    seq[12] = invalid;
    expect_rejected(">a\n" + seq + "\n", std::string("sequence 1 has character '")
                    + invalid + "' at position 12, outside the alphabet");
}

TEST_F(DBGSSHashInput, RejectsDuplicateKmers) {
    expect_rejected(">a\n" + seq_30 + "\n>b\n" + seq_30 + "\n", "contains duplicate k-mers");
}

#if ! _PROTEIN_GRAPH
TEST_F(DBGSSHashInput, RejectsReverseComplementDuplicatesCanonical) {
    std::string rc = seq_30;
    reverse_complement(rc.begin(), rc.end());
    expect_rejected(">a\n" + seq_30 + "\n>b\n" + rc + "\n", "contains duplicate k-mers",
                    21, DeBruijnGraph::CANONICAL);
}
#endif

TEST_F(DBGSSHashInput, AcceptsOneSequencePerLine) {
    write_input(">a\n" + seq_30 + "\n>b\nacgtacgtacgtacgtacgtac\n");
    DBGSSHash graph(test_input, 21, DeBruijnGraph::BASIC);
    EXPECT_EQ(10u + 2u, graph.num_nodes());
}

#if ! _PROTEIN_GRAPH
size_t minimizer_length(const DBGSSHash &graph) {
    return std::visit([](const auto &dict) -> size_t { return dict.m(); }, graph.data());
}

// random records of the given lengths
void write_random_input(const std::vector<size_t> &lengths) {
    std::mt19937 gen(42);
    std::uniform_int_distribution<int> base(0, 3);
    std::string input;
    for (size_t i = 0; i < lengths.size(); ++i) {
        input += ">" + std::to_string(i) + "\n";
        for (size_t j = 0; j < lengths[i]; ++j) {
            input += "ACGT"[base(gen)];
        }
        input += "\n";
    }
    write_input(input);
}

TEST_F(DBGSSHashInput, ChoosesMinimizerLengthFromInputSize) {
    // m = ceil(log_sigma(num_chars)) + 1, made odd: the step is at sigma^6.
#if _DNA5_GRAPH
    constexpr size_t threshold = 15'625;
    constexpr size_t large_m = 15;
#else
    constexpr size_t threshold = 4096;
    constexpr size_t large_m = 17;
#endif
    write_random_input({ 2000, threshold - 2001 });
    EXPECT_EQ(7u, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC)));

    write_random_input({ 2000, threshold - 1999 });
    EXPECT_EQ(9u, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC)));

    // A supplied input size overrides the count for the selected alphabet.
    EXPECT_EQ(large_m, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC,
                                              1'000'000'000)));
}
#endif


#if _DNA5_GRAPH
TEST_F(DBGSSHashInput, DNA5LookupTraversalAndSerialization) {
    for (auto mode : {DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL}) {
        write_input(">a\nACGNT\n>palindrome\nANT\n");
        DBGSSHash graph(test_input, 3, mode);
        std::set<std::string> expected = {"ACG", "CGN", "GNT", "ANT"};
        if (mode == DeBruijnGraph::CANONICAL) {
            for (auto word : std::vector<std::string>(expected.begin(), expected.end())) {
                reverse_complement(word.begin(), word.end());
                expected.insert(word);
            }
        }
        auto check = [&](const DBGSSHash &g) {
            EXPECT_EQ(expected.size(), g.num_nodes());
            std::set<std::string> observed;
            g.call_nodes([&](auto node) {
                EXPECT_TRUE(g.in_graph(node));
                auto word = g.get_node_sequence(node);
                EXPECT_TRUE(observed.insert(word).second);
                EXPECT_EQ(node, g.kmer_to_node(word));
                std::set<std::string> outgoing;
                g.call_outgoing_kmers(node, [&](auto next, char c) {
                    EXPECT_EQ(next, g.traverse(node, c));
                    outgoing.insert(g.get_node_sequence(next));
                });
                std::set<std::string> wanted;
                for (char c : std::string("ACGTN")) {
                    auto next = word.substr(1) + c;
                    if (expected.count(next)) wanted.insert(next);
                    EXPECT_EQ(bool(expected.count(next)), bool(g.traverse(node, c)));
                    auto prev = c + word.substr(0, 2);
                    EXPECT_EQ(bool(expected.count(prev)), bool(g.traverse_back(node, c)));
                }
                EXPECT_EQ(wanted, outgoing);
                if (mode == DeBruijnGraph::CANONICAL) {
                    auto rc = word;
                    reverse_complement(rc.begin(), rc.end());
                    EXPECT_EQ(rc, g.get_node_sequence(g.reverse_complement(node)));
                    EXPECT_EQ(node, g.reverse_complement(g.reverse_complement(node)));
                }
            });
            EXPECT_EQ(expected, observed);
            if (mode == DeBruijnGraph::CANONICAL) {
                auto palindrome = g.kmer_to_node("ANT");
                EXPECT_EQ(palindrome, g.reverse_complement(palindrome));
            }
        };
        check(graph);
        std::stringstream bytes;
        graph.serialize(bytes);
        DBGSSHash restored(3, mode);
        ASSERT_TRUE(restored.load(bytes));
        check(restored);
        graph.serialize(test_input + ".graph");
        ASSERT_TRUE(restored.load(test_input + ".graph"));
        check(restored);
    }
}

TEST_F(DBGSSHashInput, DNA5CrossesPackedWordBoundaries) {
    const std::string sequence = "ACGNTGCAAGGNTTACCATGNTACCAGTTANCGTACGNNACGTNACCGTACGNT";
    for (size_t k : {21, 22, 23, 42, 43}) {
        write_input(">a\n" + sequence + "\n");
        DBGSSHash graph(test_input, k, DeBruijnGraph::BASIC);
        std::vector<uint64_t> nodes;
        graph.map_to_nodes_sequentially(sequence, [&](auto node) { nodes.push_back(node); });
        ASSERT_EQ(sequence.size() - k + 1, nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            ASSERT_NE(0u, nodes[i]);
            EXPECT_EQ(sequence.substr(i, k), graph.get_node_sequence(nodes[i]));
        }
    }
}

TEST_F(DBGSSHashInput, DNA5EmptyGraphRoundTrip) {
    for (auto mode : {DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL}) {
        DBGSSHash graph(3, mode);
        std::stringstream bytes;
        graph.serialize(bytes);
        DBGSSHash restored(3);
        ASSERT_TRUE(restored.load(bytes));
        EXPECT_EQ(0u, restored.num_nodes());
        EXPECT_EQ(0u, restored.max_index());
    }
}

TEST_F(DBGSSHashInput, DNA5RejectsUntaggedTwoBitFiles) {
    std::stringstream bytes;
    essentials::generic_saver saver(bytes);
    uint64_t nodes = 0, k = 3;
    auto mode = DeBruijnGraph::BASIC;
    saver.visit(nodes);
    saver.visit(k);
    saver.visit(mode);
    DBGSSHash graph(3);
    EXPECT_THROW(graph.load(bytes), std::runtime_error);
}
#endif

#endif // _DNA_CASE_SENSITIVE_GRAPH

} // namespace
