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
#include "graph/representation/hash/sshash_alphabets.hpp"


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

TEST(SSHashAlphabet, DNA5EncodingAndReverseComplement) {
    auto check = [](auto word) {
        using Kmer = decltype(word);
        for (size_t k = 1; k <= Kmer::max_k; ++k) {
            std::string sequence(k, 'N');
            for (size_t i = 0; i < k; ++i)
                sequence[i] = "ACGTN"[i % 5];
            auto packed = sshash::util::string_to_uint_kmer<Kmer>(sequence.data(), k);
            auto expected = sequence;
            reverse_complement(expected.begin(), expected.end());
            packed.reverse_complement_inplace(k);
            EXPECT_EQ(expected, sshash::util::uint_kmer_to_string(packed, k));
            if (k <= Kmer::max_m) {
                const auto original = sshash::util::string_to_uint_kmer<Kmer>(sequence.data(), k);
                EXPECT_EQ(uint64_t(packed), Kmer::reverse_complement_mmer(uint64_t(original), k));
            }
            packed.reverse_complement_inplace(k);
            EXPECT_EQ(sequence, sshash::util::uint_kmer_to_string(packed, k));
        }
    };
    check(detail::DNA5Kmer<uint64_t>{});
    check(detail::DNA5Kmer<__uint128_t>{});
    check(detail::DNA5Kmer<sdsl::uint256_t>{});
    for (unsigned c = 0; c < 256; ++c) {
        EXPECT_EQ(std::string("ACGTUNacgtun").find(char(c)) != std::string::npos,
                  detail::DNA5Kmer<uint64_t>::is_valid(char(c)));
    }
}



TEST_F(DBGSSHashInput, ChoosesSmallestPackedWord) {
    constexpr auto k64 = DBGSSHash::kmer_t<uint64_t>::max_k;
    constexpr auto k128 = DBGSSHash::kmer_t<__uint128_t>::max_k;
    constexpr auto k256 = DBGSSHash::kmer_t<sdsl::uint256_t>::max_k;
    EXPECT_EQ(0u, DBGSSHash(k64).data().index());
    EXPECT_EQ(1u, DBGSSHash(k64 + 1).data().index());
    EXPECT_EQ(1u, DBGSSHash(k128).data().index());
    EXPECT_EQ(2u, DBGSSHash(k128 + 1).data().index());
    EXPECT_EQ(2u, DBGSSHash(k256).data().index());
    EXPECT_THROW(DBGSSHash(k256 + 1), std::length_error);
}

#if _PROTEIN_GRAPH
const std::string protein_sequence = "MKWVTFISRLDPQEAHNYCG";

void check_protein_nodes(const DBGSSHash &graph, const std::string &sequence,
                         const std::vector<uint64_t> &expected_nodes) {
    const size_t k = graph.get_k();
    EXPECT_EQ(sequence.size() - k + 1, graph.num_nodes());
    std::vector<uint64_t> nodes;
    graph.map_to_nodes_sequentially(sequence, [&](auto node) { nodes.push_back(node); });
    ASSERT_EQ(expected_nodes, nodes);
    for (size_t i = 0; i < nodes.size(); ++i) {
        ASSERT_NE(0u, nodes[i]);
        EXPECT_EQ(sequence.substr(i, k), graph.get_node_sequence(nodes[i]));
        EXPECT_EQ(nodes[i], graph.kmer_to_node(sequence.substr(i, k)));
        if (i + 1 < nodes.size()) {
            EXPECT_EQ(nodes[i + 1], graph.traverse(nodes[i], sequence[i + k]));
            EXPECT_EQ(nodes[i], graph.traverse_back(nodes[i + 1], sequence[i]));
        }
    }
    EXPECT_EQ(0u, graph.kmer_to_node(std::string(k, 'A')));
}

TEST_F(DBGSSHashInput, ProteinPackedWordConstructionAndRoundTrip) {
    const std::string sequence = protein_sequence + "STVQKMDF";
    for (size_t k : {11, 12, 13, 25, 26}) {
        SCOPED_TRACE(k);
        write_input(">a\n" + sequence + "\n");
        DBGSSHash graph(test_input, k);
        const size_t expected_width = k <= 12 ? 0 : k <= 25 ? 1 : 2;
        EXPECT_EQ(expected_width, graph.data().index());
        std::vector<uint64_t> nodes;
        graph.map_to_nodes_sequentially(sequence, [&](auto node) { nodes.push_back(node); });
        check_protein_nodes(graph, sequence, nodes);

        std::stringstream bytes;
        graph.serialize(bytes);
        DBGSSHash restored(3);
        ASSERT_TRUE(restored.load(bytes));
        EXPECT_EQ(k, restored.get_k());
        EXPECT_EQ(expected_width, restored.data().index());
        check_protein_nodes(restored, sequence, nodes);

        graph.serialize(test_input + ".graph");
        ASSERT_TRUE(restored.load(test_input + ".graph"));
        EXPECT_EQ(expected_width, restored.data().index());
        check_protein_nodes(restored, sequence, nodes);

        auto mapped = utils::open_ifstream(test_input + ".graph.sshashdbg", true);
        ASSERT_TRUE(restored.load(*mapped));
        EXPECT_EQ(expected_width, restored.data().index());
        check_protein_nodes(restored, sequence, nodes);
    }
}

TEST_F(DBGSSHashInput, ProteinMinimizerFitsPackedWord) {
    write_input(">a\n" + protein_sequence + "STVQKMDF\n");
    DBGSSHash graph(test_input, 25, DeBruijnGraph::BASIC, uint64_t(-1));
    EXPECT_EQ(11u, std::visit([](const auto &dict) { return dict.m(); }, graph.data()));
}

TEST_F(DBGSSHashInput, Protein12EmptyGraphRoundTrip) {
    DBGSSHash graph(12);
    std::stringstream bytes;
    graph.serialize(bytes);
    DBGSSHash restored(3);
    ASSERT_TRUE(restored.load(bytes));
    EXPECT_EQ(12u, restored.get_k());
    EXPECT_EQ(0u, restored.num_nodes());
    EXPECT_EQ(0u, restored.data().index());
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

TEST_F(DBGSSHashInput, DNA5AllThreeMers) {
    std::set<std::string> words;
    for (char a : std::string("ACGTN"))
        for (char b : std::string("ACGTN"))
            for (char c : std::string("ACGTN")) words.insert(std::string{a, b, c});
    for (auto mode : {DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY}) {
        std::string input;
        for (const auto &word : words) {
            auto rc = word;
            reverse_complement(rc.begin(), rc.end());
            if (mode == DeBruijnGraph::BASIC || word <= rc)
                input += ">" + word + "\n" + word + "\n";
        }
        write_input(input);
        DBGSSHash graph(test_input, 3, mode);
        auto check = [&](const DBGSSHash &g) {
            EXPECT_EQ(words.size(), g.num_nodes());
            std::set<std::string> observed;
            size_t self_complements = 0;
            g.call_nodes([&](auto node) {
                ASSERT_TRUE(g.in_graph(node));
                auto word = g.get_node_sequence(node);
                EXPECT_TRUE(observed.insert(word).second);
                EXPECT_EQ(node, g.kmer_to_node(word));
                for (char c : std::string("ACGTN")) {
                    EXPECT_EQ(g.kmer_to_node(word.substr(1) + c), g.traverse(node, c));
                    EXPECT_EQ(g.kmer_to_node(c + word.substr(0, 2)), g.traverse_back(node, c));
                }
                if (mode != DeBruijnGraph::BASIC) {
                    auto rc = word;
                    reverse_complement(rc.begin(), rc.end());
                    EXPECT_EQ(g.kmer_to_node(rc), g.reverse_complement(node));
                    if (rc == word) ++self_complements;
                }
            });
            EXPECT_EQ(words, observed);
            if (mode != DeBruijnGraph::BASIC) {
                EXPECT_EQ(5u, self_complements);
            }
            size_t valid_ids = 0;
            for (size_t node = 1; node <= g.max_index(); ++node) valid_ids += g.in_graph(node);
            EXPECT_EQ(g.num_nodes(), valid_ids);
        };
        check(graph);
        std::stringstream bytes;
        graph.serialize(bytes);
        DBGSSHash restored(3);
        ASSERT_TRUE(restored.load(bytes));
        check(restored);
        graph.serialize(test_input + ".graph");
        auto mapped = utils::open_ifstream(test_input + ".graph.sshashdbg", true);
        ASSERT_TRUE(restored.load(*mapped));
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
