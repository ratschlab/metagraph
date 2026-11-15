#include <fstream>
#include <random>
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
    for (auto node : map_to_nodes_sequentially(graph, sequence)) {
        nodes.push_back(node);
    }
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
        for (auto node : map_to_nodes_sequentially(graph, sequence)) {
            nodes.push_back(node);
        }
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
#if ! _PROTEIN_GRAPH
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
    // m = ceil(log_4(num_chars)) + 1, made odd, which steps from 7 to 9 at 4^6 = 4096
    write_random_input({ 2000, 2095 });
    EXPECT_EQ(7u, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC)));

    write_random_input({ 2000, 2097 });
    EXPECT_EQ(9u, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC)));

    // a given number of characters overrides the count: ceil(log_4(1e9)) + 1 = 16, made odd
    EXPECT_EQ(17u, minimizer_length(DBGSSHash(test_input, 21, DeBruijnGraph::BASIC,
                                              1'000'000'000)));
}
#endif

} // namespace
