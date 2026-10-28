#include <fstream>
#include <stdexcept>
#include <string>

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

} // namespace
