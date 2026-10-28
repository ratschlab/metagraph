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
    } catch (const std::invalid_argument &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find(reason)) << e.what();
    }
}

const std::string seq_30 = "ACGTTGCAAGGCTTACCATGGTACCAGTTA";

TEST_F(DBGSSHashInput, RejectsSequenceShorterThanK) {
    expect_rejected(">a\n" + seq_30 + "\n>b\nACGTA\n", "line 4 has a sequence of length 5");
}

TEST_F(DBGSSHashInput, RejectsWrappedFasta) {
    // a wrapped record whose tail is shorter than k, like the
    // records in tests/data/transcripts_100.fa
    expect_rejected(">a\n" + seq_30 + "\n" + seq_30 + "\nAC\n", "line 3 is not a FASTA header");
}

TEST_F(DBGSSHashInput, RejectsFastq) {
    expect_rejected("@a\n" + seq_30 + "\n+\n" + std::string(30, 'I') + "\n",
                    "line 1 is not a FASTA header");
}

TEST_F(DBGSSHashInput, RejectsInvalidCharacter) {
    std::string seq = seq_30;
    seq[10] = '*';
    expect_rejected(">a\n" + seq + "\n", "line 2 has an invalid character");
}

TEST_F(DBGSSHashInput, RejectsMissingTrailingNewline) {
    expect_rejected(">a\n" + seq_30, "line 2 is not terminated by a newline");
}

TEST_F(DBGSSHashInput, RejectsDuplicateKmers) {
    expect_rejected(">a\n" + seq_30 + "\n>b\n" + seq_30 + "\n", "duplicate k-mers");
}

#if ! _PROTEIN_GRAPH
TEST_F(DBGSSHashInput, RejectsReverseComplementDuplicatesCanonical) {
    std::string rc = seq_30;
    reverse_complement(rc.begin(), rc.end());
    expect_rejected(">a\n" + seq_30 + "\n>b\n" + rc + "\n", "duplicate k-mers",
                    21, DeBruijnGraph::CANONICAL);
}
#endif

TEST_F(DBGSSHashInput, AcceptsOneSequencePerLine) {
    write_input(">a\n" + seq_30 + "\n>b\nacgtacgtacgtacgtacgtac\n");
    DBGSSHash graph(test_input, 21, DeBruijnGraph::BASIC);
    EXPECT_EQ(10u + 2u, graph.num_nodes());
}

} // namespace
