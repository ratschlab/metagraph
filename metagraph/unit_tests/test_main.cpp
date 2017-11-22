#include <stdio.h>

#include "gtest/gtest.h"
#include "dbg_succinct_libmaus.hpp"
#include "construct.hpp"
#include "config.hpp"
#include "kmer.hpp"

int Sum(int a, int b) {
  return a + b;
}

TEST(SumTest, HandlesZeroInput) {
  // Yoda condition.
  ASSERT_EQ(3, Sum(1, 2));
  ASSERT_EQ(4, Sum(2, 2));
  printf("--- After ASSERT_EQ\n");
  EXPECT_EQ(3, Sum(1, 2));
  printf("--- After EXPECT_EQ\n");
  ASSERT_EQ(4, Sum(2, 2)) << "Something is wrong, 2+2 isn't equal to 4 !";
}

TEST(SumTest, HandlesNegativeInput) {
  ASSERT_EQ(-3, Sum(-1, -2));
}

TEST(KmerEncodeTest, Invertible) {
    int argc = 3;
    char const *argv[argc] = {"metagengraph","build","test.fa"};
    Config *config = new Config(argc, argv);
    DBG_succ *graph = new DBG_succ(config->k, config); 

    std::string test_kmer = "ATGG";
    ui256 kmer = kmer_boost::stokmer(test_kmer.c_str(), test_kmer.length(), construct::nt_lookup);
    const char *kmer_s = kmer_boost::kmertos(kmer, graph->alphabet, graph->alph_size);
    ASSERT_EQ(*kmer_s, test_kmer[test_kmer.length()-1]);
    ASSERT_EQ(std::string(kmer_s+1), test_kmer.substr(test_kmer.length()-1));
}
