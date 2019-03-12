//
//  tests.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef tests_h
#define tests_h

#include "utilities.hpp"
#include "compressed_reads.hpp"
#include "samplers.hpp"

const int test_seed = 3424;

TEST(SamplerTest,SampleNoRandom) {
    auto sampler = Sampler("AAAAAAAAA",test_seed);
    ASSERT_EQ(sampler.sample(2),"AA");
}

TEST(SamplerTest,SampleNormal) {
    auto sampler = Sampler("ADFAGADFDS",test_seed);
    ASSERT_EQ(sampler.sample(4),"ADFD");
}
TEST(SamplerTest,SampleCoverage) {
    auto sequence = "ADFAGADFDS"s;
    auto sampler = Sampler(sequence,test_seed);
    auto reads = sampler.sample_coverage(sequence.length()/2, 1);
    ASSERT_EQ(reads.size(), 2);
}

// Depends on large file -> tested and works
//TEST(CompressingReads,GetChromosomeWorks) {
//    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
//    EXPECT_EQ(chromosome.length(), 133'797'422);
//    EXPECT_EQ(chromosome.substr(0,10),"NNNNNNNNNN");
//}

TEST(CompressedReads,IdentityTest1) {
    set<string> reads = {"ATGCGATCGATATGCGAGA",
                         "ATGCGATCGAGACTACGAG",
                         "GTACGATAGACATGACGAG",
                         "ACTGACGAGACACAGATGC"};
    auto compressed_reads = CompressedReads(vector<string>(all(reads)),5);
    auto decompressed_reads = compressed_reads.get_reads();
    set<string> decompressed_read_set = set<string>(all(decompressed_reads));
    ASSERT_EQ(reads, decompressed_read_set);
}


#endif /* tests_h */
