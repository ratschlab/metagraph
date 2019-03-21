//
//  tests.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef tests_h
#define tests_h

#include <gtest/gtest.h>
#include "utilities.hpp"
#include "path_database_list_of_bifurcation_choices.hpp"
#include "samplers.hpp"
#include "path_database_baseline.hpp"
//#include "path_database.hpp"


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


vector<string> reads_for_testing_short = {"ATGCGATCGATATGCGAGA",
                                          "ATGCGATCGAGACTACGAG",
                                          "GTACGATAGACATGACGAG",
                                          "ACTGACGAGACACAGATGC"};

template <typename T>
void check_compression_decompression(PathDatabase<T>& db,vector<string>& reads) {
    auto handles = db.encode(reads);
    vector<string> decompressed_reads;
    decompressed_reads.reserve(handles.size());
    for(auto& handle : handles) {
        decompressed_reads.push_back(db.decode(handle));
    }
    ASSERT_EQ(reads,decompressed_reads);
}

TEST(PathDatabase,IdentityTest1) {
    auto db = CompressedReads(reads_for_testing_short,5);
    check_compression_decompression(db,reads_for_testing_short);
}

TEST(PathDatabase,IdentityTest2) {
    auto database = PathDatabaseBaseline(reads_for_testing_short,5);
    check_compression_decompression(database,reads_for_testing_short);
}

#if defined(__linux__) || true

TEST(PathDatabase,LongTest) {
    string reads_filename = "/cluster/home/studenyj/";

}

#endif

#endif /* tests_h */
