#include <stdio.h>

#include "gtest/gtest.h"

#include "vcf_parser.hpp"

const std::string test_data_dir = "../unit_tests/data/";
const std::string ref_file = test_data_dir + "test_vcfparse.fa";
const std::string vcf_file1 = test_data_dir + "test_vcfparse_1.vcf";
const std::string vcf_file2 = test_data_dir + "test_vcfparse_2.vcf";
const char *annots[] = {};

TEST(VcfParse, LoadVCF) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 3));
}

TEST(VcfParse, TestKmer) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 3));
    std::string annot = vcf.vcf_get_seq(annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf.seq.s, "TGCGCGC"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
}

TEST(VcfParse, TestKmerEdge) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 4));
    std::string annot = vcf.vcf_get_seq(annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf.seq.s, "ATGCGCGCG"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
}

TEST(VcfParse, TestKmerOverLeftEdge) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 5));
    std::string annot = vcf.vcf_get_seq(annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf.seq.s, "ATGCGCGCGC"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
}

TEST(VcfParse, TestKmerOverBothEdges) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 16));
    std::string annot = vcf.vcf_get_seq(annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf.seq.s, "ATGCGCGCGCGCTCTCGCGCA"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
}

TEST(VcfParse, TestKmerInfo) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file2.c_str(), 3));
    std::string annot = vcf.vcf_get_seq(annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf.seq.s, "TGCCCGC"), 0);
    EXPECT_EQ(annot, "test");
}

