#include <stdio.h>

#include "gtest/gtest.h"

#include "vcfparse.h"

const std::string working_dir = "../unit_tests";
const std::string ref_file = working_dir + "/test_vcfparse.fa";
const std::string vcf_file1 = working_dir + "/test_vcfparse_1.vcf";
const std::string vcf_file2 = working_dir + "/test_vcfparse_2.vcf";
const char *annots[] = {};

TEST(VcfParse, LoadVCF) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file1.c_str(), 3);
    ASSERT_EQ(!vcf, false);
    vcf_destroy(vcf);
}

TEST(VcfParse, TestKmer) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file1.c_str(), 3);
    ASSERT_EQ(!vcf, false);
    std::string annot = vcf_get_seq(vcf, annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf->seq.s, "TGCGCGC"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
    vcf_destroy(vcf);
}

TEST(VcfParse, TestKmerEdge) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file1.c_str(), 4);
    ASSERT_EQ(!vcf, false);
    std::string annot = vcf_get_seq(vcf, annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf->seq.s, "ATGCGCGCG"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
    vcf_destroy(vcf);
}

TEST(VcfParse, TestKmerOverLeftEdge) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file1.c_str(), 5);
    ASSERT_EQ(!vcf, false);
    std::string annot = vcf_get_seq(vcf, annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf->seq.s, "ATGCGCGCGC"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
    vcf_destroy(vcf);
}

TEST(VcfParse, TestKmerOverBothEdges) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file1.c_str(), 16);
    ASSERT_EQ(!vcf, false);
    std::string annot = vcf_get_seq(vcf, annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf->seq.s, "ATGCGCGCGCGCTCTCGCGCA"), 0);
    EXPECT_EQ(annot, "test:A:B:C");
    vcf_destroy(vcf);
}

TEST(VcfParse, TestKmerInfo) {
    vcfparse *vcf = vcf_init(ref_file.c_str(), vcf_file2.c_str(), 3);
    ASSERT_EQ(!vcf, false);
    std::string annot = vcf_get_seq(vcf, annots, sizeof(annots) / sizeof(char *));
    EXPECT_EQ(strcmp(vcf->seq.s, "TGCCCGC"), 0);
    EXPECT_EQ(annot, "test");
    vcf_destroy(vcf);
}

