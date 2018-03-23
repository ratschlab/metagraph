#include <stdio.h>

#include "gtest/gtest.h"

#include "vcf_parser.hpp"

const std::string test_data_dir = "../tests/data/";
const std::string ref_file = test_data_dir + "test_vcfparse.fa";
const std::string vcf_file1 = test_data_dir + "test_vcfparse_1.vcf";
const std::string vcf_file2 = test_data_dir + "test_vcfparse_2.vcf";

const std::vector<std::string> annots = {
    "AC_AFR",
    "AC_AMR",
    "AC_ASJ",
    "AC_EAS",
    "AC_FIN",
    "AC_NFE",
    "AC_OTH",
    "AC_SAS"
};

TEST(VcfParse, LoadVCF) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file.c_str(), vcf_file1.c_str(), 3));
}

TEST(VcfParse, TestKmer) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file, vcf_file1, 3));
    std::vector<std::string> annot;
    ASSERT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("TGCGCGC", vcf.seq);
    EXPECT_EQ(4llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    EXPECT_EQ("A", annot[1]);
    EXPECT_EQ("B", annot[2]);
    EXPECT_EQ("C", annot[3]);
    EXPECT_FALSE(vcf.get_seq(annots, &annot));
}

TEST(VcfParse, TestKmerEdge) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file, vcf_file1, 4));
    std::vector<std::string> annot;
    ASSERT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("ATGCGCGCG", vcf.seq);
    EXPECT_EQ(4llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    EXPECT_EQ("A", annot[1]);
    EXPECT_EQ("B", annot[2]);
    EXPECT_EQ("C", annot[3]);
    EXPECT_FALSE(vcf.get_seq(annots, &annot));
}

TEST(VcfParse, TestKmerOverLeftEdge) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file, vcf_file1, 5));
    std::vector<std::string> annot;
    ASSERT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("ATGCGCGCGC", vcf.seq);
    EXPECT_EQ(4llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    EXPECT_EQ("A", annot[1]);
    EXPECT_EQ("B", annot[2]);
    EXPECT_EQ("C", annot[3]);
    EXPECT_FALSE(vcf.get_seq(annots, &annot));
}

TEST(VcfParse, TestKmerOverBothEdges) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file, vcf_file1, 16));
    std::vector<std::string> annot;
    ASSERT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("ATGCGCGCGCGCTCTCGCGCA", vcf.seq);
    EXPECT_EQ(4llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    EXPECT_EQ("A", annot[1]);
    EXPECT_EQ("B", annot[2]);
    EXPECT_EQ("C", annot[3]);
    EXPECT_FALSE(vcf.get_seq(annots, &annot));
}

TEST(VcfParse, TestKmerInfoCopyNumber) {
    vcf_parser vcf;
    ASSERT_TRUE(vcf.init(ref_file, vcf_file2, 3));
    std::vector<std::string> annot;
    ASSERT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("TGCCCGC", vcf.seq);
    EXPECT_EQ(2llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    EXPECT_EQ("AC_AMR", annot[1]);
    annot.clear();

    EXPECT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("TGCTTCGC", vcf.seq);
    EXPECT_EQ(1llu, annot.size());
    EXPECT_EQ("test", annot[0]);
    annot.clear();

    EXPECT_TRUE(vcf.get_seq(annots, &annot));
    EXPECT_EQ("TGCTTTTCGC", vcf.seq);
    EXPECT_EQ(1llu, annot.size());
    EXPECT_EQ("test", annot[0]);

    EXPECT_FALSE(vcf.get_seq(annots, &annot));
}

