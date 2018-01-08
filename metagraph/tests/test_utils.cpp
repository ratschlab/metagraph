#include <stdio.h>

#include "gtest/gtest.h"

#include "utils.hpp"


TEST(get_filetype, VCF) {
    EXPECT_EQ("VCF", utils::get_filetype("file.VCF"));
    EXPECT_EQ("VCF", utils::get_filetype("file.vcf"));
    EXPECT_EQ("VCF", utils::get_filetype("file.VCF.gz"));
    EXPECT_EQ("VCF", utils::get_filetype("file.vcf.gz"));
}

TEST(get_filetype, FASTA) {
    EXPECT_EQ("FASTA", utils::get_filetype("file.FASTA"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.fasta"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.FASTA.gz"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.fasta.gz"));
}

TEST(get_filetype, FASTQ) {
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fq"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.FQ"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fastq"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fq.gz"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.FQ.gz"));
}

TEST(get_filetype, IncorrectFiletype) {
    EXPECT_EQ("", utils::get_filetype("fq"));
    EXPECT_EQ("", utils::get_filetype("fasta"));
    EXPECT_EQ("", utils::get_filetype("bz2"));
    EXPECT_EQ("", utils::get_filetype("file.gz"));
    EXPECT_EQ("", utils::get_filetype("fasta.gz"));
}

bool colexicographically_greater(const std::string &s1, const std::string &s2) {
    return utils::colexicographically_greater(s1, s2);
}


TEST(colexicographically_greater, basics) {
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAA"));

    EXPECT_TRUE(colexicographically_greater("BAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("AAAB", "AAAA"));

    EXPECT_FALSE(colexicographically_greater("AAAA", "BAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAB"));
}

TEST(colexicographically_greater, different_lengths) {
    EXPECT_TRUE(colexicographically_greater("AAAAA", "AAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAAA", "BBBB"));

    EXPECT_TRUE(colexicographically_greater("AAAAB", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("ABAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("BAAAA", "AAAA"));
}


TEST(seq_equal, basics) {
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("ABAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 1));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 2));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 3));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 4));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 5));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 6));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 7));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 100));
}
