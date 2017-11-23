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
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fq.gz"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.FQ.gz"));
}
