#include "gtest/gtest.h"

#include "seq_io/vcf_parser.hpp"

#include <set>


namespace {

using namespace mtg;
using namespace mtg::seq_io;

const std::string test_data_dir = "../tests/data/";
const std::string ref_file = test_data_dir + "test_vcfparse.fa";
const std::string vcf_file1 = test_data_dir + "test_vcfparse_1.vcf";
const std::string vcf_file1_bad = test_data_dir + "test_vcfparse_1.bad.vcf.gz";
const std::string vcf_file1_good = test_data_dir + "test_vcfparse_1.good.vcf.gz";
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

TEST(VCFParse, TestKmerNoAnnot) {
    VCFParser vcf(ref_file, vcf_file1, 3);
    std::multiset<std::string> ref { "TGCGCGC" }, obs;
    vcf.call_sequences([&](auto&& seq) { obs.emplace(std::move(seq)); });
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmer) {
    VCFParser vcf(ref_file, vcf_file1, 3);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "TGCGCGC", { "test", "A", "B", "C" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmerEdge) {
    VCFParser vcf(ref_file, vcf_file1, 4);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "ATGCGCGCG", { "test", "A", "B", "C" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmerOverLeftEdge) {
    VCFParser vcf(ref_file, vcf_file1, 5);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "ATGCGCGCGC", { "test", "A", "B", "C" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmerOverBothEdges) {
    VCFParser vcf(ref_file, vcf_file1, 16);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "ATGCGCGCGCGCTCTCGCGCA", { "test", "A", "B", "C" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmerInfoCopyNumber) {
    VCFParser vcf(ref_file, vcf_file2, 3);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "TGCCCGC", { "test", "AC_AMR" } },
        { "TGCTTCGC", { "test" } },
        { "TGCTTTTCGC", { "test" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

TEST(VCFParse, TestKmerGZIP) {
    EXPECT_THROW(VCFParser(ref_file, vcf_file1_bad, 16), std::runtime_error);
}

TEST(VCFParse, TestKmerBGZIP) {
    VCFParser vcf(ref_file, vcf_file1_good, 16);
    std::multiset<std::pair<std::string, std::multiset<std::string>>> ref {
        { "ATGCGCGCGCGCTCTCGCGCA", { "test", "A", "B", "C" } }
    }, obs;

    vcf.call_annotated_sequences(
        [&](auto&& seq, const auto &annotation) {
            obs.emplace(std::move(seq),
                        std::multiset<std::string>(annotation.begin(), annotation.end()));
        },
        annots
    );
    EXPECT_EQ(ref, obs);
}

} // namespace
