#include "gtest/gtest.h"

#include <vector>
#include <string>

#include "common/taxid_mapper.hpp"


namespace {

const std::string test_data_dir = "../tests/data";
const std::string test_accession2taxid = test_data_dir + "/nucl_gb_accession2taxid.head.gz";
const std::string test_nodes = test_data_dir + "/nodes.dmp.head";
const std::string test_dump_basename = test_data_dir + "/dump_test_taxid";


void test_gb_to_taxid(const TaxIDMapper &mapper) {
    std::vector<std::tuple<std::string, std::string, TaxIDMapper::taxid_t>> data{
        { "A00002", "A00002.1", 9913 },
        { "A00003", "A00003.1", 9913 },
        { "X17276", "X17276.1", 9646 },
        { "X60065", "X60065.1", 9913 },
        { "X55027", "X55027.1", 9913 },
        { "Z12029", "Z12029.1", 9915 },
        { "X52700", "X52700.1", 9771 },
        { "X52701", "X52701.1", 9771 },
        { "X52702", "X52702.1", 9771 },
        { "X52706", "X52706.1", 9771 },
        { "X52703", "X52703.1", 9771 },
        { "X52704", "X52704.1", 9771 },
        { "X52705", "X52705.1", 9771 },
        { "X53811", "X53811.1", 9771 },
        { "X53812", "X53812.1", 9771 },
        { "X53813", "X53813.1", 9771 },
        { "X53814", "X53814.1", 9771 },
        { "Z18633", "Z18633.1", 9770 },
        { "Z18632", "Z18632.1", 9770 },
        { "X60496", "X60496.1", 9913 },
        { "X60497", "X60497.1", 9913 }
    };

    for (const auto &triple : data) {
        EXPECT_EQ(std::get<2>(triple), mapper.gb_to_taxid(std::get<0>(triple)));
        EXPECT_EQ(std::get<2>(triple), mapper.gb_to_taxid(std::get<1>(triple)));
        EXPECT_EQ(std::get<2>(triple), mapper.gb_to_taxid("gb|" + std::get<0>(triple) + "|"));
        EXPECT_EQ(std::get<2>(triple), mapper.gb_to_taxid("gb|" + std::get<1>(triple) + "|"));
    }

    EXPECT_EQ(0u, mapper.gb_to_taxid("FAKE"));
}

void test_id_to_parent_and_label(const TaxIDMapper &mapper) {
    std::vector<std::tuple<TaxIDMapper::taxid_t,
                           TaxIDMapper::taxid_t,
                           std::string>> data{
        { 1,      1       , "no rank" },
        { 2,      131567  , "superkingdom" },
        { 6,      335928  , "genus" },
        { 7,      6       , "species" },
        { 9,      32199   , "species" },
        { 10,     1706371 , "genus" },
        { 11,     1707    , "species" },
        { 13,     203488  , "genus" },
        { 14,     13      , "species" },
        { 16,     32011   , "genus" },
        { 17,     16      , "species" },
        { 18,     213421  , "genus" },
        { 19,     18      , "species" },
        { 20,     76892   , "genus" },
        { 21,     20      , "species" },
        { 22,     267890  , "genus" },
        { 339,    338     , "species" },
        { 340,    339     , "no rank" },
        { 190485, 340     , "no rank"}
    };

    for (const auto &triple : data) {
        EXPECT_EQ(std::get<1>(triple), mapper.get_parent(std::get<0>(triple)));
        EXPECT_EQ(std::get<2>(triple), mapper.get_rank_label(std::get<0>(triple)));
    }

    EXPECT_EQ(0u, mapper.get_parent(23));
    EXPECT_EQ(std::string(), mapper.get_rank_label(23));
}

void test_id_to_ancestor_labeled(const TaxIDMapper &mapper) {
    std::vector<std::tuple<TaxIDMapper::taxid_t,
                           TaxIDMapper::taxid_t,
                           std::string>> data{
        { 1,      0      , "" },
        { 2,      2      , "superkingdom" },
        { 6,      6      , "genus" },
        { 7,      7      , "species" },
        { 9,      9      , "species" },
        { 10,     10     , "genus" },
        { 11,     11     , "species" },
        { 13,     13     , "genus" },
        { 14,     14     , "species" },
        { 16,     16     , "genus" },
        { 17,     17     , "species" },
        { 18,     18     , "genus" },
        { 19,     19     , "species" },
        { 20,     20     , "genus" },
        { 21,     21     , "species" },
        { 22,     22     , "genus" },
        { 339,    339    , "species" },
        { 340,    339    , "species" },
        { 190485, 339    , "species" }
    };

    for (const auto &triple : data) {
        std::pair<TaxIDMapper::taxid_t, std::string> ref{ std::get<1>(triple),
                                                          std::get<2>(triple) };
        EXPECT_EQ(ref, mapper.get_ancestor_labeled(std::get<0>(triple)));
    }
}

void test_id_to_ancestor_with_rank_label(const TaxIDMapper &mapper) {
    std::vector<std::tuple<TaxIDMapper::taxid_t,
                           TaxIDMapper::taxid_t,
                           std::string>> data{
        { 1,      1      , "no rank" },
        { 2,      2      , "superkingdom" },
        { 6,      6      , "genus" },
        { 7,      7      , "species" },
        { 7,      6      , "genus" },
        { 9,      9      , "species" },
        { 10,     10     , "genus" },
        { 11,     11     , "species" },
        { 13,     13     , "genus" },
        { 14,     14     , "species" },
        { 14,     13     , "genus" },
        { 16,     16     , "genus" },
        { 17,     17     , "species" },
        { 17,     16     , "genus" },
        { 18,     18     , "genus" },
        { 19,     19     , "species" },
        { 19,     18     , "genus" },
        { 20,     20     , "genus" },
        { 21,     21     , "species" },
        { 21,     20     , "genus" },
        { 22,     22     , "genus" },
        { 339,    339    , "species" },
        { 340,    339    , "species" },
        { 190485, 339    , "species"}
    };

    for (const auto &triple : data) {
        EXPECT_EQ(std::get<1>(triple),
                  mapper.get_ancestor_with_rank_label(std::get<0>(triple),
                                                      std::get<2>(triple)));
    }
}


TEST(TaxIDMapper, ReadValid) {
    TaxIDMapper mapper;
    EXPECT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    EXPECT_TRUE(mapper.parse_nodes(test_nodes));
}

TEST(TaxIDMapper, ReadInValid) {
    TaxIDMapper mapper;
    EXPECT_FALSE(mapper.parse_accession2taxid(
        test_accession2taxid.substr(0, test_accession2taxid.length() - 1)
    ));
    EXPECT_FALSE(mapper.parse_nodes(test_nodes.substr(0, test_nodes.length() - 1)));
}

TEST(TaxIDMapper, GBToTaxID) {
    TaxIDMapper mapper;
    ASSERT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    test_gb_to_taxid(mapper);
}

TEST(TaxIDMapper, IDtoParentAndLabel) {
    TaxIDMapper mapper;
    ASSERT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    ASSERT_TRUE(mapper.parse_nodes(test_nodes));
    test_id_to_parent_and_label(mapper);
}

TEST(TaxIDMapper, IDtoAncestorLabeled) {
    TaxIDMapper mapper;
    ASSERT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    ASSERT_TRUE(mapper.parse_nodes(test_nodes));
    test_id_to_ancestor_labeled(mapper);
}

TEST(TaxIDMapper, IDtoAncestorWithRankLabel) {
    TaxIDMapper mapper;
    ASSERT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    ASSERT_TRUE(mapper.parse_nodes(test_nodes));
    test_id_to_ancestor_with_rank_label(mapper);
}

TEST(TaxIDMapper, SerializeLoad) {
    TaxIDMapper mapper;
    ASSERT_TRUE(mapper.parse_accession2taxid(test_accession2taxid));
    ASSERT_TRUE(mapper.parse_nodes(test_nodes));

    std::ofstream fout(test_dump_basename);
    mapper.serialize(fout);
    fout.close();

    std::ifstream fin(test_dump_basename);
    EXPECT_TRUE(mapper.load(fin));
    fin.close();

    test_gb_to_taxid(mapper);
    test_id_to_parent_and_label(mapper);
    test_id_to_ancestor_labeled(mapper);
    test_id_to_ancestor_with_rank_label(mapper);
}

} // namespace
