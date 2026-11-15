#include <gtest/gtest.h>
#include "cli/query.hpp"

namespace {
using mtg::cli::SeqSearchResult;

TEST(SearchCutJson, ReportsCountTieWithoutClaimingMedianOrder) {
    mtg::graph::SearchCut cut{3, 1, 7, 2, 1};
    SeqSearchResult::result_type rows = SeqSearchResult::LabelCountVec{{"sample", 7}};
    SeqSearchResult result({0, "q", "ACGT"}, std::move(rows), {}, cut);
    auto json = result.to_json(false, 3);
    EXPECT_EQ(1u, json["results"].size());
    const auto &entry = json["cuts"][0];
    EXPECT_EQ(3u, entry["matched"].asUInt64());
    EXPECT_EQ(1u, entry["returned"].asUInt64());
    EXPECT_FALSE(entry["complete"].asBool());
    EXPECT_EQ(7u, entry["nth_score"].asUInt64());
    EXPECT_EQ(2u, entry["tied"].asUInt64());
    EXPECT_EQ(1u, entry["tied_returned"].asUInt64());
    EXPECT_FALSE(entry.isMember("nth_median"));
    EXPECT_TRUE(entry["graph"].isNull());
    EXPECT_TRUE(entry["pair_id"].isNull());
    EXPECT_FALSE(entry.isMember("graph_path"));
    EXPECT_FALSE(entry.isMember("annotation_path"));
}

TEST(SearchCutJson, CompleteEmptyAndUnrecordedDiffer) {
    SeqSearchResult::result_type rows = SeqSearchResult::LabelCountVec{};
    SeqSearchResult recorded({0, "q", "ACGT"}, std::move(rows), {}, mtg::graph::SearchCut{});
    auto json = recorded.to_json(false, 3);
    EXPECT_TRUE(json["cuts"][0]["complete"].asBool());
    EXPECT_TRUE(json["cuts"][0]["nth_score"].isNull());
    SeqSearchResult::result_type old_rows = SeqSearchResult::LabelCountVec{};
    SeqSearchResult legacy({0, "q", "ACGT"}, std::move(old_rows));
    EXPECT_FALSE(legacy.to_json(false, 3).isMember("cuts"));
}
} // namespace
