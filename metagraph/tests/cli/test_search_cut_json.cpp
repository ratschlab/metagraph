#include <gtest/gtest.h>

#include "cli/query.hpp"

// The /search cut's `matched` and the CLI's `cuts` (docs/DESIGN-pair-identity.md section 3.6;
// SPEC-multigraph-envelope.md section 3), adapted from master's test of the same name: the
// staging cut is SearchCut with its own members (nth_median, tied_same_median, refused), and the
// record's `cuts` are written beside to_json (cuts_json), never by it, so that the server's
// cuts are not stated twice.

namespace {

using mtg::cli::SeqSearchResult;
using mtg::graph::SearchCut;

SearchCut count_tie() {
    SearchCut cut;
    cut.matched = 3;
    cut.returned = 1;
    cut.complete = false;
    cut.nth_score = 7;
    cut.tied = 2;
    cut.tied_returned = 1;
    return cut;
}

TEST(SearchCutJson, ReportsCountTieWithMatched) {
    SeqSearchResult::result_type rows = SeqSearchResult::LabelCountVec{ { "sample", 7 } };
    SeqSearchResult result({ 0, "q", "ACGT" }, std::move(rows), {}, count_tie());
    // to_json states no cut (the server and the CLI add `cuts` beside it)
    const Json::Value json = result.to_json(false, 3);
    EXPECT_EQ(1u, json["results"].size());
    EXPECT_FALSE(json.isMember("cuts"));
    const Json::Value cuts = result.cuts_json(Json::Value(), Json::Value(), Json::Value());
    ASSERT_EQ(1u, cuts.size());
    const Json::Value &entry = cuts[0];
    EXPECT_EQ(3u, entry["matched"].asUInt64());
    EXPECT_EQ(1u, entry["returned"].asUInt64());
    EXPECT_FALSE(entry["complete"].asBool());
    EXPECT_EQ(7u, entry["nth_score"].asUInt64());
    EXPECT_EQ(2u, entry["tied"].asUInt64());
    EXPECT_EQ(1u, entry["tied_returned"].asUInt64());
    // no counts annotation: no median members
    EXPECT_FALSE(entry.isMember("nth_median"));
    EXPECT_FALSE(entry.isMember("tied_same_median"));
    // the CLI's and a single-graph server's pair tags: null, no path
    EXPECT_TRUE(entry["graph"].isNull());
    EXPECT_TRUE(entry["pair_id"].isNull());
    EXPECT_TRUE(entry["index_fp"].isNull());
    EXPECT_FALSE(entry.isMember("graph_path"));
    EXPECT_FALSE(entry.isMember("annotation_path"));
    // master's members are a subset of the staging cut's, with the same values
    for (const char *member : { "graph", "pair_id", "index_fp", "matched", "returned",
                                "complete", "nth_score", "tied", "tied_returned" }) {
        EXPECT_TRUE(entry.isMember(member)) << member;
    }
    // a multi-graph server's tags
    const Json::Value tagged = result.cuts_json("refseq", "n1-x", "3f1c");
    EXPECT_EQ("refseq", tagged[0]["graph"].asString());
    EXPECT_EQ("n1-x", tagged[0]["pair_id"].asString());
    EXPECT_EQ("3f1c", tagged[0]["index_fp"].asString());
}

TEST(SearchCutJson, CompleteAndEmpty) {
    // nothing matched: complete, matched 0
    SeqSearchResult::result_type none = SeqSearchResult::LabelCountVec{};
    SeqSearchResult empty({ 0, "q", "ACGT" }, std::move(none), {}, SearchCut());
    const Json::Value cut = empty.cuts_json(Json::Value(), Json::Value(), Json::Value())[0];
    EXPECT_TRUE(cut["complete"].asBool());
    EXPECT_TRUE(cut["nth_score"].isNull());
    EXPECT_EQ(0u, cut["matched"].asUInt64());
    EXPECT_EQ(0u, cut["returned"].asUInt64());
    // every match returned: matched equals returned
    SearchCut all;
    all.matched = 2;
    all.returned = 2;
    SeqSearchResult::result_type rows = SeqSearchResult::LabelCountVec{ { "a", 3 }, { "b", 1 } };
    SeqSearchResult complete({ 0, "q", "ACGT" }, std::move(rows), {}, all);
    const Json::Value c = complete.cuts_json(Json::Value(), Json::Value(), Json::Value())[0];
    EXPECT_TRUE(c["complete"].asBool());
    EXPECT_EQ(c["matched"], c["returned"]);
}

// A unit refused alone keeps the `matched` its search computed beside returned 0
TEST(SearchCutJson, ARefusedUnitKeepsMatched) {
    SeqSearchResult::result_type rows = SeqSearchResult::LabelCountVec{ { "sample", 7 } };
    SeqSearchResult result({ 0, "q", "ACGT" }, std::move(rows), {}, count_tie());
    result.refuse(SearchCut::Refused::VALUES, 1);
    const Json::Value cut = result.cut_to_json();
    EXPECT_EQ(3u, cut["matched"].asUInt64());
    EXPECT_EQ(0u, cut["returned"].asUInt64());
    EXPECT_FALSE(cut["complete"].asBool());
    EXPECT_TRUE(cut["nth_score"].isNull());
    EXPECT_EQ("values", cut["refused"]["account"].asString());
}

} // namespace
