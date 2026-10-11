#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "cli/stats.hpp"
#include "annotation/coord_to_header.hpp"


namespace {

using namespace mtg;

// the report of `stats --check-names` on the names given, without files
TEST(IndexNames, ReportCountsAndNamesTheFirstProblems) {
    const std::string bad_label("label\xff", 6);
    const std::string bad_header("bad\xff", 4);
    // column 0: a repeat that is not consecutive and a header that is not UTF-8;
    // column 1: one header three times
    annot::CoordToHeader mapping({ { "a", "b", "a", bad_header }, { "x", "x", "x" } },
                                 { { 1, 2, 3, 4 }, { 1, 1, 1 } });
    const Json::Value report
            = cli::index_names_report({ "ok", bad_label, "ok2" }, &mapping);

    EXPECT_EQ(3u, report["labels"].asUInt64());
    EXPECT_EQ(1u, report["labels_invalid_utf8"].asUInt64());
    EXPECT_TRUE(report["headers_checked"].asBool());
    EXPECT_EQ(7u, report["headers"].asUInt64());
    EXPECT_EQ(1u, report["headers_invalid_utf8"].asUInt64());
    EXPECT_EQ(3u, report["duplicate_record_headers"].asUInt64());
    EXPECT_FALSE(report["names_ok"].asBool());

    // the labels, then per column the headers that are not UTF-8, then its repeats
    const Json::Value &first = report["first_invalid"];
    ASSERT_EQ(5u, first.size());
    auto entry = [&](Json::ArrayIndex i) {
        const Json::Value &e = first[i];
        return std::make_tuple(e["kind"].asString(), e["problem"].asString(),
                               e["column"].asUInt64(),
                               e["record"].isNull() ? -1 : (int64_t)e["record"].asUInt64(),
                               e["hex"].asString(),
                               e.isMember("duplicate_of") ? (int64_t)e["duplicate_of"].asUInt64() : -1);
    };
    EXPECT_EQ(std::make_tuple(std::string("label"), std::string("invalid_utf8"), 1u, int64_t(-1),
                              std::string("6c6162656cff"), int64_t(-1)), entry(0));
    EXPECT_EQ(std::make_tuple(std::string("header"), std::string("invalid_utf8"), 0u, int64_t(3),
                              std::string("626164ff"), int64_t(-1)), entry(1));
    EXPECT_EQ(std::make_tuple(std::string("header"), std::string("duplicate"), 0u, int64_t(2),
                              std::string("61"), int64_t(0)), entry(2));
    EXPECT_EQ(std::make_tuple(std::string("header"), std::string("duplicate"), 1u, int64_t(1),
                              std::string("78"), int64_t(0)), entry(3));
    EXPECT_EQ(std::make_tuple(std::string("header"), std::string("duplicate"), 1u, int64_t(2),
                              std::string("78"), int64_t(0)), entry(4));
    // the raw bytes are nowhere in the report
    Json::StreamWriterBuilder builder;
    EXPECT_EQ(std::string::npos, Json::writeString(builder, report).find('\xff'));
}

TEST(IndexNames, NamesInOrderPassWithAndWithoutAMapping) {
    const annot::CoordToHeader mapping({ { "a", "b" }, { "a" } }, { { 1, 2 }, { 3 } });
    for (const annot::CoordToHeader *m : { &mapping, (const annot::CoordToHeader *)nullptr }) {
        const Json::Value report
                = cli::index_names_report({ "x", "label\xef\xbf\xbd", "" }, m);
        EXPECT_EQ(3u, report["labels"].asUInt64());
        EXPECT_EQ(0u, report["labels_invalid_utf8"].asUInt64());
        EXPECT_EQ(m != nullptr, report["headers_checked"].asBool());
        EXPECT_EQ(m ? 3u : 0u, report["headers"].asUInt64());
        EXPECT_EQ(0u, report["headers_invalid_utf8"].asUInt64());
        EXPECT_EQ(0u, report["duplicate_record_headers"].asUInt64());
        EXPECT_EQ(0u, report["first_invalid"].size());
        EXPECT_TRUE(report["names_ok"].asBool());
    }
}

TEST(IndexNames, FirstInvalidIsCappedAtTenWhileTheCountsAreNot) {
    std::vector<std::string> labels(12, std::string("\xff", 1));
    std::vector<std::string> headers(12, "same");
    annot::CoordToHeader mapping({ headers }, { std::vector<uint64_t>(12, 1) });
    const Json::Value report = cli::index_names_report(labels, &mapping);
    EXPECT_EQ(12u, report["labels_invalid_utf8"].asUInt64());
    EXPECT_EQ(11u, report["duplicate_record_headers"].asUInt64());
    EXPECT_EQ(10u, report["first_invalid"].size());
    for (const Json::Value &e : report["first_invalid"]) {
        EXPECT_EQ("label", e["kind"].asString());
    }
    EXPECT_FALSE(report["names_ok"].asBool());
}

} // namespace
