#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <json/json.h>

#include "cli/json_helpers.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "common/utils/string_utils.hpp"


// What the server reads of a request (H1, review 5 #2 and #6; SPEC-multigraph-envelope.md
// section 8, SPEC-labeled-traversal-core.md §5 and §10.3): the count of a strict-JSON body's
// values before it is parsed (parse_strict_json), the bodies of the 413 body_too_large and the
// 400 too_many_json_values (error first), request_intake, the /traverse counts checked before
// anything of a request is copied (parse_traverse_request under TraverseLimits), and make_cost's
// table by hash lookup

namespace {

using namespace mtg;
using namespace mtg::cli;
using mtg::graph::traversal::LabelChangeCost;
using mtg::graph::traversal::LabelId;

Json::Value parse(const std::string &text) {
    return parse_strict_json<RefuseRequest>(text, "");
}

std::string text(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

// whether |body| is read under a cap of |max| values
bool fits(const std::string &body, uint64_t max) {
    try {
        parse_strict_json<RefuseRequest>(body, "", max);
        return true;
    } catch (const TooManyJsonValues &e) {
        EXPECT_EQ(max, e.max_values());
        return false;
    }
}

template <class F>
std::string refusal_of(const F &f) {
    try {
        f();
    } catch (const std::exception &e) {
        return e.what();
    }
    return "accepted";
}

} // namespace

// The values counted: the root, each array element and each member's value; a member's name is
// none, an empty array or object holds none. A body of exactly the cap is read, one more value
// is refused before the parse
TEST(ServerRequestIntake, JsonValuesAreCountedExactly) {
    for (const auto &[body, values] : std::vector<std::pair<std::string, uint64_t>> {
             { "7", 1 },
             { "\"x\"", 1 },
             { "null", 1 },
             { "[]", 1 },
             { "{}", 1 },
             { "[ ]", 1 },
             { "{\n}", 1 },
             { "[1]", 2 },
             { "[1, 2, 3]", 4 },
             { "[ \n1 ,\t2\r ]", 3 },
             { "{\"a\": 1}", 2 },
             { "{\"a\": 1, \"b\": []}", 3 },
             { "{\"a\": {\"b\": {\"c\": null}}}", 4 },
             { "[[], {}, [[]]]", 5 },
             { "[[[[]]]]", 4 },
             { "[{\"a\": [1, 2]}, {\"b\": {}}]", 7 },
             // brackets, braces, commas and escaped quotes inside strings are text
             { "[\"[,{\", \"a\\\",b\", \"\\\\\", \"}]\"]", 5 },
             { "{\"[a,b]\": \"{c,d}\"}", 2 },
             { "{\"seeds\": [{\"sequence\": \"ACGT\", \"labels\": [\"a\", \"b\"]}], "
               "\"strategy\": {}}", 8 } }) {
        EXPECT_TRUE(fits(body, values)) << body;
        // (a cap of 0 is none: a body of one value is within every cap)
        if (values > 1) {
            EXPECT_FALSE(fits(body, values - 1)) << body;
            EXPECT_FALSE(fits(body, 1)) << body;
        }
        // 0: no count
        EXPECT_TRUE(fits(body, 0)) << body;
        EXPECT_EQ(text(parse(body)), text(parse_strict_json<RefuseRequest>(body, "", values)))
                << body;
    }
}

// The count is made before jsoncpp reads anything: a body over the cap is refused with its
// value count whatever follows (a comment, a syntax error, an unfinished text)
TEST(ServerRequestIntake, TheCountComesBeforeTheParse) {
    EXPECT_FALSE(fits("[1, 2, 3, x", 3));
    EXPECT_FALSE(fits("[1, 2, 3, // a comment", 3));
    EXPECT_FALSE(fits("[0, 0, 0, 0", 4));
    // within the cap, the text's own fault is refused as before
    EXPECT_EQ("P: a comment (or a '/' outside a string) at byte 9",
              refusal_of([] { parse_strict_json<RefuseRequest>("[1, 2, 3 // x", "P: ", 4); }));
    EXPECT_NE("accepted",
              refusal_of([] { parse_strict_json<RefuseRequest>("[1, 2, 3, x", "P: ", 5); }));
    // the default cap of a server, 2^18: a body of exactly that many values is read, one more
    // refused
    std::string zeros = "[0";
    for (size_t i = 2; i < 262144; ++i) {
        zeros += ",0";
    }
    zeros += "]";
    EXPECT_TRUE(fits(zeros, 262144));
    EXPECT_EQ(262143u, parse_strict_json<RefuseRequest>(zeros, "", 262144).size());
    zeros.insert(zeros.size() - 1, ",0");
    EXPECT_FALSE(fits(zeros, 262144));
    EXPECT_EQ("the request body holds more than 262144 JSON values (--max-request-json-values): "
              "it was not parsed and nothing of it ran; send fewer values (split the request)",
              refusal_of([&] { parse_strict_json<RefuseRequest>(zeros, "P: ", 262144); }));
    // a TooManyJsonValues is an invalid_argument (the CLI's catch of a request fault), and no
    // route's own refusal
    EXPECT_THROW(parse_strict_json<RefuseRequest>("[1,2]", "", 2), std::invalid_argument);
}

// The bodies of the two intake refusals: `error` first, then `code`, then the cap (the service
// recognises a refusal by its leading member; jsoncpp alone writes `code` first)
TEST(ServerRequestIntake, RefusalBodiesWriteErrorFirst) {
    const Json::Value too_large = body_too_large_json(268435457, 256);
    EXPECT_EQ(R"x({"error":"the request body (268435457 bytes) is larger than this server reads, )x"
              R"x(256 MiB (--max-request-body-mb): it was not read and nothing of it ran",)x"
              R"x("code":"body_too_large","max_request_body_mb":256})x",
              error_first_json_text(too_large));
    EXPECT_EQ(R"x({"error":"the request body (chunked) grew past what this server reads, 2 MiB )x"
              R"x((--max-request-body-mb): the rest was not read and nothing of it ran",)x"
              R"x("code":"body_too_large","max_request_body_mb":2})x",
              error_first_json_text(body_too_large_json(std::nullopt, 2)));
    EXPECT_TRUE(too_large["max_request_body_mb"].isUInt64());
    EXPECT_EQ((std::vector<std::string> { "code", "error", "max_request_body_mb" }),
              too_large.getMemberNames());

    const Json::Value values = too_many_json_values_json(262144);
    EXPECT_EQ(R"x({"error":"the request body holds more than 262144 JSON values )x"
              R"x((--max-request-json-values): it was not parsed and nothing of it ran; send )x"
              R"x(fewer values (split the request)","code":"too_many_json_values",)x"
              R"x("max_request_json_values":262144})x",
              error_first_json_text(values));
    EXPECT_EQ(std::string(TooManyJsonValues(262144).what()), values["error"].asString());
    EXPECT_EQ((std::vector<std::string> { "code", "error", "max_request_json_values" }),
              values.getMemberNames());

    // the text is the body's: the same members and values, only the order is its own
    Json::Value reread;
    std::istringstream in(error_first_json_text(values));
    in >> reread;
    EXPECT_EQ(text(values), text(reread));
    // an error with characters JSON escapes is escaped as jsoncpp escapes it
    Json::Value odd;
    odd["error"] = "a \"quoted\" \\ name\n";
    odd["code"] = "x";
    EXPECT_EQ(R"x({"error":"a \"quoted\" \\ name\n","code":"x"})x", error_first_json_text(odd));

    // an HttpError with a text is answered with that text
    const HttpError e(400, values, error_first_json_text(values));
    EXPECT_EQ(error_first_json_text(values), e.text());
    EXPECT_EQ(values["error"].asString(), std::string(e.what()));
    EXPECT_TRUE(HttpError(400, values).text().empty());
}

// GET /capabilities' request_intake: the byte cap (null when unlimited) and the 413's shape,
// the value count (null when unlimited), its routes (the strict-JSON routes and /search) and its
// 400's shape; members sorted
TEST(ServerRequestIntake, CapabilitiesStateTheIntake) {
    EXPECT_EQ(R"x({"body_too_large":{"code":"body_too_large","members":["code","error",)x"
              R"x("max_request_body_mb"],"status":413},"json_values_routes":["/pattern",)x"
              R"x("/resolve","/search","/traverse","/traverse/cancel"],"max_connections":null,)x"
              R"x("max_request_body_mb":256,)x"
              R"x("max_request_json_values":262144,"too_many_json_values":{"code":)x"
              R"x("too_many_json_values","members":["code","error","max_request_json_values"],)x"
              R"x("status":400}})x",
              text(request_intake_json(256, 262144, 0)));
    const Json::Value unlimited = request_intake_json(0, 0, 0);
    EXPECT_TRUE(unlimited["max_request_body_mb"].isNull());
    EXPECT_TRUE(unlimited["max_request_json_values"].isNull());
    EXPECT_TRUE(unlimited["max_connections"].isNull());
    EXPECT_TRUE(request_intake_json(2, 1, 0)["max_request_body_mb"].isUInt64());
    EXPECT_TRUE(request_intake_json(2, 1, 0)["max_request_json_values"].isUInt64());
    EXPECT_TRUE(request_intake_json(2, 1, 3)["max_connections"].isUInt64());
    // the members each refusal states are exactly its body's
    const Json::Value intake = request_intake_json(2, 1, 0);
    Json::Value members(Json::arrayValue);
    for (const std::string &m : body_too_large_json(1, 2).getMemberNames()) {
        members.append(m);
    }
    EXPECT_EQ(text(members), text(intake["body_too_large"]["members"]));
    members = Json::Value(Json::arrayValue);
    for (const std::string &m : too_many_json_values_json(1).getMemberNames()) {
        members.append(m);
    }
    EXPECT_EQ(text(members), text(intake["too_many_json_values"]["members"]));
}

namespace {

// a request of |seeds| seeds of |labels| explicit labels each, |extra| extra labels and a
// table of |entries| entries
Json::Value request_of(size_t seeds, size_t labels, size_t extra, size_t entries,
                       const std::string &sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT") {
    Json::Value r(Json::objectValue);
    for (size_t i = 0; i < seeds; ++i) {
        Json::Value seed(Json::objectValue);
        seed["sequence"] = sequence;
        for (size_t j = 0; j < labels; ++j) {
            seed["labels"].append("l" + std::to_string(j));
        }
        r["seeds"].append(seed);
    }
    Json::Value &l = r["strategy"]["labels"];
    l = Json::Value(Json::objectValue);
    for (size_t j = 0; j < extra; ++j) {
        l["extra"].append("x" + std::to_string(j));
    }
    if (entries) {
        l["change_cost"]["model"] = "table";
        for (size_t j = 0; j < entries; ++j) {
            Json::Value e(Json::arrayValue);
            e.append("l" + std::to_string(j % std::max<size_t>(labels, 1)));
            e.append("x" + std::to_string(j));
            e.append(1.0);
            l["change_cost"]["entries"].append(e);
        }
    } else if (extra) {
        l["change_cost"]["model"] = "constant";
        l["change_cost"]["value"] = 1;
        l["loss_budget"] = 5;
    }
    return r;
}

TraverseLimits counts(uint64_t seeds, uint64_t bp, uint64_t labels, uint64_t extra,
                      uint64_t entries) {
    TraverseLimits limits;
    limits.max_seeds = seeds;
    limits.max_seed_bp = bp;
    limits.max_explicit_labels = labels;
    limits.max_extra_labels = extra;
    limits.max_cost_entries = entries;
    return limits;
}

} // namespace

// Each count at its cap is read, one above it refused with the field, the cap and the flag;
// 0 is no cap
TEST(ServerRequestIntake, TraverseCountsAtAndAboveTheCaps) {
    const TraverseLimits caps = counts(4, 40, 3, 5, 6);
    EXPECT_NO_THROW(parse_traverse_request(request_of(4, 3, 5, 6), caps));
    EXPECT_EQ("request.seeds: at most 4 seeds per request",
              refusal_of([&] { parse_traverse_request(request_of(5, 3, 5, 6), caps); }));
    EXPECT_EQ("request.seeds: a seed is longer than the server limit of 40 bp",
              refusal_of([&] {
                  parse_traverse_request(request_of(4, 3, 5, 6, std::string(41, 'A')), caps);
              }));
    EXPECT_NO_THROW(parse_traverse_request(request_of(4, 3, 5, 6, std::string(40, 'A')), caps));
    EXPECT_EQ("request.seeds[0].labels: at most 3 labels per seed "
              "(--traverse-max-explicit-labels)",
              refusal_of([&] { parse_traverse_request(request_of(4, 4, 5, 6), caps); }));
    EXPECT_EQ("strategy.labels.extra: at most 5 labels (--traverse-max-extra-labels)",
              refusal_of([&] { parse_traverse_request(request_of(4, 3, 6, 6), caps); }));
    EXPECT_EQ("strategy.labels.change_cost.entries: at most 6 entries "
              "(--traverse-max-cost-entries)",
              refusal_of([&] { parse_traverse_request(request_of(4, 3, 5, 7), caps); }));
    // a later seed is named by its index
    Json::Value r = request_of(3, 2, 0, 0);
    r["seeds"][2]["labels"].append("a");
    r["seeds"][2]["labels"].append("b");
    EXPECT_EQ("request.seeds[2].labels: at most 3 labels per seed "
              "(--traverse-max-explicit-labels)",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    // no caps (the CLI, and parse_traverse_request without limits): all read
    EXPECT_NO_THROW(parse_traverse_request(request_of(5, 4, 6, 7, std::string(41, 'A')),
                                           TraverseLimits()));
    EXPECT_NO_THROW(parse_traverse_request(request_of(5, 4, 6, 7, std::string(41, 'A'))));
}

// The counts come before every other fault of the fields, in their order (seeds, per seed its
// length then its labels, extra, entries); a field of another type is left to the parse, which
// refuses it in its own words
TEST(ServerRequestIntake, TraverseCountsComeFirstAndLeaveTypesToTheParse) {
    const TraverseLimits caps = counts(2, 40, 2, 2, 2);
    Json::Value r = request_of(3, 3, 3, 3);
    r["unknown"] = 1;
    r["seeds"][0].removeMember("sequence");
    EXPECT_EQ("request.seeds: at most 2 seeds per request",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["seeds"].resize(2);
    EXPECT_EQ("request.seeds[0].labels: at most 2 labels per seed "
              "(--traverse-max-explicit-labels)",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["seeds"][0]["labels"].resize(2);
    r["seeds"][1]["sequence"] = std::string(41, 'A');
    EXPECT_EQ("request.seeds: a seed is longer than the server limit of 40 bp",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["seeds"][1]["sequence"] = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    EXPECT_EQ("request.seeds[1].labels: at most 2 labels per seed "
              "(--traverse-max-explicit-labels)",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["seeds"][1]["labels"].resize(2);
    EXPECT_EQ("strategy.labels.extra: at most 2 labels (--traverse-max-extra-labels)",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["strategy"]["labels"]["extra"].resize(2);
    EXPECT_EQ("strategy.labels.change_cost.entries: at most 2 entries "
              "(--traverse-max-cost-entries)",
              refusal_of([&] { parse_traverse_request(r, caps); }));
    r["strategy"]["labels"]["change_cost"]["entries"].resize(2);
    // within every count: the parse's own refusals, as before
    EXPECT_EQ("request.seeds[0].sequence is required",
              refusal_of([&] { parse_traverse_request(r, caps); }));

    // other types: the parse's words
    for (const auto &[path, value, words] :
             std::vector<std::tuple<std::vector<std::string>, Json::Value, std::string>> {
                 { { "seeds" }, Json::Value("x"), "request.seeds: expected a non-empty array" },
                 { { "strategy" }, Json::Value(7), "strategy: expected an object" },
                 { { "strategy", "labels" }, Json::Value(7),
                   "strategy.labels: expected an object" },
                 { { "strategy", "labels", "extra" }, Json::Value(7),
                   "strategy.labels.extra: expected an array of strings" },
                 { { "strategy", "labels", "change_cost" }, Json::Value("x"),
                   "strategy.labels.change_cost: expected an object" } }) {
        Json::Value q = request_of(1, 1, 0, 0);
        Json::Value *at = &q;
        for (const std::string &key : path) {
            at = &(*at)[key];
        }
        *at = value;
        EXPECT_EQ(words, refusal_of([&] { parse_traverse_request(q, caps); }))
                << path.back();
    }
    Json::Value q = request_of(1, 1, 0, 0);
    q["seeds"][0] = 5;
    EXPECT_EQ("request.seeds[0]: expected an object",
              refusal_of([&] { parse_traverse_request(q, caps); }));
    q = request_of(1, 1, 0, 0);
    q["seeds"][0]["sequence"] = 5;
    EXPECT_EQ("request.seeds[0].sequence: expected a string",
              refusal_of([&] { parse_traverse_request(q, caps); }));
    q = request_of(1, 1, 0, 0);
    q["seeds"][0]["labels"] = "a";
    EXPECT_EQ("request.seeds[0].labels: expected an array of strings",
              refusal_of([&] { parse_traverse_request(q, caps); }));
    q = Json::Value(5);
    EXPECT_EQ("request: expected an object", refusal_of([&] { parse_traverse_request(q, caps); }));
}

namespace {

// make_cost as it was: two linear scans per entry, a name's id its first position
LabelChangeCost linear_cost(const CostSpec &spec, const std::vector<std::string> &dict) {
    std::map<std::pair<LabelId, LabelId>, double> entries;
    auto id_of = [&](const std::string &name) -> std::optional<LabelId> {
        auto it = std::find(dict.begin(), dict.end(), name);
        if (it == dict.end())
            return std::nullopt;
        return static_cast<LabelId>(it - dict.begin());
    };
    for (const auto &[from, to, c] : spec.entries) {
        auto a = id_of(from), b = id_of(to);
        if (a && b)
            entries[{ *a, *b }] = c;
    }
    return LabelChangeCost::table(std::move(entries), spec.default_cost);
}

} // namespace

// The table make_cost builds by hash lookup is the one the linear scans built: ids by first
// occurrence (a name repeated in the dictionary, a seed label also in extra), entries naming a
// label outside the dictionary dropped, a later entry for the same pair replacing an earlier
TEST(ServerRequestIntake, MakeCostIsTheLinearTable) {
    std::mt19937 rng(42);
    for (size_t round = 0; round < 50; ++round) {
        std::vector<std::string> dict;
        const size_t n = 1 + rng() % 40;
        for (size_t i = 0; i < n; ++i) {
            dict.push_back("n" + std::to_string(rng() % 30));   // repeats on purpose
        }
        CostSpec spec;
        spec.model = CostSpec::TABLE;
        spec.default_cost = round % 2 ? 3.5 : graph::traversal::kInfiniteLoss;
        for (size_t e = 0; e < 60; ++e) {
            spec.entries.emplace_back("n" + std::to_string(rng() % 35),
                                      "n" + std::to_string(rng() % 35),
                                      static_cast<double>(rng() % 7));
        }
        const LabelChangeCost fast
                = make_cost(spec, std::vector<std::string_view>(dict.begin(), dict.end()));
        const LabelChangeCost slow = linear_cost(spec, dict);
        EXPECT_EQ(LabelChangeCost::TABLE, fast.model());
        EXPECT_EQ(slow.table(), fast.table()) << round;
        EXPECT_EQ(slow.default_cost(), fast.default_cost());
    }
    CostSpec constant;
    constant.model = CostSpec::CONSTANT;
    constant.value = 2;
    EXPECT_EQ(LabelChangeCost::CONSTANT, make_cost(constant, {}).model());
    EXPECT_EQ(2.0, make_cost(constant, {}).cost(0, 1));
    EXPECT_EQ(LabelChangeCost::FORBID, make_cost(CostSpec(), { "a" }).model());
}

// The reviewer's case (#6), uncapped: 100,000 labels and 50,000 entries naming absent labels
// took 10^10 string comparisons; by hash lookup it takes milliseconds. With the caps (10,000
// each) a seed's preprocessing is a fraction of that
TEST(ServerRequestIntake, MakeCostIsLinear) {
    std::vector<std::string> dict;
    for (size_t i = 0; i < 100'000; ++i) {
        dict.push_back(std::to_string(100'000 + i));
    }
    CostSpec spec;
    spec.model = CostSpec::TABLE;
    for (size_t e = 0; e < 50'000; ++e) {
        spec.entries.emplace_back("absent", "absent", 1.0);
    }
    const auto start = std::chrono::steady_clock::now();
    const LabelChangeCost cost
            = make_cost(spec, std::vector<std::string_view>(dict.begin(), dict.end()));
    const double ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
    EXPECT_TRUE(cost.table().empty());
    // generous for a loaded CI host: the quadratic scan took minutes
    EXPECT_LT(ms, 2000.0);
}

// RFC 8259 §7: a raw control character inside a string is the body's JSON fault, refused by the
// pre-scan with the route's prefix (jsoncpp kept it, and an answer quoting it wrote it as
// \u00XX, 6 bytes for 1); an escaped one, a raw byte from 0x20 up and valid UTF-8 are read as
// before
TEST(ServerRequestIntake, ARawControlCharacterInAStringIsRefused) {
    auto refusal = [](const std::string &text) {
        return refusal_of([&] { parse_strict_json<RefuseRequest>(text, "Bad json received: "); });
    };
    EXPECT_EQ("Bad json received: a control character (byte 0x01) inside a string at byte 8: "
              "RFC 8259 requires it escaped (\\u0001)",
              refusal(std::string("{\"a\": \"x\x01y\"}")));
    EXPECT_EQ("Bad json received: a control character (byte 0x09) inside a string at byte 2: "
              "RFC 8259 requires it escaped (\\u0009)",
              refusal(std::string("[\"\tx\"]")));
    EXPECT_EQ("Bad json received: a control character (byte 0x1f) inside a string at byte 2: "
              "RFC 8259 requires it escaped (\\u001f)",
              refusal(std::string("{\"\x1f\": 0}")));
    EXPECT_EQ("Bad json received: a control character (byte 0x00) inside a string at byte 2: "
              "RFC 8259 requires it escaped (\\u0000)",
              refusal(std::string("[\"\0\"]", 5)));
    // whitespace outside strings, escapes, a byte from 0x20 up, UTF-8: read
    EXPECT_EQ("x\x01y\n\t", parse("{\"a\": \"x\\u0001y\\n\\t\"}")["a"].asString());
    EXPECT_EQ(" \x7f\xc3\xa9", parse("[\" \x7f\xc3\xa9\"]")[0].asString());
    EXPECT_EQ(1, parse("{\n\t\"a\":\r\n 1}")["a"].asInt());
    // the first fault the scan meets is refused: a comment before the byte, the byte before a
    // comment, the value past the count before the byte
    EXPECT_EQ("Bad json received: a comment (or a '/' outside a string) at byte 1",
              refusal(std::string("[/\"\x01\"]")));
    EXPECT_NE(std::string::npos, refusal(std::string("[\"\x01\", /]")).find("byte 0x01"));
    EXPECT_FALSE(fits(std::string("[0, 0, \"\x01\"]"), 2));
    // /pattern's prefix: its invalid_request
    try {
        parse_strict_json<RefuseRequest>(std::string("{\"a\": \"\x02\"}"), "request: not JSON: ");
        FAIL();
    } catch (const InvalidRequest &e) {
        EXPECT_EQ(0u, std::string(e.what()).find("request: not JSON: a control character (byte "
                                                 "0x02)"));
    }
}

// RFC 8259 §8.1 (GPT review 6, finding 2): a string, member names included, is well-formed UTF-8
// (RFC 3629). jsoncpp kept any byte, and its writer wrote a byte 0xF8 to 0xFF, or a lead byte cut
// by the string's end, as \ufffd (6 bytes for 1), and let other invalid lead bytes swallow the
// next bytes, so that distinct names echoed alike. The pre-scan refuses the first byte that
// begins no well-formed sequence, naming it and its offset; valid 2-, 3- and 4-byte sequences are
// read as before
TEST(ServerRequestIntake, AStringThatIsNotUtf8IsRefused) {
    auto refusal = [](const std::string &text) {
        return refusal_of([&] { parse_strict_json<RefuseRequest>(text, "Bad json received: "); });
    };
    auto invalid_at = [](const std::string &byte, size_t at) {
        return "Bad json received: invalid UTF-8 (byte 0x" + byte + ") inside a string at byte "
                + std::to_string(at) + ": RFC 8259 requires UTF-8";
    };
    // each class: a continuation byte without a lead, C0 and C1 (overlong only), F5 to FF (no
    // lead byte), a sequence cut by the closing quote or by the end of the body, a lead byte
    // followed by something other than a continuation byte, the overlong forms of E0 and F0,
    // the surrogates ED A0 to ED BF, past U+10FFFF (F4 90)
    const std::vector<std::pair<std::string, std::string>> invalid = {
        { "\x80", "80" }, { "\xbf", "bf" }, { "\xc0\xaf", "c0" }, { "\xc1\xbf", "c1" },
        { "\xf5\x80\x80\x80", "f5" }, { "\xf8\x88\x80\x80\x80", "f8" }, { "\xfe", "fe" },
        { "\xff", "ff" }, { "\xc3", "c3" }, { "\xe2\x82", "e2" }, { "\xf0\x9f\x98", "f0" },
        { "\xc3" "A", "c3" }, { "\xe9 au", "e9" }, { "\xe2\x82" "A", "e2" },
        { "\xe0\x80\x80", "e0" }, { "\xe0\x9f\xbf", "e0" }, { "\xed\xa0\x80", "ed" },
        { "\xed\xbf\xbf", "ed" }, { "\xf0\x80\x80\x80", "f0" }, { "\xf0\x8f\xbf\xbf", "f0" },
        { "\xf4\x90\x80\x80", "f4" }, { "\xf4\xbf\xbf\xbf", "f4" },
    };
    for (const auto &[bytes, lead] : invalid) {
        EXPECT_EQ(invalid_at(lead, 2), refusal("[\"" + bytes + "\"]")) << lead;
        // in a member name, after other characters and after valid UTF-8
        EXPECT_EQ(invalid_at(lead, 2), refusal("{\"" + bytes + "\": 0}")) << lead;
        EXPECT_EQ(invalid_at(lead, 7), refusal("[\"Caf\xc3\xa9" + bytes + "x\"]")) << lead;
    }
    // cut by the end of the body (no closing quote)
    EXPECT_EQ(invalid_at("e2", 2), refusal("[\"\xe2\x82"));
    // the boundaries of each length, and the largest code points of each: read, unchanged
    for (const std::string &valid : {
            std::string("\xc2\x80"), std::string("\xdf\xbf"), std::string("\xe0\xa0\x80"),
            std::string("\xed\x9f\xbf"), std::string("\xee\x80\x80"), std::string("\xef\xbf\xbd"),
            std::string("\xef\xbf\xbf"), std::string("\xf0\x90\x80\x80"),
            std::string("\xf0\x9f\x98\x80"), std::string("\xf4\x8f\xbf\xbf"),
            std::string("Caf\xc3\xa9 \xe2\x82\xac \xf0\x9f\x98\x80 x") }) {
        EXPECT_EQ(valid, parse("[\"" + valid + "\"]")[0].asString());
        EXPECT_EQ(1, parse("{\"" + valid + "\": 1}")[valid].asInt());
    }
    // a byte from 0x80 up outside a string is jsoncpp's syntax error, with the prefix
    const std::string outside = refusal("[\xff]");
    EXPECT_EQ(0u, outside.find("Bad json received: "));
    EXPECT_EQ(std::string::npos, outside.find("invalid UTF-8"));
    // the first fault the scan meets: a comment before the byte, the byte before a comment or a
    // control character, the value past the count before the byte
    EXPECT_EQ("Bad json received: a comment (or a '/' outside a string) at byte 1",
              refusal("[/\"\xff\"]"));
    EXPECT_EQ(invalid_at("ff", 2), refusal("[\"\xff\", /]"));
    EXPECT_EQ(invalid_at("ff", 2), refusal(std::string("[\"\xff\x01\"]")));
    EXPECT_FALSE(fits("[0, 0, \"\xff\"]", 2));
    // /pattern's prefix: its invalid_request
    try {
        parse_strict_json<RefuseRequest>(std::string("{\"a\": \"\xfe\"}"), "request: not JSON: ");
        FAIL();
    } catch (const InvalidRequest &e) {
        EXPECT_EQ("request: not JSON: invalid UTF-8 (byte 0xfe) inside a string at byte 7: "
                  "RFC 8259 requires UTF-8", std::string(e.what()));
    }
    // one rule with valid_utf8 (label names, MGT strings)
    for (const auto &[bytes, lead] : invalid) {
        EXPECT_FALSE(::utils::valid_utf8(bytes)) << lead;
        EXPECT_EQ(0u, ::utils::utf8_sequence_at(bytes, 0)) << lead;
    }
    EXPECT_EQ(1u, ::utils::utf8_sequence_at("a", 0));
    EXPECT_EQ(2u, ::utils::utf8_sequence_at("\xc2\x80", 0));
    EXPECT_EQ(3u, ::utils::utf8_sequence_at("\xef\xbf\xbf", 0));
    EXPECT_EQ(4u, ::utils::utf8_sequence_at("\xf4\x8f\xbf\xbf", 0));
    EXPECT_TRUE(::utils::valid_utf8("Caf\xc3\xa9 \xe2\x82\xac \xf0\x9f\x98\x80"));
}

// A \u escape of a surrogate is half of a pair: a high one (D800-DBFF) followed at once by a low
// one (DC00-DFFF), either case of hex digit. jsoncpp refused a lone high one, decoded a lone low
// one to bytes that are no UTF-8, and joined a high one with any escape after it ("\ud800\u0041"
// read as U+10041). The pre-scan refuses each, quoting the escape as written; a pair and every
// other escape are read as before
TEST(ServerRequestIntake, AnUnpairedSurrogateEscapeIsRefused) {
    auto refusal = [](const std::string &text) {
        return refusal_of([&] { parse_strict_json<RefuseRequest>(text, "Bad json received: "); });
    };
    auto unpaired_at = [](const std::string &escape, size_t at) {
        return "Bad json received: an unpaired surrogate escape (\\u" + escape + ") inside a "
               "string at byte " + std::to_string(at) + ": a code point above U+FFFF is written as "
               "a high surrogate escape (\\ud800-\\udbff) followed by a low one (\\udc00-\\udfff)";
    };
    // a lone low one, a lone high one (at the string's end, before a character, before an
    // escape that is no low surrogate, before another high one, before an escaped character)
    EXPECT_EQ(unpaired_at("dc00", 2), refusal("[\"\\udc00\"]"));
    EXPECT_EQ(unpaired_at("DFFF", 2), refusal("[\"\\uDFFF\"]"));
    EXPECT_EQ(unpaired_at("d800", 2), refusal("[\"\\ud800\"]"));
    EXPECT_EQ(unpaired_at("d800", 2), refusal("[\"\\ud800A\"]"));
    EXPECT_EQ(unpaired_at("d800", 2), refusal("[\"\\ud800\\u0041\"]"));
    EXPECT_EQ(unpaired_at("d800", 2), refusal("[\"\\ud800\\ud800\"]"));
    EXPECT_EQ(unpaired_at("DBFF", 2), refusal("[\"\\uDBFF\\n\"]"));
    EXPECT_EQ(unpaired_at("d800", 2), refusal("[\"\\ud800"));
    // after a valid pair, after other text, in a member name
    EXPECT_EQ(unpaired_at("dc00", 14), refusal("[\"\\ud83d\\ude00\\udc00\"]"));
    EXPECT_EQ(unpaired_at("dc00", 4), refusal("[\"ab\\udc00\"]"));
    EXPECT_EQ(unpaired_at("dc00", 2), refusal("{\"\\udc00\": 0}"));
    // pairs, either case: read as their code point
    EXPECT_EQ("\xf0\x9f\x98\x80", parse("[\"\\ud83d\\ude00\"]")[0].asString());
    EXPECT_EQ("\xf0\x9f\x98\x80", parse("[\"\\uD83D\\uDE00\"]")[0].asString());
    EXPECT_EQ("\xf0\x90\x80\x80", parse("[\"\\ud800\\udc00\"]")[0].asString());
    EXPECT_EQ("\xf4\x8f\xbf\xbf", parse("[\"\\udbff\\udfff\"]")[0].asString());
    EXPECT_EQ("a\xf0\x9f\x98\x80" "b\xf0\x9f\x98\x80",
              parse("[\"a\\ud83d\\ude00b\\ud83d\\ude00\"]")[0].asString());
    // other escapes: read as before, \u0000 and the code units around the surrogates included
    EXPECT_EQ(std::string("\0", 1), parse("[\"\\u0000\"]")[0].asString());
    EXPECT_EQ("\xed\x9f\xbf\xee\x80\x80\xef\xbf\xbf\xc3\xa9",
              parse("[\"\\ud7ff\\ue000\\uffff\\u00E9\"]")[0].asString());
    // an escaped backslash before "udc00" is no escape of a surrogate
    EXPECT_EQ("\\udc00", parse("[\"\\\\udc00\"]")[0].asString());
    // a \u without four hex digits is jsoncpp's to refuse, with the prefix
    const std::string short_escape = refusal("[\"\\u12\"]");
    EXPECT_EQ(0u, short_escape.find("Bad json received: "));
    EXPECT_EQ(std::string::npos, short_escape.find("surrogate escape"));
    // /pattern's prefix: its invalid_request
    try {
        parse_strict_json<RefuseRequest>(std::string("{\"a\": \"\\udc00\"}"),
                                         "request: not JSON: ");
        FAIL();
    } catch (const InvalidRequest &e) {
        EXPECT_EQ(0u, std::string(e.what()).find("request: not JSON: an unpaired surrogate escape "
                                                 "(\\udc00) inside a string at byte 7"));
    }
}

// A backslash that begins no JSON escape (jsoncpp reads \" \\ \/ \b \f \n \r \t and \u with four
// hex digits, nothing else) ends the pre-scan: jsoncpp refuses the text at that escape, in its own
// words, so that neither a byte the escape was meant to take (a backslash before a UTF-8
// character was "invalid UTF-8" at its continuation byte) nor a later fault (a short \u before an
// unpaired surrogate escape) is named instead (train 3 review R1, N1). Past every valid escape
// the scan goes on
TEST(ServerRequestIntake, ABackslashThatBeginsNoEscapeEndsTheScan) {
    auto refusal = [](const std::string &text, uint64_t max_values = 0) {
        return refusal_of([&] {
            parse_strict_json<RefuseRequest>(text, "Bad json received: ", max_values);
        });
    };
    auto names = [](const std::string &refused, const std::string &fault) {
        return refused.find("Bad json received: ") == 0 && refused.find(fault) != std::string::npos;
    };
    const std::string bad_escape = "Bad escape sequence in string";
    const std::string bad_unicode = "Bad unicode escape sequence in string";
    // a backslash before a 2-, 3- or 4-byte character, before a byte that begins none, before a
    // letter, a digit, a space or a control character
    for (const std::string &after : { std::string("\xc3\xa9"), std::string("\xe2\x82\xac"),
                                      std::string("\xf0\x9f\x98\x80"), std::string("\xff"),
                                      std::string("q"), std::string("0"), std::string(" "),
                                      std::string("\x01") }) {
        for (const std::string &body : { "[\"\\" + after + "\"]", "{\"\\" + after + "\": 0}" }) {
            const std::string refused = refusal(body);
            EXPECT_TRUE(names(refused, bad_escape)) << refused;
            EXPECT_EQ(std::string::npos, refused.find("invalid UTF-8")) << refused;
            EXPECT_EQ(std::string::npos, refused.find("control character")) << refused;
        }
    }
    // a \u without four hex digits: jsoncpp's text, not that of a fault after it
    for (const std::string &later : { std::string("\\udc00"), std::string("\\ud800"),
                                      std::string("\xff"), std::string("\x01") }) {
        const std::string refused = refusal("[\"\\u12" + later + "\"]");
        EXPECT_TRUE(names(refused, bad_unicode)) << refused;
        EXPECT_EQ(std::string::npos, refused.find("surrogate escape (")) << refused;
        EXPECT_EQ(std::string::npos, refused.find("invalid UTF-8")) << refused;
        EXPECT_EQ(std::string::npos, refused.find("control character")) << refused;
    }
    // a bad escape before a later fault of the scan's: the escape's (an unpaired surrogate
    // escape, a comment, the value past the count)
    EXPECT_TRUE(names(refusal("[\"\\q\", \"\\udc00\"]"), bad_escape));
    const std::string comment = refusal("[\"\\q\"] /");
    EXPECT_TRUE(names(comment, bad_escape)) << comment;
    EXPECT_EQ(std::string::npos, comment.find("a comment")) << comment;
    EXPECT_TRUE(names(refusal("[\"\\q\", 0, 0, 0]", 2), bad_escape));
    // a fault before the bad escape is still named first
    EXPECT_EQ("Bad json received: invalid UTF-8 (byte 0xff) inside a string at byte 2: RFC 8259 "
              "requires UTF-8", refusal("[\"\xff\\q\"]"));
    EXPECT_EQ("Bad json received: a comment (or a '/' outside a string) at byte 1",
              refusal("[/\"\\q\"]"));
    EXPECT_FALSE(fits("[0, 0, \"\\q\"]", 2));
    // past every valid escape the scan goes on: a fault after one is the scan's, at its offset
    for (const std::string &escape : { std::string("\\\""), std::string("\\\\"),
                                       std::string("\\/"), std::string("\\b"),
                                       std::string("\\f"), std::string("\\n"),
                                       std::string("\\r"), std::string("\\t"),
                                       std::string("\\u0041"), std::string("\\ud83d\\ude00") }) {
        EXPECT_EQ("Bad json received: invalid UTF-8 (byte 0xff) inside a string at byte "
                  + std::to_string(2 + escape.size()) + ": RFC 8259 requires UTF-8",
                  refusal("[\"" + escape + "\xff\"]")) << escape;
    }
    // /pattern's prefix: its invalid_request
    try {
        parse_strict_json<RefuseRequest>(std::string("{\"a\": \"\\\xc3\xa9\"}"),
                                         "request: not JSON: ");
        FAIL();
    } catch (const InvalidRequest &e) {
        EXPECT_EQ(0u, std::string(e.what()).find("request: not JSON: "));
        EXPECT_NE(std::string::npos, std::string(e.what()).find(bad_escape));
    }
}

// An error text quotes a client's string whole up to 1 KiB, else its first 1 KiB (never inside
// a UTF-8 sequence), "..." and its size
TEST(ServerRequestIntake, QuotedClientStringsAreClipped) {
    EXPECT_EQ("", ::utils::clip_for_message(""));
    EXPECT_EQ(std::string(1024, 'n'), ::utils::clip_for_message(std::string(1024, 'n')));
    EXPECT_EQ(std::string(1024, 'n') + "... (1025 bytes)",
              ::utils::clip_for_message(std::string(1025, 'n')));
    // a two-byte character across the cut is left out whole
    std::string e;
    for (size_t i = 0; i < 1000; ++i) {
        e += "\xc3\xa9";
    }
    EXPECT_EQ("a" + e.substr(0, 1022) + "... (2001 bytes)", ::utils::clip_for_message("a" + e));
    EXPECT_EQ("abc... (6 bytes)", ::utils::clip_for_message("abcdef", 3));
    // a name that is no UTF-8 (continuation bytes only) is cut at most 3 bytes back: a clipped
    // name keeps at least 1,021 bytes, so its text is longer than any name quoted whole
    const std::string continuation(2000, '\x80');
    EXPECT_EQ(continuation.substr(0, 1021) + "... (2000 bytes)",
              ::utils::clip_for_message(continuation));
    EXPECT_LT(1024u, ::utils::clip_for_message(std::string(1025, 'n')).size());
    // the unknown field of a strict object
    Json::Value v(Json::objectValue);
    v[std::string(5000, 'k')] = 1;
    EXPECT_EQ("request: unknown field '" + std::string(1024, 'k') + "... (5000 bytes)'",
              refusal_of([&] {
                  StrictObject<RefuseRequest> o(v, "request");
                  o.finish();
              }));
}

// A log line quotes a client's text with its control characters escaped (GPT review 7, finding
// 11): a label "x\n[...] [info] FORGED" made a second log line that looked like the server's.
// Bytes below 0x20 and 0x7F become \n, \r, \t or \xNN; every other byte, a UTF-8 name's too,
// and a backslash the client wrote stay as they are
TEST(ServerRequestLog, ControlCharactersAreEscaped) {
    EXPECT_EQ("", ::utils::escape_for_log(""));
    EXPECT_EQ("Label not found: 'x\\n[2026-01-01 00:00:00.000] [info] FORGED'",
              ::utils::escape_for_log("Label not found: 'x\n[2026-01-01 00:00:00.000] [info] "
                                      "FORGED'"));
    EXPECT_EQ("a\\r\\nb\\tc\\x00d\\x1b[31me\\x7f\\x1f", ::utils::escape_for_log(
            std::string("a\r\nb\tc\0d\x1b[31me\x7f\x1f", 16)));
    EXPECT_EQ("caf\xc3\xa9 \\n kept \x80\xff",
              ::utils::escape_for_log("caf\xc3\xa9 \\n kept \x80\xff"));
    // every byte: one line, nothing below 0x20 and no 0x7F left
    std::string all;
    for (int c = 0; c < 256; ++c) {
        all.push_back(static_cast<char>(c));
    }
    const std::string escaped = ::utils::escape_for_log(all);
    EXPECT_EQ(std::string::npos, escaped.find('\n'));
    for (unsigned char c : escaped) {
        EXPECT_TRUE(c >= 0x20 && c != 0x7F) << int(c);
    }
    // 33 control bytes: \n, \r and \t two bytes each, the other 30 four
    EXPECT_EQ(256u - 33 + 3 * 2 + 30 * 4, escaped.size());
}

// jsoncpp's own messages that quote the client's text are clipped by the same rule: the quoted
// text whole up to 1,024 bytes of its own, else its first 1,024 bytes and "... (N bytes)"
TEST(ServerRequestIntake, TheReadersQuotesAreClipped) {
    auto refusal = [](const std::string &text) {
        return refusal_of([&] { parse_strict_json<RefuseRequest>(text, "P: "); });
    };
    const std::string k(3000, 'k');
    const std::string duplicate = refusal("{\"" + k + "\": 0, \"" + k + "\": 1}");
    EXPECT_NE(std::string::npos,
              duplicate.find("Duplicate key: '" + std::string(1024, 'k') + "... (3000 bytes)'"))
            << duplicate.substr(0, 200);
    EXPECT_LT(duplicate.size(), 1200u);
    // a name with a quote inside: its closing quote is the message's last
    const std::string quoted = "a'" + std::string(2000, 'q');
    EXPECT_NE(std::string::npos,
              refusal("{\"" + quoted + "\": 0, \"" + quoted + "\": 1}")
                      .find("Duplicate key: '" + quoted.substr(0, 1024) + "... (2002 bytes)'"));
    // short: whole, as before
    EXPECT_NE(std::string::npos,
              refusal("{\"" + std::string(1024, 'k') + "\": 0, \"" + std::string(1024, 'k')
                      + "\": 1}").find("Duplicate key: '" + std::string(1024, 'k') + "'\n"));
    // a number token no double holds
    const std::string number = refusal("[" + std::string(3000, '1') + "]");
    EXPECT_NE(std::string::npos,
              number.find("'" + std::string(1024, '1') + "... (3000 bytes)' is not a number."))
            << number.substr(0, 200);
    EXPECT_EQ("x", clip_reader_quotes("x"));
    EXPECT_EQ("'1e' is not a number.", clip_reader_quotes("'1e' is not a number."));
}

// A Content-Length as a byte count (RFC 9110 §8.6): digits with leading zeros, the optional
// whitespace trimmed, a list of identical counts; nothing else
TEST(ServerRequestIntake, ContentLengthValues) {
    for (const auto &[text, value] : std::vector<std::pair<std::string, uint64_t>> {
             { "3145728", 3145728 }, { " 3145728", 3145728 }, { "3145728 ", 3145728 },
             { "\t3145728\t", 3145728 }, { "000000000000003145728", 3145728 },
             { "3145728, 3145728", 3145728 }, { "3145728,3145728 ,\t3145728", 3145728 },
             { "0", 0 }, { "000", 0 }, { "18446744073709551615", 18446744073709551615ull } }) {
        EXPECT_EQ(std::optional<uint64_t>(value), content_length_value(text)) << text;
    }
    for (const std::string text : { "", " ", "-1", "+1", "3145728junk", "abc", "1 2", "1, 2",
                                    "1,", ",1", "1,,1", "0x10", "18446744073709551616",
                                    "99999999999999999999999" }) {
        EXPECT_EQ(std::nullopt, content_length_value(text)) << text;
    }
}

namespace {

// the values of |v| as the reader built it: the root, each array element and each member's value
uint64_t values_of(const Json::Value &v) {
    uint64_t n = 1;
    if (v.isArray() || v.isObject()) {
        for (const Json::Value &child : v) {
            n += values_of(child);
        }
    }
    return n;
}

// |text| as POST /search reads its body (parse_json_string: jsoncpp's defaults); none when the
// reader refuses it
std::optional<Json::Value> lenient(const std::string &text) {
    Json::CharReaderBuilder builder;
    std::unique_ptr<Json::CharReader> reader { builder.newCharReader() };
    Json::Value v;
    std::string errors;
    try {
        if (!reader->parse(text.data(), text.data() + text.size(), &v, &errors))
            return std::nullopt;
    } catch (const Json::Exception &) {
        return std::nullopt;
    }
    return v;
}

// whether |text| names a member twice in one object (jsoncpp's rejectDupKeys, else its defaults)
bool names_a_member_twice(const std::string &text) {
    Json::CharReaderBuilder builder;
    builder["rejectDupKeys"] = true;
    std::unique_ptr<Json::CharReader> reader { builder.newCharReader() };
    Json::Value v;
    std::string errors;
    try {
        return !reader->parse(text.data(), text.data() + text.size(), &v, &errors);
    } catch (const Json::Exception &) {
        return true;
    }
}

// whether |body| is read under a cap of |max| values by /search's count
bool search_fits(const std::string &body, uint64_t max) {
    try {
        count_lenient_json_values(body, max);
        return true;
    } catch (const TooManyJsonValues &e) {
        EXPECT_EQ(max, e.max_values());
        return false;
    }
}

} // namespace

// /search's body is still read by jsoncpp's lenient defaults (Fable review 2026-11-12): its
// values are counted as that reader reads them, so that the count of every body it accepts is
// the number of values of the tree it builds: a comment in each place the reader skips one
// (jsoncpp 1.9.3 takes a ' outside a string for a /, and keeps such a comment when the comment
// text it attaches to a value begins with a /), a trailing comma, the text after the root value
// (not read), strings holding brackets, commas and escapes, and the token an object takes for
// the separator after a member's value and a comment, whatever it is (a NUL byte included)
TEST(ServerRequestIntake, SearchBodiesAreCountedAsTheirReaderReadsThem) {
    for (const std::string &body : std::vector<std::string> {
             "7", "\"x\"", "null", "[]", "{}", "[ ]", "{\n}", "[1, 2, 3]",
             "{\"a\": 1, \"b\": []}", "[[], {}, [[]]]",
             "/* a, b */ [1, 2]", "// , ,\n[1, 2]", "[1 /* , */, 2]", "[1, /* , */ 2]",
             "{\"a\": /* , */ 1}", "{\"a\": 1 /* , */, \"b\": 2}", "{\"a\": 1, /* , */ \"b\": 2}",
             "{ /* , */ }", "[1, // , ,\r\n 2]", "[1, // , ,\r 2]", "{/**/\"a\": [/***/1]}",
             "[/**/'* , , */ 1]", "[//\n'/ , ,\n 1]", "/**/'* [ */ [1, 2]",
             "{\"a\": /**/'* , */ 1}",
             "[1, 2,]", "{\"a\": 1,}", "{\"a\": [1,], \"b\": {},}", "{\"a\": 1, /* c */ }",
             "[1, 2] [3, 4, 5]", "{\"a\": 1} , , ,", "7, 8, 9", "\"x\" [1, 2]", "[1]x,y,z",
             "[1] /* , , unclosed", "[1] // , ,", "true [1, 2]",
             "[\"[,{\", \"a\\\",b\", \"\\\\\", \"}]\"]", "{\"a\\\\\": \",\"}",
             "[\"/* , */\", \"'*\", \"//,\"]", "{\"FASTA\": \">q\\nACGT\", \"top_labels\": 5}",
             "{\"a\": 1 /**/ ] \"b\": 2}", "{\"a\": 1 // ,\n 7 \"b\": [2, 3]}",
             "{\"a\": 1 /**/ @ \"b\": [2]}", "{\"a\": [1] /**/ { \"b\": 2}",
             "{\"a\": 1 /**/" + std::string(1, '\0') + " \"b\": 2}",
             std::string(1000, '[') + std::string(1000, ']') }) {
        const std::optional<Json::Value> tree = lenient(body);
        ASSERT_TRUE(tree) << body;
        const uint64_t values = values_of(*tree);
        EXPECT_EQ(values, count_lenient_json_values(body, 0)) << body;
        EXPECT_TRUE(search_fits(body, values)) << body;
        if (values > 1) {
            EXPECT_FALSE(search_fits(body, values - 1)) << body;
        }
        EXPECT_TRUE(search_fits(body, 0)) << body;
    }
    // a member named twice is counted each time: the reader reads both values (it keeps the last)
    EXPECT_EQ(5u, count_lenient_json_values("{\"a\": 1, \"a\": [2, 3]}", 0));
    EXPECT_EQ(4u, values_of(*lenient("{\"a\": 1, \"a\": [2, 3]}")));
    // past its stack limit (1000 nested values) the reader throws: the count stops there too
    const std::string deep = std::string(1001, '[') + std::string(1001, ']');
    EXPECT_FALSE(lenient(deep).has_value());
    EXPECT_EQ(1001u, count_lenient_json_values(deep, 0));
    EXPECT_EQ(1001u, count_lenient_json_values(deep + deep, 0));
}

// The same on random bodies: values, arrays and objects nested up to 4 deep, with comments of
// every form wherever the reader skips one, trailing commas and text after the root value, and
// on each of them changed at a few random places by a structural byte (inserted, put in place
// of another or a byte deleted). The reader refuses some of them (a comment text that begins
// with a ' attached to a value: jsoncpp's "Comments must start with /", or a change): only the
// bodies it accepts are compared, and of the changed ones those naming no member twice
TEST(ServerRequestIntake, RandomSearchBodiesAreCountedAsTheirReaderReadsThem) {
    std::mt19937 rng(20261113);
    auto pick = [&](int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); };
    // what may stand between two tokens: whitespace, or a comment with commas and brackets
    auto gap = [&](bool comment) -> std::string {
        switch (pick(comment ? 8 : 3)) {
            case 0: return "";
            case 1: return " ";
            case 2: return "\n\t ";
            case 3: return "/* , [ { , */";
            case 4: return "// , ] } ,\n";
            case 5: return pick(2) ? "'* , , */" : "/**/'* , , */";
            case 6: return pick(2) ? "'/ , ,\r\n" : "\n/* */'/ , ,\r\n";
            default: return " /**/ ";
        }
    };
    std::function<std::string(int)> value = [&](int depth) -> std::string {
        switch (pick(depth >= 4 ? 4 : 7)) {
            case 0: return std::to_string(pick(1000));
            case 1: return pick(2) ? "\"s,[{\\\"}]\"" : "\"\"";
            case 2: return pick(2) ? "true" : "null";
            case 3: return "-1.5e3";
            case 4:
            case 5: {
                const int n = pick(5);
                if (!n)
                    return "[" + gap(false) + "]";
                std::string s = "[" + gap(true);
                for (int i = 0; i < n; ++i) {
                    if (i)
                        s += "," + gap(true);
                    s += value(depth + 1) + gap(true);
                }
                // a trailing comma, the bracket after it (a comment there is refused)
                if (!pick(3))
                    s += "," + gap(false);
                return s + "]";
            }
            default: {
                const int n = pick(5);
                std::string s = "{" + gap(true);
                for (int i = 0; i < n; ++i) {
                    if (i)
                        s += "," + gap(true);
                    s += "\"k" + std::to_string(i) + ",}\":" + gap(true) + value(depth + 1)
                            + gap(true);
                }
                if (n && !pick(3))
                    s += "," + gap(true);
                return s + "}";
            }
        }
    };
    size_t checked = 0;
    size_t changed_checked = 0;
    for (int round = 0; round < 3000; ++round) {
        std::string body = gap(true) + value(0);
        if (!pick(4))
            body += gap(true) + (pick(2) ? ", [0, 0, 0]" : "{,,} x");
        const std::optional<Json::Value> tree = lenient(body);
        if (!tree)
            continue;
        EXPECT_EQ(values_of(*tree), count_lenient_json_values(body, 0)) << body;
        ++checked;
        for (int m = 0; m < 4; ++m) {
            std::string changed = body;
            static const std::string kBytes(",[]{}\"'/*: \n\r0+-I.e\\\0", 21);
            for (int edits = 1 + pick(3); edits--; ) {
                const size_t at = pick(changed.size());
                const char byte = kBytes[pick(kBytes.size())];
                switch (pick(3)) {
                    case 0: changed[at] = byte; break;
                    case 1: changed.insert(changed.begin() + at, byte); break;
                    default: if (changed.size() > 1) changed.erase(at, 1);
                }
            }
            const std::optional<Json::Value> changed_tree = lenient(changed);
            if (!changed_tree || names_a_member_twice(changed))
                continue;
            EXPECT_EQ(values_of(*changed_tree), count_lenient_json_values(changed, 0)) << changed;
            ++changed_checked;
        }
    }
    EXPECT_GT(checked, 1000u);
    EXPECT_GT(changed_checked, 1000u);
}

// The cap: a body of exactly the server's default count is read, one more value refused before
// the parse with the strict routes' error; commas in a comment or a string are no values, nor
// are commas where the reader stops (it makes one element of "[,,," and stops at its fault)
TEST(ServerRequestIntake, SearchBodiesAreRefusedPastTheCount) {
    std::string zeros = "{\"FASTA\": \">q\\nACGT\", \"x\": [0";
    // the root, FASTA, x and its elements
    for (size_t i = 4; i < 262144; ++i) {
        zeros += ",0";
    }
    zeros += "]}";
    EXPECT_TRUE(search_fits(zeros, 262144));
    EXPECT_EQ(262144u, values_of(*lenient(zeros)));
    zeros.insert(zeros.size() - 2, ",0");
    EXPECT_FALSE(search_fits(zeros, 262144));
    EXPECT_EQ("the request body holds more than 262144 JSON values (--max-request-json-values): "
              "it was not parsed and nothing of it ran; send fewer values (split the request)",
              refusal_of([&] { count_lenient_json_values(zeros, 262144); }));
    // 0: none
    EXPECT_TRUE(search_fits(zeros, 0));
    const std::string commas(1'000'000, ',');
    EXPECT_TRUE(search_fits("{\"FASTA\": \"" + commas + "\", /* " + commas + " */ \"a\": 1} "
                            + commas, 3));
    EXPECT_TRUE(search_fits("[" + commas + "]", 2));
    EXPECT_FALSE(search_fits("[" + commas + "]", 1));
}
