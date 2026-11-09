#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>

#include "cli/query.hpp"
#include "cli/server_checks.hpp"
#include "cli/server_envelope.hpp"
#include "cli/traverse.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/walker.hpp"
#include "tests/source_path.hpp"


// The multi-graph envelope's builders (server_envelope.hpp) against the reference samples of
// SPEC-multigraph-envelope.md (tests/data/multigraph_contract): each sample's envelope built
// again from its entries and its request-level parts

namespace {

using namespace mtg::cli;

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errors;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errors)) << errors;
    return v;
}

// A reference sample, {route, status, note, request?, body}
Json::Value sample(const std::string &name) {
    const std::string path = source_path("tests/data/multigraph_contract/" + name);
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.good()) << path << " not found";
    const std::string text { std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>() };
    return parse(text);
}

// |value| as compact text, members in the writer's (sorted) order: the comparison of two
// documents, since Json::Value's == tells a parsed integer (intValue) from the same number
// written as an unsigned one (uintValue)
std::string text(const Json::Value &value) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, value);
}

// The members of |value| other than |except|, for comparing an envelope without one of them
Json::Value without(Json::Value value, const std::vector<std::string> &except) {
    for (const std::string &key : except) {
        value.removeMember(key);
    }
    return value;
}

// days since 1970-01-01 of a civil date (H. Hinnant's days_from_civil), so that the test
// reads the samples' ISO-8601 instants without the platform's timegm
int64_t days_from_civil(int64_t y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int64_t era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + static_cast<int64_t>(doe) - 719468;
}

// "2026-10-31T09:20:58.000Z" as a time point; null as none
std::optional<std::chrono::system_clock::time_point> instant(const Json::Value &iso) {
    if (iso.isNull())
        return std::nullopt;
    const std::string s = iso.asString();
    EXPECT_EQ(24u, s.size()) << s;
    const int64_t days = days_from_civil(std::stoll(s.substr(0, 4)),
                                         static_cast<unsigned>(std::stoul(s.substr(5, 2))),
                                         static_cast<unsigned>(std::stoul(s.substr(8, 2))));
    const int64_t ms = ((days * 24 + std::stoll(s.substr(11, 2))) * 60
                            + std::stoll(s.substr(14, 2))) * 60'000
                     + std::stoll(s.substr(17, 2)) * 1000 + std::stoll(s.substr(20, 3));
    return std::chrono::system_clock::time_point(std::chrono::milliseconds(ms));
}

std::optional<uint64_t> ms(const Json::Value &v) {
    return v.isNull() ? std::nullopt : std::optional<uint64_t>(v.asUInt64());
}

SelectedPair pair_of_entry(const Json::Value &entry) {
    return { entry["graph"].asString(),
             { entry["graph_path"].asString(), entry["annotation_path"].asString() } };
}

// Each entry of |body|'s answers built again from its parts by the entry builders: the
// single-graph answer (the entry without its tags) or the refusal (without its http_status)
std::vector<Json::Value> rebuilt_entries(const Json::Value &body) {
    std::vector<Json::Value> entries;
    for (const Json::Value &entry : body["answers"]) {
        const SelectedPair pair = pair_of_entry(entry);
        const std::string fp = entry["index_fp"].isNull() ? "" : entry["index_fp"].asString();
        Json::Value rebuilt;
        if (entry["outcome"] == "answered") {
            rebuilt = envelope_entry_answered(pair, fp, without(entry, { "graph", "graph_path",
                                                                         "annotation_path",
                                                                         "index_fp",
                                                                         "outcome" }));
        } else {
            EXPECT_EQ("refused", entry["outcome"].asString());
            rebuilt = envelope_entry_refused(pair, fp, entry["refusal"]["http_status"].asInt(),
                                             without(entry["refusal"], { "http_status" }));
        }
        EXPECT_EQ(text(entry), text(rebuilt));
        entries.push_back(std::move(rebuilt));
    }
    return entries;
}

// The POST envelope of |body| built again from its entries and request-level parts
Json::Value rebuilt_post_envelope(const Json::Value &body) {
    std::vector<std::string> names;
    for (const Json::Value &name : body["graphs"]) {
        names.push_back(name.asString());
    }
    EnvelopeContext context;
    if (!body["request_id"].isNull())
        context.request_id = body["request_id"].asString();
    context.server_instance = body["server_instance"].asString();
    const Json::Value &d = body["deadline"];
    context.deadline = deadline_json(ms(d["not_after_ms"]), ms(d["finish_by_ms"]),
                                     instant(d["effective_deadline_at"]),
                                     instant(d["cutoff_at"]));
    context.usage = body.get("usage", Json::Value());
    return post_envelope(body["route"].asString(), names, rebuilt_entries(body), context,
                         body["timing"]["elapsed_ms"].asDouble());
}

} // namespace


// ---------------------------------------------------------------- the envelope's builders

// Contract F4 of the release-10 review: the envelope is written as every answer of the
// server is, its members (and an entry's, and deadline's) in byte order of their names, so the
// spec states that member order carries no meaning instead of an order the server never writes
TEST(ServerEnvelope, MemberOrderIsTheWritersAndTheSpecSaysSo) {
    const Json::Value body = sample("03_resolve_envelope_cutoff.json")["body"];
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    const std::string written = Json::writeString(builder, rebuilt_post_envelope(body));
    // the top-level members in the text, in their order (a scan of the compact text: the keys
    // at depth 1, strings skipped)
    std::vector<std::string> in_text;
    int depth = 0;
    bool key_next = false;
    for (size_t i = 0; i < written.size(); ++i) {
        const char c = written[i];
        if (c == '"') {
            size_t end = i + 1;
            while (written[end] != '"') {
                end += written[end] == '\\' ? 2 : 1;
            }
            if (depth == 1 && key_next)
                in_text.push_back(written.substr(i + 1, end - i - 1));
            key_next = false;
            i = end;
        } else if (c == '{' || c == '[') {
            ++depth;
            key_next = c == '{' && depth == 1;
        } else if (c == '}' || c == ']') {
            --depth;
        } else if (c == ',' && depth == 1) {
            key_next = true;
        }
    }
    ASSERT_EQ(body.size(), in_text.size());
    std::vector<std::string> sorted = in_text;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_EQ(sorted, in_text);
    EXPECT_NE("envelope_version", in_text.front());

    std::ifstream in(source_path("docs/SPEC-multigraph-envelope.md"), std::ios::binary);
    ASSERT_TRUE(in.good()) << "the spec not found";
    const std::string spec { std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>() };
    EXPECT_EQ(std::string::npos, spec.find("Members, in this order"));
    EXPECT_NE(std::string::npos, spec.find("the order of its members, of an entry's\nand of "
                                           "`deadline`'s carries no meaning"));
}

// Sample 01: two pairs answered, nothing unreached or unfinished, no deadline given
TEST(ServerEnvelope, Sample01ResolveTwoPairs) {
    const Json::Value body = sample("01_resolve_envelope_two_pairs.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(2u, built["answered"].asUInt());
    EXPECT_EQ(0u, built["refused"].asUInt());
    EXPECT_EQ(Json::Value(Json::arrayValue), built["not_reached"]);
    EXPECT_EQ(Json::Value(Json::arrayValue), built["unfinished"]);
    for (const char *key : { "not_after_ms", "finish_by_ms", "effective_deadline_at",
                             "cutoff_at" }) {
        EXPECT_TRUE(built["deadline"][key].isNull()) << key;
    }
    EXPECT_FALSE(built.isMember("usage"));
    EXPECT_FALSE(built.isMember("pattern_contract_version"));
    // the answered pairs' accounts summed
    EXPECT_EQ(39068u, built["work"]["memory_bytes"].asUInt64());
}

// Sample 02: a named graph, one entry, no request_id
TEST(ServerEnvelope, Sample02ResolveNamed) {
    const Json::Value body = sample("02_resolve_envelope_named.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(parse(R"(["bact_other"])"), built["graphs"]);
    EXPECT_TRUE(built["request_id"].isNull());
}

// Sample 04: the discovery cut's boundary tie (section 7) passes through the entry; the
// envelope's work.memory_bytes is the answered entry's
TEST(ServerEnvelope, Sample04ResolveTruncatedTie) {
    const Json::Value body = sample("04_resolve_labels_truncated_tie.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    const Json::Value &cut = built["answers"][0]["labels_truncated"];
    EXPECT_EQ((std::vector<std::string> { "dropped_full_length", "kept", "max_dropped_kmers",
                                          "min_kept_kmers", "tied_at_boundary", "tied_kept",
                                          "total" }), cut.getMemberNames());
    EXPECT_EQ(built["answers"][0]["work"]["memory_bytes"].asUInt64(),
              built["work"]["memory_bytes"].asUInt64());
}

// Sample 09: a /traverse entry whose seed result is labels_not_in_graph (section 2); the
// envelope's work.memory_bytes is the answered entry's (the seed result's shape is the
// engine's: LabelsNotInGraph.TheSampleStatesTheSeedResultsShape)
TEST(ServerEnvelope, Sample09TraverseLabelsNotInGraph) {
    const Json::Value body = sample("09_traverse_entry_labels_not_in_graph.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ("labels_not_in_graph",
              built["answers"][0]["results"][0]["outcome"]["walks"].asString());
    EXPECT_EQ(built["answers"][0]["work"]["memory_bytes"].asUInt64(),
              built["work"]["memory_bytes"].asUInt64());
}

// Sample 03: one pair answered, one not reached, one unfinished; the lists name the pairs
// from the refusals' codes, and the deadline's instants are written as the sample states them
TEST(ServerEnvelope, Sample03ResolveCutoff) {
    const Json::Value body = sample("03_resolve_envelope_cutoff.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(1u, built["answered"].asUInt());
    EXPECT_EQ(2u, built["refused"].asUInt());
    ASSERT_EQ(1u, built["not_reached"].size());
    EXPECT_EQ("bact_pseudomonadota", built["not_reached"][0]["graph"].asString());
    ASSERT_EQ(1u, built["unfinished"].size());
    EXPECT_EQ("inv_insecta", built["unfinished"][0]["graph"].asString());
    EXPECT_EQ("2026-10-31T09:20:57.750Z", built["deadline"]["cutoff_at"].asString());
    EXPECT_EQ(1793438460000u, built["deadline"]["finish_by_ms"].asUInt64());
}

// Sample 05: a named-graph /traverse, one entry; request_id is the attempt_id and the
// attempt's usage sits at the top level; a real run (2026-11-07): the deadline's not_after_ms is
// the request's (its timeline: ServerFanout.EachSampleDeadlineTellsOneTimeline). The seed's
// final_bytes is what its result holds: its own account (at most the result's work.memory_bytes),
// without the 512 MiB budget's cache allotments, which its peak keeps (E5b)
TEST(ServerEnvelope, Sample05TraverseNamed) {
    const Json::Value s05 = sample("05_traverse_envelope_named.json");
    const Json::Value &body = s05["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ("svc-trv-91aa07-1", built["request_id"].asString());
    EXPECT_EQ(s05["request"]["attempt_id"], built["request_id"]);
    EXPECT_EQ(text(body["usage"]), text(built["usage"]));
    EXPECT_FALSE(built["answers"][0].isMember("usage"));
    EXPECT_EQ(1794089717000u, built["deadline"]["not_after_ms"].asUInt64());
    EXPECT_EQ(s05["request"]["not_after_ms"].asUInt64(), built["deadline"]["not_after_ms"].asUInt64());
    EXPECT_TRUE(built["deadline"]["finish_by_ms"].isNull());
    const Json::Value &seed = body["usage"]["per_seed"][0];
    const uint64_t own = body["answers"][0]["results"][0]["work"]["memory_bytes"].asUInt64();
    const uint64_t budget = s05["request"]["strategy"]["bounds"]["max_memory_mb"].asUInt64() << 20;
    EXPECT_EQ(639358u, seed["final_bytes"].asUInt64());
    EXPECT_LE(seed["final_bytes"].asUInt64(), own);
    EXPECT_EQ(own + mtg::graph::traversal::memory_allotments(budget),
              seed["peak_admitted_bytes"].asUInt64());
    EXPECT_EQ(budget, body["usage"]["memory"]["held_bound_bytes"].asUInt64());
}

// Sample 23 (E5b): a named /traverse whose third seed the request-wide memory bound left
// unstarted, from a real run: the envelope built again from its entry, and the unstarted
// result's statements consistent with the usage (H, the final_bytes of the seeds before it; B,
// the effective per-seed budget; R, the bound)
TEST(ServerEnvelope, Sample23TraverseRequestMemoryBound) {
    const Json::Value s23 = sample("23_traverse_request_memory_bound.json");
    const Json::Value &body = s23["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(s23["request"]["attempt_id"], built["request_id"]);
    EXPECT_EQ(text(body["usage"]), text(built["usage"]));
    ASSERT_EQ(1u, body["answers"].size());
    const Json::Value &entry = body["answers"][0];
    ASSERT_EQ(3u, entry["results"].size());
    const Json::Value &u = body["usage"];
    const uint64_t held = u["per_seed"][0]["final_bytes"].asUInt64()
                        + u["per_seed"][1]["final_bytes"].asUInt64();
    const uint64_t budget = entry["strategy"]["bounds"]["max_memory_mb"].asUInt64() << 20;
    const Json::Value &unstarted = entry["results"][2];
    const Json::Value &q = unstarted["resource_stop"];
    const uint64_t bound_mb = q["effective"].asUInt64();
    EXPECT_GT(held + budget, bound_mb << 20);
    EXPECT_LE(u["per_seed"][0]["final_bytes"].asUInt64() + budget, bound_mb << 20);
    EXPECT_FALSE(fits_request_memory(held, budget, bound_mb));
    EXPECT_TRUE(fits_request_memory(u["per_seed"][0]["final_bytes"].asUInt64(), budget, bound_mb));
    EXPECT_EQ("not started: the results already held (" + std::to_string(held) + " bytes) and "
              "this seed's memory budget (" + std::to_string(budget) + " bytes) would pass the "
              "server's request memory bound (" + std::to_string(bound_mb) + " MiB); no budget "
              "of the seed ran out: send it in a new request", unstarted["error"].asString());
    EXPECT_EQ("failed", unstarted["outcome"]["walks"].asString());
    EXPECT_EQ("request", q["scope"].asString());
    EXPECT_EQ("memory", q["resource"].asString());
    EXPECT_EQ("not_started", q["phase"].asString());
    EXPECT_EQ(bound_mb, q["requested"].asUInt64());
    EXPECT_EQ((held + (1 << 20) - 1) >> 20, q["used"].asUInt64());
    EXPECT_EQ("split_request", q["actions"][0].asString());
    EXPECT_EQ("max_request_memory_mb", unstarted["limitations"][0]["knob"].asString());
    EXPECT_EQ("memory_bound_soft", unstarted["limitations"][1]["kind"].asString());
    EXPECT_TRUE(unstarted["unknown_labels"].isArray());
    EXPECT_EQ("not_started", u["per_seed"][2]["outcome"].asString());
    EXPECT_EQ("request_memory", u["per_seed"][2]["stopped_by"].asString());
    EXPECT_EQ(unstarted["work"]["memory_bytes"], u["per_seed"][2]["final_bytes"]);
    // each walked result holds its own account: the budget's cache allotments ended with its
    // walk (its peak keeps them); E5b's run before the fix stated 1943996 bytes held here
    EXPECT_EQ(1288636u, held);
    for (Json::ArrayIndex i = 0; i < 2; ++i) {
        const uint64_t own = entry["results"][i]["work"]["memory_bytes"].asUInt64();
        EXPECT_LE(u["per_seed"][i]["final_bytes"].asUInt64(), own) << i;
        EXPECT_EQ(own + mtg::graph::traversal::memory_allotments(budget),
                  u["per_seed"][i]["peak_admitted_bytes"].asUInt64()) << i;
    }
    EXPECT_EQ(2u, u["seeds"]["started"].asUInt64());
    EXPECT_EQ(3u, u["seeds"]["requested"].asUInt64());
    for (Json::ArrayIndex i = 0; i < 2; ++i) {
        EXPECT_TRUE(entry["results"][i].isMember("arms")) << i;
        EXPECT_EQ(u["per_seed"][i]["outcome"], entry["results"][i]["outcome"]["walks"]) << i;
    }
}

// Sample 16: /pattern keeps pattern_contract_version beside the envelope's members; a pair
// not started by not_after_ms is not_reached; work.memory_bytes the answered entry's (0)
TEST(ServerEnvelope, Sample16PatternNotReached) {
    const Json::Value body = sample("16_pattern_envelope_not_reached.json")["body"];
    const Json::Value built = rebuilt_post_envelope(body);
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(0u, built["work"]["memory_bytes"].asUInt64());
    EXPECT_EQ(1, built["pattern_contract_version"].asInt());
    EXPECT_EQ(1u, built["answered"].asUInt());
    EXPECT_EQ(1u, built["refused"].asUInt());
    ASSERT_EQ(1u, built["not_reached"].size());
    EXPECT_EQ(Json::Value(Json::arrayValue), built["unfinished"]);
    EXPECT_EQ("svc-pat-0c33d1-1", built["request_id"].asString());
}

// Sample 10: the capabilities envelope over every pair: no request-level members beyond
// server_instance, the names taken from the entries
TEST(ServerEnvelope, Sample10TraverseCapabilities) {
    const Json::Value body = sample("10_traverse_capabilities_envelope.json")["body"];
    const Json::Value built = get_envelope(body["route"].asString(), rebuilt_entries(body),
                                           body["server_instance"].asString());
    EXPECT_EQ(text(body), text(built));
    for (const char *key : { "not_reached", "unfinished", "request_id", "deadline", "usage",
                             "timing", "work" }) {
        EXPECT_FALSE(built.isMember(key)) << key;
    }
}

TEST(ServerEnvelope, RequestLevelMembers) {
    const SelectedPair a { "a", { "a.dbg", "a.anno" } };
    const SelectedPair b { "b", { "b.dbg", "b.anno" } };
    Json::Value answer;
    answer["work"]["memory_bytes"] = 10;
    Json::Value unfinished;
    unfinished["error"] = "cut";
    unfinished["code"] = "unfinished";
    unfinished["reason"] = "cutoff";
    Json::Value failed;
    failed["error"] = "Label not found";
    EnvelopeContext context;
    context.server_instance = "s1";
    const Json::Value out = post_envelope(
            "/resolve", { "a", "b" },
            { envelope_entry_answered(a, "", answer),
              envelope_entry_refused(b, "fp", 503, unfinished),
              envelope_entry_refused(a, "", 400, failed) },
            context, 1.5);
    EXPECT_EQ(1, out["envelope_version"].asInt());
    EXPECT_EQ("/resolve", out["route"].asString());
    EXPECT_EQ(1u, out["answered"].asUInt());
    EXPECT_EQ(2u, out["refused"].asUInt());
    EXPECT_EQ(Json::Value(Json::arrayValue), out["not_reached"]);
    ASSERT_EQ(1u, out["unfinished"].size());
    EXPECT_EQ(envelope_pair_ref(b), out["unfinished"][0]);
    // no request_id, no deadline given: null and every deadline member null; no usage
    EXPECT_TRUE(out.isMember("request_id"));
    EXPECT_TRUE(out["request_id"].isNull());
    EXPECT_EQ(deadline_json(std::nullopt, std::nullopt, std::nullopt, std::nullopt),
              out["deadline"]);
    EXPECT_FALSE(out.isMember("usage"));
    EXPECT_EQ(10u, out["work"]["memory_bytes"].asUInt64());
    EXPECT_EQ(1.5, out["timing"]["elapsed_ms"].asDouble());
    // the entries' tags: index_fp null without a manifest
    EXPECT_TRUE(out["answers"][0]["index_fp"].isNull());
    EXPECT_EQ("fp", out["answers"][1]["index_fp"].asString());
    EXPECT_EQ(503, out["answers"][1]["refusal"]["http_status"].asInt());
    EXPECT_EQ("Label not found", out["answers"][2]["refusal"]["error"].asString());

    const Json::Value refs = envelope_pair_ref(a);
    EXPECT_EQ(parse(R"({"graph": "a", "graph_path": "a.dbg", "annotation_path": "a.anno"})"),
              refs);

    // ISO-8601 UTC with milliseconds (1793438457750 ms since the epoch is that instant)
    const auto t = std::chrono::system_clock::time_point(std::chrono::milliseconds(
            int64_t(1793438457750)));
    const Json::Value d = deadline_json(1u, std::nullopt, t, std::nullopt);
    EXPECT_EQ(1u, d["not_after_ms"].asUInt64());
    EXPECT_TRUE(d["finish_by_ms"].isNull());
    EXPECT_EQ("2026-10-31T09:20:57.750Z", d["effective_deadline_at"].asString());
    EXPECT_TRUE(d["cutoff_at"].isNull());

    // get_envelope: each name once, in the entries' order
    const Json::Value caps = get_envelope(
            "/pattern/capabilities",
            { envelope_entry_answered(a, "", Json::Value(Json::objectValue)),
              envelope_entry_answered(b, "", Json::Value(Json::objectValue)),
              envelope_entry_answered({ "b", { "b2.dbg", "b2.anno" } }, "",
                                      Json::Value(Json::objectValue)) },
            "s1");
    EXPECT_EQ(parse(R"(["a", "b"])"), caps["graphs"]);
    EXPECT_EQ(3u, caps["answered"].asUInt());
    EXPECT_EQ(0u, caps["refused"].asUInt());
}

// /pattern's entries are the shared envelope's, member for member
TEST(ServerEnvelope, PatternEntriesAreTheEnvelopes) {
    Json::Value answer = parse(R"({"patterns": [], "work": {"memory_bytes": 3}})");
    const SelectedPair p { "n", { "g.dbg", "g.anno" } };
    EXPECT_EQ(envelope_entry_answered(p, "fp", answer),
              pattern_pair_answered("n", "g.dbg", "g.anno", "fp", answer));
    EXPECT_EQ(envelope_entry_answered(p, "", answer),
              pattern_pair_answered("n", "g.dbg", "g.anno", "", answer));
    const Json::Value body = parse(R"({"error": "no", "code": "deadline"})");
    EXPECT_EQ(envelope_entry_refused(p, "fp", 503, body),
              pattern_pair_refused("n", "g.dbg", "g.anno", "fp", 503, body));
    // /pattern's envelope keeps its members until its route moves to post_envelope
    const Json::Value out = pattern_envelope(
            { "n" }, { pattern_pair_answered("n", "g.dbg", "g.anno", "fp", answer) }, 2.0);
    EXPECT_FALSE(out.isMember("envelope_version"));
    EXPECT_EQ(1u, out["answered"].asUInt());
}

// Sample 11: the pattern capabilities in the envelope, two pairs
TEST(ServerEnvelope, Sample11PatternCapabilities) {
    const Json::Value body = sample("11_pattern_capabilities_envelope.json")["body"];
    const Json::Value built = get_envelope(body["route"].asString(), rebuilt_entries(body),
                                           body["server_instance"].asString());
    EXPECT_EQ(text(body), text(built));
    EXPECT_EQ(2u, built["answered"].asUInt());
}


// ---------------------------------------------------------------- GET /capabilities (section 8)

// composite_fp: sha256 over the pairs' index_fp sorted as bytes, joined by "\n" without a
// trailing newline; null when any pair has none
TEST(ServerEnvelope, CompositeRule) {
    // the sha256 of "a\nb" (`printf 'a\nb' | sha256sum`), whichever order the pairs come in
    EXPECT_EQ("7e18f737311b2dc3b2f269dd78396b0351f14fb66efa879f768cb23181883c78",
              composite_fp_json({ "b", "a" }).asString());
    EXPECT_EQ(composite_fp_json({ "a", "b" }), composite_fp_json({ "b", "a" }));
    // one pair: its index_fp hashed alone, no separator
    EXPECT_EQ(sha256_hex("x"), composite_fp_json({ "x" }).asString());
    // sorted as bytes: "B" (0x42) before "a" (0x61)
    EXPECT_EQ("52958ed70d5d066708f7e77ec157d5df490e3f1a597fa14d0d41b513824564b6",
              composite_fp_json({ "a", "B" }).asString());
    // a pair without a manifest, or none at all: null
    EXPECT_TRUE(composite_fp_json({ "a", "" }).isNull());
    EXPECT_TRUE(composite_fp_json({}).isNull());

    // sample 12: the composite of its eleven pairs (graph_summary), with the rule beside it
    const Json::Value body = sample("12_capabilities_multigraph.json")["body"];
    std::vector<std::string> fps;
    for (const Json::Value &pair : body["graph_summary"]["pairs"]) {
        fps.push_back(pair["index_fp"].asString());
    }
    EXPECT_EQ(11u, fps.size());
    EXPECT_EQ(body["composite_fp"].asString(), composite_fp_json(fps).asString());
    EXPECT_EQ(body["composite_rule"].asString(), std::string(kCompositeRule));
    // the rule a client recomputes composite_fp by states that a pair counts once (section 8,
    // contract notice B1; contract F5 of the release-10 review): graph_summary.pairs lists a
    // pair once per name that holds it
    EXPECT_NE(std::string::npos, std::string(kCompositeRule).find(
            "each loaded (graph_path, annotation_path) pair counted once, whatever names list "
            "it"));
    std::vector<std::string> listed_twice = fps;
    listed_twice.push_back(fps.front());
    std::sort(listed_twice.begin(), listed_twice.end());
    listed_twice.erase(std::unique(listed_twice.begin(), listed_twice.end()), listed_twice.end());
    EXPECT_EQ(body["composite_fp"].asString(), composite_fp_json(listed_twice).asString());
}

// The multigraph block and routes.requests are sample 12's (its server: --max-request-body-mb
// 256, --max-request-json-values, --max-column-labels-bytes and --max-connections at their
// defaults, as staging's)
TEST(ServerEnvelope, Sample12MultigraphBlock) {
    const Json::Value body = sample("12_capabilities_multigraph.json")["body"];
    const Json::Value intake = request_intake_json(256, 262144, 0);
    EXPECT_TRUE(body["multigraph"]["request_intake"]["max_connections"].isNull());
    // GET /column_labels (H3): the cap at its default, the list of the sample's pairs
    const Json::Value &labels_total = body["multigraph"]["column_labels"]["total"];
    ASSERT_TRUE(labels_total.isUInt64());
    const Json::Value column_labels
            = column_labels_capabilities_json(16 << 20, labels_total.asUInt64());
    EXPECT_EQ(text(body["multigraph"]),
              text(multigraph_capabilities_json(0, {}, intake, column_labels)));
    // request_intake (H1): the byte cap as the top-level member states it
    EXPECT_EQ(text(body["max_request_body_mb"]),
              text(body["multigraph"]["request_intake"]["max_request_body_mb"]));
    EXPECT_EQ(body["routes"]["requests"].asString(), std::string(kRequestsRoute));
    EXPECT_EQ(body["multigraph"]["requests_route"].asString(), std::string(kRequestsRoute));
    // the request-wide memory bound of /traverse (E5b): 0 when off, as the sample's server,
    // else the flag's value, an integer
    EXPECT_TRUE(body["multigraph"]["max_request_memory_mb"].isUInt64());
    EXPECT_EQ(0u, body["multigraph"]["max_request_memory_mb"].asUInt64());
    const Json::Value set = multigraph_capabilities_json(8192, {}, intake, column_labels);
    EXPECT_TRUE(set["max_request_memory_mb"].isIntegral());
    EXPECT_EQ("8192", text(set["max_request_memory_mb"]));
    Json::Value rest = set;
    rest["max_request_memory_mb"] = 0;
    EXPECT_EQ(text(multigraph_capabilities_json(0, {}, intake, column_labels)), text(rest));
}

// Sample 10: every /traverse capabilities probe states the request-wide memory bound beside
// the per-seed maximum (E5b), 0 when off, an integer
TEST(ServerEnvelope, Sample10StatesTheRequestMemoryBound) {
    const Json::Value body = sample("10_traverse_capabilities_envelope.json")["body"];
    ASSERT_EQ(2u, body["answers"].size());
    for (const Json::Value &entry : body["answers"]) {
        EXPECT_TRUE(entry["max_memory_mb"].isUInt64());
        EXPECT_TRUE(entry["max_request_memory_mb"].isUInt64());
        EXPECT_EQ(0u, entry["max_request_memory_mb"].asUInt64());
    }
}

// Sample 25 (SPEC section 3, kmers_in_target; train 2): a coordinate query on a multi-graph server,
// a real run. Each pair's part is what the server writes for it (process_search_request): its
// rows SeqSearchResult::to_json of the coordinate result (label, k-mer count, the coordinates per
// query k-mer, the target's k-mer count or 0 without a record mapping) tagged with the pair, and
// its cut cut_to_json with the pair's tags and index_fp. A row of a pair with a .seqs states
// kmers_in_target, its record's k-mer count, every coordinate in the record below it; a pair
// without one (whole-column labels) states none. kmer_count is the query k-mers with a coordinate
TEST(ServerEnvelope, Sample25SearchCoordinateRowsAreTheWriters) {
    const Json::Value s25 = sample("25_search_coordinate_rows.json");
    ASSERT_EQ("POST /search", s25["route"].asString());
    ASSERT_EQ(200, s25["status"].asInt());
    const Json::Value &request = s25["request"];
    EXPECT_TRUE(request["query_coords"].asBool());
    const size_t n = request["top_labels"].asUInt64();
    const Json::Value &body = s25["body"];
    ASSERT_EQ(1u, body.size());
    const Json::Value &query = body[0];
    const Json::Value &rows = query["results"];
    const Json::Value &cuts = query["cuts"];
    // one cut per selected pair, in pair order (each name here one pair)
    ASSERT_EQ(request["graphs"].size(), cuts.size());
    // the coordinates of a row's kmer_coords per query k-mer: runs "position-first[-last]"
    // (positions position.., coordinates first..last), joined by ':'
    auto coordinates = [](const std::string &runs) {
        std::vector<SmallVector<uint64_t>> tuples;
        std::stringstream in(runs);
        for (std::string run; std::getline(in, run, ':'); ) {
            std::vector<uint64_t> v;
            std::stringstream parts(run);
            for (std::string x; std::getline(parts, x, '-'); ) {
                v.push_back(std::stoull(x));
            }
            EXPECT_TRUE(v.size() == 2 || v.size() == 3) << run;
            const uint64_t first = v.at(1);
            const uint64_t last = v.size() == 3 ? v[2] : first;
            for (uint64_t c = first; c <= last; ++c) {
                const size_t position = v[0] + (c - first);
                if (tuples.size() <= position)
                    tuples.resize(position + 1);
                tuples[position].push_back(c);
            }
        }
        for (auto &t : tuples) {
            std::sort(t.begin(), t.end());
        }
        return tuples;
    };
    const std::vector<std::string> tags { "graph", "graph_path", "annotation_path" };
    Json::ArrayIndex next = 0;
    size_t with_target = 0;
    size_t without_target = 0;
    for (const Json::Value &cut : cuts) {
        // this pair's rows: the next ones, tagged with it
        SeqSearchResult::LabelCountCoordsVec result;
        std::vector<Json::Value> stated;
        while (next < rows.size() && rows[next]["graph"] == cut["graph"]) {
            const Json::Value &row = rows[next++];
            for (const std::string &tag : tags) {
                EXPECT_EQ(cut[tag], row[tag]) << tag;
            }
            auto tuples = coordinates(row["kmer_coords"].asString());
            const uint64_t count = row["kmer_count"].asUInt64();
            EXPECT_EQ(count, static_cast<uint64_t>(std::count_if(tuples.begin(), tuples.end(),
                    [](const auto &t) { return !t.empty(); })));
            uint64_t in_target = 0;
            if (row.isMember("kmers_in_target")) {
                in_target = row["kmers_in_target"].asUInt64();
                EXPECT_LE(count, in_target);
                for (const auto &t : tuples) {
                    for (uint64_t c : t) {
                        EXPECT_LT(c, in_target) << "a coordinate past its record";
                    }
                }
                ++with_target;
            } else {
                ++without_target;
            }
            result.emplace_back(row["sample"].asString(), count, std::move(tuples), in_target);
            stated.push_back(without(row, tags));
        }
        ASSERT_EQ(stated.size(), cut["returned"].asUInt64()) << cut["graph"];
        mtg::graph::SearchCut c;
        c.returned = cut["returned"].asUInt64();
        c.complete = cut["complete"].asBool();
        c.nth_score = cut["nth_score"].isNull() ? 0 : cut["nth_score"].asUInt64();
        c.tied = cut["tied"].asUInt64();
        c.tied_returned = cut["tied_returned"].asUInt64();
        const SeqSearchResult written(QuerySequence { 0, query["seq_description"].asString(), "" },
                                      std::move(result), {}, c);
        const Json::Value json = written.to_json(false, 31);
        EXPECT_EQ(query["seq_description"], json["seq_description"]);
        ASSERT_EQ(stated.size(), json["results"].size());
        for (size_t i = 0; i < stated.size(); ++i) {
            EXPECT_EQ(text(stated[i]), text(json["results"][static_cast<Json::ArrayIndex>(i)])) << i;
        }
        EXPECT_EQ(text(written.cut_to_json()), text(without(cut, { "graph", "graph_path",
                                                                   "annotation_path", "index_fp" })));
        EXPECT_TRUE(cut["index_fp"].isString());
        EXPECT_LE(c.returned, n);
    }
    // every row belongs to a pair, in pair order
    EXPECT_EQ(rows.size(), next);
    // one pair of each kind: records (kmers_in_target 360 - 31 + 1), whole columns (none)
    EXPECT_EQ(2u, with_target);
    EXPECT_EQ(2u, without_target);
    EXPECT_EQ(330u, rows[0]["kmers_in_target"].asUInt64());
    EXPECT_FALSE(cuts[0]["complete"].asBool());
    EXPECT_TRUE(cuts[1]["complete"].asBool());
}
