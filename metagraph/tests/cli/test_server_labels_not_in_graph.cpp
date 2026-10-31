#include "gtest/gtest.h"

#include <fstream>
#include <iterator>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <json/json.h>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "cli/traverse.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"


// Labels a pair of a multi-graph server does not hold (SPEC-multigraph-envelope.md section 2):
// at the three sites that resolve names — /resolve's explicit labels, a traverse seed's labels
// and a traverse's extra labels — such a label is dropped and reported when the server is a
// multi-graph one (TraverseLimits::labels_not_in_graph_dropped, process_resolve_request's
// switch), and refuses the request as before on a single-graph server. The index here plays one
// pair: "Z" and "Y" stand for labels of the other pairs

namespace {

using namespace mtg;
using namespace mtg::cli;
using mtg::graph::DBGSuccinct;
using mtg::graph::DeBruijnGraph;
using mtg::graph::AnnotatedDBG;

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errors;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errors)) << errors;
    return v;
}

std::string random_seq(size_t length, std::mt19937 *rng) {
    std::string s(length, 'A');
    for (char &c : s) {
        c = "ACGT"[(*rng)() % 4];
    }
    return s;
}

// The pair: two records, labelled A and B (columns)
struct Pair {
    static constexpr size_t k = 11;
    std::vector<std::string> records;
    std::unique_ptr<AnnotatedDBG> index;

    Pair() {
        std::mt19937 rng(17);
        records = { random_seq(90, &rng), random_seq(90, &rng) };
        index = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                k, records, { "A", "B" }, DeBruijnGraph::BASIC);
    }
    // a seed inside record A, and one no record holds
    std::string seed() const { return records[0].substr(20, 30); }
};

// a traverse request of |seeds|, each {sequence, labels}, walking right
Json::Value traverse_request(const std::vector<std::pair<std::string, std::vector<std::string>>> &seeds,
                             const std::string &detail = "full") {
    Json::Value r = parse(R"({"strategy": {"direction": "right", "bounds":
                              {"max_extension_bp": 20}, "output": {"timing": false}}})");
    r["strategy"]["output"]["detail"] = detail;
    for (size_t i = 0; i < seeds.size(); ++i) {
        Json::Value s;
        s["sequence"] = seeds[i].first;
        s["seed_id"] = "s" + std::to_string(i);
        for (const std::string &label : seeds[i].second) {
            s["labels"].append(label);
        }
        r["seeds"].append(s);
    }
    return r;
}

TraverseLimits multi_graph() {
    TraverseLimits limits;
    limits.not_in_graph_per_seed = true;
    limits.labels_not_in_graph_dropped = true;
    return limits;
}

// the message of |f|'s std::invalid_argument, "" when it throws none
template <class F>
std::string refusal_of(const F &f) {
    try {
        f();
    } catch (const std::invalid_argument &e) {
        return e.what();
    }
    return "";
}

Json::Value strings(const std::vector<std::string> &values) {
    Json::Value a(Json::arrayValue);
    for (const std::string &v : values) {
        a.append(v);
    }
    return a;
}

// a walked result without what tells the requests apart: their labels, the reports of the
// absent ones, the work they did resolving them
Json::Value walk_of(Json::Value result) {
    result.removeMember("unknown_labels");
    result.removeMember("work");
    result.removeMember("annotation");
    result["seed"].removeMember("dropped_labels");
    return result;
}

} // namespace

// A seed label the pair does not hold: dropped (not_in_graph, no presence runs) and listed in
// the seed result's unknown_labels; the seed is walked under the others, as named alone; the
// dropped labels come in the request's order; on a single-graph server the request fails
TEST(LabelsNotInGraph, SeedLabelsAreDroppedAndReported) {
    const Pair pair;
    const Json::Value with_absent = traverse_request({ { pair.seed(), { "Z", "A", "Z2" } } });
    EXPECT_EQ("seed 's0': Label not found in the annotation: 'Z' (no sequence "
              "header index loaded)",
              refusal_of([&] { process_traverse_request(with_absent, *pair.index, ""); }));

    const Json::Value out = process_traverse_request(with_absent, *pair.index, "", multi_graph());
    ASSERT_EQ(1u, out["results"].size());
    const Json::Value &r = out["results"][0];
    EXPECT_EQ("complete", r["outcome"]["walks"].asString()) << r;
    EXPECT_EQ(strings({ "Z", "Z2" }), r["unknown_labels"]);
    EXPECT_EQ(strings({ "A" }), r["seed"]["labels"]);
    EXPECT_EQ(parse(R"([{"label": "Z", "reason": "not_in_graph", "runs": [],
                         "runs_kind": "presence"},
                        {"label": "Z2", "reason": "not_in_graph", "runs": [],
                         "runs_kind": "presence"}])"), r["seed"]["dropped_labels"]);
    // the walk is the one of the label the pair holds
    const Json::Value alone = process_traverse_request(
            traverse_request({ { pair.seed(), { "A" } } }), *pair.index, "", multi_graph());
    EXPECT_EQ(Json::Value(Json::arrayValue), alone["results"][0]["unknown_labels"]);
    EXPECT_EQ(walk_of(alone["results"][0]), walk_of(r));
    // no extra label: the answer's own list is empty, and present
    EXPECT_EQ(Json::Value(Json::arrayValue), out["unknown_labels"]);
    // a single-graph server's answer has neither
    const Json::Value single = process_traverse_request(
            traverse_request({ { pair.seed(), { "A" } } }), *pair.index, "");
    EXPECT_FALSE(single.isMember("unknown_labels"));
    EXPECT_FALSE(single["results"][0].isMember("unknown_labels"));

    // B is in the annotation but does not support the seed: both drops, in the request's order
    const Json::Value mixed = process_traverse_request(
            traverse_request({ { pair.seed(), { "B", "Z", "A", "Y" } } }), *pair.index, "",
            multi_graph());
    const Json::Value &dropped = mixed["results"][0]["seed"]["dropped_labels"];
    ASSERT_EQ(3u, dropped.size()) << dropped;
    EXPECT_EQ("B", dropped[0]["label"].asString());
    EXPECT_EQ("seed_unsupported", dropped[0]["reason"].asString());
    EXPECT_EQ("Z", dropped[1]["label"].asString());
    EXPECT_EQ("not_in_graph", dropped[1]["reason"].asString());
    EXPECT_EQ("Y", dropped[2]["label"].asString());
    EXPECT_EQ(strings({ "Z", "Y" }), mixed["results"][0]["unknown_labels"]);
    // a name given twice is refused as before, absent or not
    EXPECT_NE(std::string::npos, refusal_of([&] {
        process_traverse_request(traverse_request({ { pair.seed(), { "Z", "A", "Z" } } }),
                                 *pair.index, "", multi_graph());
    }).find("Duplicate seed label 'Z'"));
}

// A seed whose labels the pair holds none of is not walked: outcome.walks labels_not_in_graph,
// the names in unknown_labels and in dropped_labels, in the shape of a seed not in the graph;
// the other seeds are walked as alone. Under a memory budget it states memory_bound_soft; a
// ledger-managed request's usage states the outcome
TEST(LabelsNotInGraph, ASeedOfAbsentLabelsIsNotWalked) {
    const Pair pair;
    const Json::Value request = traverse_request({ { pair.seed(), { "Z", "Y" } },
                                                   { pair.seed(), { "A" } } });
    EXPECT_EQ("seed 's0': Label not found in the annotation: 'Z' (no sequence "
              "header index loaded)",
              refusal_of([&] { process_traverse_request(request, *pair.index, ""); }));
    const Json::Value out = process_traverse_request(request, *pair.index, "", multi_graph());
    ASSERT_EQ(2u, out["results"].size());
    const Json::Value &r = out["results"][0];
    EXPECT_EQ("labels_not_in_graph", r["outcome"]["walks"].asString()) << r;
    EXPECT_EQ("none of the seed's labels is in this graph: Z, Y", r["error"].asString());
    EXPECT_EQ(strings({ "Z", "Y" }), r["unknown_labels"]);
    EXPECT_EQ(parse(R"([{"label": "Z", "reason": "not_in_graph"},
                        {"label": "Y", "reason": "not_in_graph"}])"), r["dropped_labels"]);
    EXPECT_EQ("s0", r["seed"]["seed_id"].asString());
    EXPECT_FALSE(r["seed"]["labels_from_seed"].asBool());
    EXPECT_FALSE(r.isMember("arms"));
    EXPECT_EQ(0u, r["limitations"].size());
    EXPECT_GT(r["work"]["memory_bytes"].asUInt64(), 0u);
    // the shape of a seed not in the graph: the same members, but its own
    const Json::Value absent_seed = process_traverse_request(
            traverse_request({ { "ACGTACGTACGTACGTACGT", {} } }), *pair.index, "",
            multi_graph())["results"][0];
    ASSERT_EQ("not_in_graph", absent_seed["outcome"]["walks"].asString()) << absent_seed;
    for (const std::string &member : absent_seed.getMemberNames()) {
        if (member != "not_in_graph") {
            EXPECT_TRUE(r.isMember(member)) << member;
        }
    }
    EXPECT_EQ(absent_seed["outcome"].getMemberNames(), r["outcome"].getMemberNames());
    // the other seed is walked as alone
    const Json::Value alone = process_traverse_request(
            traverse_request({ { pair.seed(), { "A" } } }), *pair.index, "", multi_graph());
    Json::Value second = out["results"][1];
    second["seed"]["seed_id"] = "s0";
    EXPECT_EQ(walk_of(alone["results"][0]), walk_of(second));

    // under a memory budget: memory_bound_soft, as every result there; with coordinates asked
    // for, none, and why
    Json::Value budget = traverse_request({ { pair.seed(), { "Z" } } });
    budget["strategy"]["bounds"]["max_memory_mb"] = 64;
    budget["strategy"]["output"]["coordinates"] = true;
    const Json::Value budgeted = process_traverse_request(budget, *pair.index, "",
                                                          multi_graph())["results"][0];
    ASSERT_EQ(1u, budgeted["limitations"].size()) << budgeted;
    EXPECT_EQ("memory_bound_soft", budgeted["limitations"][0]["kind"].asString());
    EXPECT_TRUE(budgeted["coordinates"].isNull());
    EXPECT_TRUE(budgeted.isMember("coordinates_reason"));

    // the graphlet detail: no graphlet for it, as for a seed not in the graph
    const Json::Value g = process_traverse_request(
            traverse_request({ { pair.seed(), { "Z" } } }, "graphlet"), *pair.index, "",
            multi_graph());
    EXPECT_FALSE(g["results"][0].isMember("graphlet"));
    EXPECT_EQ("labels_not_in_graph", g["results"][0]["outcome"]["walks"].asString());

    // a ledger-managed request: its usage states the seed's outcome
    Json::Value managed = request;
    managed["attempt_id"] = "labels-1";
    const Json::Value usage = process_traverse_request(managed, *pair.index, "",
                                                       multi_graph())["usage"];
    EXPECT_EQ("labels_not_in_graph", usage["per_seed"][0]["outcome"].asString()) << usage;
    EXPECT_EQ("complete", usage["per_seed"][1]["outcome"].asString()) << usage;

    // a seed not in the graph states its absent labels too
    const Json::Value both = process_traverse_request(
            traverse_request({ { "ACGTACGTACGTACGTACGT", { "Z", "A" } } }), *pair.index, "",
            multi_graph())["results"][0];
    EXPECT_EQ("not_in_graph", both["outcome"]["walks"].asString());
    EXPECT_EQ(strings({ "Z" }), both["unknown_labels"]);
}

// The reference sample 09 (tests/data/multigraph_contract, SPEC-multigraph-envelope.md) states
// the seed result of labels_not_in_graph in the engine's shape: the members of the result, of
// its seed, outcome and limitations, under the sample request's memory budget and without
// coordinates asked for
TEST(LabelsNotInGraph, TheSampleStatesTheSeedResultsShape) {
    // run from the build directory, like every test reading ../tests/data
    std::ifstream in("../tests/data/multigraph_contract/09_traverse_entry_labels_not_in_graph.json",
                     std::ios::binary);
    ASSERT_TRUE(in.good()) << "sample 09 not found (run from the build directory)";
    const std::string text { std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>() };
    const Json::Value sample = parse(text)["body"]["answers"][0]["results"][0];
    const Pair pair;
    Json::Value request = traverse_request({ { pair.seed(), { "Z" } } });
    request["strategy"]["bounds"]["max_memory_mb"] = 512;
    const Json::Value r = process_traverse_request(request, *pair.index, "",
                                                   multi_graph())["results"][0];
    ASSERT_EQ("labels_not_in_graph", r["outcome"]["walks"].asString()) << r;
    EXPECT_EQ(sample.getMemberNames(), r.getMemberNames());
    for (const char *member : { "seed", "outcome" }) {
        EXPECT_EQ(sample[member].getMemberNames(), r[member].getMemberNames()) << member;
    }
    EXPECT_EQ(sample["outcome"], r["outcome"]);
    ASSERT_EQ(1u, r["limitations"].size());
    EXPECT_EQ(sample["limitations"][0].getMemberNames(), r["limitations"][0].getMemberNames());
    // (its effect's words depend on the annotation: a column one here, row-diff in the sample)
    EXPECT_EQ(sample["limitations"][0]["kind"], r["limitations"][0]["kind"]);
    EXPECT_EQ(sample["dropped_labels"][0].getMemberNames(),
              r["dropped_labels"][0].getMemberNames());
}

// The graphlet states a dropped label as it states every dropped seed label: an X record with
// the reason not_in_graph and no runs; the JSON beside it lists the unknown labels
TEST(LabelsNotInGraph, TheGraphletStatesTheDrop) {
    const Pair pair;
    const Json::Value out = process_traverse_request(
            traverse_request({ { pair.seed(), { "A", "Z", "B" } } }, "graphlet"), *pair.index,
            "", multi_graph());
    const Json::Value &r = out["results"][0];
    ASSERT_TRUE(r.isMember("graphlet")) << r;
    const std::string text = r["graphlet"].asString();
    EXPECT_NE(std::string::npos, text.find("\nX not_in_graph . Z\n")) << text;
    // B, unsupported, before it in no order: the request's (A, Z, B)
    const size_t z = text.find("\nX not_in_graph . Z\n");
    const size_t b = text.find("\nX seed_unsupported ");
    ASSERT_NE(std::string::npos, b) << text;
    EXPECT_LT(z, b);
    EXPECT_EQ(strings({ "Z" }), r["unknown_labels"]);
}

// strategy.labels.extra: a label the pair does not hold is no switch target here (its cost
// entries ignored) and is listed once in the answer's unknown_labels, whether or not a seed is
// walked; the walk is the one without it; a single-graph server refuses the request
TEST(LabelsNotInGraph, ExtraLabelsAreListedInTheAnswer) {
    const Pair pair;
    auto with_extra = [&](const std::vector<std::string> &extra, const std::string &seed,
                          const std::vector<std::string> &labels) {
        Json::Value r = traverse_request({ { seed, labels } });
        r["strategy"]["labels"]["extra"] = strings(extra);
        r["strategy"]["labels"]["change_cost"] = parse(R"({"model": "constant", "value": 1})");
        r["strategy"]["labels"]["loss_budget"] = 1;
        return r;
    };
    EXPECT_EQ("seed 's0': Label not found in the annotation: 'Z' (no sequence "
              "header index loaded)", refusal_of([&] {
        process_traverse_request(with_extra({ "Z", "B" }, pair.seed(), { "A" }), *pair.index, "");
    }));
    const Json::Value out = process_traverse_request(
            with_extra({ "Z", "B", "Y", "Z" }, pair.seed(), { "A" }), *pair.index, "",
            multi_graph());
    EXPECT_EQ(strings({ "Z", "Y" }), out["unknown_labels"]);
    const Json::Value &r = out["results"][0];
    EXPECT_EQ(Json::Value(Json::arrayValue), r["unknown_labels"]);
    const Json::Value without = process_traverse_request(
            with_extra({ "B" }, pair.seed(), { "A" }), *pair.index, "", multi_graph());
    EXPECT_EQ(Json::Value(Json::arrayValue), without["unknown_labels"]);
    EXPECT_EQ(walk_of(without["results"][0]), walk_of(r));
    // the strategy echoes the request's extra labels, all of them
    EXPECT_EQ(strings({ "Z", "B", "Y", "Z" }), out["strategy"]["labels"]["extra"]);

    // a table naming an absent label: its entries are ignored, the others kept
    Json::Value table = with_extra({ "Z", "B" }, pair.seed(), { "A" });
    table["strategy"]["labels"]["change_cost"] = parse(
            R"({"model": "table", "default": "forbid",
                "entries": [["A", "Z", 0.5], ["A", "B", 1]]})");
    const Json::Value tabled = process_traverse_request(table, *pair.index, "", multi_graph());
    EXPECT_EQ(strings({ "Z" }), tabled["unknown_labels"]);
    EXPECT_EQ("complete", tabled["results"][0]["outcome"]["walks"].asString()) << tabled;

    // no seed walked: the answer lists them still
    const Json::Value unwalked = process_traverse_request(
            with_extra({ "Z" }, pair.seed(), { "Y" }), *pair.index, "", multi_graph());
    EXPECT_EQ(strings({ "Z" }), unwalked["unknown_labels"]);
    EXPECT_EQ("labels_not_in_graph", unwalked["results"][0]["outcome"]["walks"].asString());
}

// /resolve's explicit labels: the absent ones are dropped and listed in unknown_labels (each
// once, in the request's order), the answer that of the others; every label absent leaves no
// label and no candidate; an explicit seed naming one is not covered by it; a single-graph
// server refuses the request and states no unknown_labels
TEST(LabelsNotInGraph, ResolveListsTheAbsentLabels) {
    const Pair pair;
    auto resolve = [&](const std::string &request, bool multi) {
        Json::Value json = parse(request);
        json["sequence"] = pair.records[0].substr(10, 50);
        Json::Value out = process_resolve_request(json, *pair.index, "", 0, nullptr, nullptr, {},
                                                  nullptr, nullptr, multi);
        out.removeMember("timing");
        out.removeMember("work");
        return out;
    };
    EXPECT_EQ("Label not found in the annotation: 'Z' (no sequence "
              "header index loaded)", refusal_of([&] {
        resolve(R"({"labels": ["Z", "A"]})", false);
    }));
    EXPECT_FALSE(resolve(R"({"labels": ["A"]})", false).isMember("unknown_labels"));

    Json::Value out = resolve(R"({"labels": ["Z", "A", "Y", "Z"]})", true);
    EXPECT_EQ(strings({ "Z", "Y" }), out["unknown_labels"]);
    Json::Value alone = resolve(R"({"labels": ["A"]})", true);
    EXPECT_EQ(Json::Value(Json::arrayValue), alone["unknown_labels"]);
    ASSERT_EQ(1u, out["labels"].size()) << out;
    EXPECT_EQ("A", out["labels"][0]["label"].asString());
    out.removeMember("unknown_labels");
    alone.removeMember("unknown_labels");
    EXPECT_EQ(alone, out);

    const Json::Value none = resolve(R"({"labels": ["Z"]})", true);
    EXPECT_EQ(strings({ "Z" }), none["unknown_labels"]);
    EXPECT_EQ(Json::Value(Json::arrayValue), none["labels"]);
    EXPECT_EQ(Json::Value(Json::arrayValue), none["candidates"]);
    EXPECT_FALSE(none["graph_runs"].empty());

    // a discovery names no label: nothing is unknown
    EXPECT_EQ(Json::Value(Json::arrayValue),
              resolve(R"({"discover": {"max_labels": 5}})", true)["unknown_labels"]);

    // an explicit seed naming an absent label: not covered by it
    const Json::Value selected = resolve(R"({"labels": ["A", "Z"], "select": {"policy":
        "explicit", "seeds": [{"kmer_interval": [0, 10], "labels": ["A", "Z"]}]}})", true);
    ASSERT_EQ(1u, selected["selection"]["seeds"].size()) << selected;
    EXPECT_EQ(strings({ "Z" }), selected["selection"]["seeds"][0]["labels_not_covering"]);
    EXPECT_NE(std::string::npos, refusal_of([&] {
        resolve(R"({"labels": ["A"], "select": {"policy": "explicit", "seeds":
            [{"kmer_interval": [0, 10], "labels": ["A", "Z"]}]}})", true);
    }).find("Explicit seed label was not profiled: 'Z'"));
}
