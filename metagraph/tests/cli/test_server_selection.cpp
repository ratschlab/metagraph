#include "gtest/gtest.h"

#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>

#include "cli/server_checks.hpp"
#include "cli/server_selection.hpp"
#include "cli/traverse.hpp"


// The selection of a multi-graph server's pairs (server_selection.hpp) on a small graph list:
// today's rules of each route, which the selectors of SPEC-multigraph-envelope.md section 6
// extend, and the body of a selector refusal

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

// A reference sample of SPEC-multigraph-envelope.md, {route, status, note, request?, body}
Json::Value sample(const std::string &name) {
    // run from the build directory, like every test reading ../tests/data
    const std::string path = "../tests/data/multigraph_contract/" + name;
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.good()) << path << " not found (run from the build directory)";
    const std::string text { std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>() };
    return parse(text);
}

} // namespace

namespace {

// a: one pair; b: two graphs, its first pair listed twice; c: one graph with two annotations
GraphIndexes small_indexes() {
    GraphIndexes indexes;
    indexes["b"] = { { "b1.dbg", "b1.anno" }, { "b2.dbg", "b2.anno" }, { "b1.dbg", "b1.anno" } };
    indexes["a"] = { { "a.dbg", "a.anno" } };
    indexes["c"] = { { "c.dbg", "c0.anno" }, { "c.dbg", "c1.anno" } };
    return indexes;
}

const SelectedPair kA { "a", { "a.dbg", "a.anno" } };
const SelectedPair kB1 { "b", { "b1.dbg", "b1.anno" } };
const SelectedPair kB2 { "b", { "b2.dbg", "b2.anno" } };
const SelectedPair kC0 { "c", { "c.dbg", "c0.anno" } };
const SelectedPair kC1 { "c", { "c.dbg", "c1.anno" } };

// the message of |f|'s std::invalid_argument (InvalidRequest is one), "" when it throws none
template <class F>
std::string refusal_of(const F &f) {
    try {
        f();
    } catch (const std::invalid_argument &e) {
        return e.what();
    }
    return "";
}

// |f|'s SelectionRefusal: its code, "" when it throws none
template <class F>
std::string selection_code_of(const F &f, Json::Value *body = nullptr) {
    try {
        f();
    } catch (const SelectionRefusal &e) {
        EXPECT_EQ(400, e.status());
        EXPECT_EQ(e.code(), e.body()["code"].asString());
        if (body)
            *body = e.body();
        return e.code();
    }
    return "";
}

std::vector<SelectedPair> select(SelectionRoute route, const std::string &request,
                                 size_t max = 10) {
    const GraphIndexes indexes = small_indexes();
    return select_pairs(route, parse(request), indexes, max);
}

std::vector<SelectedPair> probe(SelectionRoute route,
                                const std::vector<std::pair<std::string, std::string>> &query,
                                size_t max = 10) {
    const GraphIndexes indexes = small_indexes();
    return select_pairs_from_query(route, query, indexes, max);
}

using Pairs = std::vector<SelectedPair>;

} // namespace

// The answers' order: the names in byte order, each name's pairs in the list's order, a pair
// listed twice under a name once
TEST(ServerSelection, OrderedPairs) {
    const GraphIndexes indexes = small_indexes();
    EXPECT_EQ(std::vector<std::string>({ "a", "b", "c" }), sorted_names(indexes));
    const auto pairs = ordered_pairs(indexes, { "b", "c" });
    ASSERT_EQ(4u, pairs.size());
    EXPECT_EQ(kB1.pair, pairs[0].second);
    EXPECT_EQ(kB2.pair, pairs[1].second);
    EXPECT_EQ(kC0.pair, pairs[2].second);
    EXPECT_EQ("c", pairs[3].first);
    EXPECT_EQ(kC1.pair, pairs[3].second);
}

// POST /search: `graphs` deduplicated and sorted, every name of a small list without it; a
// pair listed twice once; `pairs` refused (it has no re-ask)
TEST(ServerSelection, Search) {
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(SelectionRoute::SEARCH, "{}"));
    EXPECT_EQ(Pairs({ kA, kC0, kC1 }),
              select(SelectionRoute::SEARCH, R"({"graphs": ["c", "a", "c"]})"));
    // `graph` is not /search's
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }),
              select(SelectionRoute::SEARCH, R"({"graph": "a"})"));
    EXPECT_EQ("Request with an uninitialized graph z",
              refusal_of([] { select(SelectionRoute::SEARCH, R"({"graphs": ["a", "z"]})"); }));
    EXPECT_EQ("Bad request: requests without names (no \"graphs\" field) are only supported "
              "for small indexes (<=2 names)",
              refusal_of([] { select(SelectionRoute::SEARCH, "{}", 2); }));
    EXPECT_EQ(Pairs({ kA }), select(SelectionRoute::SEARCH, R"({"graphs": ["a"]})", 2));
    Json::Value body;
    EXPECT_EQ("pairs_unsupported", selection_code_of([] {
        select(SelectionRoute::SEARCH, R"({"pairs": [{"graph": "a"}]})");
    }, &body));
    EXPECT_EQ(parse(R"({"error": "request: /search selects by graphs only (it has no re-ask): )"
                    R"(pairs is not accepted", "code": "pairs_unsupported"})"), body);
}

// POST /pattern: pattern_graph_names' rule and refusals
TEST(ServerSelection, Pattern) {
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(SelectionRoute::PATTERN, "{}"));
    EXPECT_EQ(Pairs({ kB1, kB2 }), select(SelectionRoute::PATTERN, R"({"graphs": ["b", "b"]})"));
    EXPECT_EQ("request.graphs: unknown graph 'z' (GET /capabilities lists the graphs)",
              refusal_of([] { select(SelectionRoute::PATTERN, R"({"graphs": ["z"]})"); }));
    EXPECT_EQ("request.graphs: expected a non-empty array of graph names (GET /capabilities "
              "lists them)",
              refusal_of([] { select(SelectionRoute::PATTERN, R"({"graphs": []})"); }));
    EXPECT_EQ("request.graphs: required on this server, which hosts 3 graph names (more than "
              "2; GET /capabilities lists them)",
              refusal_of([] { select(SelectionRoute::PATTERN, "{}", 2); }));
    EXPECT_EQ("pairs_unsupported", selection_code_of([] {
        select(SelectionRoute::PATTERN, R"({"pairs": []})");
    }));
}

// POST /resolve and POST /traverse: today both read exactly one pair, by the same rules and
// with the same refusals; a name is required, so the names cap does not apply
TEST(ServerSelection, ResolveAndTraverseOnePair) {
    for (SelectionRoute route : { SelectionRoute::RESOLVE, SelectionRoute::TRAVERSE }) {
        SCOPED_TRACE(route == SelectionRoute::RESOLVE ? "resolve" : "traverse");
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a"})"));
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"graphs": ["a"]})"));
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a"})", 1));
        // a pair listed twice is one: the name's two graphs, picked by graph_path
        EXPECT_EQ(Pairs({ kB1 }), select(route, R"({"graph": "b", "graph_path": "b1.dbg"})"));
        EXPECT_EQ(Pairs({ kB2 }), select(route, R"({"graphs": ["b"], "graph_path": "b2.dbg"})"));
        // graph_path is ignored where the name has one pair, as today
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a", "graph_path": "x.dbg"})"));
        EXPECT_EQ("Bad request: index 'b' spans several graphs (b1.dbg, b2.dbg); pass "
                  "'graph_path' to pick one",
                  refusal_of([&] { select(route, R"({"graph": "b"})"); }));
        EXPECT_EQ("Bad request: index 'c' lists one graph (c.dbg) with several annotations "
                  "(c0.anno, c1.anno): a traversal reads one annotation, which graph_path "
                  "cannot choose; list each annotation under a name of its own",
                  refusal_of([&] { select(route, R"({"graph": "c"})"); }));
        EXPECT_EQ("Bad request: graph_path 'c.dbg' names a graph with several annotations "
                  "under index 'c' (c0.anno, c1.anno): a traversal reads one annotation, which "
                  "graph_path cannot choose; list each annotation under a name of its own",
                  refusal_of([&] { select(route, R"({"graph": "c", "graph_path": "c.dbg"})"); }));
        EXPECT_EQ("Bad request: graph_path 'x.dbg' is not part of index 'b'",
                  refusal_of([&] { select(route, R"({"graph": "b", "graph_path": "x.dbg"})"); }));
        EXPECT_EQ("Bad request: 'graphs' names the one graph a traversal reads: expected [name]",
                  refusal_of([&] { select(route, R"({"graphs": ["a", "b"]})"); }));
        EXPECT_EQ("Bad request: give the graph as 'graph' or as 'graphs', not both",
                  refusal_of([&] { select(route, R"({"graph": "a", "graphs": ["a"]})"); }));
        EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
                  refusal_of([&] { select(route, R"({"graph": 1})"); }));
        // the graph-less request stays refused (a later release runs it)
        EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
                  refusal_of([&] { select(route, "{}"); }));
        EXPECT_EQ("Bad request: unknown graph 'z'",
                  refusal_of([&] { select(route, R"({"graph": "z"})"); }));
        EXPECT_EQ("pairs_unsupported", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "a"}]})");
        }));
    }
}

// The GET probes: ?graph= (and graph_path=) one pair, and the parameters' refusals as GET
// /traverse/capabilities and GET /pattern/capabilities state them today
TEST(ServerSelection, CapabilitiesProbes) {
    EXPECT_EQ(Pairs({ kB2 }), probe(SelectionRoute::TRAVERSE,
                                    { { "graph", "b" }, { "graph_path", "b2.dbg" } }));
    EXPECT_EQ(Pairs({ kA }), probe(SelectionRoute::PATTERN, { { "graph", "a" } }));
    // the names cap is not read by a probe of one pair
    EXPECT_EQ(Pairs({ kA }), probe(SelectionRoute::PATTERN, { { "graph", "a" } }, 1));
    // a probe without ?graph= is refused, however few names the server hosts
    EXPECT_EQ("Bad request: in multi-graph mode GET /traverse/capabilities needs "
              "?graph=<name> (and graph_path=<path> when the name spans several graphs); GET "
              "/capabilities lists the graphs",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, {}); }));
    EXPECT_EQ("Bad request: in multi-graph mode GET /pattern/capabilities needs "
              "?graph=<name> (and graph_path=<path> when the name spans several graphs); GET "
              "/capabilities lists the graphs",
              refusal_of([] { probe(SelectionRoute::PATTERN, { { "graph_path", "a.dbg" } }); }));
    EXPECT_EQ("Bad request: unknown parameter 'x' of GET /pattern/capabilities (graph, "
              "graph_path)",
              refusal_of([] { probe(SelectionRoute::PATTERN, { { "graph", "a" }, { "x", "1" } }); }));
    EXPECT_EQ("Bad request: the parameter 'graph' is given twice",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, { { "graph", "a" },
                                                                { "graph", "b" } }); }));
    EXPECT_EQ("Bad request: index 'b' spans several graphs (b1.dbg, b2.dbg); pass "
              "'graph_path' to pick one",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, { { "graph", "b" } }); }));
    EXPECT_EQ("Bad request: unknown graph 'z'",
              refusal_of([] { probe(SelectionRoute::PATTERN, { { "graph", "z" } }); }));
    // any range of (key, value): the server passes its query multimap
    const GraphIndexes indexes = small_indexes();
    const std::multimap<std::string, std::string> query { { "graph", "c" },
                                                          { "graph_path", "c.dbg" } };
    EXPECT_EQ("Bad request: graph_path 'c.dbg' names a graph with several annotations under "
              "index 'c' (c0.anno, c1.anno): a traversal reads one annotation, which "
              "graph_path cannot choose; list each annotation under a name of its own",
              refusal_of([&] {
                  select_pairs_from_query(SelectionRoute::TRAVERSE, query, indexes, 10);
              }));
    // the JSON spelling of the same probe
    EXPECT_EQ(Pairs({ kB1 }), select(SelectionRoute::CAPABILITIES,
                                     R"({"graph": "b", "graph_path": "b1.dbg"})"));
    EXPECT_EQ("Bad request: 'graph': expected a string",
              refusal_of([] { select(SelectionRoute::CAPABILITIES, R"({"graph": 1})"); }));
    EXPECT_NE("", refusal_of([] { select(SelectionRoute::CAPABILITIES, "{}"); }));

    // probe_parameters alone, as both GET handlers call it today
    EXPECT_EQ("Bad request: in multi-graph mode GET /pattern/capabilities needs ?graph=<name> "
              "(and graph_path=<path> when the name spans several graphs); GET /capabilities "
              "lists the graphs",
              refusal_of([] { probe_parameters(std::multimap<std::string, std::string>(),
                                               "/pattern/capabilities"); }));
    const auto [name, graph_path] = probe_parameters(query, "/traverse/capabilities");
    EXPECT_EQ("c", name);
    EXPECT_EQ("c.dbg", graph_path.value_or(""));
}

// A refusal of section 6: HTTP 400 {error, code, candidates?}, the candidates where a choice
// exists
TEST(ServerSelection, RefusalBody) {
    const GraphIndexes indexes = small_indexes();
    const Json::Value candidates = selection_candidates("b", indexes);
    EXPECT_EQ(parse(R"([{"graph": "b", "graph_path": "b1.dbg", "annotation_path": "b1.anno"},
                        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}])"),
              candidates);
    EXPECT_EQ(Json::Value(Json::arrayValue), selection_candidates("z", indexes));
    const SelectionRefusal refusal("request: graph 'b' holds 2 pairs", "traverse_single_pair",
                                   candidates);
    EXPECT_EQ(400, refusal.status());
    EXPECT_EQ("traverse_single_pair", refusal.code());
    EXPECT_EQ(candidates, refusal.candidates());
    EXPECT_EQ("request: graph 'b' holds 2 pairs", std::string(refusal.what()));
    EXPECT_EQ(candidates, refusal.body()["candidates"]);
    EXPECT_EQ(3u, refusal.body().size());
    // no choice: no member
    const SelectionRefusal unknown("request.graphs: unknown graph 'z'", "graph_unknown");
    EXPECT_FALSE(unknown.body().isMember("candidates"));
    EXPECT_TRUE(unknown.candidates().isNull());
    // the refusals of sample 22 are bodies of this shape
    for (const Json::Value &c : sample("22_selector_refusals.json")["body"]) {
        const Json::Value &body = c["body"];
        const SelectionRefusal r(body["error"].asString(), body["code"].asString(),
                                 body.get("candidates", Json::Value()));
        EXPECT_EQ(body, r.body()) << c["case"].asString();
    }
    // answered as a 400 with its body (answer_request): an HttpError
    const RequestAnswer answer = answer_request("", "", 0, [&](const std::string &) -> Json::Value {
        throw SelectionRefusal("request: pairs", "pairs_unsupported");
    }, true, nullptr);
    EXPECT_EQ(400, answer.status);
    EXPECT_EQ(parse(R"({"error": "request: pairs", "code": "pairs_unsupported"})"),
              parse(answer.body));
}
