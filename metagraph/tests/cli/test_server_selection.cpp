#include "gtest/gtest.h"

#include <fstream>
#include <functional>
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

// |f|'s SelectionRefusal: its code; "(none)" when it throws another std::invalid_argument (a
// refusal without a code), "" when it throws none
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
    } catch (const std::invalid_argument &) {
        return "(none)";
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
// pair listed twice once; `pairs` refused (it has no re-ask); an unknown name refused with its
// message and graph_unknown
TEST(ServerSelection, Search) {
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(SelectionRoute::SEARCH, "{}"));
    EXPECT_EQ(Pairs({ kA, kC0, kC1 }),
              select(SelectionRoute::SEARCH, R"({"graphs": ["c", "a", "c"]})"));
    // `graph` is not /search's
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }),
              select(SelectionRoute::SEARCH, R"({"graph": "a"})"));
    Json::Value body;
    EXPECT_EQ("graph_unknown", selection_code_of([] {
        select(SelectionRoute::SEARCH, R"({"graphs": ["a", "z"]})");
    }, &body));
    EXPECT_EQ(parse(R"({"error": "Request with an uninitialized graph z",
                        "code": "graph_unknown"})"), body);
    EXPECT_EQ("Bad request: requests without names (no \"graphs\" field) are only supported "
              "for small indexes (<=2 names)",
              refusal_of([] { select(SelectionRoute::SEARCH, "{}", 2); }));
    EXPECT_EQ(Pairs({ kA }), select(SelectionRoute::SEARCH, R"({"graphs": ["a"]})", 2));
    EXPECT_EQ("pairs_unsupported", selection_code_of([] {
        select(SelectionRoute::SEARCH, R"({"pairs": [{"graph": "a"}]})");
    }, &body));
    EXPECT_EQ(parse(R"({"error": "request: /search selects by graphs only (it has no re-ask): )"
                    R"(pairs is not accepted", "code": "pairs_unsupported"})"), body);
    // with graphs too: pairs is refused as /search's, not as a conflict
    EXPECT_EQ("pairs_unsupported", selection_code_of([] {
        select(SelectionRoute::SEARCH, R"({"graphs": ["a"], "pairs": [{"graph": "a"}]})");
    }));
}

// POST /pattern: `graphs` (pattern_graph_names' rule, an unknown name coded) or `pairs`; `graph`
// and `graph_path` are not its selectors (its parser refuses them)
TEST(ServerSelection, Pattern) {
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(SelectionRoute::PATTERN, "{}"));
    EXPECT_EQ(Pairs({ kB1, kB2 }), select(SelectionRoute::PATTERN, R"({"graphs": ["b", "b"]})"));
    Json::Value body;
    EXPECT_EQ("graph_unknown", selection_code_of([] {
        select(SelectionRoute::PATTERN, R"({"graphs": ["a", "z", "y"]})");
    }, &body));
    EXPECT_EQ(parse(R"x({"error": "request.graphs: unknown graph 'z' (GET /capabilities lists )x"
                    R"x(the graphs)", "code": "graph_unknown"})x"), body);
    EXPECT_EQ("request.graphs: expected a non-empty array of graph names (GET /capabilities "
              "lists them)",
              refusal_of([] { select(SelectionRoute::PATTERN, R"({"graphs": []})"); }));
    EXPECT_EQ("request.graphs: expected a non-empty array of graph names (GET /capabilities "
              "lists them)",
              refusal_of([] { select(SelectionRoute::PATTERN, R"({"graphs": ["a", 1]})"); }));
    // the names cap keeps its 400, without a code
    EXPECT_EQ("request.graphs: required on this server, which hosts 3 graph names (more than "
              "2; GET /capabilities lists them)",
              refusal_of([] { select(SelectionRoute::PATTERN, "{}", 2); }));
    EXPECT_EQ("(none)", selection_code_of([] { select(SelectionRoute::PATTERN, "{}", 2); }));
    // any number of pairs, in the index's order
    EXPECT_EQ(Pairs({ kA, kB2, kC1 }), select(SelectionRoute::PATTERN, R"({"pairs": [
        {"graph": "c", "annotation_path": "c1.anno"}, {"graph": "b", "graph_path": "b2.dbg"},
        {"graph": "a"}]})"));
    EXPECT_EQ("", refusal_of([] { select(SelectionRoute::PATTERN, R"({"graph": "zz"})"); }));
    EXPECT_EQ("selection_conflict", selection_code_of([] {
        select(SelectionRoute::PATTERN, R"({"graphs": ["a"], "pairs": [{"graph": "a"}]})");
    }));
}

// POST /resolve: any number of pairs, by graphs, by graph (with graph_path), or by pairs;
// without a selection every pair, up to the names cap
TEST(ServerSelection, Resolve) {
    const SelectionRoute route = SelectionRoute::RESOLVE;
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a"})"));
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graphs": ["a"]})"));
    EXPECT_EQ(Pairs({ kB1, kB2 }), select(route, R"({"graph": "b"})"));
    EXPECT_EQ(Pairs({ kA, kC0, kC1 }), select(route, R"({"graphs": ["c", "a"]})"));
    // a name listed twice is read once, in byte order
    EXPECT_EQ(Pairs({ kA, kB1, kB2 }), select(route, R"({"graphs": ["b", "a", "b"]})"));
    EXPECT_EQ(Pairs({ kC0, kC1 }), select(route, R"({"graph": "c"})", 1));
    // graph_path keeps the name's pairs on that graph: the two annotations of c's one graph
    EXPECT_EQ(Pairs({ kB2 }), select(route, R"({"graph": "b", "graph_path": "b2.dbg"})"));
    EXPECT_EQ(Pairs({ kC0, kC1 }), select(route, R"({"graphs": ["c"], "graph_path": "c.dbg"})"));
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(route, "{}"));
    EXPECT_EQ("request.graphs: required on this server, which hosts 3 graph names (more than "
              "2; GET /capabilities lists them)", refusal_of([&] { select(route, "{}", 2); }));
    Json::Value body;
    // a graph_path the name does not have: no pair, the name's pairs to choose from
    EXPECT_EQ("pair_not_found", selection_code_of([&] {
        select(route, R"({"graph": "b", "graph_path": "x.dbg"})");
    }, &body));
    // every refusal of a request 8fa32e553 refused keeps its text, byte for byte; the code is
    // added beside it
    EXPECT_EQ(parse(R"({"error": "Bad request: graph_path 'x.dbg' is not part of index 'b'",
                        "code": "pair_not_found", "candidates": [
                        {"graph": "b", "graph_path": "b1.dbg", "annotation_path": "b1.anno"},
                        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}]})"),
              body);
    // a name of one pair is that pair, whatever graph_path says, as the base answered it
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a", "graph_path": "x.dbg"})"));
    EXPECT_EQ("selection_conflict", selection_code_of([&] {
        select(route, R"({"graph": "a", "graphs": ["a"]})");
    }, &body));
    EXPECT_EQ("Bad request: give the graph as 'graph' or as 'graphs', not both",
              body["error"].asString());
    EXPECT_EQ("selection_conflict", selection_code_of([&] {
        select(route, R"({"graphs": ["a", "b"], "graph_path": "b1.dbg"})");
    }, &body));
    EXPECT_EQ("Bad request: 'graphs' names the one graph a traversal reads: expected [name]",
              body["error"].asString());
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([&] { select(route, R"({"graph_path": "b1.dbg"})"); }));
    EXPECT_EQ("(none)", selection_code_of([&] { select(route, R"({"graph_path": "b1.dbg"})"); }));
    EXPECT_EQ("graph_unknown", selection_code_of([&] { select(route, R"({"graph": "z"})"); },
                                                 &body));
    EXPECT_EQ(parse(R"({"error": "Bad request: unknown graph 'z'", "code": "graph_unknown"})"),
              body);
    EXPECT_EQ("graph_unknown", selection_code_of([&] {
        select(route, R"({"graphs": ["a", "z", "y"]})");
    }, &body));
    EXPECT_EQ("Bad request: unknown graph 'z'", body["error"].asString());
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([&] { select(route, R"({"graph": 1})"); }));
    for (const char *graphs : { "[]", R"("a")", R"(["a", 1])" }) {
        EXPECT_EQ("Bad request: 'graphs' names the one graph a traversal reads: expected [name]",
                  refusal_of([&] { select(route, std::string(R"({"graphs": )") + graphs + "}"); }))
                << graphs;
        EXPECT_EQ("(none)", selection_code_of([&] {
            select(route, std::string(R"({"graphs": )") + graphs + "}");
        })) << graphs;
    }
    // a graph_path that is no string is not read here: the parser refuses it, as before
    EXPECT_EQ(Pairs({ kB1, kB2 }), select(route, R"({"graph": "b", "graph_path": 2})"));
}

// POST /traverse: exactly one pair, by graph or graphs (with graph_path), or one object of
// pairs; a selection of several is refused with them as the candidates
TEST(ServerSelection, TraverseOnePair) {
    const SelectionRoute route = SelectionRoute::TRAVERSE;
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a"})"));
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graphs": ["a", "a"]})"));
    EXPECT_EQ(Pairs({ kA }), select(route, R"({"graph": "a"})", 1));
    // a pair listed twice is one: the name's two graphs, picked by graph_path
    EXPECT_EQ(Pairs({ kB1 }), select(route, R"({"graph": "b", "graph_path": "b1.dbg"})"));
    EXPECT_EQ(Pairs({ kB2 }), select(route, R"({"graphs": ["b"], "graph_path": "b2.dbg"})"));
    EXPECT_EQ(Pairs({ kC1 }), select(route, R"({"pairs": [{"graph": "c",
                                                           "annotation_path": "c1.anno"}]})"));
    // the same pair named twice is one pair
    EXPECT_EQ(Pairs({ kB2 }), select(route, R"({"pairs": [{"graph": "b", "graph_path": "b2.dbg"},
        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}]})"));

    Json::Value body;
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select(route, R"({"graph": "b"})");
    }, &body));
    EXPECT_EQ(parse(R"x({"error": "Bad request: index 'b' spans several graphs (b1.dbg, b2.dbg); )x"
                    R"x(pass 'graph_path' to pick one",
                        "code": "traverse_single_pair", "candidates": [
                        {"graph": "b", "graph_path": "b1.dbg", "annotation_path": "b1.anno"},
                        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}]})x"),
              body);
    // one graph with two annotations, by name or by graph_path
    for (const char *request : { R"({"graph": "c"})", R"({"graph": "c", "graph_path": "c.dbg"})" }) {
        EXPECT_EQ("traverse_single_pair", selection_code_of([&] { select(route, request); },
                                                            &body)) << request;
        EXPECT_EQ(selection_candidates("c", small_indexes()), body["candidates"]) << request;
    }
    selection_code_of([&] { select(route, R"({"graph": "c"})"); }, &body);
    EXPECT_EQ("Bad request: index 'c' lists one graph (c.dbg) with several annotations (c0.anno, "
              "c1.anno): a traversal reads one annotation, which graph_path cannot choose; list "
              "each annotation under a name of its own", body["error"].asString());
    selection_code_of([&] { select(route, R"({"graph": "c", "graph_path": "c.dbg"})"); }, &body);
    EXPECT_EQ("Bad request: graph_path 'c.dbg' names a graph with several annotations under "
              "index 'c' (c0.anno, c1.anno): a traversal reads one annotation, which graph_path "
              "cannot choose; list each annotation under a name of its own",
              body["error"].asString());
    // several names, several objects, no selector
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select(route, R"({"graphs": ["c", "a"]})");
    }, &body));
    EXPECT_EQ("Bad request: 'graphs' names the one graph a traversal reads: expected [name]",
              body["error"].asString());
    EXPECT_EQ(3u, body["candidates"].size());
    // several names, one of them unknown: refused as several, the known ones the candidates
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select(route, R"({"graphs": ["z", "a"]})");
    }, &body));
    EXPECT_EQ(parse(R"([{"graph": "a", "graph_path": "a.dbg", "annotation_path": "a.anno"}])"),
              body["candidates"]);
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select(route, R"({"pairs": [{"graph": "c", "annotation_path": "c1.anno"},
                                    {"graph": "a"}]})");
    }, &body));
    EXPECT_EQ("request: pairs names 2 pairs; a traverse runs on exactly one pair: list one",
              body["error"].asString());
    EXPECT_EQ(parse(R"([{"graph": "a", "graph_path": "a.dbg", "annotation_path": "a.anno"},
                        {"graph": "c", "graph_path": "c.dbg", "annotation_path": "c1.anno"}])"),
              body["candidates"]);
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] { select(route, "{}"); }, &body));
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              body["error"].asString());
    EXPECT_EQ(5u, body["candidates"].size());
    // beyond the names cap the same words, without a code (the cap's refusal has none)
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([&] { select(route, "{}", 2); }));
    EXPECT_EQ("(none)", selection_code_of([&] { select(route, "{}", 2); }));
    // a server of one pair runs a request without a selector on it
    GraphIndexes one;
    one["a"] = { { "a.dbg", "a.anno" }, { "a.dbg", "a.anno" } };
    EXPECT_EQ(Pairs({ kA }), select_pairs(route, parse("{}"), one, 10));

    EXPECT_EQ("pair_not_found", selection_code_of([&] {
        select(route, R"({"graph": "b", "graph_path": "x.dbg"})");
    }, &body));
    EXPECT_EQ("Bad request: graph_path 'x.dbg' is not part of index 'b'", body["error"].asString());
    EXPECT_EQ(selection_candidates("b", small_indexes()), body["candidates"]);
    EXPECT_EQ("selection_conflict", selection_code_of([&] {
        select(route, R"({"graph": "a", "graphs": ["a"]})");
    }, &body));
    EXPECT_EQ("Bad request: give the graph as 'graph' or as 'graphs', not both",
              body["error"].asString());
    EXPECT_EQ("graph_unknown", selection_code_of([&] { select(route, R"({"graphs": ["z"]})"); },
                                                 &body));
    EXPECT_EQ("Bad request: unknown graph 'z'", body["error"].asString());
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([&] { select(route, R"({"graph_path": "a.dbg"})"); }));
}

// The `pairs` selector (POST /resolve, /pattern, /traverse): each object names exactly one pair
// of its name; graphs, graph and graph_path beside it conflict; a malformed list is a plain 400
TEST(ServerSelection, PairsSelector) {
    for (SelectionRoute route : { SelectionRoute::RESOLVE, SelectionRoute::PATTERN,
                                  SelectionRoute::TRAVERSE }) {
        SCOPED_TRACE(route == SelectionRoute::RESOLVE ? "resolve"
                     : route == SelectionRoute::PATTERN ? "pattern" : "traverse");
        // graph_path may be omitted where the name has one graph, annotation_path where the
        // graph has one annotation
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"pairs": [{"graph": "a"}]})"));
        EXPECT_EQ(Pairs({ kA }), select(route, R"({"pairs": [{"graph": "a", "graph_path":
                                                   "a.dbg", "annotation_path": "a.anno"}]})"));
        EXPECT_EQ(Pairs({ kB1 }), select(route, R"({"pairs": [{"graph": "b",
                                                                "graph_path": "b1.dbg"}]})"));
        EXPECT_EQ(Pairs({ kB2 }), select(route, R"({"pairs": [{"graph": "b",
                                                                "annotation_path": "b2.anno"}]})"));
        EXPECT_EQ(Pairs({ kC0 }), select(route, R"({"pairs": [{"graph": "c", "graph_path":
                                                   "c.dbg", "annotation_path": "c0.anno"}]})"));

        Json::Value body;
        EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "a"}, {"graph": "b"}]})");
        }, &body));
        EXPECT_EQ(parse(R"({"error": "request.pairs[1]: 'b' has 2 pairs on 2 graphs: add )"
                        R"(graph_path", "code": "pair_ambiguous", "candidates": [
                        {"graph": "b", "graph_path": "b1.dbg", "annotation_path": "b1.anno"},
                        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}]})"),
                  body);
        for (const char *object : { R"({"graph": "c"})", R"({"graph": "c", "graph_path": "c.dbg"})" }) {
            EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
                select(route, std::string(R"({"pairs": [)") + object + "]}");
            }, &body)) << object;
            EXPECT_EQ("request.pairs[0]: 'c' has 2 pairs on graph_path 'c.dbg': add "
                      "annotation_path", body["error"].asString());
            EXPECT_EQ(selection_candidates("c", small_indexes()), body["candidates"]);
        }
        EXPECT_EQ("pair_not_found", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "c", "annotation_path": "c2.anno"}]})");
        }, &body));
        EXPECT_EQ("request.pairs[0]: no pair of 'c' has annotation_path 'c2.anno'",
                  body["error"].asString());
        EXPECT_EQ(selection_candidates("c", small_indexes()), body["candidates"]);
        EXPECT_EQ("pair_not_found", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "b", "graph_path": "b1.dbg",
                                         "annotation_path": "b2.anno"}]})");
        }, &body));
        EXPECT_EQ("request.pairs[0]: no pair of 'b' has graph_path 'b1.dbg' and annotation_path "
                  "'b2.anno'", body["error"].asString());
        // a name's graph_path does not name another name's pair
        EXPECT_EQ("pair_not_found", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "a", "graph_path": "b1.dbg"}]})");
        }));
        EXPECT_EQ("graph_unknown", selection_code_of([&] {
            select(route, R"({"pairs": [{"graph": "a"}, {"graph": "z"}]})");
        }, &body));
        EXPECT_EQ(parse(R"x({"error": "request.pairs[1]: unknown graph 'z' (GET /capabilities )x"
                        R"x(lists the graphs)", "code": "graph_unknown"})x"), body);
        // graphs beside pairs, on every route; graph and graph_path where they are selectors
        EXPECT_EQ("selection_conflict", selection_code_of([&] {
            select(route, R"({"graphs": ["a"], "pairs": [{"graph": "a"}]})");
        }, &body));
        EXPECT_EQ(parse(R"({"error": "request: graphs and pairs are exclusive: name graphs, or )"
                        R"(list pairs, not both", "code": "selection_conflict"})"), body);
        if (route != SelectionRoute::PATTERN) {
            EXPECT_EQ("selection_conflict", selection_code_of([&] {
                select(route, R"({"graph": "a", "pairs": [{"graph": "a"}]})");
            }));
            EXPECT_EQ("selection_conflict", selection_code_of([&] {
                select(route, R"({"graph_path": "a.dbg", "pairs": [{"graph": "a"}]})");
            }));
        }
        // malformed: a plain 400, before any name is looked up
        for (const char *pairs : { "[]", "{}", R"("a")", R"([{"graph": "a"}, "b"])" }) {
            EXPECT_EQ("(none)", selection_code_of([&] {
                select(route, std::string(R"({"pairs": )") + pairs + "}");
            })) << pairs;
        }
        EXPECT_EQ("request.pairs[0]: unknown member 'graphs' (graph, graph_path, annotation_path)",
                  refusal_of([&] { select(route, R"({"pairs": [{"graphs": ["a"]}]})"); }));
        EXPECT_EQ("request.pairs[0].graph: required (the pair's graph name)",
                  refusal_of([&] { select(route, R"({"pairs": [{"graph_path": "a.dbg"}]})"); }));
        EXPECT_EQ("request.pairs[0].graph_path: expected a string",
                  refusal_of([&] {
                      select(route, R"({"pairs": [{"graph": "a", "graph_path": 1}]})");
                  }));
        EXPECT_EQ("request.pairs[0].graph: expected a graph name",
                  refusal_of([&] { select(route, R"({"pairs": [{"graph": ["a"]}]})"); }));
    }
}

// The GET probes: ?graph= (with graph_path= and annotation_path=) one pair, coded as an object
// of pairs is; without ?graph= every pair, up to the names cap; the parameters' own refusals
TEST(ServerSelection, CapabilitiesProbes) {
    EXPECT_EQ(Pairs({ kB2 }), probe(SelectionRoute::TRAVERSE,
                                    { { "graph", "b" }, { "graph_path", "b2.dbg" } }));
    EXPECT_EQ(Pairs({ kA }), probe(SelectionRoute::PATTERN, { { "graph", "a" } }));
    EXPECT_EQ(Pairs({ kC1 }), probe(SelectionRoute::PATTERN, { { "graph", "c" },
                                                              { "annotation_path", "c1.anno" } }));
    EXPECT_EQ(Pairs({ kC0 }), probe(SelectionRoute::TRAVERSE, { { "annotation_path", "c0.anno" },
                                                               { "graph_path", "c.dbg" },
                                                               { "graph", "c" } }));
    // the names cap is not read by a probe of one pair
    EXPECT_EQ(Pairs({ kA }), probe(SelectionRoute::PATTERN, { { "graph", "a" } }, 1));
    // a probe without parameters describes every pair, in the envelope's order (a pair listed
    // twice under a name once), up to the names cap and refused beyond it in the words of a
    // probe without ?graph= before release 10, byte for byte
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), probe(SelectionRoute::TRAVERSE, {}));
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), probe(SelectionRoute::PATTERN, {}, 3));
    EXPECT_EQ("Bad request: in multi-graph mode GET /traverse/capabilities needs "
              "?graph=<name> (and graph_path=<path> when the name spans several graphs); GET "
              "/capabilities lists the graphs",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, {}, 2); }));
    EXPECT_EQ("(none)", selection_code_of([] { probe(SelectionRoute::TRAVERSE, {}, 2); }));
    // parameters without ?graph= are refused in the same words, however few names the server
    // hosts
    for (const char *key : { "graph_path", "annotation_path" }) {
        EXPECT_EQ("Bad request: in multi-graph mode GET /pattern/capabilities needs "
                  "?graph=<name> (and graph_path=<path> when the name spans several graphs); "
                  "GET /capabilities lists the graphs",
                  refusal_of([&] { probe(SelectionRoute::PATTERN, { { key, "a.dbg" } }); }))
                << key;
    }
    // the base's words, which named the parameters it had
    EXPECT_EQ("Bad request: unknown parameter 'x' of GET /pattern/capabilities (graph, "
              "graph_path)",
              refusal_of([] { probe(SelectionRoute::PATTERN, { { "graph", "a" }, { "x", "1" } }); }));
    EXPECT_EQ("Bad request: the parameter 'graph' is given twice",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, { { "graph", "a" },
                                                                { "graph", "b" } }); }));
    EXPECT_EQ("Bad request: the parameter 'annotation_path' is given twice",
              refusal_of([] { probe(SelectionRoute::TRAVERSE, { { "graph", "c" },
                                                                { "annotation_path", "c0.anno" },
                                                                { "annotation_path", "c1.anno" } }); }));
    Json::Value body;
    EXPECT_EQ("pair_ambiguous", selection_code_of([] {
        probe(SelectionRoute::TRAVERSE, { { "graph", "b" } });
    }, &body));
    EXPECT_EQ(parse(R"x({"error": "Bad request: index 'b' spans several graphs (b1.dbg, b2.dbg); )x"
                    R"x(pass 'graph_path' to pick one",
                        "code": "pair_ambiguous", "candidates": [
                        {"graph": "b", "graph_path": "b1.dbg", "annotation_path": "b1.anno"},
                        {"graph": "b", "graph_path": "b2.dbg", "annotation_path": "b2.anno"}]})x"),
              body);
    EXPECT_EQ("graph_unknown", selection_code_of([] {
        probe(SelectionRoute::PATTERN, { { "graph", "z" } });
    }, &body));
    EXPECT_EQ(parse(R"({"error": "Bad request: unknown graph 'z'", "code": "graph_unknown"})"),
              body);
    // with annotation_path (new in release 10) too
    EXPECT_EQ("graph_unknown", selection_code_of([] {
        probe(SelectionRoute::PATTERN, { { "graph", "z" }, { "annotation_path", "c0.anno" } });
    }, &body));
    EXPECT_EQ("Bad request: unknown graph 'z'", body["error"].asString());
    EXPECT_EQ("pair_not_found", selection_code_of([] {
        probe(SelectionRoute::PATTERN, { { "graph", "b" }, { "graph_path", "x.dbg" } });
    }, &body));
    EXPECT_EQ("Bad request: graph_path 'x.dbg' is not part of index 'b'",
              body["error"].asString());
    EXPECT_EQ(selection_candidates("b", small_indexes()), body["candidates"]);
    // a name of one pair is that pair, whatever graph_path says, as before release 10
    EXPECT_EQ(Pairs({ kA }), probe(SelectionRoute::PATTERN, { { "graph", "a" },
                                                             { "graph_path", "x.dbg" } }));
    // any range of (key, value): the server passes its query multimap
    const GraphIndexes indexes = small_indexes();
    const std::multimap<std::string, std::string> query { { "graph", "c" },
                                                          { "graph_path", "c.dbg" } };
    EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
        select_pairs_from_query(SelectionRoute::TRAVERSE, query, indexes, 10);
    }, &body));
    EXPECT_EQ("Bad request: graph_path 'c.dbg' names a graph with several annotations under "
              "index 'c' (c0.anno, c1.anno): a traversal reads one annotation, which graph_path "
              "cannot choose; list each annotation under a name of its own",
              body["error"].asString());
    EXPECT_EQ(selection_candidates("c", small_indexes()), body["candidates"]);
    // the JSON spelling of the same probe
    EXPECT_EQ(Pairs({ kB1 }), select(SelectionRoute::CAPABILITIES,
                                     R"({"graph": "b", "graph_path": "b1.dbg"})"));
    EXPECT_EQ(Pairs({ kC1 }), select(SelectionRoute::CAPABILITIES,
                                     R"({"graph": "c", "annotation_path": "c1.anno"})"));
    EXPECT_EQ("Bad request: 'graph': expected a string",
              refusal_of([] { select(SelectionRoute::CAPABILITIES, R"({"graph": 1})"); }));
    EXPECT_EQ(Pairs({ kA, kB1, kB2, kC0, kC1 }), select(SelectionRoute::CAPABILITIES, "{}"));
    EXPECT_EQ("pairs_unsupported", selection_code_of([] {
        select(SelectionRoute::CAPABILITIES, R"({"pairs": [{"graph": "a"}]})");
    }));

    // probe_parameters alone, the release-9 reading
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
    // the refusals of sample 22 are bodies of this shape (the sample is held in a variable:
    // a range-for over a member of the temporary sample() returns reads it once destroyed)
    const Json::Value refusals = sample("22_selector_refusals.json");
    for (const Json::Value &c : refusals["body"]) {
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

namespace {

// The graph lists of the deployed indexes (PLAN-multigraph-implementation.md section 6), small:
// r237, one name per pair; production RefSeq, one name over several graphs; Logan, names of
// five pairs; sra-human, one graph with several annotations under one name
const std::string kR8 = "/scratch/nvme8/refseq_r237/bact_other/";
const std::string kR10 = "/scratch/nvme10/refseq_r237/bact_pseudomonadota/";
const std::string kHuman = "/scratch/nvme0/human/graph_merged_complete_k31.primary.small.indexed.dbg";

GraphPair r237_pair(const std::string &dir) {
    return { dir + "graph.dbg", dir + "annotation.row_diff_brwt_coord.annodbg" };
}

GraphIndexes r237() {
    GraphIndexes indexes;
    indexes["bact_pseudomonadota"] = { r237_pair(kR10) };
    indexes["bact_other"] = { r237_pair(kR8) };
    indexes["viral"] = { r237_pair("/scratch/nvme9/refseq_r237/viral/") };
    return indexes;
}

GraphIndexes production_refseq() {
    GraphIndexes indexes;
    indexes["refseq"] = { r237_pair(kR8), r237_pair(kR10),
                          r237_pair("/scratch/nvme9/refseq_r237/viral/") };
    return indexes;
}

GraphIndexes logan() {
    GraphIndexes indexes;
    for (const std::string name : { "logan_b", "logan_a" }) {
        for (int i = 0; i < 5; ++i) {
            const std::string dir = "/scratch/nvme" + std::to_string(i) + "/" + name + "/";
            indexes[name].emplace_back(dir + "graph.dbg", dir + "annotation.row_diff_int_brwt.annodbg");
        }
    }
    return indexes;
}

std::string human_annotation(int nvme, int part) {
    return "/scratch/nvme" + std::to_string(nvme) + "/human/annotation_" + std::to_string(part)
         + ".relaxed.renamed.row_diff_brwt.annodbg";
}

GraphIndexes sra_human() {
    GraphIndexes indexes;
    indexes["sra-human-1/3"] = { { kHuman, human_annotation(0, 0) },
                                 { kHuman, human_annotation(3, 4) },
                                 { kHuman, human_annotation(4, 8) } };
    // the shards of every deployed index, beside it
    indexes["refseq"] = production_refseq().at("refseq");
    indexes["bact_other"] = { r237_pair(kR8) };
    return indexes;
}

std::vector<SelectedPair> select_on(const GraphIndexes &indexes, SelectionRoute route,
                                    const std::string &request, size_t max = 11) {
    return select_pairs(route, parse(request), indexes, max);
}

// |pairs| as {graph, graph_path, annotation_path} references
Json::Value refs(const std::vector<SelectedPair> &pairs) {
    Json::Value out(Json::arrayValue);
    for (const SelectedPair &pair : pairs) {
        out.append(pair_ref_json(pair));
    }
    return out;
}

} // namespace

// r237 (one name per pair): names select their one pair each, in byte order whatever the
// request's order; a pair object needs no path; a traverse names one name
TEST(ServerSelection, DeployedR237) {
    const GraphIndexes indexes = r237();
    const std::vector<SelectedPair> all = select_on(indexes, SelectionRoute::RESOLVE, "{}");
    ASSERT_EQ(3u, all.size());
    EXPECT_EQ((std::vector<std::string>{ "bact_other", "bact_pseudomonadota", "viral" }),
              (std::vector<std::string>{ all[0].name, all[1].name, all[2].name }));
    EXPECT_EQ(all, select_on(indexes, SelectionRoute::PATTERN,
                             R"({"graphs": ["viral", "bact_other", "bact_pseudomonadota"]})"));
    EXPECT_EQ(all, select_on(indexes, SelectionRoute::RESOLVE, R"({"pairs": [{"graph": "viral"},
        {"graph": "bact_pseudomonadota"}, {"graph": "bact_other"}]})"));
    EXPECT_EQ(std::vector<SelectedPair>({ all[1] }),
              select_on(indexes, SelectionRoute::TRAVERSE, R"({"graph": "bact_pseudomonadota"})"));
    // the names cap of mex (11) is not reached; one less than the names is
    EXPECT_NE("", refusal_of([&] { select_on(indexes, SelectionRoute::RESOLVE, "{}", 2); }));
    Json::Value body;
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select_on(indexes, SelectionRoute::TRAVERSE, R"({"graphs": ["viral", "bact_other"]})");
    }, &body));
    EXPECT_EQ(refs({ all[0], all[2] }), body["candidates"]);
}

// production RefSeq (one name, several graphs): /resolve and /pattern answer every pair of the
// name; a traverse and a probe name one graph; a pair object needs its graph_path
TEST(ServerSelection, DeployedProductionRefSeq) {
    const GraphIndexes indexes = production_refseq();
    const std::vector<SelectedPair> all = select_on(indexes, SelectionRoute::PATTERN,
                                                    R"({"graphs": ["refseq"]})");
    ASSERT_EQ(3u, all.size());
    EXPECT_EQ(r237_pair(kR8), all[0].pair);
    EXPECT_EQ(r237_pair(kR10), all[1].pair);
    // the list's order, not the request's
    const std::string two = R"({"pairs": [{"graph": "refseq", "graph_path": ")" + kR10
                          + R"(graph.dbg"}, {"graph": "refseq", "graph_path": ")" + kR8
                          + R"(graph.dbg"}]})";
    EXPECT_EQ(std::vector<SelectedPair>({ all[0], all[1] }),
              select_on(indexes, SelectionRoute::RESOLVE, two));
    Json::Value body;
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select_on(indexes, SelectionRoute::TRAVERSE, two);
    }, &body));
    EXPECT_EQ(refs({ all[0], all[1] }), body["candidates"]);
    EXPECT_EQ(std::vector<SelectedPair>({ all[1] }),
              select_on(indexes, SelectionRoute::TRAVERSE,
                        R"({"graph": "refseq", "graph_path": ")" + kR10 + R"(graph.dbg"})"));
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select_on(indexes, SelectionRoute::TRAVERSE, R"({"graph": "refseq"})");
    }, &body));
    EXPECT_EQ(refs(all), body["candidates"]);
    EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
        select_on(indexes, SelectionRoute::RESOLVE, R"({"pairs": [{"graph": "refseq"}]})");
    }, &body));
    EXPECT_EQ("request.pairs[0]: 'refseq' has 3 pairs on 3 graphs: add graph_path",
              body["error"].asString());
    EXPECT_EQ(refs(all), body["candidates"]);
    EXPECT_EQ(std::vector<SelectedPair>({ all[2] }),
              select_pairs_from_parameters(SelectionRoute::TRAVERSE, {
                      { "graph", "refseq" },
                      { "graph_path", "/scratch/nvme9/refseq_r237/viral/graph.dbg" } },
                      indexes, 11));
}

// Logan (names of five pairs): a re-ask names one unreached pair of a name without its four
// siblings; a whole name is five pairs, two names ten, in the list's order per name
TEST(ServerSelection, DeployedLogan) {
    const GraphIndexes indexes = logan();
    const std::vector<SelectedPair> all = select_on(indexes, SelectionRoute::PATTERN, "{}");
    ASSERT_EQ(10u, all.size());
    for (size_t i = 0; i < all.size(); ++i) {
        EXPECT_EQ(i < 5 ? "logan_a" : "logan_b", all[i].name) << i;
        EXPECT_EQ(indexes.at(all[i].name)[i % 5], all[i].pair) << i;
    }
    const std::vector<SelectedPair> one = select_on(indexes, SelectionRoute::PATTERN,
        R"({"pairs": [{"graph": "logan_b", "graph_path": "/scratch/nvme3/logan_b/graph.dbg"}]})");
    EXPECT_EQ(std::vector<SelectedPair>({ all[8] }), one);
    const std::vector<SelectedPair> mixed = select_on(indexes, SelectionRoute::RESOLVE, R"({"pairs": [
        {"graph": "logan_b", "graph_path": "/scratch/nvme0/logan_b/graph.dbg"},
        {"graph": "logan_a", "annotation_path": "/scratch/nvme4/logan_a/annotation.row_diff_int_brwt.annodbg"},
        {"graph": "logan_a", "graph_path": "/scratch/nvme1/logan_a/graph.dbg"}]})");
    EXPECT_EQ(std::vector<SelectedPair>({ all[1], all[4], all[5] }), mixed);
    // a graph path of the other name's
    Json::Value body;
    EXPECT_EQ("pair_not_found", selection_code_of([&] {
        select_on(indexes, SelectionRoute::RESOLVE, R"({"pairs": [{"graph": "logan_a",
                  "graph_path": "/scratch/nvme0/logan_b/graph.dbg"}]})");
    }, &body));
    EXPECT_EQ(5u, body["candidates"].size());
    EXPECT_EQ("logan_a", body["candidates"][0]["graph"].asString());
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select_on(indexes, SelectionRoute::TRAVERSE, R"({"graph": "logan_a"})");
    }, &body));
    EXPECT_EQ(refs({ all.begin(), all.begin() + 5 }), body["candidates"]);
    // more names than the cap: the cap's 400 without a selection, none with one
    EXPECT_NE("", refusal_of([&] { select_on(indexes, SelectionRoute::PATTERN, "{}", 1); }));
    EXPECT_EQ(5u, select_on(indexes, SelectionRoute::PATTERN, R"({"graphs": ["logan_b"]})", 1)
                          .size());
}

// sra-human (one graph, several annotations under one name): only annotation_path tells its
// pairs apart, in a pair object and in a probe; a traverse names one annotation
TEST(ServerSelection, DeployedSraHuman) {
    const GraphIndexes indexes = sra_human();
    const std::string name = "sra-human-1/3";
    const std::vector<SelectedPair> human = select_on(indexes, SelectionRoute::RESOLVE,
                                                      R"({"graph": "sra-human-1/3"})");
    ASSERT_EQ(3u, human.size());
    EXPECT_EQ(human, select_on(indexes, SelectionRoute::RESOLVE,
                               R"({"graph": "sra-human-1/3", "graph_path": ")" + kHuman + R"("})"));
    const std::string by_annotation = R"({"pairs": [{"graph": "sra-human-1/3", "annotation_path": ")"
                                    + human_annotation(3, 4) + R"("}]})";
    for (SelectionRoute route : { SelectionRoute::RESOLVE, SelectionRoute::PATTERN,
                                  SelectionRoute::TRAVERSE }) {
        EXPECT_EQ(std::vector<SelectedPair>({ human[1] }),
                  select_on(indexes, route, by_annotation));
    }
    EXPECT_EQ(std::vector<SelectedPair>({ human[2] }),
              select_pairs_from_parameters(SelectionRoute::PATTERN, {
                      { "graph", name }, { "annotation_path", human_annotation(4, 8) } },
                      indexes, 11));
    Json::Value body;
    EXPECT_EQ("traverse_single_pair", selection_code_of([&] {
        select_on(indexes, SelectionRoute::TRAVERSE,
                  R"({"graph": "sra-human-1/3", "graph_path": ")" + kHuman + R"("})");
    }, &body));
    EXPECT_EQ(refs(human), body["candidates"]);
    EXPECT_EQ("Bad request: graph_path '" + kHuman + "' names a graph with several annotations "
              "under index 'sra-human-1/3' (" + human_annotation(0, 0) + ", "
              + human_annotation(3, 4) + ", " + human_annotation(4, 8) + "): a traversal reads "
              "one annotation, which graph_path cannot choose; list each annotation under a name "
              "of its own", body["error"].asString());
    // with graph_path or without, a pair object is ambiguous until annotation_path is given
    for (const std::string &object : { std::string(R"({"graph": "sra-human-1/3"})"),
                                       R"({"graph": "sra-human-1/3", "graph_path": ")" + kHuman + R"("})" }) {
        EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
            select_on(indexes, SelectionRoute::PATTERN, R"({"pairs": [)" + object + "]}");
        }, &body)) << object;
        EXPECT_EQ(refs(human), body["candidates"]);
    }
    EXPECT_EQ("pair_ambiguous", selection_code_of([&] {
        select_pairs_from_parameters(SelectionRoute::TRAVERSE, { { "graph", name } }, indexes, 11);
    }, &body));
    EXPECT_EQ(refs(human), body["candidates"]);
}

// The refusals of sample 22, each rebuilt from its case's request on a list holding the
// sample's pairs (refseq of two): the same bodies, member for member
TEST(ServerSelection, Sample22Refusals) {
    GraphIndexes indexes = sra_human();
    indexes["refseq"] = { r237_pair(kR8), r237_pair(kR10) };
    const Json::Value cases = sample("22_selector_refusals.json")["body"];
    const std::vector<std::function<void()>> requests = {
        [&] { select_on(indexes, SelectionRoute::RESOLVE,
                        R"({"graphs": ["refseq"], "pairs": [{"graph": "refseq"}]})"); },
        [&] { select_on(indexes, SelectionRoute::PATTERN, R"({"pairs": [
                  {"graph": "bact_other"},
                  {"graph": "refseq", "graph_path": "/scratch/x/graph.dbg"}]})"); },
        [&] { select_on(indexes, SelectionRoute::RESOLVE, R"({"pairs": [
                  {"graph": "sra-human-1/3", "graph_path": ")" + kHuman + R"("}]})"); },
        [&] { select_on(indexes, SelectionRoute::TRAVERSE, R"({"graph": "refseq"})"); },
        [&] { select_on(indexes, SelectionRoute::RESOLVE, R"({"graphs": ["bact_othr"]})"); },
        [&] { select_on(indexes, SelectionRoute::SEARCH, R"({"pairs": [{"graph": "refseq"}]})"); },
        [&] { select_pairs_from_parameters(SelectionRoute::TRAVERSE,
                                           { { "graph", "sra-human-1/3" } }, indexes, 11); },
    };
    ASSERT_EQ(requests.size(), cases.size());
    for (size_t i = 0; i < cases.size(); ++i) {
        const Json::Value &c = cases[static_cast<Json::ArrayIndex>(i)];
        SCOPED_TRACE(c["case"].asString());
        EXPECT_EQ(400, c["status"].asInt());
        Json::Value body;
        EXPECT_EQ(c["body"]["code"].asString(), selection_code_of(requests[i], &body));
        EXPECT_EQ(c["body"], body);
    }
}
