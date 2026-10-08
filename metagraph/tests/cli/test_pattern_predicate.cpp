#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/pattern.hpp"
#include "cli/pattern_predicate.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"


// The predicate language of POST /pattern (src/cli/pattern_predicate.cpp, increment 5b-1;
// SPEC-pattern-search.md §19.3, §19.4, §19.9): every refusal of the parser with its path, the
// cap, the folding of unknown names into the normal form, vacuous and monotone predicates,
// the two- and three-valued evaluators, their units, the memory model and the echo's text.
// The expectations come from ORACLES written here over the request's JSON — an evaluator on
// sets of names and a folding of the table of §19.4 — never from the module.

namespace {

using namespace mtg;
using namespace mtg::cli;
using namespace mtg::cli::predicate;
using graph::pattern::Budget;
using graph::pattern::Deadline;
using Clock = Deadline::Clock;
using Names = std::set<std::string>;

Json::Value js(const std::string &text) { return parse_pattern_body(text); }

std::string compact(const Json::Value &value) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, value);
}

Budget unbounded() { return Budget(UINT64_MAX, Deadline::unbounded()); }

// the columns of a test index: |known| in order, column i the i-th
ColumnLookup columns(const std::vector<std::string> &known) {
    auto map = std::make_shared<std::map<std::string, Column>>();
    for (size_t i = 0; i < known.size(); ++i) {
        (*map)[known[i]] = i;
    }
    return [map](const std::string &name) -> std::optional<Column> {
        auto it = map->find(name);
        if (it == map->end())
            return std::nullopt;
        return it->second;
    };
}

Bound bound_json(const Json::Value &json, const std::vector<std::string> &known) {
    Budget budget = unbounded();
    Binding binding = Bound::bind(Predicate::parse(json, 1'000'000), columns(known), budget);
    EXPECT_EQ(nullptr, binding.stop);
    return std::move(*binding.bound);
}

Bound bound_of(const std::string &text, const std::vector<std::string> &known) {
    return bound_json(js(text), known);
}

// {code, message} of the refusal of |text|, or {"", ""}
std::pair<std::string, std::string> refusal(const std::string &text, uint64_t cap = 10'000,
                                            const std::string &path = "request.predicate") {
    try {
        Predicate::parse(js(text), cap, path);
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(400, e.status());
        return { e.code(), e.what() };
    }
    return { "", "" };
}


// ---------------------------------------------------------------- the oracles

const std::set<std::string> kLeafOps = { "any", "all", "none", "at_least" };

const Json::Value& leaf_list(const std::string &op, const Json::Value &arg) {
    return op == "at_least" ? arg["labels"] : arg;
}

// O-eval: the truth of a predicate (request syntax, or a folded one with constants) on the
// set |S| of the names a context carries
bool oracle_eval(const Json::Value &p, const Names &S) {
    if (p.isBool())
        return p.asBool();
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &arg = p[op];
    if (kLeafOps.count(op)) {
        const Json::Value &list = leaf_list(op, arg);
        size_t in = 0;
        for (const Json::Value &name : list) {
            in += S.count(name.asString());
        }
        if (op == "any")
            return in >= 1;
        if (op == "all")
            return in == list.size();
        if (op == "none")
            return in == 0;
        return in >= arg["n"].asUInt64();
    }
    if (op == "not")
        return !oracle_eval(arg, S);
    for (const Json::Value &operand : arg) {
        const bool v = oracle_eval(operand, S);
        if (op == "and" && !v)
            return false;
        if (op == "or" && v)
            return true;
    }
    return op == "and";
}

// O-fold: the normal form by the table of SPEC §19.4, |known| the names that are columns
Json::Value oracle_fold(const Json::Value &p, const Names &known) {
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &arg = p[op];
    if (kLeafOps.count(op)) {
        Json::Value kept(Json::arrayValue);
        bool unknown = false;
        for (const Json::Value &name : leaf_list(op, arg)) {
            if (known.count(name.asString())) {
                kept.append(name);
            } else {
                unknown = true;
            }
        }
        Json::Value out;
        if (op == "any") {
            if (kept.empty())
                return false;
            out["any"] = kept;
        } else if (op == "all") {
            if (unknown)
                return false;
            out["all"] = kept;
        } else if (op == "none") {
            if (kept.empty())
                return true;
            out["none"] = kept;
        } else {
            if (kept.size() < arg["n"].asUInt64())
                return false;
            out["at_least"]["n"] = arg["n"];
            out["at_least"]["labels"] = kept;
        }
        return out;
    }
    if (op == "not") {
        Json::Value operand = oracle_fold(arg, known);
        if (operand.isBool())
            return !operand.asBool();
        Json::Value out;
        out["not"] = operand;
        return out;
    }
    const bool absorbing = op == "or";
    std::vector<Json::Value> kept;
    for (const Json::Value &child : arg) {
        Json::Value operand = oracle_fold(child, known);
        if (operand.isBool()) {
            if (operand.asBool() == absorbing)
                return absorbing;
            continue;
        }
        kept.push_back(operand);
    }
    if (kept.empty())
        return !absorbing;
    if (kept.size() == 1)
        return kept[0];
    Json::Value out;
    out[op] = Json::Value(Json::arrayValue);
    for (Json::Value &operand : kept) {
        out[op].append(operand);
    }
    return out;
}

// the names of |p| in the order of their first appearance (document order), each once
void oracle_names(const Json::Value &p, std::vector<std::string> *out) {
    if (p.isBool())
        return;
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &arg = p[op];
    if (kLeafOps.count(op)) {
        for (const Json::Value &name : leaf_list(op, arg)) {
            if (std::find(out->begin(), out->end(), name.asString()) == out->end())
                out->push_back(name.asString());
        }
    } else if (op == "not") {
        oracle_names(arg, out);
    } else {
        for (const Json::Value &operand : arg) {
            oracle_names(operand, out);
        }
    }
}

// the operators of |p| (its nodes), the names its lists hold (a name in two lists twice), and
// whether one of them is none or not
struct Shape {
    uint64_t nodes = 0;
    uint64_t listed = 0;
    bool negative = false;
};

void oracle_shape(const Json::Value &p, Shape *shape) {
    if (p.isBool())
        return;
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &arg = p[op];
    ++shape->nodes;
    shape->negative = shape->negative || op == "none" || op == "not";
    if (kLeafOps.count(op)) {
        shape->listed += leaf_list(op, arg).size();
    } else if (op == "not") {
        oracle_shape(arg, shape);
    } else {
        for (const Json::Value &operand : arg) {
            oracle_shape(operand, shape);
        }
    }
}

// the units of deciding a context carrying |S| on the normal form |nf|: 1, and 1 per leaf a
// carried label is listed in
uint64_t oracle_units(const Json::Value &nf, const Names &S) {
    if (nf.isBool())
        return 1;
    uint64_t units = 1;
    std::function<void(const Json::Value &)> walk = [&](const Json::Value &p) {
        const std::string op = p.getMemberNames().at(0);
        const Json::Value &arg = p[op];
        if (kLeafOps.count(op)) {
            for (const Json::Value &name : leaf_list(op, arg)) {
                units += S.count(name.asString());
            }
        } else if (op == "not") {
            walk(arg);
        } else {
            for (const Json::Value &operand : arg) {
                walk(operand);
            }
        }
    };
    walk(nf);
    return units;
}

// the model of SPEC §19.9 over the normal form |nf| and the request's unknown names
uint64_t oracle_bytes(const Json::Value &nf, const std::vector<std::string> &unknown) {
    std::vector<std::string> labels;
    oracle_names(nf, &labels);
    Shape shape;
    oracle_shape(nf, &shape);
    uint64_t bytes = 64 * shape.nodes + 8 * shape.listed;
    for (const std::string &name : labels) {
        bytes += 192 + 2 * name.size();
    }
    for (const std::string &name : unknown) {
        bytes += 64 + name.size();
    }
    return bytes;
}

// the LabelIds of the names of |S| in the permitted set of |bound|
std::vector<LabelId> ids_of(const Bound &bound, const Names &S) {
    std::vector<LabelId> ids;
    for (LabelId l = 0; l < bound.labels().size(); ++l) {
        if (S.count(bound.labels()[l].name))
            ids.push_back(l);
    }
    return ids;
}


// ---------------------------------------------------------------- parsing

TEST(PatternPredicate, ParsesTheSevenOperators) {
    const Predicate p = Predicate::parse(js(R"({"or": [
            {"and": [{"any": ["562", "573"]}, {"none": ["287"]}]},
            {"not": {"all": ["546", "562"]}},
            {"at_least": {"n": 2, "labels": ["615", "546", "72407"]}}]})"), 10'000);

    // the distinct names in the order of their first appearance; the listed ones counted
    EXPECT_EQ((std::vector<std::string>{ "562", "573", "287", "546", "615", "72407" }),
              p.names);
    EXPECT_EQ(8u, p.listed);
    ASSERT_EQ(7u, p.nodes.size());
    EXPECT_EQ(Op::OR, p.nodes[p.root()].op);
    // children before their parent: every operand's id is below its node's
    for (uint32_t v = 0; v < p.nodes.size(); ++v) {
        const Predicate::Node &node = p.nodes[v];
        if (node.op == Op::AND || node.op == Op::OR || node.op == Op::NOT) {
            for (uint32_t i = node.begin; i < node.end; ++i) {
                EXPECT_LT(p.items[i], v);
            }
        }
    }
    const auto it = std::find_if(p.nodes.begin(), p.nodes.end(),
                                 [](const auto &node) { return node.op == Op::AT_LEAST; });
    ASSERT_NE(p.nodes.end(), it);
    EXPECT_EQ(2u, it->n);
    EXPECT_EQ(3u, it->end - it->begin);
}

TEST(PatternPredicate, RefusesEveryBrokenRuleWithItsPath) {
    const std::string ops = "(any, all, none, at_least, and, or, not)";
    const std::vector<std::pair<std::string, std::string>> cases = {
        // not an object, no operator, two, an unknown one
        { R"([{"any": ["562"]}])",
          "request.predicate: expected a predicate: an object with one operator " + ops },
        { R"("any")",
          "request.predicate: expected a predicate: an object with one operator " + ops },
        { R"(null)",
          "request.predicate: expected a predicate: an object with one operator " + ops },
        { R"({})", "request.predicate: expected one operator " + ops + ", got {}" },
        { R"({"any": ["562"], "all": ["573"]})",
          "request.predicate: expected exactly one operator " + ops
              + ", got 2 members ('all', 'any')" },
        { R"({"any": ["a"], "all": ["b"], "none": ["c"]})",
          "request.predicate: expected exactly one operator " + ops
              + ", got 3 members ('all', 'any', ...)" },
        { R"({"anyy": ["562"]})",
          "request.predicate: unknown operator 'anyy' (expected one of any, all, none, "
          "at_least, and, or, not)" },
        { R"({"ANY": ["562"]})",
          "request.predicate: unknown operator 'ANY' (expected one of any, all, none, "
          "at_least, and, or, not)" },
        // the lists
        { R"({"any": []})",
          "request.predicate.any: expected a list of at least one column label" },
        { R"({"none": "562"})", "request.predicate.none: expected a list of column labels" },
        { R"({"all": {"562": true}})",
          "request.predicate.all: expected a list of column labels" },
        // a number (P7: the message names the fix), another type, an empty string
        { R"({"any": [562]})",
          "request.predicate.any[0]: expected a string (a column label; write a taxid as "
          "\"562\")" },
        { R"({"and": [{"any": ["562"]}, {"none": ["287", 9606]}]})",
          "request.predicate.and[1].none[1]: expected a string (a column label; write a "
          "taxid as \"9606\")" },
        { R"({"any": [56.2]})",
          "request.predicate.any[0]: expected a string (a column label; write a taxid as "
          "\"562\")" },
        { R"({"any": [null]})", "request.predicate.any[0]: expected a string (a column label)" },
        { R"({"any": ["562", true]})",
          "request.predicate.any[1]: expected a string (a column label)" },
        { R"({"any": [["562"]]})",
          "request.predicate.any[0]: expected a string (a column label)" },
        { R"({"all": ["562", ""]})",
          "request.predicate.all[1]: expected a non-empty string (a column label)" },
        // a name twice in one list
        { R"({"any": ["562", "573", "562"]})",
          "request.predicate.any[2]: '562' listed twice in one list (first at [0])" },
        { R"({"at_least": {"n": 1, "labels": ["a", "b", "b"]}})",
          "request.predicate.at_least.labels[2]: 'b' listed twice in one list (first at [1])" },
        // at_least
        { R"({"at_least": [2, ["562", "573"]]})",
          "request.predicate.at_least: expected an object {\"n\": <integer>, \"labels\": "
          "[<column labels>]}" },
        { R"({"at_least": {"labels": ["562"]}})",
          "request.predicate.at_least.n: required: an integer from 1 to the labels listed" },
        { R"({"at_least": {"n": 1}})",
          "request.predicate.at_least.labels: required: a list of column labels" },
        { R"({"at_least": {"n": 0, "labels": ["562", "573"]}})",
          "request.predicate.at_least.n: expected an integer from 1 to 2 (the labels listed), "
          "got 0" },
        { R"({"at_least": {"n": -1, "labels": ["562"]}})",
          "request.predicate.at_least.n: expected an integer from 1 to 1 (the labels listed), "
          "got -1" },
        { R"({"at_least": {"n": 3, "labels": ["562", "573"]}})",
          "request.predicate.at_least.n: expected an integer from 1 to 2 (the labels listed), "
          "got 3" },
        { R"({"at_least": {"n": 1.5, "labels": ["562", "573"]}})",
          "request.predicate.at_least.n: expected an integer from 1 to the labels listed" },
        { R"({"at_least": {"n": "2", "labels": ["562", "573"]}})",
          "request.predicate.at_least.n: expected an integer from 1 to the labels listed" },
        { R"({"at_least": {"n": 1, "labels": []}})",
          "request.predicate.at_least.labels: expected a list of at least one column label" },
        { R"({"at_least": {"n": 1, "labels": ["562"], "of": 2}})",
          "request.predicate.at_least: unknown field 'of'" },
        // and, or, not
        { R"({"and": []})", "request.predicate.and: expected a list of at least one predicate" },
        { R"({"or": {"any": ["562"]}})", "request.predicate.or: expected a list of predicates" },
        { R"({"not": [{"any": ["562"]}]})",
          "request.predicate.not: expected a predicate: an object with one operator " + ops },
        { R"({"or": [{"all": ["a"]}, {"not": {"at_least": {"n": 2, "labels": ["a"]}}}]})",
          "request.predicate.or[1].not.at_least.n: expected an integer from 1 to 1 (the "
          "labels listed), got 2" },
        { R"({"and": [{"any": ["a"]}, {"or": [{"none": ["b"]}, {}]}]})",
          "request.predicate.and[1].or[1]: expected one operator " + ops + ", got {}" },
    };
    for (const auto &[text, message] : cases) {
        SCOPED_TRACE(text);
        EXPECT_EQ(std::make_pair(std::string("invalid_request"), message), refusal(text));
    }

    // the path is the caller's
    EXPECT_EQ("x.p.any: expected a list of column labels",
              refusal(R"({"any": 1})", 10'000, "x.p").second);
}

TEST(PatternPredicate, NestingAtMost64Deep) {
    // |depth| nested combinators around a leaf, cycling through not, and, or
    auto nested = [](unsigned depth) {
        std::string text = R"({"any": ["562"]})";
        for (unsigned d = 0; d < depth; ++d) {
            switch (d % 3) {
                case 0: text = R"({"not": )" + text + "}"; break;
                case 1: text = R"({"and": [)" + text + "]}"; break;
                case 2: text = R"({"or": [{"none": ["287"]}, )" + text + "]}"; break;
            }
        }
        return text;
    };
    EXPECT_EQ("", refusal(nested(64)).first);
    EXPECT_EQ(64u + 1 + 21, Predicate::parse(js(nested(64)), 10'000).nodes.size());

    const auto [code, message] = refusal(nested(65));
    EXPECT_EQ("invalid_request", code);
    // the outermost 64 combinators are accepted; the 65th, at the end of the path, is not
    EXPECT_NE(std::string::npos,
              message.find(": predicates nested more than 64 deep (and, or, not)"));
    EXPECT_EQ(0u, message.rfind("request.predicate.", 0));
    // 64 nots alone: accepted; 65: refused
    std::string nots = R"({"any": ["562"]})";
    for (int d = 0; d < 64; ++d) {
        nots = R"({"not": )" + nots + "}";
    }
    EXPECT_EQ("", refusal(nots).first);
    EXPECT_EQ("invalid_request", refusal(R"({"not": )" + nots + "}").first);
}

TEST(PatternPredicate, TheCapCountsTheListedNames) {
    // exactly the cap is accepted, one more is refused, naming the count and the cap
    const std::string five = R"({"and": [{"any": ["a", "b"]}, {"none": ["c", "d", "e"]}]})";
    EXPECT_EQ("", refusal(five, 5).first);
    EXPECT_EQ(5u, Predicate::parse(js(five), 5).listed);
    EXPECT_EQ(std::make_pair(std::string("predicate_too_large"),
                             std::string("request.predicate: 5 names in its lists (a name in "
                                         "two lists counts twice), above the server's cap of 4 "
                                         "(caps.max_predicate_labels): send at most 4 names; "
                                         "split the cohort")),
              refusal(five, 4));

    // a name in two lists counts twice (and is one distinct name)
    const std::string twice = R"({"or": [{"any": ["a", "b", "c"]}, {"all": ["a", "b", "c"]}]})";
    EXPECT_EQ("predicate_too_large", refusal(twice, 5).first);
    const Predicate p = Predicate::parse(js(twice), 6);
    EXPECT_EQ(6u, p.listed);
    EXPECT_EQ(3u, p.names.size());

    // the form first, then the size (SPEC §5, step 8)
    EXPECT_EQ("invalid_request", refusal(R"({"any": ["a", "b", ""]})", 1).first);
    EXPECT_EQ("invalid_request", refusal(R"({"any": ["a", "b", "c"], "x": 1})", 1).first);
    // a cap of 0 refuses every predicate by its size
    EXPECT_EQ("predicate_too_large", refusal(R"({"any": ["a"]})", 0).first);
    // past the cap a name's form is still checked, its repetition in one list no longer: the
    // walk keeps no more than the cap's names
    EXPECT_EQ("invalid_request", refusal(R"({"any": ["a", "b", "b"]})", 3).first);
    EXPECT_EQ("predicate_too_large", refusal(R"({"any": ["a", "b", "b"]})", 2).first);
    EXPECT_EQ("predicate_too_large", refusal(R"({"any": ["a", "b", "c", "c"]})", 2).first);
    EXPECT_EQ("invalid_request", refusal(R"({"any": ["a", "b", "c", 5]})", 2).first);
    EXPECT_EQ("predicate_too_large",
              refusal(R"({"or": [{"any": ["a"]}, {"at_least": {"n": 3, "labels": ["b", "c", "c"]}}]})",
                      2).first);
    EXPECT_EQ("invalid_request",
              refusal(R"({"or": [{"any": ["a"]}, {"at_least": {"n": 4, "labels": ["b", "c", "d"]}}]})",
                      2).first);

    // 10,001 names against the default cap of 10,000
    Json::Value big;
    big["any"] = Json::Value(Json::arrayValue);
    for (int i = 0; i < 10'000; ++i) {
        big["any"].append(std::to_string(i));
    }
    EXPECT_EQ(10'000u, Predicate::parse(big, 10'000).listed);
    big["any"].append("x");
    try {
        Predicate::parse(big, 10'000);
        ADD_FAILURE() << "accepted above the cap";
    } catch (const PatternRefusal &e) {
        EXPECT_EQ("predicate_too_large", e.code());
        EXPECT_EQ(0u, std::string(e.what()).rfind("request.predicate: 10001 names", 0));
    }
}


// ---------------------------------------------------------------- folding

TEST(PatternPredicate, FoldsUnknownNamesByTheTable) {
    // c0, c1, c2 are columns; u0, u1 are not
    const std::vector<std::string> known = { "c0", "c1", "c2" };
    const Names known_set(known.begin(), known.end());
    // {predicate, its normal form}: every form with no, one and every name unknown
    const std::vector<std::pair<std::string, std::string>> cases = {
        { R"({"any": ["c0", "c1"]})", R"({"any": ["c0", "c1"]})" },
        { R"({"any": ["c0", "u0"]})", R"({"any": ["c0"]})" },
        { R"({"any": ["u0", "u1"]})", "false" },
        { R"({"all": ["c0", "c1"]})", R"({"all": ["c0", "c1"]})" },
        { R"({"all": ["c0", "u0"]})", "false" },
        { R"({"all": ["u0", "u1"]})", "false" },
        { R"({"none": ["c0", "c1"]})", R"({"none": ["c0", "c1"]})" },
        { R"({"none": ["u0", "c1"]})", R"({"none": ["c1"]})" },
        { R"({"none": ["u0", "u1"]})", "true" },
        { R"({"at_least": {"n": 2, "labels": ["c0", "c1", "c2"]}})",
          R"({"at_least": {"n": 2, "labels": ["c0", "c1", "c2"]}})" },
        { R"({"at_least": {"n": 2, "labels": ["c0", "u0", "c2"]}})",
          R"({"at_least": {"n": 2, "labels": ["c0", "c2"]}})" },
        { R"({"at_least": {"n": 3, "labels": ["c0", "u0", "c2"]}})", "false" },
        { R"({"at_least": {"n": 1, "labels": ["u0", "u1"]}})", "false" },
        // and, or, not of constants
        { R"({"and": [{"any": ["u0"]}, {"any": ["c0"]}]})", "false" },
        { R"({"and": [{"none": ["u0"]}, {"any": ["c0"]}]})", R"({"any": ["c0"]})" },
        { R"({"and": [{"none": ["u0"]}, {"none": ["u1"]}]})", "true" },
        { R"({"and": [{"any": ["c0"]}, {"none": ["c1"]}]})",
          R"({"and": [{"any": ["c0"]}, {"none": ["c1"]}]})" },
        { R"({"and": [{"any": ["c0"]}]})", R"({"any": ["c0"]})" },
        { R"({"or": [{"any": ["u0"]}, {"any": ["c0"]}]})", R"({"any": ["c0"]})" },
        { R"({"or": [{"none": ["u0"]}, {"any": ["c0"]}]})", "true" },
        { R"({"or": [{"any": ["u0"]}, {"all": ["c0", "u1"]}]})", "false" },
        { R"({"or": [{"any": ["c0"]}, {"any": ["c1"]}, {"any": ["u0"]}]})",
          R"({"or": [{"any": ["c0"]}, {"any": ["c1"]}]})" },
        { R"({"not": {"any": ["u0"]}})", "true" },
        { R"({"not": {"none": ["u0"]}})", "false" },
        { R"({"not": {"any": ["c0", "u0"]}})", R"({"not": {"any": ["c0"]}})" },
        { R"({"not": {"not": {"all": ["u0"]}}})", "false" },
        // nested: what the constant absorbs is gone, a single operand stands for its node
        { R"({"or": [{"and": [{"any": ["c0"]}, {"all": ["u0", "c1"]}]}, {"none": ["c2"]}]})",
          R"({"none": ["c2"]})" },
        { R"({"and": [{"or": [{"any": ["u0"]}, {"not": {"all": ["c0", "c1"]}}]},
                      {"at_least": {"n": 1, "labels": ["u1", "c2"]}}]})",
          R"({"and": [{"not": {"all": ["c0", "c1"]}}, {"at_least": {"n": 1, "labels": ["c2"]}}]})" },
    };
    for (const auto &[text, normal] : cases) {
        SCOPED_TRACE(text);
        const Json::Value expected = js(normal);
        // the oracle's folding agrees with the table written out here
        EXPECT_EQ(compact(expected), compact(oracle_fold(js(text), known_set)));

        const Bound b = bound_of(text, known);
        EXPECT_EQ(compact(expected), compact(b.normal_form_json()));
        if (expected.isBool()) {
            ASSERT_TRUE(b.constant());
            EXPECT_EQ(expected.asBool(), *b.constant());
            EXPECT_EQ(0u, b.num_nodes());
            EXPECT_TRUE(b.labels().empty());
        } else {
            EXPECT_FALSE(b.constant());
        }
    }
}

TEST(PatternPredicate, StatesTheNamesKnownAndUnknown) {
    // SPEC §19.13 (c): a typo folds to a constant and is listed
    {
        const Bound b = bound_of(R"({"any": ["5622"]})", { "562", "573" });
        EXPECT_EQ("false", compact(b.normal_form_json()));
        EXPECT_EQ(1u, b.num_names());
        EXPECT_EQ(0u, b.num_known());
        EXPECT_EQ(R"(["5622"])", compact(b.unknown_labels_json()));
        EXPECT_FALSE(b.vacuous());
    }
    {
        const Bound b = bound_of(R"({"and": [{"any": ["562"]}, {"none": ["5622"]}]})",
                                 { "562", "573" });
        EXPECT_EQ(R"({"any":["562"]})", compact(b.normal_form_json()));
        EXPECT_EQ(R"(["5622"])", compact(b.unknown_labels_json()));
        EXPECT_EQ(2u, b.num_names());
        EXPECT_EQ(1u, b.num_known());
    }
    // known names the folding dropped count as known, are not in the permitted set; the
    // unknown ones are listed in the order of their first appearance
    {
        const Bound b = bound_of(
                R"({"or": [{"all": ["c1", "u1", "c0"]}, {"any": ["u0", "c2", "u1"]}, {"none": ["c0"]}]})",
                { "c0", "c1", "c2" });
        EXPECT_EQ(R"({"or":[{"any":["c2"]},{"none":["c0"]}]})", compact(b.normal_form_json()));
        EXPECT_EQ(5u, b.num_names());
        EXPECT_EQ(3u, b.num_known());
        EXPECT_EQ(R"(["u1","u0"])", compact(b.unknown_labels_json()));
        // the permitted set: the normal form's columns in the order of first appearance in
        // the request (c0 appears before c2 there)
        ASSERT_EQ(2u, b.labels().size());
        EXPECT_EQ("c0", b.labels()[0].name);
        EXPECT_EQ(0u, b.labels()[0].column);
        EXPECT_EQ("c2", b.labels()[1].name);
        EXPECT_EQ(2u, b.labels()[1].column);
        EXPECT_EQ(graph::traversal::LabelKind::COLUMN, b.labels()[1].kind);
    }
}

TEST(PatternPredicate, VacuousAndMonotone) {
    const std::vector<std::string> known = { "A", "B", "C" };
    // {predicate, vacuous, monotone}
    const std::vector<std::tuple<std::string, bool, bool>> cases = {
        { R"({"none": ["A"]})", true, false },
        { R"({"not": {"any": ["A"]}})", true, false },
        { R"({"or": [{"none": ["A"]}, {"any": ["B"]}]})", true, false },
        { R"({"and": [{"none": ["A"]}, {"any": ["B"]}]})", false, false },
        { R"({"at_least": {"n": 1, "labels": ["A", "B"]}})", false, true },
        { R"({"at_least": {"n": 2, "labels": ["A", "B", "C"]}})", false, true },
        { R"({"any": ["A"]})", false, true },
        { R"({"all": ["A", "B"]})", false, true },
        { R"({"or": [{"all": ["A", "B"]}, {"and": [{"any": ["C"]}, {"any": ["A"]}]}]})", false,
          true },
        { R"({"not": {"not": {"any": ["A"]}}})", false, false },
        { R"({"not": {"none": ["A"]}})", false, false },
        // folded: the none of unknowns is gone, the rest is monotone
        { R"({"and": [{"any": ["A"]}, {"none": ["X"]}]})", false, true },
        // constants: true holds on the empty set; both are monotone
        { R"({"none": ["X"]})", true, true },
        { R"({"any": ["X"]})", false, true },
    };
    for (const auto &[text, vacuous, monotone] : cases) {
        SCOPED_TRACE(text);
        const Bound b = bound_of(text, known);
        EXPECT_EQ(vacuous, b.vacuous());
        EXPECT_EQ(monotone, b.monotone());
        // vacuous: the oracle on a context carrying no label
        EXPECT_EQ(vacuous, oracle_eval(js(text), {}));
    }
}


// ---------------------------------------------------------------- evaluation

// A random predicate (request syntax) of at most |depth| nested combinators over |universe|
Json::Value random_predicate(std::mt19937 &rng, const std::vector<std::string> &universe,
                             unsigned depth) {
    std::uniform_int_distribution<int> coin(0, 9);
    Json::Value p;
    if (depth > 0 && coin(rng) < 6) {
        switch (rng() % 3) {
            case 0:
                p["not"] = random_predicate(rng, universe, depth - 1);
                break;
            default: {
                const char *op = rng() % 2 ? "and" : "or";
                p[op] = Json::Value(Json::arrayValue);
                const int operands = 1 + rng() % 3;
                for (int i = 0; i < operands; ++i) {
                    p[op].append(random_predicate(rng, universe, depth - 1));
                }
            }
        }
        return p;
    }
    std::vector<std::string> names = universe;
    std::shuffle(names.begin(), names.end(), rng);
    names.resize(1 + rng() % 4);
    Json::Value list(Json::arrayValue);
    for (const std::string &name : names) {
        list.append(name);
    }
    switch (rng() % 4) {
        case 0: p["any"] = list; break;
        case 1: p["all"] = list; break;
        case 2: p["none"] = list; break;
        default:
            p["at_least"]["n"] = static_cast<Json::UInt>(1 + rng() % names.size());
            p["at_least"]["labels"] = list;
    }
    return p;
}

// up to 20 names over a universe of 12, 8 of them columns, depth at most 6
struct RandomCase {
    Json::Value json;
    Json::Value normal;      // the oracle's
    std::vector<std::string> unknown;
};

const std::vector<std::string> kKnown = { "c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7" };
const std::vector<std::string> kUniverse
        = { "c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "u0", "u1", "u2", "u3" };

RandomCase random_case(std::mt19937 &rng) {
    const Names known(kKnown.begin(), kKnown.end());
    while (true) {
        RandomCase c;
        c.json = random_predicate(rng, kUniverse, 6);
        Shape shape;
        oracle_shape(c.json, &shape);
        if (shape.listed > 20)
            continue;
        c.normal = oracle_fold(c.json, known);
        std::vector<std::string> names;
        oracle_names(c.json, &names);
        for (const std::string &name : names) {
            if (!known.count(name))
                c.unknown.push_back(name);
        }
        return c;
    }
}

Names random_subset(std::mt19937 &rng, const std::vector<std::string> &from) {
    Names s;
    for (const std::string &name : from) {
        if (rng() % 2)
            s.insert(name);
    }
    return s;
}

TEST(PatternPredicate, EvaluatesAsTheOracle) {
    std::mt19937 rng(20261008);
    uint64_t trues = 0;
    uint64_t falses = 0;
    for (int t = 0; t < 10'000; ++t) {
        const RandomCase c = random_case(rng);
        SCOPED_TRACE(compact(c.json));
        const Bound b = bound_json(c.json, kKnown);

        // the normal form, the permitted set, the echo and the model against the oracles
        ASSERT_EQ(compact(c.normal), compact(b.normal_form_json()));
        std::vector<std::string> labels;
        oracle_names(c.normal, &labels);
        std::vector<std::string> request_order;
        oracle_names(c.json, &request_order);
        request_order.erase(std::remove_if(request_order.begin(), request_order.end(),
                                           [&](const std::string &name) {
                                               return std::find(labels.begin(), labels.end(),
                                                                name) == labels.end();
                                           }),
                            request_order.end());
        ASSERT_EQ(request_order.size(), b.labels().size());
        for (size_t l = 0; l < request_order.size(); ++l) {
            EXPECT_EQ(request_order[l], b.labels()[l].name);
        }
        Json::Value unknown(Json::arrayValue);
        for (const std::string &name : c.unknown) {
            unknown.append(name);
        }
        EXPECT_EQ(compact(unknown), compact(b.unknown_labels_json()));
        EXPECT_EQ(compact(b.normal_form_json()).size() + compact(unknown).size(),
                  b.text_bytes());
        EXPECT_EQ(oracle_bytes(c.normal, c.unknown), b.bytes());
        Shape shape;
        oracle_shape(c.normal, &shape);
        EXPECT_EQ(!shape.negative, b.monotone());
        EXPECT_EQ(oracle_eval(c.normal, {}), b.vacuous());
        EXPECT_EQ(shape.nodes, b.num_nodes());

        // the folding keeps the meaning: the module on the permitted labels of a set equals
        // the oracle on the request's predicate and the whole set (columns the folding
        // dropped included)
        for (int s = 0; s < 16; ++s) {
            const Names S = random_subset(rng, kKnown);
            const std::vector<LabelId> present = ids_of(b, S);
            uint64_t units = 0;
            const bool value = b.eval(present.data(), present.size(), &units);
            ASSERT_EQ(oracle_eval(c.json, S), value) << "S of " << S.size();
            EXPECT_EQ(oracle_units(c.normal, S), units);
            EXPECT_EQ(units, b.units(present.data(), present.size()));
            (value ? trues : falses)++;
            // every label listed twice changes nothing
            std::vector<LabelId> doubled = present;
            doubled.insert(doubled.end(), present.rbegin(), present.rend());
            EXPECT_EQ(value, b.eval(doubled.data(), doubled.size(), &units));
            EXPECT_EQ(oracle_units(c.normal, S), units);
        }

        // a monotone predicate cannot turn true when labels are taken away
        if (b.monotone()) {
            for (int s = 0; s < 4; ++s) {
                const Names T = random_subset(rng, kKnown);
                Names S;
                for (const std::string &name : T) {
                    if (rng() % 2)
                        S.insert(name);
                }
                const auto s_ids = ids_of(b, S);
                const auto t_ids = ids_of(b, T);
                EXPECT_LE(b.eval(s_ids.data(), s_ids.size()), b.eval(t_ids.data(), t_ids.size()));
            }
        }

        // the normal form parsed again is its own normal form
        if (!c.normal.isBool()) {
            EXPECT_EQ(compact(c.normal),
                      compact(bound_json(c.normal, kKnown).normal_form_json()));
        }
    }
    // both answers occur often
    EXPECT_GT(trues, 20'000u);
    EXPECT_GT(falses, 20'000u);
}

TEST(PatternPredicate, KleeneNeverContradictsACompletion) {
    std::mt19937 rng(8102026);
    uint64_t definite = 0;
    uint64_t open = 0;
    for (int t = 0; t < 10'000; ++t) {
        const RandomCase c = random_case(rng);
        SCOPED_TRACE(compact(c.json));
        const Bound b = bound_json(c.json, kKnown);
        for (int s = 0; s < 4; ++s) {
            // each permitted label sure, maybe (at most 4) or absent
            std::vector<LabelId> sure, maybe;
            Names sure_names;
            for (LabelId l = 0; l < b.labels().size(); ++l) {
                switch (rng() % 3) {
                    case 0:
                        sure.push_back(l);
                        sure_names.insert(b.labels()[l].name);
                        break;
                    case 1:
                        if (maybe.size() < 4) {
                            maybe.push_back(l);
                            break;
                        }
                        [[fallthrough]];
                    default:
                        break;
                }
            }
            uint64_t units = 0;
            const std::optional<bool> value
                    = b.eval3(sure.data(), sure.size(), maybe.data(), maybe.size(), &units);
            std::vector<LabelId> both = sure;
            both.insert(both.end(), maybe.begin(), maybe.end());
            EXPECT_EQ(b.units(both.data(), both.size()), units);

            // every completion of the maybes
            std::set<bool> outcomes;
            for (uint32_t mask = 0; mask < (1u << maybe.size()); ++mask) {
                Names S = sure_names;
                for (size_t i = 0; i < maybe.size(); ++i) {
                    if (mask >> i & 1)
                        S.insert(b.labels()[maybe[i]].name);
                }
                outcomes.insert(oracle_eval(c.json, S));
            }
            if (value) {
                ++definite;
                ASSERT_EQ((std::set<bool>{ *value }), outcomes);
            } else {
                ++open;
            }
            // a leaf alone is decided exactly; without maybes the answer is eval's
            if (!b.constant() && b.num_nodes() == 1) {
                EXPECT_EQ(outcomes.size() == 2, !value);
            }
            EXPECT_EQ(b.eval(sure.data(), sure.size()),
                      *b.eval3(sure.data(), sure.size(), nullptr, 0));
            // a label both sure and maybe is sure
            EXPECT_EQ(b.eval(sure.data(), sure.size()),
                      *b.eval3(sure.data(), sure.size(), sure.data(), sure.size()));
        }
    }
    EXPECT_GT(definite, 10'000u);
    EXPECT_GT(open, 1'000u);
}

TEST(PatternPredicate, KleeneOnCombinators) {
    const Bound b = bound_of(R"({"or": [{"any": ["A"]}, {"none": ["A"]}]})", { "A", "B" });
    const LabelId a = 0;
    // always true, but not decided with A in doubt: sound, not complete
    EXPECT_EQ(std::nullopt, b.eval3(nullptr, 0, &a, 1));
    EXPECT_EQ(true, b.eval3(&a, 1, nullptr, 0));

    const Bound m = bound_of(R"({"and": [{"any": ["A", "B"]}, {"all": ["A", "B"]}]})", { "A", "B" });
    const LabelId ids[] = { 0, 1 };
    EXPECT_EQ(false, m.eval3(nullptr, 0, nullptr, 0));
    EXPECT_EQ(std::nullopt, m.eval3(&ids[0], 1, &ids[1], 1));
    EXPECT_EQ(true, m.eval3(ids, 2, nullptr, 0));
    // the pruning question of §20.9: false once B is gone for sure
    EXPECT_EQ(false, m.eval3(&ids[0], 1, nullptr, 0));
}

TEST(PatternPredicate, UnitsAndRefusedIds) {
    // B is listed in two leaves
    const Bound b = bound_of(R"({"and": [{"any": ["A", "B"]}, {"none": ["B", "C"]}]})",
                             { "A", "B", "C" });
    ASSERT_EQ(3u, b.labels().size());
    const LabelId A = 0, B = 1, C = 2;
    EXPECT_EQ(1u, b.units(nullptr, 0));
    EXPECT_EQ(2u, b.units(&A, 1));
    EXPECT_EQ(3u, b.units(&B, 1));
    const LabelId all[] = { C, A, B, A, B };
    EXPECT_EQ(5u, b.units(all, 5));
    uint64_t units = 0;
    EXPECT_FALSE(b.eval(all, 5, &units));
    EXPECT_EQ(5u, units);
    EXPECT_TRUE(b.eval(&A, 1, &units));
    EXPECT_EQ(2u, units);

    // an id outside the permitted set is refused, and the scratch stays as it was
    const LabelId bad[] = { A, 3 };
    EXPECT_THROW(b.eval(bad, 2), std::invalid_argument);
    EXPECT_THROW(b.units(bad, 2), std::invalid_argument);
    EXPECT_THROW(b.eval3(&A, 1, &bad[1], 1), std::invalid_argument);
    EXPECT_TRUE(b.eval(&A, 1, &units));
    EXPECT_EQ(2u, units);
    EXPECT_FALSE(b.eval(nullptr, 0));

    // a constant: one unit, no label
    const Bound t = bound_of(R"({"none": ["X"]})", { "A" });
    EXPECT_TRUE(t.eval(nullptr, 0, &units));
    EXPECT_EQ(1u, units);
    EXPECT_EQ(true, t.eval3(nullptr, 0, nullptr, 0));
    EXPECT_THROW(t.eval(&A, 1), std::invalid_argument);
}

TEST(PatternPredicate, ALargePredicateEvaluatesWhatItsLabelsTouch) {
    // 5,000 leaves under one or: a context's units are its labels' leaves, not the tree
    Json::Value p;
    p["or"] = Json::Value(Json::arrayValue);
    std::vector<std::string> known;
    for (int i = 0; i < 5'000; ++i) {
        known.push_back("n" + std::to_string(i));
        Json::Value leaf;
        leaf["all"].append(known.back());
        p["or"].append(leaf);
    }
    const Bound b = bound_json(p, known);
    ASSERT_EQ(5'000u, b.labels().size());
    EXPECT_EQ(5'001u, b.num_nodes());
    uint64_t units = 0;
    const LabelId last = 4'999;
    EXPECT_TRUE(b.eval(&last, 1, &units));
    EXPECT_EQ(2u, units);
    EXPECT_FALSE(b.eval(nullptr, 0, &units));
    EXPECT_EQ(1u, units);
}


// ---------------------------------------------------------------- the answer and the model

TEST(PatternPredicate, TheModelAndTheEchoText) {
    // the bytes of SPEC §19.9, by hand: labels 562 and 287 (3 bytes each), the unknown 5622
    // (4 bytes), three nodes, three names listed
    const Bound b = bound_of(
            R"({"and": [{"any": ["562", "287"]}, {"none": ["287", "5622"]}]})", { "562", "287" });
    EXPECT_EQ(R"({"and":[{"any":["562","287"]},{"none":["287"]}]})", compact(b.normal_form_json()));
    EXPECT_EQ(2 * (192 + 2 * 3) + (64 + 4) + 3 * 64 + 3 * 8, b.bytes());
    EXPECT_EQ(compact(b.normal_form_json()).size() + compact(b.unknown_labels_json()).size(),
              b.text_bytes());

    // at_least's inner object, a constant, an empty unknown list
    const Bound a = bound_of(R"({"at_least": {"n": 12, "labels": ["a", "b", "c", "d", "e",
            "f", "g", "h", "i", "j", "k", "l"]}})",
            { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l" });
    EXPECT_EQ(compact(a.normal_form_json()).size() + 2, a.text_bytes());
    const Bound f = bound_of(R"({"any": ["X", "Y"]})", { "a" });
    EXPECT_EQ(R"(false)", compact(f.normal_form_json()));
    EXPECT_EQ(5 + compact(f.unknown_labels_json()).size(), f.text_bytes());
    EXPECT_EQ(2 * (64 + 1), f.bytes());

    // names jsoncpp escapes: exact for ASCII, an upper bound beyond it
    const std::vector<std::string> odd = { "a\"b", "c\\d", "tab\there", "nl\n", "ctl\x01",
                                           "del\x7f" };
    Json::Value list(Json::arrayValue);
    for (const std::string &name : odd) {
        list.append(name);
    }
    Json::Value p;
    p["any"] = list;
    const Bound e = bound_json(p, odd);
    EXPECT_EQ(compact(e.normal_form_json()).size() + 2, e.text_bytes());

    Json::Value q;
    q["any"].append("é");
    q["any"].append("Escherichia coli 🦠");
    const Bound u = bound_json(q, { "é" });
    const uint64_t written = compact(u.normal_form_json()).size()
            + compact(u.unknown_labels_json()).size();
    EXPECT_GE(u.text_bytes(), written);
    EXPECT_LE(u.text_bytes(), 6 * written);
}

TEST(PatternPredicate, TheEchoChecksItsDeliveryEvery4096Names) {
    Json::Value p;
    std::vector<std::string> known;
    for (int i = 0; i < 10'000; ++i) {
        known.push_back(std::to_string(i));
        p["none"].append(known.back());
    }
    for (int i = 0; i < 5'000; ++i) {
        p["none"].append("x" + std::to_string(i));
    }
    const Bound b = bound_json(p, known);
    // one node and 10,000 names written: checked before the first and every 4,096 after
    int checks = 0;
    const Json::Value nf = b.normal_form_json([&]() { ++checks; });
    EXPECT_EQ(10'000u, nf["none"].size());
    EXPECT_EQ(3, checks);
    checks = 0;
    EXPECT_EQ(5'000u, b.unknown_labels_json([&]() { ++checks; }).size());
    EXPECT_EQ(2, checks);
    // a check that throws ends the echo
    EXPECT_THROW(b.normal_form_json([]() { throw std::runtime_error("deadline"); }),
                 std::runtime_error);

    // 4,097 writes (a node and 4,096 names): before the first and before the 4,097th
    Json::Value q;
    for (int i = 0; i < 4'096; ++i) {
        q["any"].append(known[i]);
    }
    q["any"].append("y");
    const Bound c = bound_json(q, known);
    checks = 0;
    c.normal_form_json([&]() { ++checks; });
    EXPECT_EQ(2, checks);
    // a single unknown name: checked once, before it
    checks = 0;
    c.unknown_labels_json([&]() { ++checks; });
    EXPECT_EQ(1, checks);
}


// ---------------------------------------------------------------- binding

// A clock that reads |start| for the first |good| readings and 10 s after it from then on
std::function<Clock::time_point()> expiring_clock(Clock::time_point start, int good,
                                                  std::shared_ptr<int> readings) {
    return [=]() {
        return (*readings)++ < good ? start : start + std::chrono::seconds(10);
    };
}

TEST(PatternPredicate, BindingReadsTheClockEvery4096Names) {
    Json::Value p;
    std::vector<std::string> known;
    for (int i = 0; i < 6'000; ++i) {
        known.push_back(std::to_string(i));
        p["any"].append(known.back());
    }
    const Predicate predicate = Predicate::parse(p, 10'000);
    const ColumnLookup lookup = columns(known);

    // the work time passes after the first reading: 4,096 names resolved, then the stop
    {
        const Clock::time_point start = Clock::now();
        auto readings = std::make_shared<int>(0);
        Budget budget(UINT64_MAX, Deadline(start, 1'000, 0, expiring_clock(start, 1, readings)));
        uint64_t lookups = 0;
        Binding binding = Bound::bind(predicate, [&](const std::string &name) {
            ++lookups;
            return lookup(name);
        }, budget);
        EXPECT_STREQ("time", binding.stop);
        EXPECT_FALSE(binding.bound);
        EXPECT_EQ(0u, binding.bytes);
        EXPECT_EQ(0u, binding.admitted);
        EXPECT_EQ(4'096u, lookups);
        EXPECT_EQ(2, *readings);
        EXPECT_EQ(graph::pattern::StopReason::TIME, budget.stopped());
    }
    // every reading before the work: each 4,096 units (names, nodes, items) of binding
    {
        const Clock::time_point start = Clock::now();
        auto readings = std::make_shared<int>(0);
        Budget budget(UINT64_MAX, Deadline(start, 1'000, 0, expiring_clock(start, 1'000,
                                                                            readings)));
        Binding binding = Bound::bind(predicate, lookup, budget);
        ASSERT_TRUE(binding.bound);
        // 6,000 names resolved, 6,001 nodes and items folded, 6,001 reached, 6,000 names
        // priced, 6,000 copied, 6,001 nodes and items built, the 6,000 names of the leaf
        // counted and placed in the labels' lists of leaves, one base value: 48,004 units, a
        // reading before the first and before each 4,096 after
        EXPECT_EQ(12, *readings);
    }
}

TEST(PatternPredicate, BindingAdmitsItsBytesBeforeItCopies) {
    const Json::Value p = js(R"({"and": [{"any": ["562", "287"]}, {"none": ["287", "5622"]}]})");
    const Predicate predicate = Predicate::parse(p, 10);
    const uint64_t expected = oracle_bytes(oracle_fold(p, { "562", "287" }), { "5622" });
    {
        Budget budget = unbounded();
        std::vector<uint64_t> asked;
        Binding binding = Bound::bind(predicate, columns({ "562", "287" }), budget,
                               [&](uint64_t bytes) { asked.push_back(bytes); return false; });
        EXPECT_STREQ("max_memory", binding.stop);
        EXPECT_FALSE(binding.bound);
        EXPECT_EQ(expected, binding.bytes);
        EXPECT_EQ(0u, binding.admitted);
        EXPECT_EQ(std::vector<uint64_t>{ expected }, asked);
    }
    {
        Budget budget = unbounded();
        Binding binding = Bound::bind(predicate, columns({ "562", "287" }), budget,
                               [&](uint64_t) { return true; });
        EXPECT_EQ(nullptr, binding.stop);
        ASSERT_TRUE(binding.bound);
        EXPECT_EQ(expected, binding.bytes);
        EXPECT_EQ(expected, binding.admitted);
        EXPECT_EQ(expected, binding.bound->bytes());
    }
    // the time passing while it builds (after the admission): stopped, what was admitted
    // stated for the caller to release
    {
        Json::Value big;
        std::vector<std::string> known;
        for (int i = 0; i < 1'000; ++i) {
            known.push_back(std::to_string(i));
            big["any"].append(known.back());
        }
        // 1,000 names resolved, 1,001 folded, 1,001 reached, 1,000 priced: 4,002 units, one
        // reading; the copy of the names crosses the next
        const Clock::time_point start = Clock::now();
        auto readings = std::make_shared<int>(0);
        Budget budget(UINT64_MAX, Deadline(start, 1'000, 0, expiring_clock(start, 1, readings)));
        bool admitted = false;
        Binding binding = Bound::bind(Predicate::parse(big, 10'000), columns(known), budget,
                               [&](uint64_t) { admitted = true; return true; });
        EXPECT_TRUE(admitted);
        EXPECT_STREQ("time", binding.stop);
        EXPECT_FALSE(binding.bound);
        EXPECT_GT(binding.bytes, 0u);
        EXPECT_EQ(binding.bytes, binding.admitted);
    }
    // nothing to admit to: nothing admitted, nothing for the caller to release
    {
        Budget budget = unbounded();
        Binding binding = Bound::bind(predicate, columns({ "562", "287" }), budget);
        ASSERT_TRUE(binding.bound);
        EXPECT_GT(binding.bytes, 0u);
        EXPECT_EQ(0u, binding.admitted);
    }
}

TEST(PatternPredicate, BindsToTheColumnsOfAnIndexNotItsHeaders) {
    // three records, one per column, with a record mapping: a header is a label of
    // /traverse, never a predicate term (P20)
    const size_t k = 5;
    const std::vector<std::string> seqs = { "ACGTTGCAAGT", "TTGACCATGGA", "GGCATCCATTA" };
    const std::vector<std::string> cols = { "562", "573", "287" };
    auto anno = test::build_anno_graph<graph::DBGSuccinct, annot::RowDiffColumnAnnotator>(
            k, seqs, cols, graph::DeBruijnGraph::BASIC, true, { 0, 0, 0 });
    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(encoder.size());
    std::vector<std::vector<uint64_t>> num_kmers(encoder.size());
    for (size_t i = 0; i < seqs.size(); ++i) {
        headers[encoder.encode(cols[i])].push_back("NZ_CP00" + std::to_string(i) + ".1");
        num_kmers[encoder.encode(cols[i])].push_back(seqs[i].size() - k + 1);
    }
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    graph::traversal::LabelOracle oracle(*anno, &cth);
    ASSERT_TRUE(oracle.find_header("NZ_CP000.1"));

    const Predicate predicate = Predicate::parse(
            js(R"({"and": [{"any": ["562", "NZ_CP000.1"]}, {"none": ["5622", "287"]}]})"), 10);
    Budget budget = unbounded();
    Binding binding = Bound::bind(predicate, oracle, budget);
    ASSERT_TRUE(binding.bound);
    const Bound &b = *binding.bound;
    EXPECT_EQ(R"({"and":[{"any":["562"]},{"none":["287"]}]})", compact(b.normal_form_json()));
    EXPECT_EQ(R"(["NZ_CP000.1","5622"])", compact(b.unknown_labels_json()));
    EXPECT_EQ(4u, b.num_names());
    EXPECT_EQ(2u, b.num_known());
    ASSERT_EQ(2u, b.labels().size());
    for (const LabelRef &label : b.labels()) {
        EXPECT_EQ(graph::traversal::LabelKind::COLUMN, label.kind);
        EXPECT_EQ(oracle.find_column(label.name), std::optional<Column>(label.column));
    }
    EXPECT_EQ("562", b.labels()[0].name);
    EXPECT_EQ("287", b.labels()[1].name);
    // the names are case-sensitive and not trimmed
    Budget again = unbounded();
    Binding other = Bound::bind(Predicate::parse(js(R"({"any": [" 562", "562 "]})"), 10), oracle, again);
    ASSERT_TRUE(other.bound);
    EXPECT_EQ("false", compact(other.bound->normal_form_json()));
}

} // namespace
