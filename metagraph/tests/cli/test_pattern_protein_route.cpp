/**
 * The /pattern route and its capabilities on a protein graph (docs/DESIGN-pattern-protein.md,
 * SPEC-pattern-search.md §29, E16): what the Protein build's capabilities state, the entry of a
 * protein pattern read in residues, the refusals (the nucleotide kinds and '*' in their slots,
 * genetic_code and strands "reverse" as the request's 400), substitutions, labels, a predicate
 * (its strands evaluated and echoed as "context") and the paths of a pattern longer than k.
 *
 * Alphabet: the Protein build only (_PROTEIN_GRAPH).
 */
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

#if _PROTEIN_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;

const size_t kK = 6;
// three proteins: p1 and p2 share MKVLAAG, p1 and p3 share RDNPQE then L or I; p2 holds the
// stored ambiguity letters B, X, Z; p3 holds U and O
const std::vector<std::string> kRecords = {
    "MKVLAAGIVGLLLAQSEWHTRDNPQELKYSTAMKV",
    "MKVLAAGPVGRWHEDNQSTAEEKLIBXZQRSTW",
    "GGSTARDNPQEIKYSTAMUOPWQ",
};
const std::vector<std::string> kLabels = { "p1", "p2", "p3" };

std::unique_ptr<AnnotatedDBG> protein_index() {
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, kRecords, kLabels, DeBruijnGraph::BASIC);
}

PatternLimits limits_of() {
    PatternLimits limits;
    limits.min_information_bits = 4;
    return limits;
}

Json::Value answer(const AnnotatedDBG &index, const std::string &body) {
    return process_pattern_request(parse_pattern_body(body), index, limits_of(), "rel");
}

std::vector<std::string> strings(const Json::Value &list) {
    std::vector<std::string> out;
    for (const Json::Value &v : list) {
        out.push_back(v.asString());
    }
    return out;
}

// the k-mers of the records holding |motif| (exact, residues) at some offset: the oracle of an
// exact pattern of at most k residues
std::set<std::pair<std::string, uint32_t>> contexts_of(const std::string &motif) {
    std::set<std::pair<std::string, uint32_t>> out;
    for (const std::string &r : kRecords) {
        for (size_t i = 0; i + kK <= r.size(); ++i) {
            const std::string kmer = r.substr(i, kK);
            for (size_t p = 0; p + motif.size() <= kK; ++p) {
                if (kmer.compare(p, motif.size(), motif) == 0)
                    out.emplace(kmer, p);
            }
        }
    }
    return out;
}

TEST(PatternProteinRoute, Capabilities) {
    auto index = protein_index();
    const pattern::GraphSupport support = route_support(index->get_graph());
    EXPECT_TRUE(support.supported) << support.reason;
    const Json::Value caps = pattern_capabilities_json(index.get(), limits_of(), false);
    EXPECT_TRUE(caps["available"].asBool());
    EXPECT_TRUE(caps["unavailable_reason"].isNull());
    EXPECT_EQ("$ABCDEFGHIJKLMNOPQRSTUVWYZX", caps["alphabet"].asString());
    EXPECT_EQ("basic", caps["graph_mode"].asString());
    EXPECT_EQ(kK, caps["k"].asUInt64());
    EXPECT_EQ(std::vector<std::string> { "protein" }, strings(caps["kinds"]));
    EXPECT_EQ("residues", caps["protein_reading"].asString());
    EXPECT_EQ((std::vector<std::string> { "[...]", "[^...]" }),
              strings(caps["protein_classes"]));
    std::string residues;
    for (const std::string &r : strings(caps["protein_residues"])) {
        residues += r;
    }
    EXPECT_EQ("ACDEFGHIKLMNPQRSTVWYUOXBZJ", residues);
    EXPECT_TRUE(caps["genetic_codes"].isArray());
    EXPECT_EQ(0u, caps["genetic_codes"].size());
    EXPECT_TRUE(caps["default_genetic_code"].isNull());
    EXPECT_EQ("SPEC-pattern-search.md section 29", caps["protein_rule"].asString());
    EXPECT_EQ((std::vector<std::string> { "both", "forward" }), strings(caps["strands"]));
    EXPECT_EQ("both", caps["default_strands"].asString());
    EXPECT_EQ(std::vector<std::string> { "basic" }, caps["scopes_by_graph_mode"].getMemberNames());
    EXPECT_EQ(std::vector<std::string> { "protein" }, caps["mismatches"]["max"].getMemberNames());
    EXPECT_EQ(3u, caps["mismatches"]["max"]["protein"].asUInt64());
    EXPECT_EQ(std::vector<std::string> { "protein" }, caps["mismatches"]["unit"].getMemberNames());
    EXPECT_EQ("residue", caps["mismatches"]["unit"]["protein"].asString());
    EXPECT_EQ("residues", caps["mismatches"]["protein_substitutes"].asString());
    EXPECT_TRUE(caps["strand_stated"].asBool());
    // the gate block of GET /traverse/capabilities states what the full block states of its keys
    const Json::Value gate = pattern_traverse_block(caps);
    EXPECT_EQ(caps["kinds"], gate["kinds"]);
    EXPECT_EQ(caps["genetic_codes"], gate["genetic_codes"]);
    EXPECT_EQ(caps["protein_residues"], gate["protein_residues"]);
}

TEST(PatternProteinRoute, EntriesAndSlotRefusals) {
    auto index = protein_index();
    const Json::Value a = answer(*index, R"({"patterns": [
        {"protein": "dnpqe", "id": "exact"},
        {"protein": "DNPQ[EL]"},
        {"protein": "[^P]NPQE"},
        {"protein": "RDNPQEJ"},
        {"protein": "XNPQE"},
        {"protein": "EDNQ"},
        {"dna": "ACGTAC"},
        {"iupac": "ACGTN"},
        {"protein": "MK*"},
        {"protein": "MK[]"}
    ]})");
    EXPECT_EQ("$ABCDEFGHIJKLMNOPQRSTUVWYZX", a["index"]["alphabet"].asString());
    const Json::Value &e = a["patterns"];
    ASSERT_EQ(10u, e.size());
    const std::vector<std::string> answered = { "DNPQE", "DNPQ[EL]", "[^P]NPQE", "RDNPQEJ",
                                                "XNPQE", "EDNQ" };
    const std::vector<size_t> lengths = { 5, 5, 5, 7, 5, 4 };
    for (size_t i = 0; i < answered.size(); ++i) {
        SCOPED_TRACE(answered[i]);
        const Json::Value &p = e[static_cast<Json::ArrayIndex>(i)];
        ASSERT_FALSE(p.isMember("error")) << p["error"];
        EXPECT_EQ("protein", p["kind"].asString());
        EXPECT_EQ(answered[i], p["pattern"].asString());
        EXPECT_EQ(lengths[i], p["length"].asUInt64());
        EXPECT_EQ(lengths[i], p["residues"].asUInt64());
        EXPECT_TRUE(p.isMember("genetic_code"));
        EXPECT_TRUE(p["genetic_code"].isNull());
        EXPECT_EQ(std::vector<std::string> { "+" }, strings(p["strands"]));
        EXPECT_FALSE(p["palindromic"].asBool());
        if (lengths[i] <= kK) {
            EXPECT_EQ(std::vector<std::string> { "+" },
                      p["counts"]["contexts"]["by_strand"].getMemberNames());
        }
    }
    EXPECT_EQ("exact", e[0]["id"].asString());
    // DNPQE: RDNPQE at 1, DNPQEL and DNPQEI at 0; [EL], [^P] and X admit the same here
    const auto dnpqe = contexts_of("DNPQE");
    EXPECT_EQ(3u, dnpqe.size());
    for (Json::ArrayIndex i : { 0u, 1u, 2u, 4u }) {
        EXPECT_EQ(dnpqe.size(), e[i]["counts"]["contexts"]["value"].asUInt64()) << i;
        std::set<std::pair<std::string, uint32_t>> released;
        for (const Json::Value &r : e[i]["results"]) {
            released.emplace(r["kmer"].asString(), r["offset"].asUInt());
            EXPECT_EQ("+", r["strand"].asString());
        }
        EXPECT_EQ(dnpqe, released) << i;
    }
    // J (I or L) at the 7th residue of a pattern longer than k: both anchors' paths
    EXPECT_EQ("long", e[3]["scope"].asString());
    EXPECT_EQ(1u, e[3]["counts"]["anchors"]["value"].asUInt64());
    // EDNQ: p2's EDNQ only
    EXPECT_EQ(contexts_of("EDNQ").size(), e[5]["counts"]["contexts"]["value"].asUInt64());
    // the nucleotide kinds, the stop and an empty class: bad_alphabet in their slots
    for (Json::ArrayIndex i : { 6u, 7u, 8u, 9u }) {
        EXPECT_EQ("bad_alphabet", e[i]["error"]["code"].asString()) << i;
    }
    EXPECT_NE(std::string::npos, e[6]["error"]["message"].asString().find("protein graph"));
    EXPECT_EQ("dna", e[6]["kind"].asString());
    EXPECT_EQ("iupac", e[7]["kind"].asString());
}

TEST(PatternProteinRoute, StoredAmbiguityIsMatchedByNothing) {
    auto index = protein_index();
    // p2 holds LIBXZQ: no pattern position matches its B, X or Z, X and a negation included,
    // and no substitution reaches them
    const Json::Value a = answer(*index, R"({"patterns": [
        {"protein": "LIBXZQ"}, {"protein": "LIXXXQ"}, {"protein": "LI[^P][^P][^P]Q"},
        {"protein": "KLI"}
    ], "max_mismatches": 0})");
    for (Json::ArrayIndex i : { 0u, 1u, 2u }) {
        EXPECT_EQ(0u, a["patterns"][i]["counts"]["contexts"]["value"].asUInt64()) << i;
        EXPECT_EQ("exact", a["patterns"][i]["counts"]["contexts"]["relation"].asString()) << i;
    }
    EXPECT_LT(0u, a["patterns"][3]["counts"]["contexts"]["value"].asUInt64());
    // within 3 substitutions of LIAAAQ, LIBXZQ would be one: its stored letters are none
    const Json::Value m = answer(*index, R"({"patterns": [{"protein": "LIAAAQ"}],
                                             "max_mismatches": 3})");
    for (const Json::Value &r : m["patterns"][0]["results"]) {
        const std::string instance = r["kmer"].asString().substr(r["offset"].asUInt(), 6);
        EXPECT_EQ(std::string::npos, instance.find_first_of("BJZX")) << instance;
    }
}

TEST(PatternProteinRoute, RequestRefusals) {
    auto index = protein_index();
    for (const char *body : {
            R"({"patterns": [{"protein": "DNPQE"}], "genetic_code": 1})",
            R"({"patterns": [{"protein": "DNPQE"}], "genetic_code": 11})",
            R"({"patterns": [{"protein": "DNPQE"}], "strands": "reverse"})" }) {
        try {
            answer(*index, body);
            ADD_FAILURE() << "answered " << body;
        } catch (const PatternRefusal &e) {
            EXPECT_EQ(400, e.status()) << body;
            EXPECT_EQ("invalid_request", e.code()) << body;
        }
    }
    // a pattern naming no kind: the message names the kind a protein graph reads (the review's S2)
    try {
        answer(*index, R"({"patterns": [{"id": "a"}]})");
        ADD_FAILURE() << "answered a pattern without a kind";
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(400, e.status());
        EXPECT_NE(std::string::npos, std::string(e.what()).find("(a protein graph reads 'protein')"))
                << e.what();
    }
    for (const char *strands : { "both", "forward" }) {
        const Json::Value a = answer(*index, std::string(R"({"patterns": [{"protein": "DNPQE"}],
                                                             "strands": ")") + strands + "\"}");
        EXPECT_EQ(3u, a["patterns"][0]["counts"]["contexts"]["value"].asUInt64()) << strands;
    }
}

TEST(PatternProteinRoute, SubstitutionsAreResidues) {
    auto index = protein_index();
    // RDNPQE[ST]: within 1 substitution, RDNPQE then L (p1) or I (p3) substitutes the class
    const Json::Value a = answer(*index, R"({"patterns": [{"protein": "DNPQE[ST]"}],
                                             "max_mismatches": 1})");
    const Json::Value &e = a["patterns"][0];
    ASSERT_FALSE(e.isMember("error")) << e["error"];
    EXPECT_EQ(1u, e["max_mismatches"].asUInt64());
    EXPECT_EQ(2u, e["counts"]["contexts"]["by_mismatches"].size());
    EXPECT_EQ(0u, e["counts"]["contexts"]["by_mismatches"]["0"]["value"].asUInt64());
    EXPECT_EQ(2u, e["counts"]["contexts"]["by_mismatches"]["1"]["value"].asUInt64());
    ASSERT_EQ(2u, e["results"].size());
    std::set<std::string> found;
    for (const Json::Value &r : e["results"]) {
        EXPECT_EQ(1u, r["mismatches"].asUInt());
        ASSERT_EQ(1u, r["substitutions"].size());
        const Json::Value &s = r["substitutions"][0];
        EXPECT_EQ(5u, s["position"].asUInt());
        EXPECT_EQ("[ST]", s["expected"].asString());
        found.insert(s["found"].asString());
        EXPECT_FALSE(s.isMember("codon"));
    }
    EXPECT_EQ((std::set<std::string> { "I", "L" }), found);
    // a pattern longer than k with substitutions: its slot (cap 0)
    const Json::Value l = answer(*index, R"({"patterns": [{"protein": "RDNPQELKY"}],
                                             "max_mismatches": 1})");
    EXPECT_EQ("mismatches_above_cap", l["patterns"][0]["error"]["code"].asString());
    EXPECT_EQ(0u, l["patterns"][0]["max_mismatches_cap"].asUInt64());
    // its message counts residues, not bases (the review's S2)
    const std::string message = l["patterns"][0]["error"]["message"].asString();
    EXPECT_NE(std::string::npos, message.find("a pattern of 9 residues, longer than k")) << message;
    EXPECT_NE(std::string::npos, message.find("windows of at most k residues")) << message;
    EXPECT_EQ(std::string::npos, message.find("bases")) << message;
}

TEST(PatternProteinRoute, LabelsAndAPredicate) {
    auto index = protein_index();
    // MKVLAAG: p1 and p2 (one k-mer each side of the branch, MKVLAA shared)
    // (a column annotation has no budget-aware decode: its reads are asked for explicitly)
    const Json::Value a = answer(*index, R"({"patterns": [{"protein": "KVLAAG"}],
                                             "output": {"labels": "all"},
                                             "allow_unbudgeted_annotation": true})");
    const Json::Value &e = a["patterns"][0];
    ASSERT_EQ(1u, e["results"].size());
    std::vector<std::string> labels;
    for (const Json::Value &l : e["results"][0]["labels"]) {
        labels.push_back(l["column"].asString());
    }
    std::sort(labels.begin(), labels.end());
    EXPECT_EQ((std::vector<std::string> { "p1", "p2" }), labels);

    // a predicate, asked with predicate_strands "either": evaluated and echoed as "context"
    // (one strand), in the answer's predicate and limits and in every selected context
    const Json::Value p = answer(*index, R"({"patterns": [{"protein": "DNPQE"}, {"protein": "AAG"}],
                                             "predicate": {"any": ["p1"]},
                                             "predicate_strands": "either",
                                             "allow_unbudgeted_annotation": true,
                                             "output": {"labels": "predicate_only"}})");
    EXPECT_EQ("context", p["predicate"]["strands"].asString());
    EXPECT_EQ("context", p["limits"]["predicate_strands"].asString());
    // DNPQE in p1: RDNPQE (shared with p3) and DNPQEL; AAG in p1: LAAGIV, AAGIVG, VLAAGI,
    // KVLAAG (KVLAAG and MKVLAA... carry p1 and p2)
    const Json::Value &d = p["patterns"][0];
    ASSERT_FALSE(d.isMember("error")) << d["error"];
    EXPECT_EQ(3u, d["counts"]["contexts"]["value"].asUInt64());
    EXPECT_EQ(2u, d["counts"]["selected"]["value"].asUInt64());
    std::set<std::string> kmers;
    for (const Json::Value &r : d["results"]) {
        kmers.insert(r["kmer"].asString());
        EXPECT_EQ(std::vector<std::string> { "p1" }, strings(r["selection_labels"]));
        for (const std::string &s : strings(r["selection_strands"])) {
            EXPECT_EQ("context", s);
        }
    }
    EXPECT_EQ((std::set<std::string> { "RDNPQE", "DNPQEL" }), kmers);
    const Json::Value &g = p["patterns"][1];
    for (const Json::Value &r : g["results"]) {
        EXPECT_NE(std::string::npos, kRecords[0].find(r["kmer"].asString())) << r["kmer"];
    }
}

TEST(PatternProteinRoute, PathsOfALongPattern) {
    auto index = protein_index();
    const Json::Value a = answer(*index, R"({"patterns": [
        {"protein": "MKVLAAG[IP]VG"}, {"protein": "RDNPQEJKYSTAM"}
    ], "long_search": "paths"})");
    const Json::Value &e = a["patterns"][0];
    ASSERT_FALSE(e.isMember("error")) << e["error"];
    EXPECT_EQ("long", e["scope"].asString());
    EXPECT_EQ(10u, e["length"].asUInt64());
    EXPECT_EQ(1u, e["counts"]["anchors"]["value"].asUInt64());
    EXPECT_EQ(2u, e["counts"]["paths"]["value"].asUInt64());
    std::set<std::string> sequences;
    for (const Json::Value &r : e["results"]) {
        sequences.insert(r["sequence"].asString());
        EXPECT_EQ(10u - kK + 1, r["nodes"].size());
    }
    EXPECT_EQ((std::set<std::string> { "MKVLAAGIVG", "MKVLAAGPVG" }), sequences);
    const Json::Value &j = a["patterns"][1];
    EXPECT_EQ(2u, j["counts"]["paths"]["value"].asUInt64());
}

// more than four allowed children at a node (the review of 2026-11-14, [P1], its example: k 6,
// CDEFGK followed by A, C, D, F and H; and a node followed by every symbol of the alphabet):
// long_search "paths" and "supported_paths" answer every walk, each with its record's label,
// and the counts are exact
TEST(PatternProteinRoute, PathsBranchingBeyondFour) {
    std::vector<std::string> records;
    std::vector<std::string> labels;
    for (char r : std::string("ACDFH")) {
        records.push_back(std::string("CDEFGK") + r);
        labels.push_back(std::string("c") + r);
    }
    // the 22 residues read as themselves and the stored ambiguity letters B, J, Z, X
    const std::string all = "ACDEFGHIKLMNPQRSTVWYUOBJZX";
    for (char r : all) {
        records.push_back(std::string("MNPQRS") + r);
        labels.push_back(std::string("m") + r);
    }
    auto index = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, records, labels, DeBruijnGraph::BASIC);
    // (pattern, the residues after its first k that its walks hold)
    const std::vector<std::pair<std::string, std::string>> cases = {
        { "CDEFGKX", "ACDFH" }, { "CDEFGK[^C]", "ADFH" }, { "MNPQRSX", all.substr(0, 22) },
        { "MNPQRS[^AW]", all.substr(1, 17) + all.substr(19, 3) },
    };
    for (const char *mode : { "paths", "supported_paths" }) {
        // (a column annotation has no budget-aware decode: its reads are asked for explicitly)
        const Json::Value a = answer(*index, std::string(R"({"patterns": [
            {"protein": "CDEFGKX"}, {"protein": "CDEFGK[^C]"}, {"protein": "MNPQRSX"},
            {"protein": "MNPQRS[^AW]"}
        ], "long_search": ")") + mode + R"(", "output": {"labels": "all"},
            "allow_unbudgeted_annotation": true})");
        ASSERT_EQ(cases.size(), a["patterns"].size());
        for (size_t i = 0; i < cases.size(); ++i) {
            const auto &[text, after] = cases[i];
            SCOPED_TRACE(std::string(mode) + " " + text);
            const Json::Value &e = a["patterns"][static_cast<Json::ArrayIndex>(i)];
            ASSERT_FALSE(e.isMember("error")) << e["error"];
            EXPECT_EQ("long", e["scope"].asString());
            const std::string prefix = text.substr(0, kK);
            std::set<std::string> expected;
            for (char r : after) {
                expected.insert(prefix + r);
            }
            const Json::Value &count = e["counts"][std::string(mode) == "paths"
                                                       ? "paths" : "supported_paths"];
            EXPECT_EQ("exact", count["relation"].asString());
            EXPECT_EQ(expected.size(), count["value"].asUInt64());
            EXPECT_EQ(expected.size(), e["counts"]["paths"]["value"].asUInt64());
            EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
            ASSERT_EQ(expected.size(), e["results"].size());
            std::set<std::string> sequences;
            for (const Json::Value &r : e["results"]) {
                const std::string sequence = r["sequence"].asString();
                sequences.insert(sequence);
                EXPECT_EQ(2u, r["nodes"].size());
                // the record holding the walk: its label, and that one only
                ASSERT_EQ(1u, r["labels"].size()) << r["labels"];
                EXPECT_EQ(std::string(1, prefix == "CDEFGK" ? 'c' : 'm') + sequence.back(),
                          r["labels"][0]["column"].asString());
            }
            EXPECT_EQ(expected, sequences);
        }
    }
}

#endif // _PROTEIN_GRAPH

} // namespace
