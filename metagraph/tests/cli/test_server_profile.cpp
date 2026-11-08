#include "gtest/gtest.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <json/json.h>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotation.hpp"
#include "cli/server_profile.hpp"
#include "cli/traverse.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "tests/annotation/coordinate_test_graph.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/source_path.hpp"


// The capability profile of a (graph, annotation) pair (E17, SPEC-multigraph-envelope.md
// section 8): GET /capabilities' graph_summary.pairs[] entry and a single-graph server's
// pair_summary, and support_fp over its capability fields

namespace {

using namespace mtg;
using namespace mtg::cli;
using mtg::graph::AnnotatedDBG;
using mtg::graph::DeBruijnGraph;
using mtg::graph::traversal::LabelOracle;

std::vector<std::string> random_sequences(size_t n, size_t length, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<std::string> out(n);
    for (std::string &s : out) {
        for (size_t i = 0; i < length; ++i) {
            s.push_back("ACGT"[rng() % 4]);
        }
    }
    return out;
}

std::vector<std::string> names(const std::string &prefix, size_t n) {
    std::vector<std::string> out;
    for (size_t i = 0; i < n; ++i) {
        out.push_back(prefix + std::to_string(i));
    }
    return out;
}

// the entry a multi-graph server states for |index| under the name |name|, its meta
// fingerprint computed as the server computes it at start-up
Json::Value profile(const AnnotatedDBG &index, const std::string &name = "g") {
    IndexIdentity id;
    id.meta_fp = index_meta_fingerprint(LabelOracle(index));
    return pair_summary_json(index, id, name, "/data/" + name + "/graph.dbg",
                             "/data/" + name + "/annotation.annodbg");
}

// The rule as the SPEC states it, implemented apart from the server's: FNV-1a-64 (offset
// basis 14695981039346656037, prime 1099511628211) over the version line and one line per
// field, "<name>=<value>", each value's compact JSON, as 16 lowercase hex digits
std::string json_scalar(const Json::Value &v) {
    if (v.isNull())
        return "null";
    if (v.isBool())
        return v.asBool() ? "true" : "false";
    if (v.isUInt64())
        return std::to_string(v.asUInt64());
    EXPECT_TRUE(v.isString());
    const std::string s = v.asString();
    for (unsigned char c : s) {
        // the values are ASCII without characters JSON escapes (names and an alphabet)
        EXPECT_TRUE(c >= 0x20 && c < 0x7f && c != '"' && c != '\\') << s;
    }
    return "\"" + s + "\"";
}

std::string recomputed(const Json::Value &pair, const std::string &version,
                       const std::vector<std::string> &fields) {
    std::string text = version + "\n";
    for (const std::string &name : fields) {
        const size_t dot = name.find('.');
        const Json::Value &v = dot == std::string::npos
                ? pair[name] : pair[name.substr(0, dot)][name.substr(dot + 1)];
        text += name + "=" + json_scalar(v) + "\n";
    }
    uint64_t h = 14695981039346656037ULL;
    for (unsigned char c : text) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return buf;
}

// The fenced block of SPEC-multigraph-envelope.md section 8 that states the text support_fp
// hashes: its version line and its field names (the lines "<name>=<value>")
std::vector<std::string> spec_support_block() {
    std::ifstream in(source_path("docs/SPEC-multigraph-envelope.md"));
    EXPECT_TRUE(in.good());
    std::vector<std::string> block;
    std::string line;
    bool inside = false;
    while (std::getline(in, line)) {
        if (!inside && line.rfind("mgt-support ", 0) == 0) {
            inside = true;
        } else if (inside && line.rfind("```", 0) == 0) {
            break;
        }
        if (inside)
            block.push_back(line);
    }
    return block;
}

TEST(MultiGraphProfile, TheSpecStatesTheFieldsAndTheirOrder) {
    const std::vector<std::string> block = spec_support_block();
    ASSERT_FALSE(block.empty()) << "no `mgt-support` block in section 8";
    EXPECT_EQ(std::string(kSupportVersion), block.front());
    std::vector<std::string> fields;
    for (size_t i = 1; i < block.size(); ++i) {
        const std::string suffix = "=<value>";
        ASSERT_GT(block[i].size(), suffix.size()) << block[i];
        ASSERT_EQ(suffix, block[i].substr(block[i].size() - suffix.size())) << block[i];
        fields.push_back(block[i].substr(0, block[i].size() - suffix.size()));
    }
    EXPECT_EQ(support_fields(), fields);
    // in byte order of the names, each once
    std::vector<std::string> sorted = fields;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    EXPECT_EQ(sorted, fields);
}

// Every annotation type the CLI knows: annotation_form is its --anno-type name, and
// decode_charged the oracle's predicate on the matrix the loader constructs for it — true on
// the row-diff forms over an index matrix with a single-row read and on plain brwt and
// column (E14), false on the others
TEST(MultiGraphProfile, AnnotationFormAndDecodeChargedOnEveryForm) {
    const std::set<std::string> charged {
        "column", "brwt", "row_diff", "row_diff_coord", "row_diff_brwt", "row_diff_brwt_coord",
        "row_diff_disk", "row_diff_sparse", "row_diff_flat", "row_diff_int_brwt",
    };
    const std::set<std::string> not_charged {
        "row", "bin_rel_wt", "flat", "rbfish", "rb_brwt", "row_sparse", "int_brwt",
        "row_diff_int_disk", "column_coord", "brwt_coord", "row_diff_disk_coord",
    };
    std::set<std::string> seen;
    for (int t = Config::ColumnCompressed; t <= Config::RowDiffDiskCoord; ++t) {
        const auto type = static_cast<Config::AnnotationType>(t);
        auto annotation = initialize_annotation(type);
        ASSERT_TRUE(annotation) << t;
        const std::string form = Config::annotype_to_string(type);
        EXPECT_EQ(form, annotation_form(*annotation));
        ASSERT_TRUE(charged.count(form) || not_charged.count(form)) << form;
        EXPECT_EQ(charged.count(form) > 0,
                  LabelOracle::decode_charged(annotation->get_matrix())) << form;
        seen.insert(form);
    }
    EXPECT_EQ(charged.size() + not_charged.size(), seen.size());
}

// On indexes built with data, the entry states what the oracle over the index says
TEST(MultiGraphProfile, DecodeChargedIsTheOraclesOnBuiltIndexes) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    const auto seqs = random_sequences(3, 80, 17);
    const auto labels = names("L", 3);
    std::vector<std::pair<std::string, std::unique_ptr<AnnotatedDBG>>> indexes;
    indexes.emplace_back("column", test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(
            21, seqs, labels));
    indexes.emplace_back("column_coord", test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(
            21, seqs, labels, DeBruijnGraph::BASIC, true));
    indexes.emplace_back("row_diff", test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, seqs, labels));
    indexes.emplace_back("row_diff_coord",
            test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
                    21, seqs, labels, DeBruijnGraph::BASIC, true));
    indexes.emplace_back("row_diff_disk", test::build_anno_graph<DBGSuccinct, RowDiffDiskAnnotator>(
            21, seqs, labels));
    indexes.emplace_back("flat", test::build_anno_graph<DBGSuccinct, RowFlatAnnotator>(
            21, seqs, labels));
    for (const auto &[form, index] : indexes) {
        const Json::Value p = profile(*index);
        const LabelOracle oracle(*index);
        EXPECT_EQ(form, p["annotation_form"].asString());
        EXPECT_EQ(oracle.decode_charged(), p["decode_charged"].asBool()) << form;
        EXPECT_EQ(LabelOracle::decode_charged(index->get_annotator().get_matrix()),
                  p["decode_charged"].asBool()) << form;
        EXPECT_EQ(form != "flat" && form != "column_coord", p["decode_charged"].asBool()) << form;
        // what /traverse states
        EXPECT_EQ(oracle.graph().alphabet(), p["alphabet"].asString());
        EXPECT_EQ(oracle.has_coordinates(), p["traversal"]["has_coordinates"].asBool());
        EXPECT_TRUE(p["traversal"]["counts"].isNull()) << form;
        EXPECT_EQ("succinct", p["graph_representation"].asString());
        EXPECT_EQ(index_meta_fingerprint(oracle), p["index_meta_fp"].asString());
    }
}

// support_fp is recomputed from the entry's own fields by the rule the SPEC states: its
// version line, then one line per field in the SPEC's order
TEST(MultiGraphProfile, SupportFpIsRecomputableFromTheStatedFields) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    const std::vector<std::string> block = spec_support_block();
    ASSERT_FALSE(block.empty());
    std::vector<std::string> fields;
    for (size_t i = 1; i < block.size(); ++i) {
        fields.push_back(block[i].substr(0, block[i].find('=')));
    }
    const auto seqs = random_sequences(3, 80, 23);
    const auto labels = names("L", 3);
    for (bool coordinates : { false, true }) {
        auto index = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
                21, seqs, labels, DeBruijnGraph::BASIC, coordinates);
        const Json::Value p = profile(*index);
        EXPECT_EQ(recomputed(p, block.front(), fields), p["support_fp"].asString());
        EXPECT_EQ(support_fingerprint(p), p["support_fp"].asString());
        EXPECT_EQ(16u, p["support_fp"].asString().size());
        // the text names each field once, with its value
        const std::string text = support_text(p);
        EXPECT_EQ(0u, text.find(std::string(kSupportVersion) + "\n"));
        EXPECT_NE(std::string::npos, text.find("\nannotation_form=\"" + std::string(
                coordinates ? "row_diff_coord" : "row_diff") + "\"\n")) << text;
        EXPECT_NE(std::string::npos, text.find("\ntraversal.counts=null\n")) << text;
        EXPECT_NE(std::string::npos, text.find("\nk=21\n")) << text;
    }
}

// Two pairs of one build (the same graph representation, k, mode and annotation form, other
// sequences and columns: two Logan chunks) share support_fp, while their index_meta_fp and
// num_labels differ; a pair of another annotation form or without coordinates does not
TEST(MultiGraphProfile, TwoPairsOfOneBuildShareTheirSupportFp) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    auto chunk1 = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, random_sequences(3, 90, 1), names("SRR1", 3), DeBruijnGraph::BASIC, true);
    auto chunk2 = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, random_sequences(5, 70, 2), names("SRR2", 5), DeBruijnGraph::BASIC, true);
    const Json::Value p1 = profile(*chunk1, "chunk1");
    const Json::Value p2 = profile(*chunk2, "chunk2");
    EXPECT_EQ(p1["support_fp"], p2["support_fp"]);
    EXPECT_NE(p1["index_meta_fp"], p2["index_meta_fp"]);
    EXPECT_NE(p1["traversal"]["num_labels"], p2["traversal"]["num_labels"]);
    EXPECT_NE(p1["graph_path"], p2["graph_path"]);

    auto plain = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, random_sequences(3, 90, 1), names("SRR1", 3), DeBruijnGraph::BASIC, false);
    auto column = test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(
            21, random_sequences(3, 90, 1), names("SRR1", 3), DeBruijnGraph::BASIC, true);
    auto other_k = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            25, random_sequences(3, 90, 1), names("SRR1", 3), DeBruijnGraph::BASIC, true);
    std::set<std::string> fps { p1["support_fp"].asString() };
    for (const AnnotatedDBG *index : { plain.get(), column.get(), other_k.get() }) {
        EXPECT_TRUE(fps.insert(profile(*index)["support_fp"].asString()).second);
    }
}

// Each capability field changes support_fp alone; the identity, the size, a member outside
// the list and support_fp itself do not
TEST(MultiGraphProfile, ExactlyTheCapabilityFieldsEnterSupportFp) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    auto index = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, random_sequences(3, 80, 5), names("L", 3), DeBruijnGraph::BASIC, true);
    const Json::Value p = profile(*index);
    const std::string fp = support_fingerprint(p);
    ASSERT_EQ(p["support_fp"].asString(), fp);
    auto member = [](Json::Value &pair, const std::string &name) -> Json::Value& {
        const size_t dot = name.find('.');
        return dot == std::string::npos ? pair[name]
                                        : pair[name.substr(0, dot)][name.substr(dot + 1)];
    };
    for (const std::string &name : support_fields()) {
        Json::Value changed = p;
        Json::Value &v = member(changed, name);
        if (v.isBool()) {
            v = !v.asBool();
        } else if (v.isUInt64()) {
            v = v.asUInt64() + 1;
        } else if (v.isNull()) {
            v = "x";
        } else {
            v = v.asString() + "x";
        }
        EXPECT_NE(fp, support_fingerprint(changed)) << name;
        // and null differs from absent nowhere: an absent member is stated as null
        Json::Value nulled = p;
        member(nulled, name) = Json::Value();
        Json::Value removed = p;
        if (name.find('.') == std::string::npos) {
            removed.removeMember(name);
        } else {
            removed["traversal"].removeMember(name.substr(name.find('.') + 1));
        }
        EXPECT_EQ(support_fingerprint(nulled), support_fingerprint(removed)) << name;
    }
    for (const char *name : { "graph", "graph_path", "annotation_path", "index_ns", "index_fp",
                              "index_meta_fp", "support_fp", "traversal.num_labels",
                              "a_member_added_later" }) {
        Json::Value changed = p;
        member(changed, name) = "other";
        EXPECT_EQ(fp, support_fingerprint(changed)) << name;
    }
    // the version line is part of the text: another version, another fingerprint
    EXPECT_EQ(0u, support_text(p).find("mgt-support 2\n"));
}

// seqs_order (H4, SPEC-multigraph-envelope.md section 9.1): what the load verified of the
// order of the record mapping's columns, a capability field (hashed): "verified" once its
// stated labels matched the annotation's, "count_only" for a mapping that states none, null
// without a mapping; pairs that differ in it have different profiles
TEST(MultiGraphProfile, SeqsOrderStatesWhatTheLoadVerified) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    const std::vector<test::CoordinateRecord> records {
        { "ACGTTGCATGACCTAGGATCCGATTACGGA", "A", "a1" },
        { "GGCATTACCAGTTCAGGA", "B", "b1" },
    };
    auto count_only = test::coordinate_test_graph(5, records, true, true);
    auto verified = test::coordinate_test_graph(5, records, true, true, true);
    auto without = test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, random_sequences(3, 80, 5), names("L", 3), DeBruijnGraph::BASIC, true);
    const Json::Value c = profile(*count_only);
    const Json::Value v = profile(*verified);
    const Json::Value w = profile(*without);
    EXPECT_EQ("count_only", c["seqs_order"].asString());
    EXPECT_EQ("verified", v["seqs_order"].asString());
    ASSERT_TRUE(w.isMember("seqs_order"));
    EXPECT_TRUE(w["seqs_order"].isNull());
    EXPECT_FALSE(w["traversal"]["has_coord_to_header"].asBool());
    EXPECT_TRUE(c["traversal"]["has_coord_to_header"].asBool());
    // the one difference of the two pairs, and support_fp tells them apart
    Json::Value c_but_verified = c;
    c_but_verified["seqs_order"] = "verified";
    c_but_verified.removeMember("support_fp");
    Json::Value v_without_fp = v;
    v_without_fp.removeMember("support_fp");
    EXPECT_EQ(v_without_fp, c_but_verified);
    EXPECT_NE(c["support_fp"], v["support_fp"]);
    // its line after mask's, before traversal's (the names' byte order)
    const std::string text = support_text(v);
    const size_t line = text.find("\nseqs_order=\"verified\"\ntraversal.counts=");
    ASSERT_NE(std::string::npos, line) << text;
    EXPECT_LT(text.find("\nmask="), line) << text;
}

// Every member of the entry, and of its `traversal`, is either a field support_fp hashes or
// one it leaves out on purpose (SPEC-multigraph-envelope.md section 8.1: the identity, the
// size, support_fp itself). A member added later must go into one of the two lists: left out
// of support_fields() unnoticed, it would merge two pairs of different capabilities under one
// support_fp (a_member_added_later above shows that a member outside the list leaves the hash
// alone). Both lists name members the entry has, so neither keeps a stale name.
TEST(MultiGraphProfile, EveryMemberIsHashedOrExcluded) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    const std::set<std::string> excluded {
        "graph", "graph_path", "annotation_path", "index_ns", "index_fp", "index_meta_fp",
        "traversal.num_labels", "support_fp",
    };
    const std::set<std::string> hashed(support_fields().begin(), support_fields().end());
    for (const std::string &name : excluded) {
        EXPECT_FALSE(hashed.count(name)) << name << ": both hashed and excluded";
    }
    const auto seqs = random_sequences(3, 80, 31);
    const auto labels = names("L", 3);
    std::vector<std::unique_ptr<AnnotatedDBG>> indexes;
    indexes.push_back(test::build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(
            21, seqs, labels, DeBruijnGraph::BASIC, true));
    indexes.push_back(test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(21, seqs, labels));
    indexes.push_back(test::build_anno_graph<DBGSuccinct, RowFlatAnnotator>(21, seqs, labels));
    indexes.push_back(test::build_anno_graph<DBGHashOrdered, ColumnCompressed<>>(
            21, seqs, labels));
    for (const auto &index : indexes) {
        // a multi-graph server's entry, and a single-graph server's pair_summary (tags null)
        IndexIdentity id;
        id.meta_fp = index_meta_fingerprint(LabelOracle(*index));
        const Json::Value single = pair_summary_json(*index, id, Json::Value(), Json::Value(),
                                                     Json::Value());
        for (const Json::Value &p : { profile(*index), single }) {
            ASSERT_TRUE(p.isObject());
            ASSERT_TRUE(p["traversal"].isObject());
            std::set<std::string> members;
            for (const std::string &name : p.getMemberNames()) {
                if (name == "traversal") {
                    for (const std::string &member : p[name].getMemberNames()) {
                        members.insert("traversal." + member);
                    }
                    continue;
                }
                // another object or array would need its members listed as traversal's are
                EXPECT_FALSE(p[name].isObject() || p[name].isArray()) << name;
                members.insert(name);
            }
            for (const std::string &name : members) {
                EXPECT_TRUE(hashed.count(name) || excluded.count(name))
                        << name << ": a member of the entry neither in support_fields() nor"
                                   " excluded from support_fp";
            }
            for (const std::string &name : hashed) {
                EXPECT_TRUE(members.count(name)) << name << ": hashed, not a member";
            }
            for (const std::string &name : excluded) {
                EXPECT_TRUE(members.count(name)) << name << ": excluded, not a member";
            }
        }
    }
}

// graph_representation: the loaded graph's --graph name
TEST(MultiGraphProfile, GraphRepresentationIsTheGraphsType) {
    using namespace mtg::annot;
    using namespace mtg::graph;
    const auto seqs = random_sequences(2, 60, 9);
    const auto labels = names("L", 2);
    EXPECT_EQ("succinct", profile(*test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
    EXPECT_EQ("hash", profile(*test::build_anno_graph<DBGHashOrdered, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
    EXPECT_EQ("hashfast", profile(*test::build_anno_graph<DBGHashFast, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
    EXPECT_EQ("hashstr", profile(*test::build_anno_graph<DBGHashString, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
#if ! _PROTEIN_GRAPH
    EXPECT_EQ("bitmap", profile(*test::build_anno_graph<DBGBitmap, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
    EXPECT_EQ("sshash", profile(*test::build_anno_graph<DBGSSHash, ColumnCompressed<>>(
            11, seqs, labels))["graph_representation"].asString());
    // a PRIMARY graph wrapped in CanonicalDBG (DNA: a reverse complement): the wrapped graph's
    const Json::Value primary = profile(*test::build_anno_graph<DBGSuccinct, ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::PRIMARY));
    EXPECT_EQ("succinct", primary["graph_representation"].asString());
    EXPECT_EQ("primary", primary["traversal"]["regime"].asString());
#endif
    // each name is one the CLI's --graph takes
    for (const char *name : { "succinct", "hash", "hashfast", "hashstr", "bitmap", "sshash" }) {
        EXPECT_NE(Config::GraphType::INVALID, Config::string_to_graphtype(name)) << name;
    }
}

} // namespace
