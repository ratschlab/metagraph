#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "tests/graph/traversal/walker_test_helpers.hpp"
#include "tests/graph/traversal/test_trie_reference.hpp"

#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/traverse.hpp"


/**
 * End-to-end tests against a real index in the format of the public `refseq33m` index:
 * basic-mode DBGSuccinct at k=31 (dummy k-mers not masked), RowDiff<BRWT> annotation with
 * k-mer coordinates, columns = NCBI taxids, and a CoordToHeader (.seqs) sidecar, so that
 * labels are accession headers.
 *
 * The fixture is the published copy of the mini index (42 blaNDM-carrying RefSeq records,
 * ~9.3 Mbp; pre-release test-data-mini-refseq-v1 of this repository), fetched by
 * scripts/traversal/fetch_mini_refseq.sh into metagraph/build/mini_refseq. The tests skip
 * when it is absent, so they never block a normal unit-test run:
 *
 *     cd metagraph/build
 *     ../scripts/traversal/fetch_mini_refseq.sh
 *     ./unit_tests --gtest_filter='MiniRefSeq*'
 */

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

const std::string kIndexDir = "mini_refseq";
const std::string kGraph = kIndexDir + "/graph_k31.dbg";
const std::string kAnnoBase = kIndexDir + "/annotation.relaxed.relabeled";
const std::string kQueryFasta = kIndexDir + "/sanity/ndm1.fa";

bool index_available() {
    return std::filesystem::exists(kGraph)
        && std::filesystem::exists(kAnnoBase + annot::RowDiffBRWTCoordAnnotator::kExtension)
        && std::filesystem::exists(kAnnoBase + annot::CoordToHeader::kExtension)
        && std::filesystem::exists(kQueryFasta);
}

// the records of the FASTA files of |dir|, in file order: (header up to the first blank, sequence)
std::vector<std::pair<std::string, std::string>> read_records(const std::string &dir) {
    std::vector<std::filesystem::path> files;
    for (const auto &entry : std::filesystem::directory_iterator(dir)) {
        files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    std::vector<std::pair<std::string, std::string>> records;
    for (const auto &path : files) {
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            if (!line.empty() && line[0] == '>') {
                records.emplace_back(line.substr(1, line.find_first_of(" \t") - 1), "");
            } else if (!records.empty()) {
                records.back().second += line;
            }
        }
    }
    return records;
}

std::string read_query() {
    std::ifstream in(kQueryFasta);
    std::string line, seq;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '>')
            seq += line;
    }
    return seq;
}

class MiniRefSeq : public ::testing::Test {
  protected:
    void SetUp() override {
        if (!index_available()) {
            GTEST_SKIP() << "mini_refseq index not found in " << std::filesystem::current_path()
                         << "/" << kIndexDir << " (fetch it with scripts/traversal/fetch_mini_refseq.sh)";
        }
        auto graph = std::make_shared<DBGSuccinct>(2);
        ASSERT_TRUE(graph->load(kGraph));

        auto annotation = std::make_unique<annot::RowDiffBRWTCoordAnnotator>();
        ASSERT_TRUE(annotation->load(kAnnoBase + annot::RowDiffBRWTCoordAnnotator::kExtension));
        using Matrix = annot::RowDiffBRWTCoordAnnotator::binary_matrix_type;
        const_cast<Matrix&>(annotation->get_matrix()).set_graph(graph.get());

        auto cth = std::make_unique<annot::CoordToHeader>();
        ASSERT_TRUE(cth->load(kAnnoBase));

        anno_graph_ = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation),
                                                     false, std::move(cth));
        ASSERT_TRUE(anno_graph_->check_compatibility());
        oracle_ = std::make_unique<LabelOracle>(*anno_graph_);
        query_ = read_query();
        ASSERT_FALSE(query_.empty());
    }

    SupportProfile discover_headers(const std::string &query) const {
        ResolveOptions options;
        options.discover = true;
        options.discover_kind = LabelKind::HEADER;
        options.discover_max_labels = 100;
        return resolve_support(*oracle_, query, options);
    }

    // the longest seed of the blaNDM query: the whole gene, carried by 19 accessions
    FrozenSeed blandm_seed() const {
        SelectionPolicy policy;
        policy.policy = SelectionPolicy::LONGEST_FIRST;
        policy.max_seeds = 1;
        auto selection = select_seeds(discover_headers(query_), query_, policy, false);
        EXPECT_FALSE(selection.seeds.empty());
        return selection.seeds.empty() ? FrozenSeed() : selection.seeds[0];
    }

    std::unique_ptr<AnnotatedDBG> anno_graph_;
    std::unique_ptr<LabelOracle> oracle_;
    std::string query_;
};


TEST_F(MiniRefSeq, IndexShape) {
    EXPECT_EQ(Regime::BASIC, oracle_->regime());
    EXPECT_EQ(31u, oracle_->get_k());
    EXPECT_EQ(9u, oracle_->num_columns());          // one column per taxid
    EXPECT_TRUE(oracle_->has_coordinates());
    ASSERT_NE(nullptr, oracle_->coord_to_header());
    uint64_t sequences = 0;
    for (Column c = 0; c < oracle_->coord_to_header()->num_columns(); ++c) {
        sequences += oracle_->coord_to_header()->num_sequences(c);
    }
    EXPECT_EQ(42u, sequences);                       // 42 indexed accessions

    EXPECT_EQ(LabelKind::COLUMN, oracle_->resolve_label("573").kind);   // Klebsiella pneumoniae
    auto accession = oracle_->resolve_label("NZ_CM008882.1");
    EXPECT_EQ(LabelKind::HEADER, accession.kind);
    EXPECT_EQ("NZ_CM008882.1", oracle_->header_name(accession.column, accession.seq_id));
    EXPECT_THROW(oracle_->resolve_label("NZ_DOES_NOT_EXIST.1"), std::invalid_argument);
}

TEST_F(MiniRefSeq, DiscoverBlaNDMCarriers) {
    auto profile = discover_headers(query_);
    EXPECT_EQ(31u, profile.k);
    EXPECT_EQ(813u - 31 + 1, profile.num_kmers);     // 783 k-mers
    ASSERT_EQ(1u, profile.graph_runs.size());
    EXPECT_EQ((KmerInterval{ 0, 783 }), profile.graph_runs[0]);

    // 25 accessions carry at least one k-mer; 19 carry the whole gene
    EXPECT_EQ(25u, profile.labels.size());
    size_t full = 0;
    for (const auto &label : profile.labels) {
        EXPECT_EQ(LabelKind::HEADER, label.label.kind);
        full += label.kmers_supported == 783;
    }
    EXPECT_EQ(19u, full);
    for (size_t i = 1; i < profile.labels.size(); ++i) {
        EXPECT_GE(profile.labels[i - 1].kmers_supported, profile.labels[i].kmers_supported);
    }
    auto it = std::find_if(profile.labels.begin(), profile.labels.end(),
                           [](const LabelProfile &l) { return l.label.name == "NZ_CM008882.1"; });
    ASSERT_NE(profile.labels.end(), it);
    EXPECT_EQ("x783", encode_runs(it->runs, profile.num_kmers));

    // the same labels named explicitly give the same profiles
    ResolveOptions named;
    for (const auto &l : profile.labels) {
        named.labels.push_back(l.label.name);
    }
    auto explicit_ = resolve_support(*oracle_, query_, named);
    ASSERT_EQ(profile.labels.size(), explicit_.labels.size());
    for (size_t i = 0; i < profile.labels.size(); ++i) {
        EXPECT_EQ(profile.labels[i].runs, explicit_.labels[i].runs) << profile.labels[i].label.name;
    }
}

// the hits of header labels are the coordinates in their sequences, mapped one by one
TEST_F(MiniRefSeq, HeaderHitsFilterRequestedRanges) {
    const SupportProfile discovered = discover_headers(query_);
    ASSERT_GE(discovered.labels.size(), 12u);
    std::vector<LabelRef> refs;
    for (size_t i = 0; i < 12; ++i) {
        refs.push_back(discovered.labels[i].label);
    }
    std::vector<node_index> keys;
    std::vector<Row> rows;
    for (node_index key : oracle_->keys_of_sequence(query_)) {
        keys.push_back(key);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    const auto tuples = oracle_->get_row_tuples(rows);
    std::vector<LabelQuery::NodeHits> expected(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        std::set<LabelId> labels;
        for (const auto &[c, coords] : tuples[i]) {
            for (uint64_t coord : coords) {
                const uint64_t seq_id = oracle_->map_coord(c, coord).first;
                for (LabelId l = 0; l < refs.size(); ++l) {
                    if (refs[l].column == c && refs[l].seq_id == seq_id)
                        labels.insert(l);
                }
            }
        }
        expected[i].assign(labels.begin(), labels.end());
    }
    EXPECT_EQ(expected, LabelQuery(*oracle_, refs).fetch(keys));
}

TEST_F(MiniRefSeq, PartialCarrierRunsMatchPresenceMask) {
    ResolveOptions options;
    options.labels = { "NZ_MOIC01000010.1" };        // x231o31x167o31x323 on refseq33m
    auto profile = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, profile.labels.size());
    EXPECT_EQ("x231o31x167o31x323", encode_runs(profile.labels[0].runs, profile.num_kmers));
    EXPECT_EQ(721u, profile.labels[0].kmers_supported);

    options.labels = { "562" };                      // Escherichia coli
    auto column_profile = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, column_profile.labels.size());
    EXPECT_EQ(LabelKind::COLUMN, column_profile.labels[0].label.kind);
    EXPECT_GE(column_profile.labels[0].kmers_supported, 721u);
}

// a basic graph: the reverse complement of the gene hits other records
TEST_F(MiniRefSeq, ForwardStrandOnly) {
    std::string rc = query_;
    ::reverse_complement(rc);
    auto forward = discover_headers(query_);
    auto reverse = discover_headers(rc);
    EXPECT_EQ(25u, forward.labels.size());
    EXPECT_EQ(17u, reverse.labels.size());
    std::set<std::string> fw, rv;
    for (const auto &l : forward.labels) {
        fw.insert(l.label.name);
    }
    for (const auto &l : reverse.labels) {
        rv.insert(l.label.name);
    }
    std::vector<std::string> both;
    std::set_intersection(fw.begin(), fw.end(), rv.begin(), rv.end(), std::back_inserter(both));
    EXPECT_EQ(3u, both.size());
}

TEST_F(MiniRefSeq, SelectFullLengthSeed) {
    auto profile = discover_headers(query_);
    SelectionPolicy policy;
    policy.max_seeds = 1;
    auto by_support = select_seeds(profile, query_, policy, false);
    ASSERT_EQ(1u, by_support.seeds.size());
    EXPECT_EQ((KmerInterval{ 0, 231 }), by_support.seeds[0].kmers);
    EXPECT_EQ(25u, by_support.seeds[0].labels.size());

    policy.policy = SelectionPolicy::LONGEST_FIRST;
    auto selection = select_seeds(profile, query_, policy, false);
    ASSERT_FALSE(selection.seeds.empty());
    const auto &seed = selection.seeds[0];
    EXPECT_EQ((KmerInterval{ 0, 783 }), seed.kmers);
    EXPECT_EQ(query_, seed.sequence);
    EXPECT_EQ(19u, seed.labels.size());
    EXPECT_EQ(0u, seed.population.dropped_count);

    policy.max_labels_per_seed = 5;
    auto capped = select_seeds(profile, query_, policy, false);
    ASSERT_FALSE(capped.seeds.empty());
    EXPECT_EQ(5u, capped.seeds[0].labels.size());
    EXPECT_EQ(14u, capped.seeds[0].population.dropped_count);
    EXPECT_FALSE(capped.seeds[0].population.dropped_digest.empty());
}

// the labels derived from the gene are its 19 carriers: the selected seed's labels, with its id
TEST_F(MiniRefSeq, DeriveCarriersFromSeed) {
    Seed seed;
    seed.sequence = query_;
    Strategy strategy;
    strategy.max_extension_bp = 3000;
    auto derived = traverse_seed(*oracle_, seed, strategy);
    EXPECT_TRUE(derived.labels_from_seed);
    EXPECT_EQ(19u, derived.label_dict.size());
    EXPECT_EQ(19u, derived.labels_supporting_total);
    std::vector<std::string> names;
    for (size_t i = 0; i < derived.label_dict.size(); ++i) {
        const LabelRef &ref = derived.label_dict[i];
        EXPECT_EQ(LabelKind::HEADER, ref.kind);
        if (i) {
            const LabelRef &prev = derived.label_dict[i - 1];
            EXPECT_LT(std::make_pair(prev.column, prev.seq_id), std::make_pair(ref.column, ref.seq_id));
        }
        names.push_back(ref.name);
    }
    const FrozenSeed frozen = blandm_seed();
    EXPECT_EQ(std::set<std::string>(frozen.labels.begin(), frozen.labels.end()),
              std::set<std::string>(names.begin(), names.end()));
    EXPECT_EQ(frozen.seed_id, derived.validated_seed_id);

    Seed named = seed;
    named.labels = names;
    auto explicitly = traverse_seed(*oracle_, named, strategy);
    EXPECT_FALSE(explicitly.labels_from_seed);
    for (size_t a : { kLeft, kRight }) {
        EXPECT_EQ(derived.arms[a].steps, explicitly.arms[a].steps);
        EXPECT_EQ(derived.arms[a].paths.size(), explicitly.arms[a].paths.size());
    }

    strategy.max_seed_labels = 5;
    auto capped = traverse_seed(*oracle_, seed, strategy);
    EXPECT_EQ(5u, capped.label_dict.size());
    EXPECT_EQ(19u, capped.labels_supporting_total);
    EXPECT_EQ(14u, capped.labels_dropped);
    uint64_t digest = 0xcbf29ce484222325ULL;
    for (size_t i = 5; i < names.size(); ++i) {
        digest = fnv1a64(names[i] + "\n", digest);
    }
    EXPECT_EQ(hex64(digest), capped.labels_dropped_digest);
}

// every label end of a walk from the gene: the seed with its flank occurs in that record
TEST_F(MiniRefSeq, ExtendBlaNDMAndValidateAgainstSourceRecords) {
    std::map<std::string, std::string> records;
    for (auto &[name, seq] : read_records(kIndexDir + "/fasta")) {
        records[name] = std::move(seq);
    }
    ASSERT_EQ(42u, records.size());

    const FrozenSeed frozen = blandm_seed();
    Seed seed;
    seed.seed_id = frozen.seed_id;
    seed.sequence = frozen.sequence;
    seed.labels = frozen.labels;
    Strategy strategy;
    strategy.max_extension_bp = 3000;
    strategy.max_label_branches = 2;
    auto result = traverse_seed(*oracle_, seed, strategy);
    EXPECT_FALSE(result.seed_id_mismatch);
    EXPECT_EQ(19u, result.label_dict.size());

    size_t checked = 0;
    for (size_t a : { kLeft, kRight }) {
        const ArmResult &arm = result.arms[a];
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << "arm " << a;
        ASSERT_FALSE(arm.paths.empty());
        for (const auto &path : arm.paths) {
            const std::string flank = spell_path(arm, path);
            const std::string full = a == kRight ? seed.sequence + flank : flank + seed.sequence;
            for (const auto &end : path.end_labels) {
                const std::string &accession = result.label_dict[end.label].name;
                ASSERT_TRUE(records.count(accession)) << accession;
                // a walk of k-mers of the record: contiguous unless the record repeats them
                if (end.branches == 0) {
                    EXPECT_NE(std::string::npos, records.at(accession).find(full))
                        << "arm " << a << " path " << path.id << " (" << full.size() << " bp)"
                        << " is attributed to " << accession << " but does not occur in it";
                    checked++;
                }
            }
        }
    }
    EXPECT_GE(checked, 19u);
}

// the walker against the string model of the 42 source records (test_trie_reference.hpp)
TEST_F(MiniRefSeq, WalkerAgainstTheSourceRecords) {
    std::vector<std::string> seqs, labels;
    for (auto &[name, seq] : read_records(kIndexDir + "/fasta")) {
        labels.push_back(name);
        seqs.push_back(std::move(seq));
    }
    ASSERT_EQ(42u, seqs.size());
    trie::StringIndex ix(31, DeBruijnGraph::BASIC, seqs, labels);

    std::vector<std::pair<std::string, uint64_t>> cases { { blandm_seed().sequence, 1000 } };
    for (size_t i : { 0u, 17u, 33u }) {
        ASSERT_GT(seqs[i].size(), 1200u) << labels[i];
        cases.emplace_back(seqs[i].substr(1000, 150), 300);
    }
    for (const auto &[seed_seq, radius] : cases) {
        const std::set<std::string> carriers = ix.carriers(seed_seq);
        ASSERT_FALSE(carriers.empty());
        Seed seed;
        seed.sequence = seed_seq;
        Strategy st;
        st.max_label_branches = Strategy::kUnlimited;
        st.max_splits_per_path = Strategy::kUnlimited;
        st.max_extension_bp = radius;
        st.max_live_paths = 1'000'000;
        st.max_paths = 1'000'000;
        const SeedResult res = traverse_seed(*oracle_, seed, st);
        std::set<std::string> names;
        for (const LabelRef &l : res.label_dict) {
            names.insert(l.name);
        }
        EXPECT_EQ(carriers, names);
        for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
            const ArmResult &a = res.arms[static_cast<size_t>(arm)];
            EXPECT_EQ(ArmResult::COMPLETE, a.status);
            trie::RefClaims got;
            for (const Segment &seg : a.segments) {
                for (const LabelEnd &e : seg.label_ends) {
                    // the outward walk of the segment's chain up to the label end
                    std::vector<size_t> chain;
                    for (size_t s = seg.id; ; s = a.segments[s].parents[0]) {
                        chain.push_back(s);
                        if (a.segments[s].parents.empty())
                            break;
                    }
                    std::string walk;
                    for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
                        std::string piece = a.segments[*it].sequence;
                        if (arm == Arm::LEFT)
                            std::reverse(piece.begin(), piece.end());
                        walk += piece;
                    }
                    got[{ walk.substr(0, e.at_bp), res.label_dict[e.label].name }] = e.reason;
                }
            }
            const trie::RefClaims expected = trie::claims(ix, carriers, seed_seq, arm, radius);
            EXPECT_FALSE(expected.empty());
            EXPECT_EQ(expected, got) << "arm " << to_string(arm) << " radius " << radius;
        }
    }
}

// the request layer on the real index: a /resolve and a /traverse of its selected seed
TEST_F(MiniRefSeq, ResolveThenTraverseRequests) {
    cli::TraverseLimits limits;
    Json::Value resolve;
    resolve["sequence"] = query_;
    resolve["discover"]["kind"] = "header";
    resolve["select"]["policy"] = "longest_first";
    resolve["select"]["max_seeds"] = 1;
    const Json::Value profile = cli::process_resolve_request(resolve, *anno_graph_, limits);
    ASSERT_EQ(1u, profile["selection"]["seeds"].size());
    const Json::Value &frozen = profile["selection"]["seeds"][0];
    EXPECT_EQ(19u, frozen["labels"].size());

    Json::Value traverse;
    traverse["seeds"][0]["seed_id"] = frozen["seed_id"];
    traverse["seeds"][0]["sequence"] = frozen["sequence"];
    traverse["seeds"][0]["labels"] = frozen["labels"];
    traverse["strategy"]["bounds"]["max_extension_bp"] = 1000;
    const Json::Value out = cli::process_traverse_request(traverse, *anno_graph_, limits);
    const Json::Value &result = out["results"][0];
    EXPECT_FALSE(result["seed"]["seed_id_mismatch"].asBool());
    EXPECT_EQ(19u, result["label_dict"].size());
    for (const char *arm : { "left", "right" }) {
        EXPECT_EQ("complete", result["arms"][arm]["status"].asString()) << arm;
        EXPECT_GE(result["arms"][arm]["paths"].size(), 1u) << arm;
    }
    // too many seeds for the server's limit
    limits.max_seeds = 0;
    EXPECT_THROW(cli::process_traverse_request(traverse, *anno_graph_, limits), cli::InvalidRequest);
}

} // namespace
