#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "tests/graph/traversal/test_trie_oracle.hpp"
#include "tests/graph/traversal/test_trie_checks.hpp"

#include "graph/traversal/resolve.hpp"
#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/unix_tools.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
namespace trie = mtg::test::trie;

/**
 * Opt-in tests against a real index in the format of the public `refseq33m`
 * index: basic-mode DBGSuccinct at k=31 (forward strand only, dummy k-mers not
 * masked), RowDiff<BRWT> annotation with k-mer coordinates, columns = NCBI
 * taxids, and a CoordToHeader (.seqs) sidecar so labels are accession headers.
 *
 * The fixture is built by scripts/traversal/build_mini_refseq.sh (42 real
 * blaNDM-carrying RefSeq records, ~9.3 Mbp, about 25 s). These tests skip when
 * it is absent, so they never block a normal unit-test run. Run them with:
 *
 *     cd metagraph/build
 *     ../scripts/traversal/build_mini_refseq.sh ./mini_refseq ./metagraph \
 *         ../scripts/traversal/mini_refseq_accessions.tsv
 *     ./unit_tests --gtest_filter='MiniRefSeq*'
 */
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

std::string read_query() {
    std::ifstream in(kQueryFasta);
    std::string line, seq;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '>')
            seq += line;
    }
    return seq;
}

// A traversal result reduced to what it claims: the seed's verdict (ids, kept and dropped
// labels, per-label summary, which annotation accessor was used) and, per arm, its shape
// and every leaf's flank with the accessions reaching it. Everything a derived run must
// reproduce from an explicit one is in here; the per-request annotation COUNTERS are not,
// because a derived set is obtained by reading whole rows and an explicit one by reading
// its own columns, so only some of them can be equal (checked separately, see
// DeriveCarriersFromSeed).
std::string reduce(const SeedResult &r) {
    std::ostringstream os;
    os << r.validated_seed_id << ' ' << r.num_seed_labels << ' ' << r.length_bp
       << ' ' << r.num_kmers << ' ' << r.access_path;
    for (const auto &d : r.dropped_labels) {
        os << "\nDROP " << d.name << ' ' << d.reason;
        for (const auto &[from, to] : d.runs) os << ' ' << from << '-' << to;
    }
    for (size_t i = 0; i < r.label_summary.size(); ++i) {
        os << "\nSUM " << r.label_dict[i].name << ' ' << to_string(r.label_dict[i].kind);
        for (size_t a : { static_cast<size_t>(Arm::LEFT), static_cast<size_t>(Arm::RIGHT) }) {
            const LabelArmSummary &s = r.label_summary[i][a];
            os << ' ' << s.direct_bp << ':' << s.reach_bp << ':' << s.reentries
               << ':' << s.runs.size();
        }
    }
    for (size_t a : { static_cast<size_t>(Arm::LEFT), static_cast<size_t>(Arm::RIGHT) }) {
        const ArmResult &arm = r.arms[a];
        os << "\nARM " << a << ' ' << arm.requested << ' ' << arm.status << ' ' << arm.steps
           << ' ' << arm.output_bp << ' ' << arm.segments.size() << ' ' << arm.splits.size()
           << ' ' << arm.runs.size() << ' ' << arm.growth.size() << ' '
           << arm.branch_events_total;
        for (const auto &run : arm.runs) {
            os << "\n R " << r.label_dict[run.label].name << ' ' << run.from_bp << ' '
               << run.to_bp << ' ' << run.entered_by_switch << ' ' << run.route_bp << ' '
               << run.ended << ' ' << (run.ended ? to_string(run.end_reason) : "-");
        }
        std::vector<std::string> leaves;
        for (const auto &path : arm.paths) {
            std::string line = std::to_string(path.length_bp) + " " + spell_path(arm, path);
            std::vector<std::string> ends;
            for (const auto &end : path.end_labels) {
                ends.push_back(r.label_dict[end.label].name + ":" + std::to_string(end.loss)
                               + ":" + std::to_string(end.route_bp));
            }
            std::sort(ends.begin(), ends.end());
            for (const auto &e : ends) line += " " + e;
            leaves.push_back(line);
        }
        std::sort(leaves.begin(), leaves.end());
        for (const auto &l : leaves) os << "\n " << l;
    }
    return os.str();
}

class MiniRefSeq : public ::testing::Test {
  protected:
    void SetUp() override {
        if (!index_available()) {
            GTEST_SKIP() << "mini_refseq index not found in " << std::filesystem::current_path()
                         << "/" << kIndexDir << " (build it with scripts/traversal/build_mini_refseq.sh)";
        }
        auto graph = std::make_shared<DBGSuccinct>(2);
        ASSERT_TRUE(graph->load(kGraph));

        auto annotation = std::make_unique<annot::RowDiffBRWTCoordAnnotator>();
        ASSERT_TRUE(annotation->load(kAnnoBase + annot::RowDiffBRWTCoordAnnotator::kExtension));
        // RowDiff needs the base graph to reconstruct rows (its anchors and fork
        // successors are serialized inside the .annodbg for this representation)
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

    std::unique_ptr<AnnotatedDBG> anno_graph_;
    std::unique_ptr<LabelOracle> oracle_;
    std::string query_;
};


// the oracle recognizes the refseq33m-style index and its label semantics
TEST_F(MiniRefSeq, IndexShape) {
    EXPECT_EQ(Regime::BASIC, oracle_->regime());
    EXPECT_EQ(31u, oracle_->get_k());
    EXPECT_EQ(9u, oracle_->num_columns());          // one column per taxid
    EXPECT_TRUE(oracle_->has_coordinates());
    ASSERT_NE(nullptr, oracle_->coord_to_header());
    EXPECT_EQ(9u, oracle_->coord_to_header()->num_columns());
    uint64_t sequences = 0;
    for (Column c = 0; c < oracle_->coord_to_header()->num_columns(); ++c) {
        sequences += oracle_->coord_to_header()->num_sequences(c);
    }
    EXPECT_EQ(42u, sequences);                       // 42 indexed accessions

    // a taxid is a column label; an accession is a header label
    auto taxid = oracle_->resolve_label("573");      // Klebsiella pneumoniae
    EXPECT_EQ(LabelKind::COLUMN, taxid.kind);
    auto accession = oracle_->resolve_label("NZ_CM008882.1");
    EXPECT_EQ(LabelKind::HEADER, accession.kind);
    EXPECT_EQ("NZ_CM008882.1", oracle_->header_name(accession.column, accession.seq_id));
    EXPECT_THROW(oracle_->resolve_label("NZ_DOES_NOT_EXIST.1"), std::invalid_argument);

    // coordinate-aware annotation: the header's k-mer count is known
    EXPECT_GT(oracle_->num_kmers_in_sequence(accession.column, accession.seq_id), 0u);
}

// discovery over a real gene query reproduces the per-accession support that
// `metagraph query` reports on the same index
TEST_F(MiniRefSeq, DiscoverBlaNDMCarriers) {
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;

    auto profile = resolve_support(*oracle_, query_, options);
    EXPECT_EQ(31u, profile.k);
    EXPECT_EQ(813u - 31 + 1, profile.num_kmers);     // 783 k-mers
    // the query is present end to end in this index
    ASSERT_EQ(1u, profile.graph_runs.size());
    EXPECT_EQ((KmerInterval{ 0, 783 }), profile.graph_runs[0]);

    // 25 accessions carry at least one k-mer; 19 carry the whole gene
    EXPECT_EQ(25u, profile.labels.size());
    size_t full = 0;
    for (const auto &label : profile.labels) {
        EXPECT_EQ(LabelKind::HEADER, label.label.kind);
        EXPECT_GT(label.kmers_supported, 0u);
        if (label.kmers_supported == 783)
            full++;
    }
    EXPECT_EQ(19u, full);
    // discovery is ordered by supported k-mers
    for (size_t i = 1; i < profile.labels.size(); ++i) {
        EXPECT_GE(profile.labels[i - 1].kmers_supported, profile.labels[i].kmers_supported);
    }

    // one of the known full-length carriers, with a single support run
    auto it = std::find_if(profile.labels.begin(), profile.labels.end(),
                           [](const LabelProfile &l) { return l.label.name == "NZ_CM008882.1"; });
    ASSERT_NE(profile.labels.end(), it);
    EXPECT_EQ(783u, it->kmers_supported);
    EXPECT_EQ("x783", encode_runs(it->runs, profile.num_kmers));
}

// explicit labels: a partial carrier's support runs match the presence mask that
// `metagraph query --query-mode signature` reports (SNPs break runs of k-mers)
TEST_F(MiniRefSeq, PartialCarrierRunsMatchPresenceMask) {
    ResolveOptions options;
    options.labels = { "NZ_MOIC01000010.1" };        // x231o31x167o31x323 on refseq33m
    auto profile = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, profile.labels.size());
    EXPECT_EQ("x231o31x167o31x323", encode_runs(profile.labels[0].runs, profile.num_kmers));
    EXPECT_EQ(721u, profile.labels[0].kmers_supported);
    // two SNPs, each disrupting the 31 k-mers that overlap it, split the support
    // into three runs: this is the gap pattern an agent would investigate
    EXPECT_EQ(3u, profile.labels[0].runs.size());

    // the same accession under its taxid column still covers the whole query,
    // because other records of that taxid carry the missing k-mers
    options.labels = { "562" };                      // Escherichia coli
    auto column_profile = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, column_profile.labels.size());
    EXPECT_EQ(LabelKind::COLUMN, column_profile.labels[0].label.kind);
    EXPECT_GE(column_profile.labels[0].kmers_supported, 721u);
}

// the forward-strand-only index does not report the reverse complement of a
// query under the same accession (unless the record carries both orientations)
TEST_F(MiniRefSeq, ForwardStrandOnly) {
    std::string rc = query_;
    std::reverse(rc.begin(), rc.end());
    for (char &c : rc) {
        c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';
    }
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;

    auto forward = resolve_support(*oracle_, query_, options);
    auto reverse = resolve_support(*oracle_, rc, options);
    // 25 forward carriers, 17 reverse carriers, and the two sets are nearly disjoint
    EXPECT_EQ(25u, forward.labels.size());
    EXPECT_EQ(17u, reverse.labels.size());
    std::set<std::string> fw, rv;
    for (const auto &l : forward.labels) fw.insert(l.label.name);
    for (const auto &l : reverse.labels) rv.insert(l.label.name);
    std::vector<std::string> both;
    std::set_intersection(fw.begin(), fw.end(), rv.begin(), rv.end(), std::back_inserter(both));
    // only records carrying two copies in opposite orientations are in both sets
    EXPECT_EQ(3u, both.size()) << "forward and RC hit sets should be nearly disjoint";
    // NZ_CP034369.1 carries two blaNDM copies in opposite orientations
    EXPECT_EQ(1u, fw.count("NZ_CP034369.1"));
    EXPECT_EQ(1u, rv.count("NZ_CP034369.1"));
}

// selection freezes a seed that is carried in full by many accessions: this is
// the input the traversal takes, and it must be reproducible
TEST_F(MiniRefSeq, SelectFullLengthSeed) {
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;
    auto profile = resolve_support(*oracle_, query_, options);

    // max_support maximizes the NUMBER of labels carrying the whole block, so it
    // prefers the breadth-optimal prefix (every carrier agrees on the first 231
    // k-mers; the SNPs that break 6 of them come later) over the full gene
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.max_seeds = 1;
    policy.release_id = "mini_refseq";
    auto by_support = select_seeds(profile, query_, policy, /* canonical */ false);
    ASSERT_EQ(1u, by_support.seeds.size());
    EXPECT_EQ((KmerInterval{ 0, 231 }), by_support.seeds[0].kmers);
    EXPECT_EQ(25u, by_support.seeds[0].labels.size());

    // longest_first is the policy for "the longest possible seed": the whole gene,
    // carried in full by the 19 full-length accessions
    policy.policy = SelectionPolicy::LONGEST_FIRST;
    auto selection = select_seeds(profile, query_, policy, false);
    ASSERT_FALSE(selection.seeds.empty());
    const auto &seed = selection.seeds[0];
    EXPECT_EQ((KmerInterval{ 0, 783 }), seed.kmers);
    EXPECT_EQ(query_, seed.sequence);
    EXPECT_EQ(19u, seed.labels.size());
    EXPECT_EQ(19u, seed.population.supporting_total);
    EXPECT_EQ(0u, seed.population.dropped_count);
    // reproducible identity
    auto again = select_seeds(profile, query_, policy, false);
    EXPECT_EQ(seed.seed_id, again.seeds[0].seed_id);

    // capping the label set reports what was dropped so the complement can be run
    policy.max_labels_per_seed = 5;
    auto capped = select_seeds(profile, query_, policy, false);
    ASSERT_FALSE(capped.seeds.empty());
    EXPECT_EQ(5u, capped.seeds[0].labels.size());
    EXPECT_EQ(14u, capped.seeds[0].population.dropped_count);
    EXPECT_FALSE(capped.seeds[0].population.dropped_digest.empty());
}

// The design note's `permit: per_hit` default on the real index: seeding the whole
// blaNDM gene with NO labels must derive exactly the full-length carriers — the same 19
// accessions the explicit longest_first selection freezes — and traverse identically.
TEST_F(MiniRefSeq, DeriveCarriersFromSeed) {
    Seed seed;
    seed.sequence = query_;          // no labels: derive them from the seed
    Strategy strategy;
    strategy.max_extension_bp = 3000;
    strategy.profile_bin_bp = 500;
    auto derived = traverse_seed(*oracle_, seed, strategy, LabelChangeCost::forbid());

    EXPECT_TRUE(derived.labels_from_seed);
    EXPECT_EQ(19u, derived.num_seed_labels);
    EXPECT_EQ(19u, derived.labels_supporting_total);
    EXPECT_EQ(0u, derived.labels_dropped);
    EXPECT_TRUE(derived.labels_dropped_digest.empty());
    EXPECT_TRUE(derived.dropped_labels.empty());
    // the index has a CoordToHeader, so the derived labels are accessions, ascending
    std::vector<std::string> names;
    for (size_t i = 0; i < derived.label_dict.size(); ++i) {
        const LabelRef &ref = derived.label_dict[i];
        EXPECT_EQ(LabelKind::HEADER, ref.kind);
        EXPECT_EQ(ref.name, oracle_->header_name(ref.column, ref.seq_id));
        if (i) {
            const LabelRef &prev = derived.label_dict[i - 1];
            EXPECT_LT(std::make_pair(prev.column, prev.seq_id),
                      std::make_pair(ref.column, ref.seq_id));
        }
        names.push_back(ref.name);
    }

    // exactly the set `resolve` + `select` freeze as the longest seed's carriers
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;
    auto profile = resolve_support(*oracle_, query_, options);
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::LONGEST_FIRST;
    policy.max_seeds = 1;
    auto selection = select_seeds(profile, query_, policy, /* canonical */ false);
    ASSERT_FALSE(selection.seeds.empty());
    const FrozenSeed &frozen = selection.seeds[0];
    ASSERT_EQ(19u, frozen.labels.size());
    EXPECT_EQ(std::set<std::string>(frozen.labels.begin(), frozen.labels.end()),
              std::set<std::string>(names.begin(), names.end()));
    // seed ids are computed over the sorted labels, so an equal set means an equal id
    EXPECT_EQ(frozen.seed_id, derived.validated_seed_id);

    // ... and the traversal is the one naming those accessions explicitly: the derived
    // list is resubmittable VERBATIM and reproduces the whole result but the three
    // derivation-only fields
    Seed named;
    named.sequence = query_;
    named.labels = names;
    // a fresh oracle: SeedResult::annotation_counters are the oracle's running totals,
    // and |oracle_| has meanwhile also served the resolve/select above
    LabelOracle fresh(*anno_graph_);
    auto explicitly = traverse_seed(fresh, named, strategy, LabelChangeCost::forbid());
    EXPECT_FALSE(explicitly.labels_from_seed);
    EXPECT_EQ(reduce(explicitly), reduce(derived));
    EXPECT_EQ(19u, explicitly.labels_supporting_total);   // its own total, never 0
    EXPECT_EQ(0u, explicitly.labels_dropped);
    EXPECT_TRUE(explicitly.labels_dropped_digest.empty());
    // Deriving adds no annotation reads: it intersects the very rows the per-label seed
    // validation reads anyway, one per seed k-mer. What it does differ in is how much of
    // each row it looks at — whole rows versus the named columns — so the coordinate
    // mappings and the fetch time may differ while the row count may not.
    EXPECT_EQ(explicitly.annotation_counters.rows_requested,
              derived.annotation_counters.rows_requested);
    EXPECT_EQ(explicitly.annotation_counters.keys_mapped,
              derived.annotation_counters.keys_mapped);

    // the cap bounds the traversal state, not the discovery: the labels above it are
    // found either way and reported
    strategy.max_seed_labels = 5;
    auto capped = traverse_seed(*oracle_, seed, strategy, LabelChangeCost::forbid());
    EXPECT_EQ(5u, capped.num_seed_labels);
    EXPECT_EQ(19u, capped.labels_supporting_total);
    EXPECT_EQ(14u, capped.labels_dropped);
    EXPECT_EQ(16u, capped.labels_dropped_digest.size());
    uint64_t digest = kFnvOffsetBasis;
    for (size_t i = 5; i < names.size(); ++i) {
        digest = fnv1a64(names[i] + "\n", digest);
    }
    EXPECT_EQ(hex64(digest), capped.labels_dropped_digest);
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ(names[i], capped.label_dict[i].name);
    }
}

// trace support distinguishes a repeat occurrence from mere k-mer presence
TEST_F(MiniRefSeq, TraceSupportOnRealCoordinates) {
    ResolveOptions options;
    options.labels = { "NZ_CP034369.1" };            // two blaNDM copies
    options.support = Support::KMER;
    auto by_kmer = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, by_kmer.labels.size());
    EXPECT_EQ(783u, by_kmer.labels[0].kmers_supported);

    options.support = Support::TRACE;
    auto by_trace = resolve_support(*oracle_, query_, options);
    ASSERT_EQ(1u, by_trace.labels.size());
    // the gene occurs contiguously, so the trace is unbroken
    EXPECT_EQ(783u, by_trace.labels[0].kmers_supported);
    EXPECT_TRUE(by_trace.labels[0].trace_breaks.empty());
    EXPECT_EQ(1u, by_trace.labels[0].runs.size());

    // the annotation access path used for header labels is the coordinate path
    EXPECT_GT(oracle_->counters().tuple_rows_fetched, 0u);
    EXPECT_GT(oracle_->counters().coords_mapped, 0u);
}

// The decisive end-to-end check: extend the blaNDM gene through a real index and
// verify every recovered flank against the source record it is attributed to.
// A label reaching a leaf with route_bp == 0 travelled the path's own spelled bases,
// so seed+flank must occur verbatim in that accession. Labels merged in at a
// reconvergence (route_bp > 0) support the leaf node but not the spelled bases, and
// are explicitly excluded — that distinction is the point of the field.
TEST_F(MiniRefSeq, ExtendBlaNDMAndValidateAgainstSourceRecords) {
    // the source records, by accession
    std::map<std::string, std::string> records;
    for (const auto &entry : std::filesystem::directory_iterator(kIndexDir + "/fasta")) {
        std::ifstream in(entry.path());
        std::string line, name, seq;
        auto flush = [&]() {
            if (!name.empty())
                records[name] = seq;
        };
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '>') {
                flush();
                name = line.substr(1, line.find_first_of(" \t\r") - 1);
                seq.clear();
            } else {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                seq += line;
            }
        }
        flush();
    }
    ASSERT_EQ(42u, records.size());

    // freeze the longest seed: the whole gene, carried by its full-length carriers
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;
    auto profile = resolve_support(*oracle_, query_, options);
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::LONGEST_FIRST;
    policy.max_seeds = 1;
    auto selection = select_seeds(profile, query_, policy, false);
    ASSERT_FALSE(selection.seeds.empty());
    const FrozenSeed &frozen = selection.seeds[0];
    ASSERT_EQ(19u, frozen.labels.size());

    Seed seed;
    seed.seed_id = frozen.seed_id;
    seed.sequence = frozen.sequence;
    seed.labels = frozen.labels;

    Strategy strategy;
    strategy.max_extension_bp = 3000;
    strategy.profile_bin_bp = 500;
    auto result = traverse_seed(*oracle_, seed, strategy, LabelChangeCost::forbid());
    EXPECT_EQ(19u, result.num_seed_labels);
    EXPECT_TRUE(result.dropped_labels.empty());

    size_t checked = 0, merged_in = 0;
    for (size_t a : { static_cast<size_t>(Arm::LEFT), static_cast<size_t>(Arm::RIGHT) }) {
        const ArmResult &arm = result.arms[a];
        ASSERT_TRUE(arm.requested);
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << "arm " << a;
        // the label set can only shrink along an arm, so leaves stay bounded by |S|
        EXPECT_LE(arm.paths.size(), 19u);
        ASSERT_FALSE(arm.paths.empty());

        for (const auto &path : arm.paths) {
            const std::string flank = spell_path(arm, path);
            const std::string full = a == static_cast<size_t>(Arm::RIGHT)
                                         ? seed.sequence + flank
                                         : flank + seed.sequence;
            for (const auto &end : path.end_labels) {
                const std::string &accession = result.label_dict[end.label].name;
                ASSERT_TRUE(records.count(accession)) << accession;
                if (end.route_bp) {
                    merged_in++;
                    continue;   // supports the leaf node, not these spelled bases
                }
                EXPECT_EQ(0.0, end.loss);   // forbid: no switching
                EXPECT_NE(std::string::npos, records.at(accession).find(full))
                    << "arm " << a << " path " << path.id << " (" << full.size() << " bp)"
                    << " is attributed to " << accession << " but does not occur in it";
                checked++;
            }
        }
    }
    // the fixture really does exercise both cases
    EXPECT_GE(checked, 20u);
    EXPECT_GT(merged_in, 0u) << "expected at least one label merged in at a reconvergence";

    // the left context of a mobile element fragments faster than the right one:
    // this is the tuning evidence an agent reads off the growth profile
    const ArmResult &left = result.arms[static_cast<size_t>(Arm::LEFT)];
    const ArmResult &right = result.arms[static_cast<size_t>(Arm::RIGHT)];
    ASSERT_FALSE(left.growth.empty());
    ASSERT_FALSE(right.growth.empty());
    EXPECT_EQ(19u, left.growth.front().max_live_labels);
    EXPECT_EQ(19u, right.growth.front().max_live_labels);
    EXPECT_LT(left.growth.back().max_live_labels, right.growth.back().max_live_labels);
}

// The verification contract (spec §6.9) on the real index. T is the structural trie
// around the whole blaNDM gene with the labels recorded, E its walks filtered by the
// 19 derived carriers with test-side code (tests/graph/traversal/test_trie_oracle.hpp),
// A the label-constrained exhaustive trie over the same carriers: leaves(A) == leaves(E)
// at radius 100 and 300 on both arms, every run complete. The structural trie's size
// against the constrained one is how much work the label constraint does at this locus;
// it is recorded as test properties and printed.
TEST_F(MiniRefSeq, StructuralOracleMatchesTheConstrainedTrie) {
    for (uint64_t radius : { 100u, 300u }) {
        const std::string tag = "r" + std::to_string(radius);
        Strategy sc;
        sc.exhaustive = true;
        sc.max_label_branches = Strategy::kUnlimited;
        sc.max_splits_per_path = Strategy::kUnlimited;
        sc.merge_reconverge = false;
        sc.max_extension_bp = radius;
        Seed seed;
        seed.sequence = query_;          // no labels: the carriers are derived
        auto A = traverse_seed(*oracle_, seed, sc, LabelChangeCost::forbid());
        ASSERT_EQ(19u, A.num_seed_labels);
        std::set<std::string> carriers;
        for (size_t i = 0; i < A.num_seed_labels; ++i) carriers.insert(A.label_dict[i].name);

        Strategy st = sc;
        st.label_mode = LabelMode::ANNOTATE;
        st.max_labels_per_node = 100;    // 42 accessions in the index: nothing is cut
        auto T = traverse_seed(*oracle_, seed, st, LabelChangeCost::forbid());
        EXPECT_EQ(0u, T.num_seed_labels);

        for (size_t a : { static_cast<size_t>(Arm::LEFT), static_cast<size_t>(Arm::RIGHT) }) {
            const ArmResult &ta = T.arms[a], &aa = A.arms[a];
            const std::string what = tag + " " + to_string(ta.arm);
            ASSERT_EQ(ArmResult::COMPLETE, ta.status) << what << " structural";
            ASSERT_EQ(ArmResult::COMPLETE, aa.status) << what << " constrained";
            EXPECT_EQ(radius, ta.complete_to_bp);
            EXPECT_EQ(radius, aa.complete_to_bp);
            EXPECT_EQ(0u, ta.nodes_labels_truncated) << what;
            EXPECT_LE(ta.max_labels_at_node, 42u);

            auto v = trie::verify(T, A, a, carriers);
            EXPECT_EQ(radius, v.depth);
            EXPECT_FALSE(v.cut);
            trie::expect_equal(v.expected, v.actual, what);
            EXPECT_FALSE(v.actual.empty()) << what;

            // the contrast: structural vs label-constrained trie
            const std::string key = tag + "_" + to_string(ta.arm) + "_";
            RecordProperty(key + "structural_leaves", static_cast<int>(ta.paths.size()));
            RecordProperty(key + "structural_segments", static_cast<int>(ta.segments.size()));
            RecordProperty(key + "structural_splits", static_cast<int>(ta.splits.size()));
            RecordProperty(key + "structural_bp", static_cast<int>(ta.output_bp));
            RecordProperty(key + "structural_max_labels_at_node", static_cast<int>(ta.max_labels_at_node));
            RecordProperty(key + "constrained_leaves", static_cast<int>(aa.paths.size()));
            RecordProperty(key + "constrained_segments", static_cast<int>(aa.segments.size()));
            RecordProperty(key + "constrained_splits", static_cast<int>(aa.splits.size()));
            RecordProperty(key + "constrained_bp", static_cast<int>(aa.output_bp));
            std::cout << "[  METRIC  ] blaNDM radius " << radius << ' ' << to_string(ta.arm)
                      << " arm: structural trie " << ta.paths.size() << " leaves / "
                      << ta.segments.size() << " segments / " << ta.splits.size()
                      << " splits / " << ta.output_bp << " bp (up to "
                      << ta.max_labels_at_node << " labels at a node) vs label-constrained "
                      << aa.paths.size() << " leaves / " << aa.segments.size()
                      << " segments / " << aa.splits.size() << " splits / " << aa.output_bp
                      << " bp; oracle leaves " << v.expected.size() << std::endl;
        }
    }
}


/*
 * The trie contract (spec §6.9) on the real index, with the 42 source records as the
 * third party: the string model of test_trie_reference.hpp is built from the very
 * FASTA the index was made from, so on this graph — k = 31, basic, dummy k-mers not
 * masked, RowDiff<BRWT> with coordinates, header labels — the structural trie, the
 * label-constrained walker and the records can all be compared (test_trie_checks.hpp),
 * and support: trace against the records themselves.
 */
std::pair<std::vector<std::string>, std::vector<std::string>> load_source_records() {
    std::vector<std::filesystem::path> files;
    for (const auto &entry : std::filesystem::directory_iterator(kIndexDir + "/fasta"))
        files.push_back(entry.path());
    std::sort(files.begin(), files.end());
    std::vector<std::string> seqs, labels;
    for (const auto &path : files) {
        std::ifstream in(path);
        std::string line, name, seq;
        auto flush = [&]() {
            if (!name.empty()) {
                labels.push_back(name);
                seqs.push_back(seq);
            }
        };
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '>') {
                flush();
                name = line.substr(1, line.find_first_of(" \t\r") - 1);
                seq.clear();
            } else {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                seq += line;
            }
        }
        flush();
    }
    return { seqs, labels };
}

// the whole blaNDM gene as the longest frozen seed, as the extension test picks it
std::string blandm_seed(LabelOracle &oracle, const std::string &query) {
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 100;
    auto profile = resolve_support(oracle, query, options);
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::LONGEST_FIRST;
    policy.max_seeds = 1;
    auto selection = select_seeds(profile, query, policy, false);
    return selection.seeds.empty() ? "" : selection.seeds[0].sequence;
}

TEST_F(MiniRefSeq, TrieContractAgainstTheSourceRecords) {
    auto [seqs, labels] = load_source_records();
    ASSERT_EQ(42u, seqs.size());
    const size_t k = oracle_->get_k();
    ASSERT_EQ(31u, k);
    trie::PhaseTimer timer;
    trie::StringIndex ix(k, DeBruijnGraph::BASIC, seqs, labels);
    timer.lap("string model of 42 records");

    struct SeedCase { std::string name, seq; };
    std::vector<SeedCase> cases;
    const std::string blandm = blandm_seed(*oracle_, query_);
    ASSERT_FALSE(blandm.empty());
    cases.push_back({ "blaNDM", blandm });
    // and 150 bp windows cut from three records at 1 kb
    for (size_t i : { 0u, 17u, 33u }) {
        ASSERT_GT(seqs[i].size(), 1200u) << labels[i];
        cases.push_back({ "window:" + labels[i], seqs[i].substr(1000, 150) });
    }

    const uint64_t radius = 300;
    size_t structural_leaves = 0;
    for (const SeedCase &sc : cases) {
        const std::set<std::string> carriers = ix.carriers(sc.seq);
        ASSERT_FALSE(carriers.empty()) << sc.name;
        std::vector<std::set<std::string>> permitted { carriers, { *carriers.begin() } };
        if (carriers.size() > 1)
            permitted.push_back({ *std::prev(carriers.end()) });
        trie::CaseSpec c(sc.name, k, seqs, labels, sc.seq, radius, permitted);
        trie::CaseResult r = trie::check_case_on(*anno_graph_, ix, DeBruijnGraph::BASIC, c);
        EXPECT_EQ(carriers, r.carriers) << sc.name;
        for (size_t a : { trie::kLeft, trie::kRight })
            structural_leaves += r.T.arms[a].paths.size();
        // the index is unmasked: no '$' may ever reach the output
        for (size_t a : { trie::kLeft, trie::kRight }) {
            for (const auto &seg : r.T.arms[a].segments)
                EXPECT_EQ(std::string::npos, seg.sequence.find('$')) << sc.name;
        }
        timer.lap(sc.name + " contract");
        SeedResult trace;
        trie::check_trace(*anno_graph_, ix, c, carriers, r.A.at(carriers), &trace, sc.name);
        timer.lap(sc.name + " trace");
        if (sc.name == "blaNDM") {
            EXPECT_EQ(19u, carriers.size());
            // the shape of this locus at 300 bp, as first measured when the trie oracle
            // was built: 24 structural walks on the left, 2 on the right
            EXPECT_EQ(24u, r.T.arms[trie::kLeft].paths.size());
            EXPECT_EQ(2u, r.T.arms[trie::kRight].paths.size());
        }
    }
    EXPECT_GT(structural_leaves, 2 * cases.size());

    // deeper, label-constrained only (the structural trie is a small-radius tool): the
    // claims of the 19 carriers at 1000 bp under k-mer and under trace support, against
    // the records
    const std::set<std::string> carriers = ix.carriers(blandm);
    for (Support support : { Support::KMER, Support::TRACE }) {
        Strategy st = trie::exhaustive(LabelMode::CONSTRAIN, 1000);
        st.support = support;
        SeedResult A = trie::run(*anno_graph_, blandm, trie::as_list(carriers), st);
        EXPECT_TRUE(A.dropped_labels.empty());
        for (size_t a : { trie::kLeft, trie::kRight }) {
            const Arm arm = A.arms[a].arm;
            const std::string what = std::string("blaNDM 1000 bp ") + to_string(support)
                                     + " arm " + to_string(arm);
            EXPECT_EQ(ArmResult::COMPLETE, A.arms[a].status) << what;
            const trie::RefClaims ref = support == Support::KMER
                    ? trie::claims(ix, carriers, blandm, arm, 1000)
                    : trie::trace_claims(ix, carriers, blandm, arm, 1000);
            trie::expect_equal(ref.leaves, trie::constrained_claims(A, a, 1000), what);
            EXPECT_FALSE(ref.leaves.empty()) << what;
        }
        timer.lap(std::string("blaNDM 1000 bp ") + to_string(support));
    }
}

// Stage 3 of DESIGN-traverse-graphlet.md §14.1 on the real index (RowDiff<BRWT> with
// coordinates, whose reads are budget-aware): budgets swept across their stops — memory
// stops in both phases (traversal, annotation_decode), work stops — give the same result
// whatever annotation.batch_kmers is, and every read is charged with its dependency rows
TEST_F(MiniRefSeq, BudgetAwareReadsDoNotDependOnBatchKmers) {
    ASSERT_TRUE(oracle_->decode_charged());
    const std::string seed = query_.substr(0, 120);
    size_t decode_stops = 0, compared = 0;
    for (LabelMode mode : { LabelMode::CONSTRAIN, LabelMode::ANNOTATE }) {
        std::vector<std::pair<uint64_t, uint64_t>> budgets;     // memory bytes, work units
        for (uint64_t m = uint64_t(1) << 15; m <= (uint64_t(16) << 20); m = m * 3 / 2)
            budgets.emplace_back(m, 0);
        for (uint64_t w : { 2000, 20000, 200000, 2000000 })
            budgets.emplace_back(0, w);
        for (const auto &[memory, work] : budgets) {
            std::string reference;
            for (size_t batch : { 1, 3, 64, 1000 }) {
                Strategy st;
                st.label_mode = mode;
                st.max_extension_bp = 150;
                st.max_labels_per_node = mode == LabelMode::ANNOTATE ? 64 : st.max_labels_per_node;
                st.max_memory_bytes = memory;
                st.max_work_units = work;
                st.batch_kmers = batch;
                Seed s;
                s.sequence = seed;
                std::ostringstream os;
                try {
                    // an oracle per run, as per request: its counters accumulate
                    LabelOracle oracle(*anno_graph_);
                    const SeedResult r = traverse_seed(oracle, s, st, LabelChangeCost::forbid());
                    EXPECT_TRUE(r.account.decode_charged);
                    os << reduce(r) << "\nwork " << r.account.work_used << " soft "
                       << r.account.soft_overshoot << " keys " << r.annotation_counters.keys_mapped
                       << " requested " << r.annotation_counters.rows_requested;
                    for (const ArmResult &arm : r.arms) {
                        os << "\narm " << arm.complete_to_bp << " " << arm.work_units;
                    }
                    if (r.resource_stop) {
                        const ResourceStop &q = *r.resource_stop;
                        os << "\nstop " << q.resource << " " << q.phase << " " << q.at_bp
                           << " " << q.used << " " << q.demand;
                        EXPECT_TRUE(std::string(q.phase) == "traversal"
                                    || (std::string(q.phase) == "annotation_decode"
                                        && q.resource == ResourceStop::MEMORY));
                        decode_stops += batch == 64 && std::string(q.phase) == "annotation_decode";
                    }
                } catch (const SeedBudgetError &e) {
                    os << "failed " << e.what() << " " << e.stop().phase << " "
                       << e.account().work_used << " " << e.account().soft_overshoot;
                    decode_stops += batch == 64 && std::string(e.stop().phase) == "annotation_decode";
                }
                if (reference.empty()) {
                    reference = os.str();
                } else {
                    EXPECT_EQ(reference, os.str()) << to_string(mode) << " memory " << memory
                                                   << " work " << work << " batch_kmers " << batch;
                }
                compared++;
            }
        }
    }
    EXPECT_GT(compared, 80u);
    EXPECT_GT(decode_stops, 0u);
}


// Pass 5, the chunked deadlines: every read of a request decoded in chunks of one row (with
// the deadline checked between them) gives the response of the unchunked request, byte for
// byte, when no deadline is reached — constrain and annotate, derived and named labels, no
// budget, memory budgets near and at their stops (the budget-aware reads, refusals included)
// and work budgets
TEST_F(MiniRefSeq, PacedReadsAreByteIdentical) {
    const std::string seed = query_.substr(0, 120);
    size_t compared = 0, stopped = 0;
    for (const char *mode : { "constrain", "annotate" }) {
        for (int budget = 0; budget < 6; ++budget) {
            Json::Value r;
            r["seeds"][0]["sequence"] = seed;
            Json::Value &st = r["strategy"];
            st["labels"]["mode"] = mode;
            st["bounds"]["max_extension_bp"] = 150;
            st["bounds"]["time_budget_ms"] = 600000;
            st["output"]["detail"] = budget % 2 ? "graphlet" : "full";
            st["output"]["timing"] = false;
            if (std::string(mode) == "annotate")
                st["labels"]["max_labels_per_node"] = 64;
            switch (budget) {
                case 1: st["bounds"]["max_memory_mb"] = 1; break;
                case 2: st["bounds"]["max_memory_mb"] = 3; break;
                case 3: st["bounds"]["max_memory_mb"] = 16; break;
                case 4: st["bounds"]["max_work_units"] = 20000; break;
                case 5: st["bounds"]["max_work_units"] = 2000000; break;
                default: break;
            }
            mtg::cli::TraverseLimits whole, paced;
            paced.chunk_target_ms = 1e-9;
            const std::string a = mtg::cli::json_text(
                    mtg::cli::process_traverse_request(r, *anno_graph_, "", whole), true);
            const std::string b = mtg::cli::json_text(
                    mtg::cli::process_traverse_request(r, *anno_graph_, "", paced), true);
            EXPECT_EQ(a, b) << mode << " budget " << budget;
            stopped += a.find("\"resource_stop\"") != std::string::npos;
            compared++;
        }
    }
    EXPECT_EQ(12u, compared);
    EXPECT_GT(stopped, 0u);
}

// The efficiency pass, the row-diff path cache: the rows a request's reads reconstruct are
// kept so that later reads stop their row-diff paths at them. The response is byte for byte
// the one without the cache — constrain and annotate, derived labels, no budget, memory
// budgets near and at their stops (the cache is then each seed's, within what the label
// cache leaves of its allotment), work budgets, batch_kmers 1 and 64, a cache that keeps
// everything and one that evicts all the time, reads in one piece and in one-row chunks —
// and the cache is used (its hits) and emptied after a seed under a memory budget
TEST_F(MiniRefSeq, PathCacheKeepsTheResponse) {
    const std::string seed = query_.substr(0, 120);
    size_t compared = 0, stopped = 0;
    for (const char *mode : { "constrain", "annotate" }) {
        for (int budget = 0; budget < 6; ++budget) {
            for (int batch : { 1, 64 }) {
                Json::Value r;
                r["seeds"][0]["sequence"] = seed;
                Json::Value &st = r["strategy"];
                st["labels"]["mode"] = mode;
                st["bounds"]["max_extension_bp"] = 150;
                st["bounds"]["time_budget_ms"] = 600000;
                st["annotation"]["batch_kmers"] = batch;
                st["output"]["detail"] = budget % 2 ? "graphlet" : "full";
                st["output"]["timing"] = false;
                if (std::string(mode) == "annotate")
                    st["labels"]["max_labels_per_node"] = 64;
                switch (budget) {
                    case 1: st["bounds"]["max_memory_mb"] = 1; break;
                    case 2: st["bounds"]["max_memory_mb"] = 3; break;
                    case 3: st["bounds"]["max_memory_mb"] = 16; break;
                    case 4: st["bounds"]["max_work_units"] = 20000; break;
                    case 5: st["bounds"]["max_work_units"] = 2000000; break;
                    default: break;
                }
                mtg::cli::TraverseLimits off;
                const std::string reference = mtg::cli::json_text(
                        mtg::cli::process_traverse_request(r, *anno_graph_, "", off), true);
                stopped += reference.find("\"resource_stop\"") != std::string::npos;
                for (uint64_t bytes : { uint64_t(128) << 20, uint64_t(64) << 10 }) {
                    for (double chunk : { 0.0, 1e-9 }) {
                        mtg::cli::TraverseLimits on;
                        on.path_cache_bytes = bytes;
                        on.chunk_target_ms = chunk;
                        EXPECT_EQ(reference, mtg::cli::json_text(
                                mtg::cli::process_traverse_request(r, *anno_graph_, "", on), true))
                            << mode << " budget " << budget << " batch " << batch << " cache "
                            << bytes << " chunk " << chunk;
                        compared++;
                    }
                }
            }
        }
    }
    EXPECT_EQ(96u, compared);
    EXPECT_GT(stopped, 0u);

    // the cache is used, kept from seed to seed without a memory budget, and a seed's own
    // (emptied when it ends) under one
    for (uint64_t memory : { uint64_t(0), uint64_t(16) << 20 }) {
        LabelOracle oracle(*anno_graph_);
        oracle.set_path_cache_max(uint64_t(128) << 20);
        Strategy st;
        st.max_extension_bp = 150;
        st.max_memory_bytes = memory;
        st.batch_kmers = 1;
        Seed s;
        s.sequence = seed;
        traverse_seed(oracle, s, st, LabelChangeCost::forbid());
        EXPECT_GT(oracle.path_cache().hits(), 0u) << memory;
        if (memory) {
            EXPECT_EQ(0u, oracle.path_cache().bytes());
        } else {
            EXPECT_GT(oracle.path_cache().bytes(), 0u);
        }
        EXPECT_EQ(uint64_t(128) << 20, oracle.path_cache().limit());
    }
}

// Review of the efficiency pass, finding 1: under a memory budget the path cache is off until
// the depth-0 state is admitted. A refused annotate root states whether its read alone was
// refused (a lower bound) or its standalone demand did not fit; with the cache on while the
// roots were read, the left root's row-diff path held the right root's row, whose read then
// completed and was refused by its demand ("the row's standalone demand ... is 28400 bytes")
// where the seed without the cache states "the row, read alone ..., needs more than the
// 18850 bytes". The seed_id, charged with the depth-0 state, moves the account across the
// right root's refusal: every such response is the same with the cache on and off
TEST_F(MiniRefSeq, PathCacheKeepsRootRefusals) {
    // 32 bp of a mini refseq record (mini_batch3): the left root's path holds the right root
    const std::string seed = "AAGCGGGGACATTCTTCTCGGCTGACTCAGTC";
    auto request = [&](size_t id_bytes) {
        Json::Value r;
        r["seeds"][0]["sequence"] = seed;
        r["seeds"][0]["seed_id"] = std::string(id_bytes, 'x');
        Json::Value &st = r["strategy"];
        st["labels"]["mode"] = "annotate";
        st["labels"]["max_labels_per_node"] = 1;
        st["bounds"]["max_extension_bp"] = 33;
        st["bounds"]["max_memory_mb"] = 1;
        st["output"]["detail"] = "summary";
        st["output"]["timing"] = false;
        return r;
    };
    mtg::cli::TraverseLimits off, on;
    on.path_cache_bytes = uint64_t(128) << 20;
    auto text = [&](const Json::Value &r, const mtg::cli::TraverseLimits &limits) {
        return mtg::cli::json_text(mtg::cli::process_traverse_request(r, *anno_graph_, "", limits),
                                   true);
    };
    auto refused = [](const std::string &out, const char *arm) {
        return out.find(std::string("annotation row of the ") + arm + " arm's root")
                != std::string::npos;
    };
    // where the right root's read is refused: between the walks and the left root's refusals
    std::vector<size_t> window;
    for (size_t bytes = 0; bytes < (size_t(1) << 20); bytes += 8192) {
        const std::string out = text(request(bytes), off);
        if (refused(out, "right"))
            window.push_back(bytes);
        if (refused(out, "left"))
            break;
    }
    ASSERT_FALSE(window.empty());
    size_t compared = 0, read_alone = 0;
    for (size_t bytes = window.front() > 8192 ? window.front() - 8192 : 0;
            bytes < window.back() + 8192; bytes += 256) {
        const Json::Value r = request(bytes);
        const std::string a = text(r, off);
        if (!refused(a, "right"))
            continue;
        compared++;
        read_alone += a.find("the row, read alone with its row-diff dependency rows, needs more")
                        != std::string::npos;
        EXPECT_EQ(a, text(r, on)) << "seed_id of " << bytes << " bytes";
    }
    EXPECT_GT(compared, 10u);
    // the refusals the cache changed: a read alone refused (its lower bound stated)
    EXPECT_GT(read_alone, 0u);
}

// R16 (the owner's decision; feature level 4): the server's maxima of the budgets. Set, a
// larger budget is lowered to the maximum and an omitted one set to it, both echoed in
// strategy.clamped (requested "unlimited" for an omitted one) and in strategy.bounds; a smaller
// budget is kept. A seed its clamped budget stopped states it as server_clamp, and its stop
// what the request asked for ("unlimited" when it gave none: representable in a graphlet's
// K and Q records, which null is not). Unset (0), nothing changes.
TEST_F(MiniRefSeq, ServerBudgetMaxima) {
    const std::string seed = query_.substr(0, 120);
    auto request = [&](Json::Value bounds) {
        Json::Value r;
        r["seeds"][0]["sequence"] = seed;
        bounds["max_extension_bp"] = 150;
        r["strategy"]["bounds"] = bounds;
        r["strategy"]["output"]["detail"] = "summary";
        r["strategy"]["output"]["timing"] = false;
        return r;
    };
    auto clamped = [](const Json::Value &out, const std::string &field) {
        for (const Json::Value &c : out["strategy"]["clamped"]) {
            if (c["field"].asString() == field)
                return c;
        }
        return Json::Value();
    };
    mtg::cli::TraverseLimits limits;
    limits.max_memory_mb = 16;
    limits.max_work_units = 20000;
    // omitted: set to the maxima
    Json::Value out = mtg::cli::process_traverse_request(request(Json::Value()), *anno_graph_,
                                                         "", limits);
    Json::Value m = clamped(out, "bounds.max_memory_mb");
    ASSERT_FALSE(m.isNull());
    EXPECT_EQ("unlimited", m["requested"].asString());
    EXPECT_EQ(16u, m["effective"].asUInt64());
    Json::Value w = clamped(out, "bounds.max_work_units");
    ASSERT_FALSE(w.isNull());
    EXPECT_EQ("unlimited", w["requested"].asString());
    EXPECT_EQ(20000u, w["effective"].asUInt64());
    EXPECT_EQ(16u, out["strategy"]["bounds"]["max_memory_mb"].asUInt64());
    EXPECT_EQ(20000u, out["strategy"]["bounds"]["max_work_units"].asUInt64());
    // the work maximum stops this walk: stated as the server's clamp
    const Json::Value &result = out["results"][0];
    ASSERT_TRUE(result.isMember("resource_stop")) << mtg::cli::json_text(result, true);
    EXPECT_EQ("work", result["resource_stop"]["resource"].asString());
    EXPECT_EQ("unlimited", result["resource_stop"]["requested"].asString());
    bool stated = false;
    for (const Json::Value &l : result["limitations"]) {
        stated |= l["kind"].asString() == "server_clamp"
                && l["knob"].asString() == "bounds.max_work_units";
    }
    EXPECT_TRUE(stated) << mtg::cli::json_text(result["limitations"], true);
    // a request cannot raise it: no action raising it, the walk_domain states the maximum
    // (review of the efficiency pass, finding 2)
    for (const Json::Value &a : result["resource_stop"]["actions"]) {
        EXPECT_NE("raise_work_budget", a.asString());
    }
    size_t at_max = 0;
    for (const char *arm : { "left", "right" }) {
        for (const Json::Value &l : result["arms"][arm]["limitations"]) {
            if (l["kind"].asString() != "walk_domain"
                    || l["knob"].asString() != "bounds.max_work_units")
                continue;
            at_max++;
            EXPECT_EQ(20000u, l["server_limit"].asUInt64());
            EXPECT_EQ(std::string::npos, l["effect"].asString().find("raise the knob"));
            EXPECT_NE(std::string::npos, l["effect"].asString().find("the server's maximum"));
        }
    }
    EXPECT_GT(at_max, 0u);
    // representable as a graphlet (its K and Q records)
    Json::Value as_graphlet = request(Json::Value());
    as_graphlet["strategy"]["output"]["detail"] = "graphlet";
    out = mtg::cli::process_traverse_request(as_graphlet, *anno_graph_, "", limits);
    ASSERT_TRUE(out["results"][0].isMember("graphlet")) << mtg::cli::json_text(out, true);
    EXPECT_NE(std::string::npos, out["results"][0]["graphlet"].asString()
                                     .find("server_clamp bounds.max_work_units i:20000 u "));
    // larger: lowered (the requested value echoed); smaller: kept, nothing clamped
    Json::Value larger;
    larger["max_memory_mb"] = 64;
    larger["max_work_units"] = 1000000;
    out = mtg::cli::process_traverse_request(request(larger), *anno_graph_, "", limits);
    EXPECT_EQ(64u, clamped(out, "bounds.max_memory_mb")["requested"].asUInt64());
    EXPECT_EQ(16u, clamped(out, "bounds.max_memory_mb")["effective"].asUInt64());
    EXPECT_EQ(1000000u, clamped(out, "bounds.max_work_units")["requested"].asUInt64());
    EXPECT_EQ(20000u, out["strategy"]["bounds"]["max_work_units"].asUInt64());
    Json::Value smaller;
    smaller["max_memory_mb"] = 8;
    smaller["max_work_units"] = 10000;
    out = mtg::cli::process_traverse_request(request(smaller), *anno_graph_, "", limits);
    EXPECT_TRUE(clamped(out, "bounds.max_memory_mb").isNull());
    EXPECT_TRUE(clamped(out, "bounds.max_work_units").isNull());
    EXPECT_EQ(8u, out["strategy"]["bounds"]["max_memory_mb"].asUInt64());
    EXPECT_EQ(10000u, out["strategy"]["bounds"]["max_work_units"].asUInt64());
    // unset: the request as it was
    const std::string a = mtg::cli::json_text(
            mtg::cli::process_traverse_request(request(Json::Value()), *anno_graph_, "", {}), true);
    EXPECT_EQ(std::string::npos, a.find("max_memory_mb"));
    EXPECT_EQ(std::string::npos, a.find("max_work_units"));

    // Review of the efficiency pass, finding 2: a seed the server's maximum fails before any
    // walk states it as a walked seed does — the server_clamp, the stop's `requested` (what the
    // request asked for), no action raising the budget and a walk_domain with server_limit that
    // says a request cannot raise it — whether it failed in its seed phase (work: the derived
    // set's rows) or at its depth-0 state (memory: an annotate root's read)
    auto failed_at_max = [&](const Json::Value &r, const mtg::cli::TraverseLimits &limits,
                             const std::string &field, const char *resource,
                             const Json::Value &requested, uint64_t maximum) {
        // a knob's value: an integer, or "unlimited"
        auto same = [](const Json::Value &a, const Json::Value &b) {
            return a.isString() ? b.isString() && a.asString() == b.asString()
                                : b.isNumeric() && a.asUInt64() == b.asUInt64();
        };
        const Json::Value out = mtg::cli::process_traverse_request(r, *anno_graph_, "", limits);
        const Json::Value &res = out["results"][0];
        const std::string where = mtg::cli::json_text(res, true);
        ASSERT_TRUE(res.isMember("error")) << where;
        EXPECT_EQ("failed", res["outcome"]["walks"].asString()) << where;
        ASSERT_TRUE(res.isMember("resource_stop")) << where;
        EXPECT_EQ(resource, res["resource_stop"]["resource"].asString()) << where;
        EXPECT_TRUE(same(requested, res["resource_stop"]["requested"])) << where;
        EXPECT_EQ(maximum, res["resource_stop"]["effective"].asUInt64()) << where;
        for (const Json::Value &action : res["resource_stop"]["actions"]) {
            EXPECT_EQ(std::string::npos, action.asString().find("raise_")) << where;
        }
        size_t clamps = 0, domains = 0;
        for (const Json::Value &l : res["limitations"]) {
            if (l["knob"].asString() != field)
                continue;
            if (l["kind"].asString() == "server_clamp") {
                clamps++;
                EXPECT_EQ(maximum, l["limit"].asUInt64());
                EXPECT_TRUE(same(requested, l["observed"])) << where;
            } else if (l["kind"].asString() == "walk_domain") {
                domains++;
                EXPECT_EQ(maximum, l["server_limit"].asUInt64());
                EXPECT_EQ(std::string::npos, l["effect"].asString().find("raise the knob"))
                    << where;
                EXPECT_NE(std::string::npos, l["effect"].asString().find("the server's maximum"));
            }
        }
        EXPECT_EQ(1u, clamps) << where;
        EXPECT_EQ(1u, domains) << where;
    };
    // the seed phase: the derived set's rows of an 813 bp seed under 200 units, the request
    // without a budget and with a larger one
    Json::Value derive;
    derive["seeds"][0]["sequence"] = query_;
    derive["strategy"]["bounds"]["max_extension_bp"] = 2000;
    derive["strategy"]["output"]["detail"] = "full";
    derive["strategy"]["output"]["timing"] = false;
    mtg::cli::TraverseLimits work;
    work.max_work_units = 200;
    failed_at_max(derive, work, "bounds.max_work_units", "work", Json::Value("unlimited"), 200);
    derive["strategy"]["bounds"]["max_work_units"] = 5000;
    failed_at_max(derive, work, "bounds.max_work_units", "work", Json::Value(5000), 200);
    // the same seed failed by the request's own budget keeps the lever to raise it
    work.max_work_units = 0;
    derive["strategy"]["bounds"]["max_work_units"] = 200;
    const Json::Value own = mtg::cli::process_traverse_request(derive, *anno_graph_, "", work);
    EXPECT_EQ("raise_work_budget", own["results"][0]["resource_stop"]["actions"][0].asString());
    for (const Json::Value &l : own["results"][0]["limitations"]) {
        EXPECT_NE("server_clamp", l["kind"].asString());
        EXPECT_FALSE(l.isMember("server_limit"));
    }
    // the depth-0 state: the right annotate root's read under the server's 1 MiB (the seed_id,
    // charged with the depth-0 state, scanned to where that read is refused)
    mtg::cli::TraverseLimits memory;
    memory.max_memory_mb = 1;
    bool root_refused = false;
    for (size_t bytes = 0; bytes < (size_t(1) << 20) && !root_refused; bytes += 8192) {
        Json::Value r;
        r["seeds"][0]["sequence"] = "AAGCGGGGACATTCTTCTCGGCTGACTCAGTC";
        r["seeds"][0]["seed_id"] = std::string(bytes, 'x');
        r["strategy"]["labels"]["mode"] = "annotate";
        r["strategy"]["labels"]["max_labels_per_node"] = 1;
        r["strategy"]["bounds"]["max_extension_bp"] = 33;
        r["strategy"]["output"]["detail"] = "summary";
        r["strategy"]["output"]["timing"] = false;
        const Json::Value out = mtg::cli::process_traverse_request(r, *anno_graph_, "", memory);
        if (out["results"][0]["error"].asString().find("arm's root") == std::string::npos)
            continue;
        root_refused = true;
        failed_at_max(r, memory, "bounds.max_memory_mb", "memory", Json::Value("unlimited"), 1);
    }
    EXPECT_TRUE(root_refused);
}

// The efficiency pass, C4: the coordinate mapping by runs (LabelOracle::map_coords) against
// one map_single_coord per coordinate, on the tuple rows of the blaNDM query (same results,
// checked; the time printed). A measurement, not a check of speed: run it with
// --gtest_also_run_disabled_tests
TEST_F(MiniRefSeq, DISABLED_CoordMappingBenchmark) {
    std::vector<Row> rows;
    for (node_index key : oracle_->keys_of_sequence(query_)) {
        if (key != npos)
            rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    const auto tuples = oracle_->get_row_tuples(rows);
    uint64_t coords = 0;
    for (const auto &row : tuples) {
        for (const auto &entry : row) {
            coords += entry.second.size();
        }
    }
    ASSERT_GT(coords, 0u);
    const int reps = 200;
    uint64_t sum_a = 0, sum_b = 0;
    Timer timer;
    for (int rep = 0; rep < reps; ++rep) {
        for (const auto &row : tuples) {
            for (const auto &[c, cs] : row) {
                for (Coord coord : cs) {
                    const auto [seq, local] = oracle_->map_coord(c, coord);
                    sum_a += seq * 31 + local;
                }
            }
        }
    }
    const double per_coord_a = timer.elapsed() * 1e9 / (double(coords) * reps);
    timer.reset();
    for (int rep = 0; rep < reps; ++rep) {
        for (const auto &row : tuples) {
            for (const auto &[c, cs] : row) {
                oracle_->map_coords(c, cs.data(), cs.size(), [&](Coord, uint64_t seq, Coord local) {
                    sum_b += seq * 31 + local;
                });
            }
        }
    }
    const double per_coord_b = timer.elapsed() * 1e9 / (double(coords) * reps);
    EXPECT_EQ(sum_a, sum_b);
    std::cerr << rows.size() << " rows, " << coords << " coordinates: map_single_coord "
              << per_coord_a << " ns a coordinate, map_coords " << per_coord_b
              << " ns (x" << per_coord_a / per_coord_b << ")" << std::endl;
}

} // namespace
