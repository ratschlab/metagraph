#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <thread>
#include <tuple>

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
#include "cli/traverse_attempts.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/unix_tools.hpp"
#include "common/seq_tools/reverse_complement.hpp"


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

// R6 (review of pass 5) and finding 7 of the review of its fixes: a /resolve decodes the
// present rows in batches of at most ResolveOptions::batch_rows rows sized to about
// batch_bytes (the first kResolveFirstBatchRows, each at most twice the one before), and holds
// one batch at a time — a discovery accumulates its labels' support while it reads (no profile
// pass decoding the rows again) and keeps a repeated k-mer's row for its later occurrences
// within kept_bytes, an explicit profile primes its query with each distinct row once. The
// profile is the same however the rows are batched or kept, and a discovery's profiles are
// exactly those of the same labels given explicitly: column and header labels, presence and
// trace, on the query and on a query repeating most of it (its k-mers occurring twice)
TEST_F(MiniRefSeq, ResolveDecodesEachRowOnceInBoundedBatches) {
    auto text = [](const SupportProfile &p) {
        std::ostringstream os;
        for (const auto &l : p.labels) {
            os << l.label.name << ' ' << l.kmers_supported << ' ' << encode_runs(l.runs, p.num_kmers);
            for (uint64_t b : l.trace_breaks) {
                os << ' ' << b;
            }
            os << '\n';
        }
        return os.str();
    };
    const std::string repeated = query_ + query_.substr(0, 2 * query_.size() / 3);
    size_t compared = 0;
    for (const std::string &query : { query_, repeated }) {
        const std::vector<node_index> keys = oracle_->keys_of_sequence(query);
        std::vector<node_index> present_keys;
        std::copy_if(keys.begin(), keys.end(), std::back_inserter(present_keys),
                     [](node_index key) { return key != npos; });
        const uint64_t present = present_keys.size();
        std::sort(present_keys.begin(), present_keys.end());
        const uint64_t distinct = std::unique(present_keys.begin(), present_keys.end())
                                - present_keys.begin();
        ASSERT_GT(distinct, 2 * kResolveFirstBatchRows);
        if (&query == &repeated) {
            ASSERT_GT(present, distinct + 2 * kResolveFirstBatchRows);
        }
        for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
            for (Support support : { Support::KMER, Support::TRACE }) {
                ResolveOptions options;
                options.discover = true;
                options.discover_kind = kind;
                options.support = support;
                options.discover_max_labels = 50;
                LabelOracle reference_oracle(*anno_graph_);
                const SupportProfile reference = resolve_support(reference_oracle, query, options);
                ASSERT_FALSE(reference.labels.empty());
                // the same labels given explicitly (their profile pass, LabelQuery's hits)
                ResolveOptions explicit_options;
                explicit_options.support = support;
                for (const auto &l : reference.labels) {
                    explicit_options.labels.push_back(l.label.name);
                }
                {
                    LabelOracle oracle(*anno_graph_);
                    EXPECT_EQ(text(reference), text(resolve_support(oracle, query, explicit_options)))
                        << int(kind) << ' ' << int(support);
                }
                // the doubling clamp (batch 7), one row a read (bytes 1), the re-reads without
                // kept rows (kept 0, batch 1) and the defaults
                for (size_t batch : { size_t(1), size_t(7), kResolveBatchRows }) {
                    for (uint64_t bytes : { uint64_t(1), kResolveBatchBytes }) {
                        for (uint64_t kept : { uint64_t(0), kResolveKeptBytes }) {
                            for (ResolveOptions *o : { &options, &explicit_options }) {
                                if (o == &explicit_options && kept != kResolveKeptBytes)
                                    continue;   // an explicit profile keeps no rows
                                o->batch_rows = batch;
                                o->batch_bytes = bytes;
                                o->kept_bytes = kept;
                                LabelOracle oracle(*anno_graph_);
                                std::vector<size_t> reads;
                                oracle.test_read_hook = [&](size_t rows) { reads.push_back(rows); };
                                EXPECT_EQ(text(reference), text(resolve_support(oracle, query, *o)))
                                    << batch << ' ' << bytes << ' ' << kept;
                                const uint64_t fetched = oracle.counters().rows_fetched
                                                       + oracle.counters().tuple_rows_fetched;
                                const bool direct = oracle.counters().direct_reads > 0;
                                compared++;
                                if (direct)
                                    continue;
                                // no read beyond batch_rows rows, the first at most
                                // kResolveFirstBatchRows, each at most twice the one before, and
                                // with a byte target below one row's bytes one row a read after
                                // the first
                                ASSERT_FALSE(reads.empty());
                                EXPECT_EQ(std::min(batch, kResolveFirstBatchRows), reads[0]);
                                for (size_t i = 0; i < reads.size(); ++i) {
                                    EXPECT_LE(reads[i], batch) << i;
                                    if (i > 0) {
                                        EXPECT_LE(reads[i], 2 * reads[i - 1]) << i;
                                    }
                                    if (bytes == 1 && i > 0) {
                                        EXPECT_EQ(1u, reads[i]) << i;
                                    }
                                }
                                // every distinct row once when the repeated ones can be kept
                                // (an explicit profile always), at most once per k-mer else
                                if (o == &explicit_options || kept == kResolveKeptBytes) {
                                    EXPECT_EQ(distinct, fetched) << batch << ' ' << bytes;
                                } else {
                                    EXPECT_LE(fetched, present) << batch << ' ' << bytes;
                                    EXPECT_GE(fetched, distinct) << batch << ' ' << bytes;
                                    if (kept == 0 && batch == 1 && &query == &repeated) {
                                        EXPECT_EQ(present, fetched);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    EXPECT_EQ(2u * 4 * (3 * 2 * (2 + 1)), compared);
}

// R6 (permitted-range filtering): a query of header labels looks a coordinate up among the
// ranges of its own sequences instead of mapping it to its sequence (rank and select); its
// hits are those of mapping every coordinate (map_coord) and keeping the requested sequences'
// — with and without coordinates, by the default and the budget-aware reads
TEST_F(MiniRefSeq, HeaderHitsFilterRequestedRanges) {
    ResolveOptions options;
    options.discover = true;
    options.discover_kind = LabelKind::HEADER;
    options.discover_max_labels = 12;
    const SupportProfile discovered = resolve_support(*oracle_, query_, options);
    ASSERT_GE(discovered.labels.size(), 6u);
    std::vector<LabelRef> refs;
    for (const auto &l : discovered.labels) {
        refs.push_back(l.label);
    }
    const std::vector<node_index> keys = oracle_->keys_of_sequence(query_);
    std::vector<Row> rows;
    std::vector<node_index> present;
    for (node_index key : keys) {
        if (key == npos)
            continue;
        present.push_back(key);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    const auto tuples = oracle_->get_row_tuples(rows);
    for (bool with_coords : { false, true }) {
        // the reference: every coordinate mapped, the requested sequences kept
        std::vector<LabelQuery::NodeHits> expected(rows.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            std::map<LabelId, std::vector<Coord>> by_label;
            for (const auto &[c, coords] : tuples[i]) {
                for (Coord coord : coords) {
                    const auto [seq_id, local] = oracle_->map_coord(c, coord);
                    for (LabelId l = 0; l < refs.size(); ++l) {
                        if (refs[l].column == c && refs[l].seq_id == seq_id)
                            by_label[l].push_back(local);
                    }
                }
            }
            for (auto &[l, locals] : by_label) {
                std::sort(locals.begin(), locals.end());
                LabelQuery::Hit hit{ l, {} };
                if (with_coords)
                    hit.coords.assign(locals.begin(), locals.end());
                expected[i].push_back(hit);
            }
        }
        LabelQuery query(*oracle_, refs, with_coords);
        EXPECT_EQ(expected, query.fetch(present)) << with_coords;
        if (oracle_->decode_charged()) {
            LabelQuery budgeted(*oracle_, refs, with_coords);
            annot::matrix::DecodeBudget budget;
            std::vector<LabelQuery::NodeHits> out;
            std::vector<KeyCost> costs;
            out.reserve(present.size());
            costs.reserve(present.size());
            size_t refused_at = 0;
            ASSERT_TRUE(budgeted.fetch(present.data(), present.size(), budget, &out, &costs,
                                       &refused_at, nullptr));
            EXPECT_EQ(expected, out) << with_coords;
        }
    }
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


// R8 (staging at feature level 3: a warm tuple walk with a 5,000 ms budget ran 6,136 ms, and
// nothing said which piece overran): under output.timing every seed states its deadline record
// — its longest uninterruptible piece (kind, rows, coordinates mapped to headers in it) and,
// when a stop ended its walk, what stopped it and how long after its deadline — and nowhere
// else. A read maps its rows' coordinates inside its piece, so its measured rate includes the
// mapping: with slow reads the longest piece is a read of header rows, with its coordinates
TEST_F(MiniRefSeq, DeadlineRecordNamesTheLongestPiece) {
    const std::string seed = query_.substr(0, 120);
    // header labels: the reads map coordinates; slow reads (2 ms a row) make them the longest
    std::vector<std::string> headers;
    {
        Json::Value probe;
        probe["sequence"] = seed;
        probe["discover"]["kind"] = "header";
        probe["discover"]["max_labels"] = 3;
        const Json::Value out = mtg::cli::process_resolve_request(probe, *anno_graph_, "", 0,
                                                                   nullptr, nullptr);
        for (const Json::Value &l : out["labels"]) {
            headers.push_back(l["label"].asString());
        }
    }
    ASSERT_FALSE(headers.empty());
    for (bool stopped : { false, true }) {
        LabelOracle oracle(*anno_graph_);
        oracle.pacer().target_ms = stopped ? 1 : 50;
        oracle.test_read_hook = [](size_t rows) {
            std::this_thread::sleep_for(std::chrono::microseconds(2000 * rows));
        };
        Strategy st;
        st.max_extension_bp = 150;
        st.time_budget_ms = stopped ? 50 : 600000;
        Seed s;
        s.sequence = seed;
        s.labels = { headers[0] };
        const SeedResult r = traverse_seed(oracle, s, st, LabelChangeCost::forbid());
        const DeadlineRecord &d = r.deadline;
        EXPECT_GT(d.longest.ms, 0);
        EXPECT_GT(d.longest.rows, 0u) << d.longest.kind;
        const std::string kind = d.longest.kind;
        EXPECT_TRUE(kind == "read" || kind == "chunk" || kind == "rest") << kind;
        EXPECT_GT(d.longest.coordinates, 0u) << kind;
        if (stopped) {
            ASSERT_TRUE(r.resource_stop);
            EXPECT_EQ(ResourceStop::TIME, r.resource_stop->resource);
            EXPECT_EQ("time_budget", d.stopped_by);
            ASSERT_TRUE(d.after_deadline_ms);
            EXPECT_GE(*d.after_deadline_ms, 0);
        } else {
            EXPECT_FALSE(r.resource_stop);
            EXPECT_EQ("", d.stopped_by);
            EXPECT_FALSE(d.after_deadline_ms);
        }
    }
    // in the response's timing only
    for (bool timing : { false, true }) {
        Json::Value r;
        r["seeds"][0]["sequence"] = seed;
        r["seeds"][0]["labels"][0] = headers[0];
        r["strategy"]["bounds"]["max_extension_bp"] = 50;
        r["strategy"]["output"]["timing"] = timing;
        mtg::cli::TraverseLimits limits;
        limits.chunk_target_ms = 50;
        const Json::Value out = mtg::cli::process_traverse_request(r, *anno_graph_, "", limits);
        const std::string text = mtg::cli::json_text(out, true);
        EXPECT_EQ(timing, text.find("\"longest_piece\"") != std::string::npos);
        if (timing) {
            const Json::Value &d = out["results"][0]["timing"]["deadline"];
            for (const char *field : { "ms", "kind", "rows", "coordinates" }) {
                EXPECT_TRUE(d["longest_piece"].isMember(field)) << field;
            }
            EXPECT_FALSE(d.isMember("stopped_by"));
            EXPECT_TRUE(out["results"][0]["timing"].isMember("seed_phase_ms"));
        }
    }
}

// The efficiency pass, the row-diff path cache, and the chunked deadlines: the rows a request's
// reads reconstruct are kept so that later reads stop their row-diff paths at them, and every
// read may be decoded in chunks of one row (with the deadline checked between them). When no
// deadline is reached the response is byte for byte the one without the cache, in one piece —
// constrain and annotate, derived labels, no budget, memory budgets near and at their stops (the
// cache is then each seed's, within what the label cache leaves of its allotment; the
// budget-aware reads' refusals included), work budgets, batch_kmers 1 and 64, no cache, a cache
// that keeps everything and one that evicts all the time, reads in one piece and in one-row
// chunks — and the cache is used (its hits) and emptied after a seed under a memory budget
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
                for (uint64_t bytes : { uint64_t(0), uint64_t(128) << 20, uint64_t(64) << 10 }) {
                    for (double chunk : { 0.0, 1e-9 }) {
                        if (!bytes && !chunk)
                            continue;   // the reference itself
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
                // R10: the retention rule with no row narrow (mini refseq's rows all are), as
                // on a wide index — checkpoints of 4 and of 16, and only the requested rows
                for (const auto &rule : { std::make_tuple(4u, 1u, uint64_t(0)),
                                          std::make_tuple(16u, 8u, uint64_t(0)),
                                          std::make_tuple(1u << 30, 0u, uint64_t(0)) }) {
                    mtg::cli::TraverseLimits on;
                    on.path_cache_bytes = uint64_t(128) << 20;
                    on.path_cache_retention = rule;
                    EXPECT_EQ(reference, mtg::cli::json_text(
                            mtg::cli::process_traverse_request(r, *anno_graph_, "", on), true))
                        << mode << " budget " << budget << " batch " << batch << " checkpoint "
                        << std::get<0>(rule);
                    compared++;
                }
            }
        }
    }
    EXPECT_EQ(192u, compared);
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

// R6 (review of pass 5, path reuse): a seed walked with the path cache warm — the rows an
// earlier seed of the request kept (without a memory budget the cache is the request's) — gives
// the result it gives cold, alone in its request, and the result without the cache: the same
// walks, the same stops under memory and work budgets, at batch_kmers 1 and 64, in one piece
// and in one-row chunks, with a default, a tiny and no cache
TEST_F(MiniRefSeq, PathCacheWarmSeedIsTheColdSeed) {
    const std::string seed = query_.substr(0, 120);
    size_t compared = 0;
    for (int budget = 0; budget < 4; ++budget) {
        for (int batch : { 1, 64 }) {
            auto request = [&](size_t copies) {
                Json::Value r;
                for (size_t i = 0; i < copies; ++i) {
                    r["seeds"][Json::ArrayIndex(i)]["sequence"] = seed;
                }
                Json::Value &st = r["strategy"];
                st["bounds"]["max_extension_bp"] = 150;
                st["annotation"]["batch_kmers"] = batch;
                st["output"]["detail"] = "full";
                st["output"]["timing"] = false;
                switch (budget) {
                    case 1: st["bounds"]["max_memory_mb"] = 3; break;
                    case 2: st["bounds"]["max_work_units"] = 20000; break;
                    case 3: st["bounds"]["max_work_units"] = 2000000; break;
                    default: break;
                }
                return r;
            };
            mtg::cli::TraverseLimits off;
            Json::Value cold = mtg::cli::process_traverse_request(request(1), *anno_graph_, "",
                                                                  off)["results"][0];
            const std::string reference = mtg::cli::json_text(cold, true);
            for (uint64_t bytes : { uint64_t(0), uint64_t(64) << 10, uint64_t(128) << 20 }) {
                for (double chunk : { 0.0, 1e-9 }) {
                    mtg::cli::TraverseLimits on;
                    on.path_cache_bytes = bytes;
                    on.chunk_target_ms = chunk;
                    const Json::Value both = mtg::cli::process_traverse_request(
                            request(2), *anno_graph_, "", on)["results"];
                    ASSERT_EQ(2u, both.size());
                    EXPECT_EQ(reference, mtg::cli::json_text(both[0], true));
                    Json::Value warm = both[1];
                    EXPECT_TRUE(warm["duplicate"].asBool());
                    warm.removeMember("duplicate");
                    EXPECT_EQ(reference, mtg::cli::json_text(warm, true))
                        << "budget " << budget << " batch " << batch << " cache " << bytes
                        << " chunk " << chunk;
                    compared++;
                }
            }
        }
    }
    EXPECT_EQ(48u, compared);
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

// The server's maxima of the budgets (SPEC §10.3) on a seed they fail in its seed phase -- the
// derived set's rows of an 813 bp seed under the server's 200 work units -- which the
// integration test (TestTraverseWideIndex.test_server_budget_maxima: the clamps, their echo,
// server_clamp, the K and Q records, a depth-0 failure) leaves to this index. Unset (0),
// nothing changes.
TEST_F(MiniRefSeq, ServerBudgetMaxima) {
    // a seed the server's maximum fails before any walk states it as a walked seed does: the
    // server_clamp, the stop's `requested` (what the request asked for), no action raising the
    // budget and a walk_domain with server_limit that says a request cannot raise it
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
    // unset: the request as it was
    const std::string unset = mtg::cli::json_text(
            mtg::cli::process_traverse_request(derive, *anno_graph_, "", {}), true);
    EXPECT_EQ(std::string::npos, unset.find("server_clamp"));
    EXPECT_EQ(std::string::npos, unset.find("\"clamped\":[{"));
}

// ---------------------------------------------------------------- record coordinates (§18)

namespace {

// The source records: by accession (header labels: positions within the record) and, per
// column (a taxid: its FASTA file's name), in file order (column labels: record i's k-mer j has
// the column coordinate offset_i + j, offset_i the k-mers of the records before it)
struct SourceRecords {
    std::map<std::string, std::string> by_accession;
    std::map<std::string, std::vector<std::pair<std::string, uint64_t>>> by_column;
};

SourceRecords read_records(size_t k) {
    SourceRecords out;
    for (const auto &entry : std::filesystem::directory_iterator(kIndexDir + "/fasta")) {
        if (entry.path().extension() != ".fa")
            continue;
        const std::string column = entry.path().stem().string();
        std::ifstream in(entry.path());
        std::string line, name, seq;
        std::vector<std::pair<std::string, std::string>> records;
        auto flush = [&]() {
            if (!name.empty())
                records.emplace_back(name, seq);
        };
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            if (!line.empty() && line[0] == '>') {
                flush();
                name = line.substr(1, line.find_first_of(" \t") - 1);
                seq.clear();
            } else {
                for (char &c : line) c = std::toupper(static_cast<unsigned char>(c));
                seq += line;
            }
        }
        flush();
        uint64_t offset = 0;
        for (const auto &[acc, text] : records) {
            out.by_accession[acc] = text;
            out.by_column[column].emplace_back(text, offset);
            offset += text.size() >= k ? text.size() - k + 1 : 0;
        }
    }
    return out;
}

size_t overlapping_count(const std::string &hay, const std::string &needle) {
    size_t n = 0;
    for (size_t i = hay.find(needle); i != std::string::npos; i = hay.find(needle, i + 1)) {
        n++;
    }
    return n;
}

struct PositionalCheck {
    size_t runs = 0, occurrences = 0, seed_lists = 0, lower_bound = 0, strictly_lower = 0,
           switch_runs = 0, cut = 0, crossing = 0, record_end = 0;
};

// The positional oracle (DESIGN §18.3): every reported interval's bases in the source record are
// the run's own spelled bases (the seed for a seed occurrence), in natural orientation on both
// arms; the true count is the number of occurrences of the chain's string — the seed and the
// walk for a seed-entered run, the switch node's k-mer to the last node for a switch-entered one —
// exactly, and at least where the run is marked a lower bound
void check_positions(const SeedResult &r, const std::string &seed, const SourceRecords &src,
                     PositionalCheck *c, const std::string &what) {
    ASSERT_TRUE(r.coordinates_recorded) << what;
    const size_t k = r.k;
    auto records_of = [&](LabelId l) {
        const LabelRef &ref = r.label_dict[l];
        if (ref.kind == LabelKind::HEADER)
            return std::vector<std::pair<std::string, uint64_t>>{ { src.by_accession.at(ref.name), 0 } };
        return src.by_column.at(ref.name);
    };
    // The bases at [s, e) of every record whose numbering holds the interval: one, except for
    // a column label's interval in a record's last k - 1 bases, whose numbers are the next
    // record's first k - 1 positions too (review of W1, finding 6; the capabilities' rule says
    // so): an occurrence holds if one of them is the run's string. |second|: it was not the
    // first record holding the numbers (the ambiguity occurred)
    auto matches = [&](LabelId l, uint64_t s, uint64_t e, const std::string &want,
                       bool *second) -> int {
        int found = -1;
        size_t n = 0;
        for (const auto &[text, start] : records_of(l)) {
            if (s >= start && e - start <= text.size()) {
                if (found < 0 && text.compare(s - start, e - s, want) == 0) {
                    found = 1;
                    *second = n > 0;
                }
                n++;
            }
        }
        return n ? std::max(found, 0) : -1;      // -1: no record holds [s, e)
    };
    auto count = [&](LabelId l, const std::string &needle) {
        size_t n = 0;
        for (const auto &rec : records_of(l)) n += overlapping_count(rec.first, needle);
        return n;
    };
    for (size_t l = 0; l < r.seed_coordinates.size(); ++l) {
        const SeedCoordinates &sc = r.seed_coordinates[l];
        EXPECT_EQ(count(l, seed), sc.total) << what << " seed list of " << r.label_dict[l].name;
        for (Coord s : sc.starts) {
            bool second = false;
            const int m = matches(l, s, s + seed.size(), seed, &second);
            if (m < 0) {
                c->crossing++;
                continue;
            }
            EXPECT_EQ(1, m) << what << " seed [" << s << ", " << s + seed.size() << ")";
            c->record_end += second;
        }
        c->seed_lists++;
    }
    for (const ArmResult &arm : r.arms) {
        if (!arm.requested)
            continue;
        ASSERT_EQ(arm.runs.size(), arm.run_coordinates.size()) << what;
        for (size_t i = 0; i < arm.runs.size(); ++i) {
            const LabelRun &run = arm.runs[i];
            const RunCoordinates &rc = arm.run_coordinates[i];
            const std::string at = what + " " + to_string(arm.arm) + " run " + std::to_string(i);
            std::vector<size_t> chain;
            for (size_t s = run.segment; ; s = arm.segments[s].parents[0]) {
                chain.push_back(s);
                if (arm.segments[s].parents.empty())
                    break;
            }
            std::string flank;
            if (arm.arm == Arm::RIGHT) {
                for (auto it = chain.rbegin(); it != chain.rend(); ++it) flank += arm.segments[*it].sequence;
            } else {
                for (size_t s : chain) flank += arm.segments[s].sequence;
            }
            const uint64_t length = run.to_bp - run.from_bp, n = flank.size();
            std::string own, string;
            if (arm.arm == Arm::RIGHT) {
                own = flank.substr(run.from_bp, length);
                const std::string full = seed + flank;
                string = run.entered_by_switch ? full.substr(seed.size() + run.from_bp + 1 - k, length + k - 1)
                                               : full.substr(0, seed.size() + run.to_bp);
            } else {
                own = flank.substr(n - run.to_bp, length);
                const std::string full = flank + seed;
                string = run.entered_by_switch ? full.substr(n - run.to_bp, length + k - 1)
                                               : full.substr(n - run.to_bp);
            }
            for (Coord e : rc.ends) {
                const uint64_t s = arm.arm == Arm::RIGHT ? e + k - length : e;
                bool second = false;
                const int m = matches(run.label, s, s + length, own, &second);
                if (m < 0) {
                    c->crossing++;
                    continue;
                }
                EXPECT_EQ(1, m) << at << " [" << s << ", " << s + length << ")";
                c->record_end += second;
                c->occurrences++;
            }
            const size_t truth = count(run.label, string);
            if (rc.lower_bound) {
                EXPECT_LE(rc.total, truth) << at;
                c->lower_bound++;
                c->strictly_lower += rc.total < truth;
            } else {
                EXPECT_EQ(truth, rc.total) << at << (run.entered_by_switch ? " switch" : " seed");
            }
            c->runs++;
            c->switch_runs += run.entered_by_switch;
            c->cut += rc.total > rc.ends.size();
        }
    }
}

} // namespace

// The positional oracle on the real index (C2, C4; plan revision 1's switch cells): the whole
// blaNDM gene to 1,000 bp and a six-copy repeat window of NZ_CP030345.1 (150 bp), under the
// branch limit 0, 2 and the exhaustive preset, switch cells (a constant cost of 0.5 and 1 within
// a loss budget of 2, limits 0 and 2, to 3,000 bp as recorded) on these and on the seeds of the
// recorded switch requests (next/coords-plan/swreq_*: blaNDM reverse-complemented, three 200-bp
// windows of its carriers, a second repeat), caps 1 and "unlimited" (no list here holds more than
// 16 chains, asserted, so a cap of 16 is "unlimited"), header labels and column (taxid)
// labels, the latter also switching (where review of W1's finding 6 shows: an interval in a
// record's last k - 1 bases has the next record's first numbers)
TEST_F(MiniRefSeq, CoordinatesAgainstTheSourceRecords) {
    const SourceRecords src = read_records(31);
    ASSERT_EQ(42u, src.by_accession.size());
    const std::string repeat = "CAAAGTTAGCGATGAGGCAGCCTTTTGTCTTATTCAAAGGCCTTACATTTCAAAAACTCTGCTTACC"
                               "AGGCGCATTTCGCCCAGGGGATCACCATAATAAAATGCTGAGGCCTGGCCTTTGCGTAGTGCACGCAT"
                               "CACCTCAATACCTTT";
    ASSERT_EQ(150u, repeat.size());
    ASSERT_NE(std::string::npos, src.by_accession.at("NZ_CP030345.1").find(repeat));
    struct Cell {
        std::string name;
        std::string seed;
        std::function<void(Strategy*, LabelChangeCost*)> set;
    };
    std::vector<Cell> cells;
    for (const auto &[name, seed] : { std::make_pair(std::string("ndm1"), query_),
                                      std::make_pair(std::string("repeat"), repeat) }) {
        cells.push_back({ name + " limit 0", seed, [](Strategy*, LabelChangeCost*) {} });
        cells.push_back({ name + " limit 2", seed, [](Strategy *st, LabelChangeCost*) {
            st->max_label_branches = 2;
        } });
        cells.push_back({ name + " exhaustive", seed, [](Strategy *st, LabelChangeCost*) {
            st->exhaustive = true;
            st->max_label_branches = Strategy::kUnlimited;
            st->max_splits_per_path = Strategy::kUnlimited;
            st->max_extension_bp = 600;
        } });
        // the seeds of the recorded switch requests swreq_mini_ndm1_* and swreq_rep0_*: to
        // their radius of 3,000 bp, like the recorded ones below (review of W2: these four
        // cells a seed ran to the default 1,000 bp)
        for (double c : { 0.5, 1.0 }) {
            for (size_t limit : { size_t(0), size_t(2) }) {
                cells.push_back({ name + " switch " + std::to_string(c) + " limit " + std::to_string(limit),
                                  seed, [c, limit](Strategy *st, LabelChangeCost *cost) {
                    *cost = LabelChangeCost::constant(c);
                    st->loss_budget = 2;
                    st->max_label_branches = limit;
                    st->max_extension_bp = 3000;
                } });
            }
        }
    }
    cells.push_back({ "repeat column labels", repeat, [](Strategy *st, LabelChangeCost*) {
        st->seed_label_kind = LabelKind::COLUMN;
        st->max_label_branches = 2;
    } });
    // the seeds of the recorded switch requests (swreq_*): each under the four switch cells, to
    // the recorded requests' radius of 3,000 bp
    std::string ndm1_rc = query_;
    reverse_complement(ndm1_rc.begin(), ndm1_rc.end());
    const std::vector<std::pair<std::string, std::string>> recorded {
        { "ndm1_rc", ndm1_rc },
        { "win200_00", "ACCCGACCAAGGTCACCCGCACCGCGCTGCAGAACGCCGCGTCGATCGCGGGCCTGATGATCACCACCGAAGC"
                       "GATGGTGGCCGAGGCCCCGAAGAAGGACGAGCCGGCGATGCCGGCCGGCGGCGGCATGGGCGGCATGGGCGG"
                       "CATGGATTTCTAAGCCCCGCGATCCATCAAGCAAGACCACAAAGCCCGGCCTCGT" },
        { "win200_03", "GCGATCCTTCCAACTCGTCGCAAAGCCCAGCTTCGCATAAAACGCCTCTGTCACATCGAAATCGCGCGATGG"
                       "CAGATTGGGGGTGACGTGGTCAGCCATGGCTCAGCGCAGCTTGTCGGCCATGCGGGCCGTATGAGTGATTGC"
                       "GGCGCGGCTATCGGGGGCGGAATGGCTCATCACGATCATGCTGGCCTTGGGGAACG" },
        { "win200_05", "CGCCCCGTGCGGTTACGTCGAATGTCGCGGGCGCTTTGACATCGCGCGCAGCTGGCCAGATCGCCATGGTCG"
                       "GTTTGTTCGTCGATGCGGATGATGCTGTCATCGCCGACGCACTGGTGGCAGCCAAGCTGAACGCGCTGCAGC"
                       "TGCACGGTTCGGAATCGCCCGAACGCGTGGCCCAGTTGCGCGCGCGGTTTGGCAAG" },
        { "rep3", "AGCGGTAAATCGTGGAGTGATCGACATTCACTCCGCGTTCAGCCAGCATCTCCTGCAGCTCACGGTAACTGATG"
                  "CCGTATTTGCAGTACCAGCGTACGGCCCACAGAATGATGTCACGCTGAAAATGCCGGCCTTTGAATGGGTTCATGT" },
    };
    for (const auto &[name, seed] : recorded) {
        ASSERT_TRUE(name == "ndm1_rc" || seed.size() >= 150) << name;
        for (double c : { 0.5, 1.0 }) {
            for (size_t limit : { size_t(0), size_t(2) }) {
                cells.push_back({ name + " switch " + std::to_string(c) + " limit " + std::to_string(limit),
                                  seed, [c, limit](Strategy *st, LabelChangeCost *cost) {
                    *cost = LabelChangeCost::constant(c);
                    st->loss_budget = 2;
                    st->max_label_branches = limit;
                    st->max_extension_bp = 3000;
                } });
            }
        }
    }
    // column labels switching (the reviewer's col_sw cells, to the recorded requests' 3,000 bp)
    for (const auto &[name, seed] : { std::make_pair(std::string("ndm1"), query_),
                                      std::make_pair(std::string("repeat"), repeat),
                                      std::make_pair(std::string("ndm1_rc"), recorded[0].second),
                                      std::make_pair(std::string("win200_03"), recorded[2].second) }) {
        cells.push_back({ name + " column labels switch 0.5 limit 2", seed,
                          [](Strategy *st, LabelChangeCost *cost) {
            st->seed_label_kind = LabelKind::COLUMN;
            *cost = LabelChangeCost::constant(0.5);
            st->loss_budget = 2;
            st->max_label_branches = 2;
            st->max_extension_bp = 3000;
        } });
    }
    PositionalCheck total;
    // per seed, the deepest end of a run in its header-label switch cells: past the default
    // radius of 1,000 bp, so that the positions of what the recorded switch requests reach
    // beyond it are checked too (the review of W2 found the ndm1 and repeat seeds' switch cells
    // walked to 1,000 bp only); some single cells end sooner (a limit-0 walk from the repeat:
    // 36 bp)
    std::map<std::string, uint64_t> switch_deepest;
    for (const Cell &cell : cells) {
        for (size_t cap : { size_t(1), Strategy::kUnlimited }) {
            Strategy st;
            st.support = Support::TRACE;
            st.merge_reconverge = false;
            st.max_extension_bp = 1000;
            st.coordinates = true;
            st.max_coordinate_occurrences = cap;
            LabelChangeCost cost = LabelChangeCost::forbid();
            cell.set(&st, &cost);
            Seed seed;
            seed.sequence = cell.seed;
            const SeedResult r = traverse_seed(*oracle_, seed, st, cost);
            check_positions(r, cell.seed, src, &total, cell.name + " cap " + std::to_string(cap));
            if (cap == Strategy::kUnlimited) {
                for (const SeedCoordinates &sc : r.seed_coordinates) {
                    EXPECT_LE(sc.total, 16u) << cell.name;
                }
                for (const ArmResult &arm : r.arms) {
                    for (const RunCoordinates &rc : arm.run_coordinates) {
                        EXPECT_LE(rc.total, 16u) << cell.name;
                    }
                }
            }
            if (cell.name.find(" switch ") != std::string::npos
                    && cell.name.find("column") == std::string::npos) {
                uint64_t &deepest = switch_deepest[cell.name.substr(0, cell.name.find(' '))];
                for (const ArmResult &arm : r.arms) {
                    for (const LabelRun &run : arm.runs) {
                        deepest = std::max(deepest, run.to_bp);
                    }
                }
            }
        }
    }
    EXPECT_EQ(7u, switch_deepest.size());
    for (const auto &[name, deepest] : switch_deepest) {
        EXPECT_GT(deepest, 1000u) << name;
    }
    std::cerr << cells.size() << " cells x 2 caps: " << total.runs << " runs ("
              << total.switch_runs << " switch-entered, " << total.lower_bound
              << " lower bounds, " << total.strictly_lower << " strictly), " << total.occurrences
              << " occurrences, " << total.seed_lists << " seed lists, " << total.cut
              << " lists cut, " << total.crossing << " across records, " << total.record_end
              << " in a record's last k - 1 bases read as the next record's" << std::endl;
    EXPECT_GT(total.runs, 1000u);
    EXPECT_GT(total.switch_runs, 0u);
    EXPECT_GT(total.lower_bound, 0u);       // the switch cells mark some
    EXPECT_GT(total.cut, 0u);               // the repeat's six copies at cap 1
    EXPECT_EQ(0u, total.crossing);
    EXPECT_GT(total.record_end, 0u);        // the column record-end numbering is exercised
}

// Opt-in responses (trace, coordinates) are the same whatever the row-diff path cache holds,
// the reads' chunking and annotation.batch_kmers, unbudgeted and at their budget stops (memory
// and work): coordinates are read from the rows the walk reads, which none of these change
TEST_F(MiniRefSeq, CoordinatesKeepTheResponseUnderThePathCache) {
    const std::string seed = query_.substr(0, 120);
    size_t compared = 0, stopped = 0;
    for (int budget = 0; budget < 6; ++budget) {
        std::string reference;
        for (int batch : { 64, 1 }) {
            Json::Value r;
            r["seeds"][0]["sequence"] = seed;
            Json::Value &st = r["strategy"];
            st["support"] = "trace";
            st["branching"]["on_reconverge"] = "keep";
            st["branching"]["max_label_branches"] = 2;
            st["bounds"]["max_extension_bp"] = 400;
            st["bounds"]["time_budget_ms"] = 600000;
            st["annotation"]["batch_kmers"] = batch;
            st["output"]["detail"] = budget % 2 ? "graphlet" : "full";
            st["output"]["timing"] = false;
            st["output"]["coordinates"] = true;
            st["output"]["max_coordinate_occurrences"] = budget % 3 == 0 ? Json::Value(1) : Json::Value(16);
            switch (budget) {
                case 1: st["bounds"]["max_memory_mb"] = 2; break;
                case 2: st["bounds"]["max_memory_mb"] = 4; break;
                case 3: st["bounds"]["max_memory_mb"] = 16; break;
                case 4: st["bounds"]["max_work_units"] = 20000; break;
                case 5: st["bounds"]["max_work_units"] = 200000; break;
                default: break;
            }
            for (uint64_t bytes : { uint64_t(0), uint64_t(64) << 10, uint64_t(128) << 20 }) {
                for (double chunk : { 0.0, 1e-9, 1.0, 50.0 }) {
                    mtg::cli::TraverseLimits limits;
                    limits.path_cache_bytes = bytes;
                    limits.chunk_target_ms = chunk;
                    Json::Value out = mtg::cli::process_traverse_request(r, *anno_graph_, "", limits);
                    // the echo states the batch: compared without it
                    out["strategy"]["annotation"].removeMember("batch_kmers");
                    const std::string text = mtg::cli::json_text(out, true);
                    ASSERT_NE(std::string::npos, text.find("\"coordinates\":{")) << budget;
                    if (reference.empty()) {
                        reference = text;
                        stopped += text.find("\"resource_stop\"") != std::string::npos;
                    } else {
                        EXPECT_EQ(reference, text) << "budget " << budget << " batch " << batch
                                                   << " cache " << bytes << " chunk " << chunk;
                    }
                    compared++;
                }
            }
        }
    }
    EXPECT_EQ(144u, compared);
    EXPECT_GT(stopped, 0u);
}

// Plan revision 3 on the real index: a seed's opt-in result text less coordinates_text_bytes is
// its opt-out text — the delivery reserve's ratio sample of a request with coordinates is the
// one the request without them gives — in every detail, at caps 1 (the repeat's lists cut, a K
// record in a graphlet) and "unlimited" (16 is "unlimited" here: CoordinatesAgainstTheSourceRecords),
// for header and column labels, switch cells and the
// null form (support kmer); unbudgeted, so the walks are the same. compact_json_size is the
// writer's length of each whole response
TEST_F(MiniRefSeq, CoordinateShareIsExact) {
    const std::string repeat = "CAAAGTTAGCGATGAGGCAGCCTTTTGTCTTATTCAAAGGCCTTACATTTCAAAAACTCTGCTTACC"
                               "AGGCGCATTTCGCCCAGGGGATCACCATAATAAAATGCTGAGGCCTGGCCTTTGCGTAGTGCACGCAT"
                               "CACCTCAATACCTTT";
    const char *const trace = R"({"support": "trace", "branching": {"on_reconverge": "keep",
                                  "max_label_branches": 2}, "bounds": {"max_extension_bp": 1000}})";
    const std::vector<std::pair<std::string, std::string>> cells {
        { query_, trace },
        { repeat, trace },
        { repeat, R"({"support": "trace", "labels": {"seed_label_kind": "column"},
                      "branching": {"on_reconverge": "keep", "max_label_branches": 2},
                      "bounds": {"max_extension_bp": 1000}})" },
        { repeat, R"({"support": "trace", "labels": {"change_cost": {"model": "constant", "value": 0.5},
                      "loss_budget": 2}, "branching": {"on_reconverge": "keep", "max_label_branches": 2},
                      "bounds": {"max_extension_bp": 1000}})" },
        { query_, R"({"branching": {"on_reconverge": "keep"}, "bounds": {"max_extension_bp": 500}})" },
    };
    auto parse = [](const std::string &text) {
        Json::Value v;
        std::string errors;
        std::unique_ptr<Json::CharReader> reader(Json::CharReaderBuilder().newCharReader());
        EXPECT_TRUE(reader->parse(text.data(), text.data() + text.size(), &v, &errors)) << errors;
        return v;
    };
    size_t results = 0, cut = 0, graphlet_cut = 0;
    for (const auto &[seed, strategy] : cells) {
        for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
            Json::Value plain;
            plain["seeds"][0]["sequence"] = seed;
            plain["strategy"] = parse(strategy);
            plain["strategy"]["output"]["detail"] = detail;
            plain["strategy"]["output"]["timing"] = false;
            const Json::Value off = mtg::cli::process_traverse_request(plain, *anno_graph_, "");
            EXPECT_EQ(mtg::cli::json_text(off, true).size(), mtg::cli::compact_json_size(off));
            for (const Json::Value &cap : { Json::Value(1), Json::Value("unlimited") }) {
                Json::Value r = plain;
                r["strategy"]["output"]["coordinates"] = true;
                r["strategy"]["output"]["max_coordinate_occurrences"] = cap;
                const Json::Value on = mtg::cli::process_traverse_request(r, *anno_graph_, "");
                EXPECT_EQ(mtg::cli::json_text(on, true).size(), mtg::cli::compact_json_size(on));
                ASSERT_EQ(off["results"].size(), on["results"].size());
                for (Json::ArrayIndex i = 0; i < on["results"].size(); ++i) {
                    const Json::Value &a = on["results"][i], &b = off["results"][i];
                    const uint64_t coordinates = mtg::cli::coordinates_text_bytes(a);
                    EXPECT_GT(coordinates, 0u);
                    EXPECT_EQ(0u, mtg::cli::coordinates_text_bytes(b));
                    EXPECT_EQ(mtg::cli::json_text(b, true).size(),
                              mtg::cli::json_text(a, true).size() - coordinates)
                        << detail << " cap " << mtg::cli::json_text(cap, true) << " " << strategy;
                    const bool is_cut = a["coordinates"].isObject() && !a["coordinates"]["complete"].asBool();
                    cut += is_cut;
                    graphlet_cut += is_cut && a.isMember("graphlet")
                                  && a["graphlet"].asString().find("\nK * coordinates ") != std::string::npos;
                    results++;
                }
            }
        }
    }
    std::cerr << results << " results, " << cut << " incomplete, " << graphlet_cut
              << " with a coordinates K record" << std::endl;
    EXPECT_EQ(40u, results);
    EXPECT_GT(graphlet_cut, 0u);
}

// R21 (4) in annotate mode compares the parents by their OWN segment (the fewest labels present
// at a node of it, true counts; DESIGN §26.6), not by the walks they display. A parent that is
// itself a merged segment holds the union of its parents' labels, so after nested merges the
// rule can display a walk carried by fewer labels upstream. Pinned here as the rule stands, so
// that the library's check of it (derive.carried_labels) follows a known rule: on
// mini_win200_02 (the real cache's annotate_merge cell, right arm) the merge at 288 bp takes the
// 5-bp merged segment from 283 (6 labels) before the 32-bp segment from 256 (5 labels), though
// the walk through the former goes on through a 7-bp segment of 4 labels (276) — its path's
// continuation is carried by label 3 only, where level 5's (the 256 parent's) was by 0, 1, 2, 4
TEST_F(MiniRefSeq, AnnotateMergeRanksParentsByTheirOwnSegment) {
    Json::Value r;
    r["seeds"][0]["sequence"] = "TTAGCTTGGCGTGAGATTACCAATGTGTGACGGTTCGGTAGAGGCTTGCCGATAGACTCAAAGGTCTTTC"
                                "GCCCCATGACAACGACTTTTCCCTCAGTGAGTCTGCGAAAAATCTTCTGCTCACCCGGAATTTTCCAGG"
                                "GGATATTAGGACCATTGCCAATAACCCGATTGGCTCCCATCGCAGCAACGAGATAAATGCG";
    r["strategy"]["labels"]["mode"] = "annotate";
    r["strategy"]["branching"]["on_reconverge"] = "merge";
    r["strategy"]["bounds"]["max_extension_bp"] = 300;
    r["strategy"]["bounds"]["max_live_paths"] = 300;
    r["strategy"]["output"]["detail"] = "full";
    r["strategy"]["output"]["timing"] = false;
    const Json::Value out = mtg::cli::process_traverse_request(r, *anno_graph_, "");
    const Json::Value &arm = out["results"][0]["arms"]["right"];
    std::map<uint64_t, const Json::Value*> segments;
    for (const Json::Value &s : arm["segments"]) {
        segments[s["id"].asUInt64()] = &s;
    }
    auto fewest = [](const Json::Value &s) {
        uint64_t f = s["labels_total"].asUInt64();
        for (const Json::Value &run : s["label_sets"]) {
            f = std::min(f, run["labels_total"].asUInt64());
        }
        return f;
    };
    const Json::Value *merge = nullptr;
    for (const Json::Value &s : arm["segments"]) {
        if (s["from_bp"].asUInt64() == 288 && s["parents"].size() == 2)
            merge = &s;
    }
    ASSERT_NE(nullptr, merge);
    const Json::Value &first = *segments.at((*merge)["parents"][0].asUInt64());
    const Json::Value &second = *segments.at((*merge)["parents"][1].asUInt64());
    EXPECT_EQ(283u, first["from_bp"].asUInt64());
    EXPECT_EQ(5u, first["length_bp"].asUInt64());
    EXPECT_EQ(2u, first["parents"].size());            // itself a merge
    EXPECT_EQ(6u, fewest(first));
    EXPECT_EQ(256u, second["from_bp"].asUInt64());
    EXPECT_EQ(32u, second["length_bp"].asUInt64());
    EXPECT_EQ(5u, fewest(second));
    // the displayed walk through the first goes on through its own first parent: 4 labels
    const Json::Value &upstream = *segments.at(first["parents"][0].asUInt64());
    EXPECT_EQ(276u, upstream["from_bp"].asUInt64());
    EXPECT_EQ(4u, fewest(upstream));
    size_t through = 0;
    for (const Json::Value &path : arm["paths"]) {
        bool via = false;
        for (const Json::Value &id : path["segments"]) {
            via |= id.asUInt64() == (*merge)["id"].asUInt64();
        }
        if (!via)
            continue;
        through++;
        ASSERT_EQ(1u, path["continuation"]["labels"].size());
        EXPECT_EQ(3u, path["continuation"]["labels"][0].asUInt64());
    }
    EXPECT_EQ(1u, through);
}

// The review of 2026-10-06, U05-01 (W11), at the walker: under a memory budget the budget-aware
// lookahead decoded every run of a warm larger than its cache and evicted its own earlier runs,
// and the walk decoded the near rows again. On the review's walk (the first 70 bp of
// NZ_LPPQ01000025.1 to the right, 10,000 bp, constrain) 16 MiB decoded 31-39% more tuple rows
// than a work budget alone (whose caches have no byte bound) at batch_kmers 2,048-8,192, the
// walk the same. The warms stop before evicting their own runs now: +11-20% here, held below
// +25%. What remains is stated, not fixed (SPEC §6.8): a later warm's first run still evicts the
// cache wholesale, rows earlier warms read ahead that the walk had not reached among them.
// Evicting the oldest rows first instead kept the label cache full and starved the row-diff
// path cache that shares its allotment (2.5 times the stored rows read on mini_refseq)
TEST_F(MiniRefSeq, LookaheadKeepsItsRunsUnderAMemoryBudget) {
    const std::string fasta = kIndexDir + "/fasta/1296536.fa";
    std::ifstream in(fasta);
    if (!in)
        GTEST_SKIP() << fasta << " not found (scripts/traversal/build_mini_refseq.sh writes it)";
    std::string line, seed;
    while (seed.size() < 70 && std::getline(in, line)) {
        if (line.empty() || line[0] == '>') {
            if (!seed.empty())
                break;      // the first record only
            continue;
        }
        seed += line;
    }
    ASSERT_GE(seed.size(), 70u);
    seed.resize(70);
    // what the work-only walk and the memory-budgeted one must share: everything but the
    // physical work (timing) and the work units the budgets count differently
    std::function<void(Json::Value&)> walked = [&](Json::Value &v) {
        if (v.isObject()) {
            v.removeMember("timing");
            v.removeMember("work_units");
            for (const std::string &name : v.getMemberNames()) {
                walked(v[name]);
            }
        } else if (v.isArray()) {
            for (Json::Value &e : v) {
                walked(e);
            }
        }
    };
    auto walk = [&](size_t batch, bool memory) {
        Json::Value r;
        r["seeds"][0]["sequence"] = seed;
        Json::Value &st = r["strategy"];
        st["direction"] = "right";
        st["labels"]["mode"] = "constrain";
        st["bounds"]["max_extension_bp"] = 10000;
        if (memory) {
            st["bounds"]["max_memory_mb"] = 16;
        } else {
            st["bounds"]["max_work_units"] = Json::UInt64(1'000'000'000'000);
        }
        st["output"]["detail"] = "summary";
        st["output"]["timing"] = true;
        st["annotation"]["batch_kmers"] = Json::UInt64(batch);
        return mtg::cli::process_traverse_request(r, *anno_graph_, "")["results"][0];
    };
    for (size_t batch : { 2048, 4096, 8192 }) {
        Json::Value work = walk(batch, false);
        Json::Value memory = walk(batch, true);
        const uint64_t work_rows = work["timing"]["tuple_rows_fetched"].asUInt64();
        const uint64_t memory_rows = memory["timing"]["tuple_rows_fetched"].asUInt64();
        ASSERT_GT(work_rows, 20000u) << batch;
        EXPECT_FALSE(memory.isMember("resource_stop")) << batch;
        EXPECT_EQ(10000u, memory["arms"]["right"]["complete_to_bp"].asUInt64()) << batch;
        // at ea285c2e: 36,275, 36,838 and 38,441 against 27,642 (+31%, +33%, +39%)
        EXPECT_LE(memory_rows * 4, work_rows * 5)
            << "batch_kmers " << batch << ": " << memory_rows << " tuple rows under 16 MiB, "
            << work_rows << " under the work budget alone";
        walked(work);
        walked(memory);
        EXPECT_EQ(work["arms"], memory["arms"]) << batch;
    }
}

} // namespace
