#include <gtest/gtest.h>
#include <set>

#include "all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "graph/representation/canonical_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"

#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_match.hpp"
#include "graph/alignment_redone/aligner_config.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include "graph/alignment/dbg_aligner.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align_redone;
using namespace mtg::test;
using namespace mtg::kmer;


template <typename Graph>
class DBGAlignerRedoneTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGHashFast, DBGSuccinctCached, DBGSSHash> AlignGraphTypes;
TYPED_TEST_SUITE(DBGAlignerRedoneTest, AlignGraphTypes);

TEST(AlignmentQuery, ForwardAndReverseComplement) {
    DBGHashFast graph(3);
    Query empty(graph, "");
    EXPECT_TRUE(empty.get_query().empty());
    EXPECT_TRUE(empty.get_query(true).empty());

    Query short_query(graph, "AACG");
    EXPECT_EQ("AACG", short_query.get_query());
    EXPECT_EQ("CGTT", short_query.get_query(true));

    Query long_query(graph, std::string(1000, 'A'));
    EXPECT_EQ(std::string(1000, 'A'), long_query.get_query());
    EXPECT_EQ(std::string(1000, 'T'), long_query.get_query(true));
}


TEST(AlignmentCigarRedone, EmptyAndMalformedLengths) {
    Alignment empty;
    EXPECT_EQ(0u, empty.get_clipping());
    EXPECT_EQ(0u, empty.get_end_clipping());
    EXPECT_TRUE(empty.generate_cigar().is_valid({}, {}));
    EXPECT_EQ("2G3=", Cigar("2G3=").to_string());
    EXPECT_TRUE(Cigar("2G3=").is_valid("ACG", "ACG"));
    for (const auto &cigar : { "4294967295=", "4294967295X", "4294967295I",
                              "4294967295D", "4294967295S", "1=1I1D" }) {
        EXPECT_FALSE(Cigar(cigar).is_valid("A", "AC")) << cigar;
        EXPECT_FALSE(Cigar(cigar).is_valid({}, {})) << cigar;
    }
}

TEST(AlignmentCigarRedone, StrictParsingAndClipping) {
    EXPECT_TRUE(Cigar("").empty());
    for (const auto &text : { "3", "2=1", "1Q", "1a=", "1=xyz", "=1",
                              "-1=", "+1=", "1 =", "0=", "1=0I" }) {
        EXPECT_THROW(Cigar{ text }, std::invalid_argument) << text;
    }
    EXPECT_THROW(Cigar{ "4294967296=" }, std::out_of_range);
    EXPECT_THROW(Cigar{ "999999999999999999999999=" }, std::out_of_range);

    Cigar empty;
    empty.extend_clipping(0);
    EXPECT_TRUE(empty.empty());
    Cigar exact("2=");
    exact.extend_clipping(0);
    EXPECT_EQ("2=", exact.to_string());
    EXPECT_TRUE(exact.is_valid("AA", "AA"));
    EXPECT_FALSE(Cigar("1=1S1=").is_valid("AA", "AAA"));
    EXPECT_TRUE(Cigar("1S1=1S").is_valid("A", "AAA"));
}

TEST(AlignmentCigarRedone, CountersUseSizeType) {
    const Cigar cigar("2147483648=1I2D1X");
    EXPECT_EQ(size_t{2147483648}, cigar.get_num_matches());
    EXPECT_EQ(size_t{2147483650}, cigar.get_num_query());
    EXPECT_EQ(size_t{2147483651}, cigar.get_num_ref());
}

TEST(AlignmentScoringRedone, DefaultsAndInvalidBytes) {
    DBGAlignerConfig config;
    config.set_scoring_matrix();
    EXPECT_FALSE(config.alignment_edit_distance);
#if _PROTEIN_GRAPH
    EXPECT_EQ(4, config.match_score("A"));
#else
    EXPECT_EQ(2, config.match_score("A"));
#endif
    for (unsigned int byte = 128; byte < 256; ++byte) {
        const std::string invalid(1, static_cast<char>(byte));
        EXPECT_EQ(config.score_matrix[0]['A'], config.score_sequences(invalid, "A"));
        EXPECT_EQ(config.score_matrix['A'][0], config.score_sequences("A", invalid));
        DBGHashFast graph(3);
        EXPECT_THROW(Query(graph, invalid), std::invalid_argument);
    }
    EXPECT_THROW(config.score_sequences("AA", "A"), std::invalid_argument);
#if _DNA_CASE_SENSITIVE_GRAPH
    EXPECT_EQ(8, config.match_score("acgt"));
    EXPECT_EQ(-3, config.score_sequences("a", "g"));
    EXPECT_EQ(-3, config.score_sequences("A", "a"));
    config.alignment_edit_distance = true;
    config.set_scoring_matrix();
    EXPECT_EQ(5, config.match_score("acgtN"));
    EXPECT_EQ(-1, config.score_sequences("A", "a"));
#endif
}

#if _DNA5_GRAPH
TEST(LegacyQueryDNA5, AmbiguousMatches) {
    const std::string reference = "AACNNNACGNNNTTA";
    auto graph = build_graph_batch<DBGHashFast>(3, { reference });
    for (bool edit_distance : { false, true }) {
        mtg::graph::align::DBGAlignerConfig config{};
        config.alignment_edit_distance = edit_distance;
        config.alignment_match_score = 2;
        config.alignment_mm_transition_score = config.alignment_mm_transversion_score = 3;
        config.min_seed_length = 3;
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.forward_and_reverse_complement = false;
        config.seed_complexity_filter = false;
        config.set_scoring_matrix();
        const auto match_score = edit_distance ? 1 : 2;
        EXPECT_EQ(match_score, config.match_score("N"));
        EXPECT_EQ(mtg::graph::align::Cigar::MATCH, mtg::graph::align::kCharToOp['N']['N']);
        mtg::graph::align::DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(reference);
        ASSERT_EQ(1u, paths.size());
        EXPECT_EQ("15=", paths[0].get_cigar().to_string());
        EXPECT_EQ(15 * match_score, paths[0].get_score());
        EXPECT_TRUE(paths[0].get_cigar().is_valid(reference, reference));
    }
}
#endif

#if ! _PROTEIN_GRAPH
TEST(LegacyAlignment, ReverseOrientationWithoutTwoSidedExtension) {
    const std::string reference = "GTAGTGCTAGCTTGTAGTCGTGCTGATGCA";
    auto graph = build_graph_batch<DBGHashFast>(11, {reference});
    std::string reverse = reference;
    reverse_complement(reverse.begin(), reverse.end());
    mtg::graph::align::DBGAlignerConfig config{};
    config.min_seed_length = 11;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.two_sided = false;
    config.seed_complexity_filter = false;
    config.score_matrix = mtg::graph::align::DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    auto results = mtg::graph::align::DBGAligner<>(*graph, config).align(reverse);
    ASSERT_EQ(1u, results.size());
    EXPECT_TRUE(results[0].get_orientation());
    EXPECT_EQ("30=", results[0].get_cigar().to_string());
    EXPECT_TRUE(results[0].is_valid(*graph, &config));
    reverse[4] = reverse[4] == 'A' ? 'C' : 'A';
    results = mtg::graph::align::DBGAligner<>(*graph, config).align(reverse);
    ASSERT_EQ(1u, results.size());
    EXPECT_TRUE(results[0].get_orientation());
    EXPECT_EQ("25=1X4=", results[0].get_cigar().to_string());
    EXPECT_TRUE(results[0].is_valid(*graph, &config));
}
#endif

TEST(BossRanges, EachLabelIsCalledOnce) {
    auto graph = build_graph_batch<DBGSuccinct>(3, { "AAACAAAGAAATAACACCAGCAT" });
    const auto &boss = dynamic_cast<const DBGSuccinct&>(*graph).get_boss();
    bool tested_repeated_label = false;
    for (size_t left = 1; left < boss.num_edges(); ++left) {
        for (size_t length = 1; length < static_cast<size_t>(boss.alph_size - 1) && left + length - 1 <= boss.num_edges(); ++length) {
            const size_t right = left + length - 1;
            std::set<unsigned int> expected, actual;
            for (unsigned int symbol = 1; symbol < boss.alph_size; ++symbol) {
                boss::BOSS::edge_index start = left, end = right;
                if (boss.tighten_range(&start, &end, symbol)) {
                    expected.insert(symbol);
                    size_t occurrences = 0;
                    for (size_t i = left; i <= right; ++i)
                        occurrences += static_cast<unsigned int>(boss.get_W(i) % boss.alph_size) == symbol;
                    tested_repeated_label |= occurrences > 1;
                }
            }
            size_t callbacks = 0;
            boss.call_tighter_ranges(left, right, [&](auto, auto, auto symbol) {
                ++callbacks;
                actual.insert(symbol);
            });
            EXPECT_EQ(expected, actual);
            EXPECT_EQ(actual.size(), callbacks) << left << "," << right;
        }
    }
    EXPECT_TRUE(tested_repeated_label);
}

#if ! _PROTEIN_GRAPH
TEST(CanonicalBounds, InvalidNodesAreAbsent) {
    auto base = std::make_shared<DBGHashFast>(3, DeBruijnGraph::PRIMARY);
    base->add_sequence("AACGT");
    CanonicalDBG graph(base);
    EXPECT_FALSE(graph.in_graph(0));
    EXPECT_FALSE(graph.in_graph(graph.max_index() + 1));
    EXPECT_FALSE(graph.in_graph(std::numeric_limits<DeBruijnGraph::node_index>::max()));
    auto node = graph.kmer_to_node("AAC");
    ASSERT_NE(DeBruijnGraph::npos, node);
    EXPECT_TRUE(graph.in_graph(node));
}
#endif

void run_alignment(const DeBruijnGraph &graph,
                   DBGAlignerConfig config,
                   std::string_view query,
                   const std::vector<std::string> &reference,
                   const std::vector<std::string> &cigar_str,
                   size_t end_trim = 0,
                   bool needs_extension = false,
                   bool needs_extension_long_seed = false) {
    config.chain_alignments = true;
    config.two_sided = true;
    size_t k = graph.get_k();
    if (config.min_seed_length == 0)
        config.min_seed_length = k;

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        config.max_seed_length = std::max(mx, config.min_seed_length);
        bool check_chaining = mx == std::numeric_limits<size_t>::max()
            ? !needs_extension_long_seed
            : !needs_extension;

        Query aln_query(graph, query);
        ExactSeeder seeder(aln_query, config);
        std::vector<Alignment> paths;
        Extender extender(aln_query, config);
        std::vector<Alignment> paths_no_extend = seeder.get_inexact_anchors();
        std::vector<Alignment> paths_no_align = seeder.get_inexact_anchors(false);
        for (const auto &base_path : paths_no_extend) {
            extender.extend(base_path, [&](Alignment&& path) {
                paths.emplace_back(std::move(path));
            });
        }

        auto aln_sort = [&config](const auto &a, const auto &b) {
            return std::make_pair(score_match(a, config), b.get_orientation())
                 > std::make_pair(score_match(b, config), a.get_orientation());
        };

        std::sort(paths_no_extend.begin(), paths_no_extend.end(), aln_sort);
        std::sort(paths_no_align.begin(), paths_no_align.end(), aln_sort);
        std::sort(paths.begin(), paths.end(), aln_sort);

        ASSERT_LE(reference.size(), paths.size()) << mx;
        paths.resize(reference.size());

        for (size_t i = 0; i < reference.size(); ++i) {
            auto check_ref = [&](const Alignment &path, const std::string &reference, const std::string &type) {
                if (reference.size()) {
                    EXPECT_EQ(reference.size() - k + 1 + end_trim, path.get_path().size()) << mx << "\t" << type;
                    EXPECT_EQ(reference, path.get_spelling()) << mx << "\t" << type;
                }
            };

            auto check_aln = [&](const Alignment &path, const Cigar &cigar, const std::string &reference, const std::string &type) {
                ASSERT_EQ(cigar.to_string(), path.get_cigar().to_string()) << mx << "\t" << type;
                EXPECT_EQ(end_trim, path.get_end_trim()) << mx << "\t" << type;
                if (reference.size()) {
                    ASSERT_TRUE(cigar.is_valid(reference, path.get_query()));
                    EXPECT_EQ(config.score_cigar(reference, path.get_query(), cigar),
                              score_match(path, config)) << mx << "\t" << type;
                }

                EXPECT_EQ(cigar.get_clipping(), path.get_clipping()) << mx << "\t" << type;
                EXPECT_EQ(cigar.get_end_clipping(), path.get_end_clipping()) << mx << "\t" << type;
            };

            std::string cur_reference = reference[i];
            if (paths[i].get_orientation() && !paths[i].get_end_trim())
                ::reverse_complement(cur_reference.begin(), cur_reference.end());

            check_ref(paths[i], cur_reference, "extend");

            if (i < cigar_str.size() && cigar_str[i].size()) {
                std::string cur_cigar = cigar_str[i];
                Cigar cigar(cur_cigar);
                if (paths[i].get_orientation() && !paths[i].get_end_trim()) {
                    std::reverse(cigar.data().begin(), cigar.data().end());
                    cur_cigar = cigar.to_string();
                }
                check_aln(paths[i], cigar, cur_reference, "extend");

                if (check_chaining) {
                    // this alignment should work with chaining alone
                    ASSERT_LT(i, paths_no_extend.size());
                    check_ref(paths_no_extend[i], cur_reference, "chain");
                    check_aln(paths_no_extend[i], cigar, cur_reference, "chain");

                    ASSERT_LT(i, paths_no_align.size());
                    EXPECT_EQ(cigar.to_string(), paths_no_align[i].get_cigar().to_string()) << mx << "\tnoalign";
                }
            }
        }
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_empty) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, NonuniformSubstitutionScores) {
    const std::string reference = "AACGTTGCAATGCCGATTCAGGCTAAC";
    auto graph = build_graph_batch<TypeParam>(5, { reference });
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -3);
    config.score_matrix['C']['C'] = 4;
    config.score_matrix['G']['G'] = 6;
    config.score_matrix['T']['T'] = 8;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.forward_and_reverse_complement = false;
    run_alignment(*graph, config, "AACGTTGC" "G" "ATGCCGATTCAGGCTAAC",
                  { reference }, { "8=1X18=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, LongBranchingExtension) {
    uint32_t state = 7;
    std::string reference;
    for (size_t i = 0; i < 400; ++i) {
        state = state * 1664525u + 1013904223u;
        reference += "ACGT"[state >> 30];
    }
    auto change = [](char c) { return c == 'A' ? 'C' : 'A'; };
    std::string branch = reference;
    branch.at(100) = change(branch.at(100));
    auto graph = build_graph_batch<TypeParam>(17, { reference, branch });
    std::string query = reference;
    for (size_t i : { 8, 199, 391 })
        query.at(i) = change(query.at(i));
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.left_end_bonus = config.right_end_bonus = 5;
    config.forward_and_reverse_complement = false;
    run_alignment(*graph, config, query, { reference },
                  { "8=1X190=1X191=1X8=" }, 0, true, true);
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, ExactShortcutComparesBothOrientations) {
    auto graph = build_graph_batch<TypeParam>(3, { "AAACCC", "GGGTTT" });
    Query query(*graph, "AAACCC");
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -3);
    config.score_matrix['C']['C'] = 4;
    config.score_matrix['G']['G'] = 6;
    config.score_matrix['T']['T'] = 8;
    config.min_seed_length = 3;
    ExactSeeder seeder(query, config);
    Extender extender(query, config);
    size_t count = 0;
    align_query(query, seeder, extender, [&](Alignment&& aln) {
        ++count;
        EXPECT_TRUE(aln.get_orientation());
        EXPECT_EQ(42, score_match(aln, config));
        EXPECT_EQ("6=", aln.get_cigar().to_string());
    });
    EXPECT_EQ(1u, count);
}
#endif

#if _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, ProteinBlosumExtension) {
    const std::string reference = "MKWVTFISLLFLFSSAYS";
    auto graph = build_graph_batch<TypeParam>(5, { reference });
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::score_matrix_blosum62;
    config.max_nodes_per_seq_char = 5;
    config.max_ram_per_alignment = 200;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.forward_and_reverse_complement = false;
    run_alignment(*graph, config, "MKWVTFISLL" "Y" "LFSSAYS",
                  { reference }, { "10=1X7=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, ProteinHighDegreeExtension) {
    const std::string reference = "MKWVTFISLLFLFSSAYS";
    std::vector<std::string> references{ reference };
    for (char c : std::string("ACGTW")) {
        auto branch = reference;
        branch[5] = c;
        references.push_back(branch);
    }
    auto graph = build_graph_batch<TypeParam>(5, references);
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::score_matrix_blosum62;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.forward_and_reverse_complement = false;
    run_alignment(*graph, config, "MKWVTFISLL" "Y" "LFSSAYS",
                  { reference }, { "10=1X7=" }, 0, true, true);
}
#endif

#if _DNA5_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, Dna5RetainsAmbiguousReferenceBases) {
    if constexpr(!std::is_same_v<DBGSSHash, TypeParam>) {
        const std::string reference = "AACGTNACGTNACCT";
        auto graph = build_graph_batch<TypeParam>(5, { reference });
        DBGAlignerConfig config{};
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -3);
        config.left_end_bonus = config.right_end_bonus = 5;
        config.forward_and_reverse_complement = false;
        run_alignment(*graph, config, reference, { reference }, { "15=" });
    }
}
#endif

TYPED_TEST(DBGAlignerRedoneTest, PipelineExtendsIsolatedSeedsAndLimitsResults) {
    const std::string reference = "TAGCTTCGAGGCCAA";
    auto graph = build_graph_batch<TypeParam>(5, { reference });
    // The five-base seed before N must survive selection, and extension must
    // stop exact traversal at N before continuing with a mismatch.
    Query query(*graph, "TAGCTNCGAGGCCAA");
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = 5;
    config.max_seed_length = 5;
    config.forward_and_reverse_complement = false;
    for (size_t limit : { 0, 1, 2 }) {
        config.num_alternative_paths = limit;
        ExactSeeder seeder(query, config);
        Extender extender(query, config);
        std::vector<Alignment> results;
        align_query(query, seeder, extender, [&](Alignment&& aln) {
            results.push_back(std::move(aln));
        });
        ASSERT_LE(results.size(), limit);
        if (!limit) continue;
        ASSERT_FALSE(results.empty());
        EXPECT_EQ(reference, results[0].get_spelling());
        EXPECT_EQ("5=1X9=", results[0].get_cigar().to_string());
        for (const auto &aln : results) EXPECT_FALSE(aln.get_orientation());
        if (results.size() == 2) {
            EXPECT_NE(results[0].get_cigar().to_string(), results[1].get_cigar().to_string());
        }
    }
    config.min_exact_match = 1.0;
    ExactSeeder seeder(query, config);
    Extender extender(query, config);
    align_query(query, seeder, extender, [&](Alignment&&) {
        ADD_FAILURE() << "A mismatching alignment cannot satisfy min_exact_match=1";
    });
}

TYPED_TEST(DBGAlignerRedoneTest, SearchLimitsRetainValidPartialAlignments) {
    const std::string reference = "TAGCTTCGAGGCCAA";
    auto graph = build_graph_batch<TypeParam>(5, { reference });
    Query query(*graph, "TAGCTNCGAGGCCAA");
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = config.max_seed_length = 5;
    config.forward_and_reverse_complement = false;
    for (bool two_sided : {false, true}) {
        SCOPED_TRACE(two_sided);
        config.two_sided = two_sided;
        for (double nodes : { 0.0, 0.1, 1.0, 1000.0 }) {
            for (double ram : { 0.0, 0.0001, 0.001, 0.01, 1.0 }) {
                SCOPED_TRACE(nodes);
                SCOPED_TRACE(ram);
                config.max_nodes_per_seq_char = nodes;
                config.max_ram_per_alignment = ram;
                ExactSeeder seeder(query, config);
                Extender extender(query, config);
                std::vector<Alignment> results;
                align_query(query, seeder, extender, [&](Alignment&& aln) {
                    ASSERT_TRUE(aln.get_cigar().is_valid(aln.get_spelling(), aln.get_query()));
                    results.push_back(std::move(aln));
                });
                ASSERT_FALSE(results.empty());
                if (nodes == 1000 && ram == 1) {
                    EXPECT_EQ("5=1X9=", results[0].get_cigar().to_string());
                } else if (!nodes || !ram) {
                    EXPECT_LT(results[0].get_seed().size(), query.get_query().size());
                }
            }
        }
    }
}

TYPED_TEST(DBGAlignerRedoneTest, SeedLimitRejectsAmbiguousLoci) {
    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        size_t ambiguous = 0;
        for (const std::string prefix : { "AAAAAA", "AACGTT", "CGTACG", "GACTGA" }) {
            auto graph = build_graph_batch<TypeParam>(7,
                { prefix + "A", prefix + "C", prefix + "G", prefix + "T" });
            Query query(*graph, prefix);
            DBGAlignerConfig config{};
            config.min_seed_length = 6;
            config.max_seed_length = 7;
            config.forward_and_reverse_complement = false;
            auto unlimited = ExactSeeder(query, config).get_anchors();
            // SSHash can seed a sub-k query only when it contains a stored
            // minimizer; exercise the loci that have multiple candidates.
            if (unlimited.size() <= 1) continue;
            ++ambiguous;
            config.max_num_seeds_per_locus = 1;
            EXPECT_TRUE(ExactSeeder(query, config).get_anchors().empty());
            config.max_num_seeds_per_locus = unlimited.size();
            EXPECT_EQ(unlimited.size(), ExactSeeder(query, config).get_anchors().size());
        }
        EXPECT_GT(ambiguous, 0u);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_sequence_much_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, align_big_self_loop) {
    size_t k = 3;
    std::string reference = "AAAA";
    std::string query =     "AAAAAAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "3=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "14=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight_with_N) {
    size_t k = 5;
    std::string reference = "TAGCTTCGAGGCCAA";
    std::string query     = "TAGCTNCGAGGCCAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

    run_alignment(*graph, config, query, { reference }, { "5=1X9=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X9=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    //                                 X

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_branch) {
    size_t k = 6;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "17=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_branch_with_cycle) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "17=" });
}

TYPED_TEST(DBGAlignerRedoneTest, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =     "AGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "6=" });
}

TYPED_TEST(DBGAlignerRedoneTest, variation) {
    size_t k = 5;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X6=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X6=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAACTCGAAA";
    std::string reference_2 = "TTAAGCAAGTCGAAA";
    std::string query =       "TTAAGCAATGGGAAA";
    //                                 XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -2;

    run_alignment(*graph, config, query, { "" }, { "8=3X4=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { "" }, { "8=3X4=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, multiple_variations) {
    size_t k = 6;
    std::string reference = "TACGCAACTCTCTGAACTTGT";
    std::string query =     "TACGCAATTCTCTGTATTTGT";
    //                              X      X X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "7=1X6=1X1=1X4=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 4;
        run_alignment(*graph, config, query, { reference }, { "7=1X6=1X1=1X4=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_noise_in_branching_point) {
    size_t k = 5;
    std::string reference_1 = "AAAAACTTTTTT";
    std::string reference_2 = "AAAAATTGGGGG";
    std::string query =       "AAAAATTTTTTT";
    //                              D

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);

    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { "AAAAACTTTTTTT" }, { "5=1D7=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { "AAAAACTTTTTTT" }, { "5=1D7=" });
    }

    // TODO: too many possibilities here
    // config.gap_opening_penalty = -5;
    // config.gap_extension_penalty = -2;
    // if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
    //     config.min_seed_length = 2;
    //     run_alignment(*graph, config, query, { "AAAACTTTTTT" }, { "4=1X6=" });
    // } else {
    //     run_alignment(*graph, config, query, { "AAAACTTTTTT" }, { "4=1X6=" }, 0, true, true);
    // }
}

TYPED_TEST(DBGAlignerRedoneTest, alternative_path_basic) {
    size_t k = 5;
    std::vector<std::string> references = {
                        "AGACAATTTTTTAAA",
                        "AGACAATTTTTGAAA",
                        "AGACAAGTTTTTAAA",
                        "AGACAAGTTTTGAAA" };
    std::string query = "AGACAACTTTTCAAA";
    //                         X    X

    auto graph = build_graph_batch<TypeParam>(k, references);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

    run_alignment(*graph, config, query, { "" }, { "6=1X4=1X3=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_multiple_misalignment) {
    size_t k = 5;
    std::string reference = "AAAAGCGGACCCTTTCCGTTAT";
    std::string query =     "AAAAGGGGACCCTTTTCGTTAT";
    //                            X         X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X9=1X6=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X9=1X6=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCATTGTT";
    //                            I

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "5=1I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=1I5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_multi) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAATTGTT";
    //                            II

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "5=2I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=2I5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAAAAAAAAATTGTT";
    //                            IIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference }, { "5=9I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=9I5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long_offset) {
    size_t k = 5;
    std::string reference = "TTTCCGGTTGTTA";
    std::string query =     "TTTCCGCAAAAAAAAATTGTTA";
    //                             XIIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference }, { "6=1X9I6=" }, 0, true, true);
    // if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
    //     config.min_seed_length = 3;
    //     run_alignment(*graph, config, query, { reference }, { "6=1X9I6=" });
    // } else {
    // }
}

TYPED_TEST(DBGAlignerRedoneTest, align_delete) {
    size_t k = 5;
    //                               D
    std::string reference = "TTCGAT""A""GGCCT";
    std::string query =     "TTCGAT"   "GGCCT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "6=1D5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "6=1D5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_gap) {
    size_t k = 10;
    //                                  DDDD
    std::string reference = "TTTTCTGTATACCTTGGCGCTCTCC";
    std::string query =     "TTTTCTGTATA"  "GGCGCTCTCC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "11=4D10=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference }, { "11=4D10=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_gap_after_seed) {
    size_t k = 4;
    //                           DDDD
    std::string reference = "TTTCCCAAGGCGCTCTC";
    std::string query =     "TTTC"  "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -2;
    config.left_end_bonus = 4;
    config.right_end_bonus = 4;
    run_alignment(*graph, config, query, { reference }, { "4=4D9=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_loop_deletion) {
    size_t k = 4;
    std::string reference = "AAAATTTTCGAGGCCAA";
    std::string query =     "AAAACGAGGCCAA";
    //                           DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::unit_scoring_matrix(
        2, alphabets::kAlphabetDNA, alphabets::kCharToDNA
    );
    config.gap_opening_penalty = -2;
    config.gap_extension_penalty = -2;
    run_alignment(*graph, config, query, { "AAAATTTCGAGGCCAA" }, { "4=3D9=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight_long_xdrop) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAGGCCAAGCCTGACTGATCGATGCATGCTAGCTAGTCAGTCAGCGTGAGCTAGCAT";
    std::string reference_2 = "AGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    std::string query = reference_1;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 30;
    config.rel_score_cutoff = 0.8;
    run_alignment(*graph, config, query, { query }, { "63=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_drop_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTGGCGCTCTC";
    std::string query =     "TTTCCGGGGCGCTCTC";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -6, -6);
    config.gap_opening_penalty = -10;
    config.gap_extension_penalty = -4;
    config.xdrop = 6;
    run_alignment(*graph, config, query, { reference.substr(7) }, { "7S9=" });
}

// TODO: double-check
// TYPED_TEST(DBGAlignerRedoneTest, align_long_gap_after_seed) {
//     size_t k = 4;
//     std::string reference = "TTTCCCTTAAGGCGCTCTC";
//     std::string query =     "TTTC"    "GGCGCTCTC";
//     //                       SSSS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -1;
//     run_alignment(*graph, config, query, { reference.substr(10) }, { "4S9=" });
// }

TYPED_TEST(DBGAlignerRedoneTest, align_repeat_sequence_no_delete_after_insert) {
    size_t k = 27;
    std::string reference = "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                                                    IIIIIII        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "45=7I8=1X39=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 25;
        run_alignment(*graph, config, query, { reference }, { "45=7I8=1X39=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_clipping1) {
    size_t k = 4;
    std::string reference = "GGCCTGTTTG";
    std::string query =     "ACCCTGTTTG";
    //                       SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(2) }, { "2S8=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_clipping2) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAA";
    std::string query =      "TTAGCTTCGAGGCCAA";
    //                        SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(3) }, { "2S14=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_long_clipping) {
    size_t k = 4;
    std::string reference = "TTTTTTTAAAAGCTTCGAGGCCAA";
    std::string query =     "CCCCCCCAAAAGCTTCGAGGCCAA";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(7) }, { "7S17=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_end_clipping) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAATTTTTTT";
    std::string query =     "AAAAGCTTCGAGGCCAACCCCCCC";
    //                                        SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(0, 17) }, { "17=7S" });
}

// TODO: implement min cell score
// TYPED_TEST(DBGAlignerTest, align_clipping_min_cell_score) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =        "ACCTTTCGAGGCCAA";
//     //                          SS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(7u, path.size());
//     EXPECT_EQ(reference.substr(5), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
//     EXPECT_EQ("2S13=", path.get_cigar().to_string());
//     EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(2u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity) {
    size_t k = 27;
    std::string reference = "CTAGAACTTAAAGTATAATAATACTAATAATAAAATAAAATACA";
    std::string query =     "CTAGAACTTAAAGTATAATAATACTAATAA""AAGTACAATACA";
    //                                                     DD  X  X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, { reference }, { "30=2D2=1X2=1X6=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity2) {
    size_t k = 27;
    std::string reference = "GCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                         IIIIIII        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);

    run_alignment(*graph, config, query, { reference.substr(3) }, { "10S4=1X9=1X8=1X39=" }, 0, true, true);

    config.gap_opening_penalty = -3;
    run_alignment(*graph, config, query, { reference },           {       "18=7I8=1X39=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity2_del) {
    size_t k = 27;
    std::string query =         "GCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string reference =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                             DDDDDDD        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);

    config.gap_opening_penalty = -5;
    run_alignment(*graph, config, query, { reference.substr(10) }, { "3S4=1X9=1X8=1X39=" }, 0, true, true);

    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference },           {       "18=7D8=1X39=" }, 0, true, true);
}

// TODO: implement seed complexity filter
// TYPED_TEST(DBGAlignerTest, align_low_similarity3) {
//     size_t k = 27;
//     std::string reference = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCTGGGATTATAGGTGTGAACCACCACACCTGGCTAATTTTTTTTGTGTGTGTGTGTGTTTTTTC";
//     std::string query =     "AAAAAAAAAAAAAAAAAAAAAAAAAAACGCCAAAAAGGGGGAATAGGGGGGGGGGAACCCCAACACCGGTATGTTTTTTTGTGTGTGGGGGATTTTTTTC";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     for (bool seed_complexity_filter : { false, true }) {
//         DBGAlignerConfig config;
//         config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//         config.seed_complexity_filter = seed_complexity_filter;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
// #if ! _PROTEIN_GRAPH
//         EXPECT_EQ(seed_complexity_filter, paths.empty());
// #else
//         EXPECT_FALSE(paths.empty());
// #endif
//     }
// }

// TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4) {
//     size_t k = 6;
//     std::vector<std::string> seqs;
//     mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
//                                           [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
//     auto graph = build_graph_batch<TypeParam>(k, std::move(seqs));

//     std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGA";

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.left_end_bonus = 5;
//     config.right_end_bonus = 5;
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     // config.min_seed_length = 6;
//     run_alignment(*graph, config, query, { "" }, { "" });
// }

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity5) {
    size_t k = 31;
    std::string reference = "GTCGTCAGATCGGAAGAGCGTCGTGTAGGGAAAG" "GTCTTC""GCCTGTGTAGATCTC";
    std::string query =        "GTCAGATCGGAAGAGCGTCGTGTAGGGAAAG""AGTGTTCCTGGTGGTGTAGATCTC";
    //                                                           I  X   II XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 11;
        run_alignment(*graph, config, query, { reference.substr(3, 52) }, { "31=1I2=1X3=2I1=3X11=" }, 0, true, true);
    }
}

// TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =        "ACCTTTCGAGGCCAA";
//     //                          SS

//     auto graph = std::make_shared<DBGSuccinct>(k);
//     graph->add_sequence(reference);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.min_seed_length = 2;
//     config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
//     config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(7u, path.size());
//     EXPECT_EQ(reference.substr(5), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
//     EXPECT_EQ("2S13=", path.get_cigar().to_string());
//     EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(2u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

#if ! _PROTEIN_GRAPH
// TEST(DBGAlignerTest, align_suffix_seed_snp_canonical) {
//     size_t k = 18;
//     std::string reference = "AAAAACTTTCGAGGCCAA";
//     std::string query =     "GGGGGCTTTCGAGGCCAA";
//     //                       SSSSS

//     std::string reference_rc = "TTGGCCTCGAAAGTTTTT";
//     std::string query_rc     = "TTGGCCTCGAAAGCCCCC";

//     for (auto mode : { DeBruijnGraph::PRIMARY, DeBruijnGraph::CANONICAL }) {
//         auto dbg_succ = std::make_shared<DBGSuccinct>(k, mode);
//         dbg_succ->add_sequence(reference_rc);

//         std::shared_ptr<DeBruijnGraph> graph = dbg_succ;
//         if (mode == DeBruijnGraph::PRIMARY)
//             graph = std::make_shared<CanonicalDBG>(graph);

//         DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//         config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
//         config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//         config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//         config.min_seed_length = 13;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
//         ASSERT_EQ(1ull, paths.size());
//         auto path = paths[0];

//         EXPECT_EQ(1u, path.size()); // includes dummy k-mers
//         if (path.get_sequence() == reference.substr(5)) {
//             EXPECT_EQ(config.score_sequences(query.substr(5), reference.substr(5)),
//                       path.get_score());
//             EXPECT_EQ("5S13=", path.get_cigar().to_string());
//             EXPECT_EQ(5u, path.get_clipping());
//             EXPECT_EQ(0u, path.get_end_clipping());
//         } else {
//             ASSERT_EQ(reference_rc.substr(0, 13), path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query_rc.substr(0, 13),
//                                              reference_rc.substr(0, 13)),
//                       path.get_score());
//             EXPECT_EQ("13=5S", path.get_cigar().to_string());
//             EXPECT_EQ(0u, path.get_clipping());
//             EXPECT_EQ(5u, path.get_end_clipping());
//         }
//         EXPECT_EQ(5u, path.get_offset());
//         EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//         EXPECT_FALSE(is_exact_match(path));
//         EXPECT_TRUE(path.is_valid(*graph, &config));
//         check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//         // TODO: sub-k seeds which are sink tips in the underlying graph currently can't be found
//         if (mode != DeBruijnGraph::PRIMARY)
//             check_extend(graph, aligner.get_config(), paths, query);
//     }
// }

// TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4_rep_primary) {
//     size_t k = 6;
//     std::vector<std::string> seqs;
//     mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
//                                           [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
//     auto graph = build_graph_batch<TypeParam>(k, std::move(seqs), DeBruijnGraph::PRIMARY);

//     std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGA";

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.left_end_bonus = 5;
//     config.right_end_bonus = 5;
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     // config.min_seed_length = 6;
//     run_alignment(*graph, config, query, { "" }, { "" });
// }
#endif

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, align_both_directions) {
    size_t k = 7;
    std::string reference =    "AAAAGCTTTCGAGGCCAA";
    std::string query =        "AAAAGTTTTCGAGGCCAA";
    //                               X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::CANONICAL);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X12=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 5;
        if constexpr(std::is_same_v<DBGSSHash, TypeParam>) {
            const auto &sshash = static_cast<const DBGSSHash&>(*graph);
            std::visit([&](const auto &dict) {
                using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
                auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(reference.data(), k);
                auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                // we can only find the first seed of it's also the minimizer
                if (sshash::util::compute_minimizer<kmer_t>(first_kmer_ref, k, dict.m(), dict.hasher()).pos_in_seq
                        || sshash::util::compute_minimizer<kmer_t>(first_kmer_query, k, dict.m(), dict.hasher()).pos_in_seq) {
                    run_alignment(*graph, config, query, { reference }, { "5=1X12=" }, 0, true, true);
                } else {
                    run_alignment(*graph, config, query, { reference }, { "5=1X12=" });
                }
            }, sshash.data());
        } else {
            run_alignment(*graph, config, query, { reference }, { "5=1X12=" });
        }
    }
}

#endif

TEST(AlignmentTraceback, RetainsReplacedPredecessors) {
    // Two branches with repeats replace a DP predecessor after a deletion
    // successor has been queued. Traceback must use the displaced distance.
    const std::string first =
        "AGCCACCACGCCTGGTCCAAGAACCCAACTTTTAGATCTAGAGTGATGTCAGCATGACATTGATTTCCTGAGGCCCAGGGGTGAAGGAGC"
        "TGAGGACAGCAGAGGGGTGAAGGAAGTCAGCTACAGACAGCAGCAGCTGATGCACAGGCCTCCCAGTGCCTGAAGTCACCCGGAATTGGG"
        "AAGTGC";
    const std::string second =
        "GGCCTCCCAAAGTGCTGGGATTGCAGGCGTGAGACAC";
    const std::string sequence =
        "GGCCTCCCAAAGTGCTGGGATTCCACGTGGAGCCACCACGCCTGGTCCAAGAACCCAACTTTAAAATCAAGAGTGATGTCAGCATGACAT"
        "TGATTTCCTGAGGCCCAGGGGTGAAGGAGCTGAGGACAGCAGAGGGGTGATAGGAAGTCAGTTACGACCAGAGCAGCTGCCCACAGGCCC"
        "TCCCATGTGCCTGAAGTCACCCGGA";
    auto graph = build_graph_batch<DBGSuccinctCached>(31, {first, second}, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    config.min_seed_length = 19;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.two_sided = true;
    config.xdrop = 27;
    config.min_exact_match = 0.7;
    config.gap_opening_penalty = -6;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    Query query(*graph, sequence);
    ExactSeeder seeder(query, config);
    Extender extender(query, config);
    size_t count = 0;
    align_query(query, seeder, extender, [&](Alignment&& a) {
        ++count;
        EXPECT_TRUE(a.is_spelling_valid(*graph));
        EXPECT_TRUE(a.get_cigar().is_valid(a.get_spelling(), a.get_query()));
        EXPECT_GE(score_match(a, config), 275);
    });
    EXPECT_EQ(1u, count);
}

TYPED_TEST(DBGAlignerRedoneTest, ChainingAndTwoSidedSearchAreIndependentOptions) {
    const std::string reference = "GTAGTGCTAGCTTGTAGTCGTGCTGATGCA";
    const std::string sequence = "GTAGTGCTACCTTGTAGTCGTGGTGATGCA";
    auto graph = build_graph_batch<TypeParam>(11, {reference}, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    EXPECT_FALSE(config.chain_alignments);
    EXPECT_FALSE(config.two_sided);
    config.min_seed_length = config.max_seed_length = 11;
    config.forward_and_reverse_complement = false;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.xdrop = 27;
    for (bool chain : {false, true}) {
        config.chain_alignments = chain;
        for (bool two_sided : {false, true}) {
            config.two_sided = two_sided;
            Query query(*graph, sequence);
            ExactSeeder seeder(query, config);
            Extender extender(query, config);
            std::vector<Alignment> results;
            align_query(query, seeder, extender, [&](Alignment&& a) {
                EXPECT_FALSE(a.get_orientation());
                EXPECT_TRUE(a.is_spelling_valid(*graph));
                EXPECT_TRUE(a.get_cigar().is_valid(a.get_spelling(), a.get_query()));
                results.push_back(std::move(a));
            });
            ASSERT_EQ(1u, results.size()) << chain << " " << two_sided;
            if (two_sided) {
                EXPECT_EQ("9=1X12=1X7=", results[0].get_cigar().to_string());
            } else {
                EXPECT_GT(results[0].get_clipping(), 0u);
            }
        }
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_both_directions2) {
    size_t k = 11;
    std::string reference =    "GTAGTGCTAGCTTGTAGTCGTGCTGATGCA";
    std::string query =        "GTAGTGCTACCTTGTAGTCGTGGTGATGCA";
    //                                   X            X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(0, query.size()) }, { "9=1X12=1X7=" }, 0, true, true);
}

// TYPED_TEST(DBGAlignerTest, align_nodummy) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =     "AAAAGTTTTCGAGGCCAA";
//     //                            X

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

//     for (bool both_directions : { false, true }) {
// #if _PROTEIN_GRAPH
//         if (both_directions)
//             continue;
// #endif
//         config.forward_and_reverse_complement = both_directions;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
//         ASSERT_EQ(1ull, paths.size());
//         auto path = paths[0];

//         if (both_directions) {
//             EXPECT_EQ(12u, path.size());
//             EXPECT_EQ(reference, path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
//             EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
//             EXPECT_EQ(17u, path.get_cigar().get_num_matches());
//         } else {
//             EXPECT_EQ(6u, path.size());
//             EXPECT_EQ(reference.substr(6), path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query.substr(6), reference.substr(6)), path.get_score());
//             EXPECT_EQ("6S12=", path.get_cigar().to_string());
//             EXPECT_EQ(12u, path.get_cigar().get_num_matches());
//             EXPECT_EQ(6u, path.get_clipping());
//         }
//         EXPECT_FALSE(is_exact_match(path));
//         EXPECT_EQ(0u, path.get_end_clipping());
//         EXPECT_EQ(0u, path.get_offset());
//         EXPECT_TRUE(path.is_valid(*graph, &config));
//         check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//         check_extend(graph, aligner.get_config(), paths, query);
//     }
// }

TYPED_TEST(DBGAlignerRedoneTest, align_seed_to_end) {
    size_t k = 5;
    //                       SSSSSSSSSSSSSSSSSSSSSS
    std::string reference =                  "ATCCCTTTTAAAA";
    std::string query =     "ATCCCGGGGGGGGGGGGGGGGGTTTTAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(5) }, { "22S8=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_bfs_vs_dfs_xdrop) {
    size_t k = 31;
    std::string reference_1 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGGATCAGAGAAGAAAGGAATGGCAAATACTCCATAGAAGAGTACCTTCAGGACAGGGACAGATACTATGAGGAGGTGGCCAT";
    std::string reference_2 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGAATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";
    std::string query =       "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGATGAGTACGATGAGTACAAGAGCATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 27;
    run_alignment(*graph, config, query, { "" }, { "48=1X20=1X80=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 19;
        run_alignment(*graph, config, query, { "" }, { "48=1X20=1X80=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy_short) {
    size_t k = 7;
    std::string reference = "GAAAAGCT";
    std::string query =     "GAAAAGTT";
    //                             X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference }, { "6=1X1=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy_short2) {
    size_t k = 31;
    std::string reference = "GCTAGTCTATAGCGTACGATCTAGCTAGTCT";
    std::string query =     "GCTAGTCTATAGCGTACGATCTTTTTAGTC";
    //                                             XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 19;
        run_alignment(*graph, config, query, { reference.substr(0, 30) }, { "22=3X5=" }, 1, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy) {
    size_t k = 7;
    std::string reference = "AAAAGCTTTCGATT";
    std::string query =     "AAAAGTTTTCGAA";
    //                            X      X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.left_end_bonus = 2;
    config.right_end_bonus = 2;
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference.substr(0, query.size()) }, { "5=1X6=1X" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_extended_insert_after_match) {
    size_t k = 27;
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG" "GCCCAGGCCCAGGC""CCA" "GGCCCAGGC""CC" "AGGCCCAGGCCCAGGCCCAAGCC";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG" "CCCCAGGCCCAGGC""CCA" "GGCCCAGGC""CC" "AGGCCCAGGCCCAGGCCCAAGCC";
    std::string query =       "CGTGGCCCAGGCCCAGGCCCAGTGGGCGTTGGCCCAGGCGGCCA""CGG" "TGGCTGCG""CAGGCCCGCCTGGCACAAGCCACGCTG";
    //                                               III  XXX         II     I  DDDX   II X  I

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 15;
        run_alignment(*graph, config, query,
                      { "CGTGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC" },
                      { "" }, 0, true, true);
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, align_suffix_seed_no_full_seeds) {
    size_t k = 31;
    std::string reference =     "CTGCTGCGCCATCGCAACCCACGGTTGCTTTTTGAGTCGCTGCTCACGTTAGCCATCACACTGACGTTAAGCTGGCTTTCGATGCTGTATCTTTTTTTT";
    std::string query     = "CTTACTGCTGCGCTCTTCGCAAACCCCACGGTTTCTTGTTTTGAGCTCGCCTGCTCACGATACCCATACACACTGACGTTCAAGCTGGCTTTCGATGTTGTATC";
    //                       SSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::PRIMARY);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, {}, {});
    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 9;
        run_alignment(*graph, config, query,
                      { reference.substr(0, 91) },
                      { "4S9=1I1=1X6=2I9=1X2=2I7=1I3=1I9=1X2=1X4=1I12=1I16=1X6=" }, 0, true, true);
    }
}
#endif

} // namespace
