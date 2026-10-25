#include <gtest/gtest.h>
#include <fstream>
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

    Query mixed_case(graph, "AaCg");
#if _DNA_CASE_SENSITIVE_GRAPH
    EXPECT_EQ("AaCg", mixed_case.get_query());
    EXPECT_EQ("cGtT", mixed_case.get_query(true));
#else
    EXPECT_EQ("AACG", mixed_case.get_query());
    EXPECT_EQ("CGTT", mixed_case.get_query(true));
#endif
}

#if ! _DNA_CASE_SENSITIVE_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, LowercaseQueriesMatchUppercaseInAllModes) {
#if _PROTEIN_GRAPH
    const std::string reference = "MKWVTFISLLFLFSSAYS";
    const std::string sequence = "MKWVTFISLLYLFSSAYS";
#else
    const std::string reference = "AACGTCAGTACCTGATGCATGCA";
    const std::string sequence = "AACGTCAGTACCTGATTCATGCA";
#endif
    auto graph = build_graph_batch<TypeParam>(5, { reference });
    std::string lowercase = sequence;
    std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    DBGAlignerConfig config{};
    config.set_scoring_matrix();
    config.min_seed_length = 5;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.xdrop = 27;
    using Result = std::tuple<int32_t, std::string, std::string, bool>;
    auto run = [&](const std::string &input) {
        Query query(*graph, input);
        EXPECT_EQ(sequence, query.get_query());
        ExactSeeder seeder(query, config);
        Extender extender(query, config);
        std::vector<Result> results;
        align_query(query, seeder, extender, [&](Alignment&& alignment) {
            EXPECT_TRUE(alignment.get_cigar().is_valid(alignment.get_spelling(), query.get_query()));
            results.emplace_back(score_match(alignment, config), alignment.get_spelling(),
                                 alignment.get_cigar().to_string(), alignment.get_orientation());
        });
        return results;
    };
    for (bool chain : { false, true }) {
        for (bool two_sided : { false, true }) {
            config.chain_alignments = chain;
            config.two_sided = two_sided;
            const auto expected = run(sequence);
            ASSERT_FALSE(expected.empty());
            EXPECT_EQ(expected, run(lowercase)) << chain << "," << two_sided;
        }
    }
}
#endif


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
    EXPECT_TRUE(mtg::graph::align::Cigar("").empty());
    for (const auto &text : { "3", "2=1", "1Q", "1a=", "1=xyz", "=1",
                              "-1=", "+1=", "1 =", "0=", "1=0I" }) {
        EXPECT_THROW(Cigar{ text }, std::invalid_argument) << text;
        EXPECT_THROW(mtg::graph::align::Cigar{ text }, std::invalid_argument) << text;
    }
    for (const auto *text : { "4294967296=", "999999999999999999999999=" }) {
        EXPECT_THROW(Cigar{ text }, std::out_of_range) << text;
        EXPECT_THROW(mtg::graph::align::Cigar{ text }, std::out_of_range) << text;
    }

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
    const char *text = "2147483648=1I2147483648=2D1X";
    const Cigar cigar(text);
    EXPECT_EQ(size_t{4294967296}, cigar.get_num_matches());
    EXPECT_EQ(size_t{4294967296}, mtg::graph::align::Cigar(text).get_num_matches());
    EXPECT_EQ(size_t{4294967298}, cigar.get_num_query());
    EXPECT_EQ(size_t{4294967299}, cigar.get_num_ref());
}

TEST(AlignmentCigar, ConversionPreservesStoredOperations) {
    // Conversion must not parse or normalize the stored representation.
    const std::vector<Cigar::value_type> operations{
        {Cigar::CLIPPED, 1}, {Cigar::MATCH, 2}, {Cigar::MATCH, 3},
        {Cigar::INSERTION, 0}, {Cigar::NODE_INSERTION, 2},
        {Cigar::DELETION, 1}, {Cigar::CLIPPED, 2}
    };
    const mtg::graph::align::Cigar legacy(operations);
    const Cigar converted(legacy.data());
    EXPECT_EQ(operations, converted.data());
    EXPECT_EQ(operations, mtg::graph::align::Cigar(converted.data()).data());
}

TEST(AlignmentCigar, ValidatorsKeepTheirQueryViewContracts) {
    const Cigar full_query("1S2=1S");
    const mtg::graph::align::Cigar aligned_slice(full_query.data());
    EXPECT_TRUE(full_query.is_valid("AC", "TACG"));
    EXPECT_TRUE(aligned_slice.is_valid("AC", "AC"));
    EXPECT_FALSE(full_query.is_valid("AC", "AC"));
    EXPECT_FALSE(aligned_slice.is_valid("AC", "TACG"));
}

TEST(AlignmentCigarRedone, PostChainGraphJumpsAllowInternalClipping) {
    EXPECT_TRUE(Cigar("1=1S1G1D1=").is_valid("A$C", "ATC"));
    EXPECT_TRUE(Cigar("1=1S1D1=").is_valid("A$C", "ATC"));
    EXPECT_TRUE(Cigar("1=1D1G1S1=").is_valid("G$T", "GAT"));
    EXPECT_TRUE(Cigar("1=1D1S1=").is_valid("G$T", "GAT"));
    EXPECT_TRUE(Cigar("1=1S3D1=").is_valid("AC$GT", "ATT"));
    EXPECT_TRUE(Cigar("1=3D1S1=").is_valid("AC$GT", "AAT"));
    EXPECT_FALSE(Cigar("1=1S1=").is_valid("AC", "ATC"));
    EXPECT_FALSE(Cigar("1=1S1D1=").is_valid("AGC", "ATC"));
    EXPECT_FALSE(Cigar("1=1D1S1=").is_valid("GGT", "GAT"));
    EXPECT_FALSE(Cigar("1=1S3D1=").is_valid("ACCGT", "ATT"));
    EXPECT_FALSE(Cigar("1=3D1S1=").is_valid("ACCGT", "AAT"));
    EXPECT_FALSE(Cigar("1=1S9D1=").is_valid("AC$GT", "ATT"));
    EXPECT_FALSE(Cigar("1=9D1S1=").is_valid("AC$GT", "AAT"));
    EXPECT_FALSE(Cigar("1=9S1G1D1=").is_valid("A$C", "ATC"));
    EXPECT_FALSE(Cigar("1=1D1G9S1=").is_valid("G$T", "GAT"));
}

TEST(AlignmentCigar, GraphJumpsRoundTrip) {
    for (const auto *text : { "10=1G13=", "2G3=1G4X", "3=2G" }) {
        EXPECT_EQ(text, Cigar(text).to_string());
        EXPECT_EQ(text, mtg::graph::align::Cigar(text).to_string());
    }
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

TEST(CanonicalSubKSeeds, IncludesEveryMatchingSuffixRangeNode) {
    for (bool mask_dummy : { false, true }) {
        SCOPED_TRACE(mask_dummy);
        auto base = std::make_shared<DBGSuccinct>(5, DeBruijnGraph::PRIMARY);
        base->add_sequence("AACGT");
        base->add_sequence("CCCGT");
        if (mask_dummy)
            base->mask_dummy_kmers(1, false);
        CanonicalDBG graph(base);
        Query query(graph, "ACG");
        DBGAlignerConfig config;
        config.min_seed_length = 3;
        config.max_seed_length = 5;
        config.forward_and_reverse_complement = false;
        config.seed_complexity_filter = false;
        std::set<std::string> spellings;
        for (const auto &anchor : ExactSeeder(query, config).get_anchors()) {
            EXPECT_FALSE(anchor.get_orientation());
            EXPECT_EQ("ACG", anchor.get_seed());
            EXPECT_TRUE(anchor.is_spelling_valid(graph));
            spellings.insert(anchor.get_path_spelling());
        }
        EXPECT_EQ((std::set<std::string>{ "ACGGG", "ACGTT" }), spellings);
    }
}

TEST(CanonicalAlternatives, MirrorsDoNotDisplaceDistinctPaths) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        SCOPED_TRACE(mode);
        DBGSuccinct graph(5, mode);
        for (const auto *sequence : { "AACGTAC", "AACGTTC", "GTACGTT", "GAACGTT" })
            graph.add_sequence(sequence);
        graph.mask_dummy_kmers(1, false);
        DBGAlignerConfig config;
        config.set_scoring_matrix();
        config.min_seed_length = config.max_seed_length = 5;
        config.num_alternative_paths = 2;
        config.rel_score_cutoff = 0;
        config.seed_complexity_filter = false;
        Query query(graph, "AACGTAC");
        ExactSeeder seeder(query, config);
        Extender extender(query, config);
        using Result = std::pair<bool, std::string>;
        std::set<Result> results;
        align_query(query, seeder, extender, [&](Alignment&& alignment) {
            results.emplace(alignment.get_orientation(), alignment.get_spelling());
        });
        const std::set<Result> expected = mode == DeBruijnGraph::BASIC
            ? std::set<Result>{ { false, "AACGTAC" }, { true, "GTACGTT" } }
            : std::set<Result>{ { false, "AACGTAC" }, { false, "AACGTTC" } };
        EXPECT_EQ(expected, results);
    }
}

TEST(CanonicalAlternatives, AsymmetricScoresKeepTheirOrientation) {
    for (bool end_bonuses : { false, true }) {
        SCOPED_TRACE(end_bonuses);
        DBGAlignerConfig config;
        config.set_scoring_matrix();
        config.min_seed_length = config.max_seed_length = 5;
        config.num_alternative_paths = 2;
        config.rel_score_cutoff = 0;
        config.seed_complexity_filter = false;
        if (end_bonuses) {
            config.left_end_bonus = 10;
            config.right_end_bonus = 0;
        } else {
            config.score_matrix['T']['T'] = 9;
        }
        using Result = std::tuple<bool, int32_t, std::string, std::string>;
        auto run = [&](DeBruijnGraph::Mode mode) {
            DBGSuccinct graph(5, mode);
            // The BASIC control stores both strands explicitly, with no mirror
            // normalization; its scoring and available paths are otherwise equal.
            graph.add_sequence("AACGTAC");
            graph.add_sequence("GTACGTT");
            graph.mask_dummy_kmers(1, false);
            Query query(graph, end_bonuses ? "AACGTACNNNNNNNNNNNNNNNNNNNN" : "AACGTAC");
            ExactSeeder seeder(query, config);
            Extender extender(query, config);
            std::set<Result> results;
            align_query(query, seeder, extender, [&](Alignment&& alignment) {
                results.emplace(alignment.get_orientation(), score_match(alignment, config),
                                alignment.get_cigar().to_string(), alignment.get_spelling());
            });
            return results;
        };
        const auto expected = run(DeBruijnGraph::BASIC);
        ASSERT_EQ(2u, expected.size());
        if (!end_bonuses) {
            ASSERT_TRUE(std::any_of(expected.begin(), expected.end(),
                                   [](const auto &result) { return std::get<0>(result); }));
        }
        EXPECT_EQ(expected, run(DeBruijnGraph::CANONICAL));
    }
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

        for (const auto *candidates : { &paths_no_extend, &paths_no_align, &paths }) {
            for (const auto &path : *candidates)
                ASSERT_TRUE(path.get_cigar().is_valid(path.get_spelling(), path.get_query()))
                    << path;
        }

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

TYPED_TEST(DBGAlignerRedoneTest, ZeroSeedLengthUsesGraphK) {
    const std::string reference = "CGTACGATGGCATTCAGGCTAACGATCGATCG";
    auto graph = build_graph_batch<TypeParam>(11, {reference}, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    config.min_seed_length = 0;
    config.max_seed_length = 11;
    config.seed_complexity_filter = false;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    Query short_query(*graph, reference.substr(0, 10));
    EXPECT_TRUE(ExactSeeder(short_query, config).get_anchors().empty());
    Query query(*graph, reference.substr(0, 15) + "A" + reference.substr(16));
    auto defaults = ExactSeeder(query, config).get_anchors();
    config.min_seed_length = graph->get_k();
    auto explicit_k = ExactSeeder(query, config).get_anchors();
    EXPECT_EQ(explicit_k, defaults);
    ASSERT_FALSE(defaults.empty());
    for (const auto &anchor : defaults) {
        EXPECT_EQ(graph->get_k(), anchor.get_seed().size());
        EXPECT_FALSE(anchor.empty());
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, SubKSeedsHonorComplexityFilter) {
    if constexpr (std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        for (auto mode : {DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL}) {
            SCOPED_TRACE(mode);
            auto graph = build_graph_batch<TypeParam>(11, {std::string(31, 'A')}, mode);
            Query query(*graph, std::string(9, 'A'));
            DBGAlignerConfig config;
            config.min_seed_length = 9;
            config.max_seed_length = std::numeric_limits<size_t>::max();
            config.forward_and_reverse_complement = false;
            EXPECT_TRUE(ExactSeeder(query, config).get_anchors().empty());
            config.seed_complexity_filter = false;
            EXPECT_FALSE(ExactSeeder(query, config).get_anchors().empty());

            // The column seeder keeps full MEMs, even when their individual
            // k-mers have low complexity. Preserve that behavior here too.
            Query full_query(*graph, std::string(20, 'A'));
            config.seed_complexity_filter = true;
            EXPECT_FALSE(ExactSeeder(full_query, config).get_anchors().empty());
            for (size_t max_seed_length : {size_t(0), graph->get_k()}) {
                SCOPED_TRACE(max_seed_length);
                config.max_seed_length = max_seed_length;
                config.seed_complexity_filter = true;
                EXPECT_TRUE(ExactSeeder(full_query, config).get_anchors().empty());
                config.seed_complexity_filter = false;
                EXPECT_FALSE(ExactSeeder(full_query, config).get_anchors().empty());
            }
        }
    }
}
#endif

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

#if ! _PROTEIN_GRAPH
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

#if ! _PROTEIN_GRAPH
TEST(SSHashSubKSeeds, RepeatedMinimizerUsesEveryMatchingOccurrence) {
    // The sampled m-mer at position 10 also occurs at position 4. Aligning
    // only to that first occurrence misses the exact query beginning at 7.
    const std::string reference = "CTAAAATCTAAATCTCCA";
    const std::string sequence = "CTAAATCTC";
    const auto filename = test_dump_dir() + "/repeated_minimizer_seed.fa";
    {
        std::ofstream out(filename);
        out << ">reference\n" << reference << '\n';
    }
    DBGSSHash graph(filename, 11, DeBruijnGraph::BASIC, 0, test_dump_dir());
    DBGAlignerConfig config;
    config.min_seed_length = 9;
    config.max_seed_length = 11;
    config.min_exact_match = 1;
    config.forward_and_reverse_complement = false;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    Query query(graph, sequence);
    ExactSeeder seeder(query, config);
    Extender extender(query, config);
    std::vector<Alignment> results;
    align_query(query, seeder, extender, [&](Alignment&& a) {
        results.push_back(std::move(a));
    });
    ASSERT_EQ(1u, results.size());
    EXPECT_FALSE(results[0].get_orientation());
    EXPECT_EQ(sequence, results[0].get_spelling());
    EXPECT_EQ("9=", results[0].get_cigar().to_string());
}

TEST(SSHashSubKSeeds, CanonicalReverseSeedIncludesContigEndpoint) {
    // The sampled five-base minimizer is the final m-mer of this k-mer.
    const std::string reference = "GAGCCATGGGC";
    const std::string sequence = "GCCCATGGCT";
    const auto filename = test_dump_dir() + "/reverse_endpoint_seed.fa";
    {
        std::ofstream out(filename);
        out << ">reference\n" << reference << '\n';
    }
    DBGSSHash graph(filename, 11, DeBruijnGraph::CANONICAL, 0, test_dump_dir());
    DBGAlignerConfig config;
    config.min_seed_length = 9;
    config.max_seed_length = 11;
    config.min_exact_match = 1;
    config.forward_and_reverse_complement = false;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    Query query(graph, sequence);
    ExactSeeder seeder(query, config);
    Extender extender(query, config);
    std::vector<Alignment> results;
    align_query(query, seeder, extender, [&](Alignment&& a) {
        results.push_back(std::move(a));
    });
    ASSERT_EQ(1u, results.size());
    EXPECT_FALSE(results[0].get_orientation());
    EXPECT_EQ(sequence, results[0].get_spelling());
    EXPECT_EQ("10=", results[0].get_cigar().to_string());
}

TYPED_TEST(DBGAlignerRedoneTest, BasicSubKSeedsSearchBothQueryOrientations) {
    if constexpr (std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        const std::string reference =
            "AAAGGATTCTAGCGTTCCAGCCAGACTATTCAAGTCACACCTGCGTTAGAAGACAGCCTTGGATGTAAAT";
        const std::string forward = "AAGGATTCT";
        auto reverse = forward;
        reverse_complement(reverse.begin(), reverse.end());
        auto graph = build_graph_batch<TypeParam>(11, {reference}, DeBruijnGraph::BASIC);
        DBGAlignerConfig config;
        config.min_seed_length = 9;
        config.max_seed_length = 11;
        config.min_exact_match = 1;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
        for (bool reverse_query : {false, true}) {
            for (bool both_orientations : {false, true}) {
                SCOPED_TRACE(::testing::Message() << "reverse query=" << reverse_query
                                                  << " both orientations=" << both_orientations);
                config.forward_and_reverse_complement = both_orientations;
                Query query(*graph, reverse_query ? reverse : forward);
                ExactSeeder seeder(query, config);
                Extender extender(query, config);
                std::vector<Alignment> results;
                align_query(query, seeder, extender, [&](Alignment&& a) {
                    results.push_back(std::move(a));
                });
                if (reverse_query && !both_orientations) {
                    EXPECT_TRUE(results.empty());
                    continue;
                }
                ASSERT_EQ(1u, results.size());
                EXPECT_EQ(reverse_query, results[0].get_orientation());
                EXPECT_EQ(forward, results[0].get_spelling());
                EXPECT_EQ("9=", results[0].get_cigar().to_string());
            }
        }
    }
}
#endif

TYPED_TEST(DBGAlignerRedoneTest, ChainingWithZeroMatchScore) {
    const std::string reference =
        "AAAGGATTCTAGCGTTCCAGCCAGACTATTCAAGTCACACCTGCGTTAGAAGACAGCCTTGGATGTAAAT";
    auto sequence = reference;
    sequence.at(20) = 'A';
    auto graph = build_graph_batch<TypeParam>(11, {reference}, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    config.min_seed_length = 11;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.chain_alignments = true;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(0, -3, -3);
    Query query(*graph, sequence);
    ExactSeeder seeder(query, config);
    auto anchors = seeder.get_inexact_anchors();
    ASSERT_FALSE(anchors.empty());
    for (const auto &anchor : anchors) {
        EXPECT_TRUE(anchor.is_spelling_valid(*graph));
        EXPECT_TRUE(anchor.get_cigar().is_valid(anchor.get_spelling(), anchor.get_query()));
    }
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

TEST(AlignmentCigar, BothImplementationsCheckEveryMismatch) {
    struct Case {
        const char *cigar;
        const char *reference;
        const char *query;
        bool valid;
    };
    for (const auto &[cigar, reference, query, valid] : {
            Case{"2X", "AA", "AT", false},
            Case{"2X", "AA", "TA", false},
            Case{"2X", "AA", "AA", false},
            Case{"2X", "AA", "TC", true},
            Case{"2=", "AA", "AT", false},
            Case{"2=", "AA", "AA", true},
            Case{"1S2X", "AA", "CTC", true},
            Case{"1S2X", "AA", "CAT", false},
            Case{"2X", "AAA", "TC", false},
            Case{"2X", "AA", "TCC", false}}) {
        SCOPED_TRACE(std::string(cigar) + " " + reference + " " + query);
        EXPECT_EQ(valid, Cigar(cigar).is_valid(reference, query));
        const mtg::graph::align::Cigar legacy(cigar);
        // The legacy API receives only the aligned query interval.
        auto legacy_query = std::string_view(query);
        legacy_query.remove_prefix(legacy.get_clipping());
        legacy_query.remove_suffix(legacy.get_end_clipping());
        EXPECT_EQ(valid, legacy.is_valid(reference, legacy_query));
    }
}

TEST(AlignmentTrimming, KeepsLabelClassesWithTrimmedPath) {
    const std::string reference = "AACGTT";
    auto graph = build_graph_batch<DBGHashFast>(4, { reference });
    for (size_t length : {3, 4, 5}) {
        const std::string query = reference.substr(0, length);
        Alignment alignment(*graph, query, false, map_to_nodes_sequentially(*graph, reference),
                            Cigar(Cigar::MATCH, length), std::string(reference),
                            reference.size() - length, 1, {{17}});
        const auto cigar = alignment.get_cigar();
        const auto coords = alignment.get_coords();
        const auto removed = std::min(reference.size() - length, alignment.size() - 1);
        EXPECT_EQ(removed, alignment.trim_end());
        EXPECT_EQ(alignment.size(), alignment.get_label_classes().size());
        EXPECT_EQ(std::vector<Anchor::label_class_t>(alignment.size(), 1),
                  alignment.get_label_classes());
        EXPECT_EQ(reference.size() - length - removed, alignment.get_end_trim());
        EXPECT_EQ(query, alignment.get_spelling());
        EXPECT_EQ(cigar, alignment.get_cigar());
        EXPECT_EQ(coords, alignment.get_coords());
        EXPECT_TRUE(alignment.is_spelling_valid(*graph));
        EXPECT_EQ(0u, alignment.trim_end());
    }
}

TEST(AlignmentExtension, KeepsBetterPathsAfterACompleteQuery) {
    const std::string query = "ACGTTAGGACCTATGCGATTCAGACT";
    const std::string mismatches = "ACGTTAGTACATATGCGATTCAGACT";
    const std::string deletion = "ACGTTAGCCCGACCTATGCGATTCAGACT";
    constexpr size_t k = 7;
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 100;
    config.min_seed_length = config.max_seed_length = k;
    config.forward_and_reverse_complement = false;

    for (bool backward : {false, true}) {
        for (bool add_alternative : {false, true}) {
            SCOPED_TRACE(backward);
            SCOPED_TRACE(add_alternative);
            auto sequence = query;
            auto reference = deletion;
            auto alternative = mismatches;
            if (backward) {
                std::reverse(sequence.begin(), sequence.end());
                std::reverse(reference.begin(), reference.end());
                std::reverse(alternative.begin(), alternative.end());
            }
            auto graph = build_graph_batch<DBGHashFast>(k,
                add_alternative ? std::vector<std::string>{reference, alternative}
                                : std::vector<std::string>{reference});
            Query aln_query(*graph, sequence);
            const size_t begin = backward ? sequence.size() - k : 0;
            Anchor seed(sequence, begin, begin + k, false,
                        {graph->kmer_to_node(sequence.substr(begin, k))});
            Extender extender(aln_query, config);
            size_t found = 0;
            extender.extend(Alignment(seed), [&](Alignment&& alignment) {
                ++found;
                EXPECT_EQ(43, score_match(alignment, config));
                EXPECT_EQ(backward ? "19=3D7=" : "7=3D19=",
                          alignment.get_cigar().to_string());
                EXPECT_TRUE(alignment.is_spelling_valid(*graph));
                EXPECT_TRUE(alignment.get_cigar().is_valid(alignment.get_spelling(),
                                                         alignment.get_query()));
            }, !backward, backward);
            EXPECT_GT(found, 0u);
        }
    }
}

TEST(AlignmentChaining, RetainsEarlierSeedsAfterConnectionLimit) {
    const std::string reference =
        "CGGCGCGAACTACCGCTTGATTCTTTAATCCATAAGATATAACTCGATACTCAGACCCAT";
    std::string sequence = reference;
    sequence.at(17) = 'A';
    sequence.at(38) = 'C';
    auto graph = build_graph_batch<DBGSSHash>(7, {reference});
    Query query(*graph, sequence);
    DBGAlignerConfig config{};
    config.min_seed_length = 7;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.seed_complexity_filter = false;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    const auto seeds = ExactSeeder(query, config).get_inexact_anchors();
    ASSERT_EQ(3u, seeds.size());

    config.chain_alignments = true;
    for (bool limit_ram : {false, true}) {
        SCOPED_TRACE(limit_ram);
        config.max_ram_per_alignment = limit_ram ? 0 : 100;
        config.max_nodes_per_seq_char = limit_ram ? 10000 : 0;
        auto results = ExactSeeder(query, config).get_inexact_anchors();
        ASSERT_EQ(seeds.size(), results.size());
        std::sort(results.begin(), results.end(), [](const auto &a, const auto &b) {
            return a.get_clipping() < b.get_clipping();
        });
        for (size_t i = 0; i < seeds.size(); ++i) {
            EXPECT_EQ(seeds[i].get_cigar(), results[i].get_cigar());
            EXPECT_EQ(seeds[i].get_spelling(), results[i].get_spelling());
            EXPECT_TRUE(results[i].is_spelling_valid(*graph));
        }
    }
}

} // namespace
