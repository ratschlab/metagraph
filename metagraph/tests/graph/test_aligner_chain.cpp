#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "test_aligner_helpers.hpp"

#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_methods.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;

typedef DBGAligner<>::score_t score_t;

template <typename Graph>
class DBGAlignerChainTest : public DeBruijnGraphTest<Graph> {};

TYPED_TEST_SUITE(DBGAlignerChainTest, FewGraphTypes);

inline void check_chain(const DBGAligner<>::DBGQueryAlignment &paths,
                        const DeBruijnGraph &graph,
                        const DBGAlignerConfig &config,
                        size_t expected_size) {
    for (const auto &path : paths) {
        EXPECT_TRUE(path.is_valid(graph, &config)) << path;
        check_json_dump_load(graph, path, paths.get_query(), paths.get_query(true));
    }

    if (!std::is_sorted(paths.begin(), paths.end(),
                        [](const auto &a, const auto &b) {
                            return a.get_query().end() < b.get_query().begin();
                        })
            || paths.size() != expected_size) {
        for (const auto &path : paths) {
            TEST_COUT << uintptr_t(path.get_query().data()) << " "
                      << uintptr_t(path.get_query().data() + path.get_query().size()) << " "
                      << path << "\n";
        }
    }
}

TYPED_TEST(DBGAlignerChainTest, align_chain_swap) {
    size_t k = 5;
    std::string reference = "ATGATATGA" "TGACCCCGG";
    std::string query     = "TGACCCCGG" "ATGATATGA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.chain_alignments = true;
    graph->add_sequence(reference);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    check_chain(paths, *graph, config, 1);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(query, paths[0].get_sequence());
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerChainTest, align_chain_overlap_2) {
    size_t k = 5;
    std::string reference1 = "TGAGGATCAG";
    std::string reference2 =        "CAGCTAGCTAGCTAGC";
    std::string query      = "TGAGGATCAGCTAGCTAGCTAGC";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.chain_alignments = true;
    graph->add_sequence(reference1);
    graph->add_sequence(reference2);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    check_chain(paths, *graph, config, 1);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(query, paths[0].get_sequence());
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerChainTest, align_chain_overlap_3) {
    size_t k = 5;
    std::string reference1 = "TGAGGATCAG";
    std::string reference2 =        "CAGCTAGCT";
    std::string reference3 =              "GCTTGCTAGC";
    std::string query      = "TGAGGATCAGCTAGCTTGCTAGC";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.chain_alignments = true;
    graph->add_sequence(reference1);
    graph->add_sequence(reference2);
    graph->add_sequence(reference3);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    check_chain(paths, *graph, config, 1);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(query, paths[0].get_sequence());
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerChainTest, align_chain_insert_no_chain_if_full_coverage) {
    size_t k = 10;
    std::string reference = "TGAGGATCAGTTCTAGCTTGCTAGC";
    std::string query     = "TGAGGATCAG""CTAGCTTGCTAGC";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.chain_alignments = true;
    graph->add_sequence(reference);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    check_chain(paths, *graph, config, 1);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(reference, paths[0].get_sequence());
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerChainTest, align_chain_insert1) {
    size_t k = 10;
    std::string reference1 = "TGAGGATCAGTTCTAGCTTG";
    std::string reference2 =             "CTAGCTTGCTAGC";
    std::string query      = "TGAGGATCAG""CTAGCTTGCTAGC";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.chain_alignments = true;
    graph->add_sequence(reference1);
    graph->add_sequence(reference2);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    check_chain(paths, *graph, config, 1);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ("TGAGGATCAGTTCTAGCTTGCTAGC", paths[0].get_sequence());
    check_extend(graph, aligner.get_config(), paths, query);
}

} // namespace
