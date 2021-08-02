#include "gtest/gtest.h"

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#include "seq_io/sequence_io.hpp"
#include "../test_annotated_dbg_helpers.hpp"

#define private public
#define protected public

#include "annotation/taxonomy/tax_classifier.hpp"

namespace {

TEST(TaxonomyTest, ClsAnno_DfsStatistics) {
    std::unique_ptr<mtg::annot::TaxonomyClsAnno> tax = std::make_unique<mtg::annot::TaxonomyClsAnno>();
    tsl::hopscotch_map<uint32_t, std::vector<uint32_t>> tree {
        {0, {1, 2, 3}},      // node 0 -> root
        {1, {4, 5}},         // node 1
        {2, {}},             // node 2
        {3, {6}},            // node 3
        {4, {7, 8}},         // node 4
        {5, {}},
        {6, {}},
        {7, {}},
        {8, {}},
    };

    std::vector<uint32_t> expected_linearization = {
        0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0
    };
    tsl::hopscotch_map<uint32_t, uint32_t> expected_node_depths = {
        {0, 4},
        {1, 3},
        {2, 1},
        {3, 2},
        {4, 2},
        {5, 1},
        {6, 1},
        {7, 1},
        {8, 1},
    };
    tsl::hopscotch_map<uint32_t, uint32_t> expected_node_to_linearization_idx = {
        {0, 0},
        {1, 1},
        {2, 11},
        {3, 13},
        {4, 2},
        {5, 8},
        {6, 14},
        {7, 3},
        {8, 5},
    };

    std::vector<uint32_t> tree_linearization;
    tax->dfs_statistics(0, tree, &tree_linearization);
    EXPECT_EQ(expected_linearization, tree_linearization);
    EXPECT_EQ(expected_node_depths, tax->node_depth_);
    EXPECT_EQ(expected_node_to_linearization_idx, tax->node_to_linearization_idx_);
}

TEST(TaxonomyTest, ClsAnno_RmqPreprocessing) {
    std::unique_ptr<mtg::annot::TaxonomyClsAnno> tax = std::make_unique<mtg::annot::TaxonomyClsAnno>();

    tax->node_depth_ = {
        {0, 4},
        {1, 3},
        {2, 1},
        {3, 2},
        {4, 2},
        {5, 1},
        {6, 1},
        {7, 1},
        {8, 1},
    };

    std::vector<uint32_t> linearization = {
        0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0
    };
    std::vector<std::vector<uint32_t>> expected_rmq = {
        {0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0},
        {0, 1, 4, 4, 4, 4, 1, 1, 1, 0, 0, 0, 0, 3, 3, 0, 0},
        {0, 1, 4, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };

    tax->rmq_preprocessing(linearization);
    EXPECT_EQ(expected_rmq, tax->rmq_data_);
}

TEST (TaxonomyTest, ClsAnno_FindLca) {
    mtg::annot::TaxonomyClsAnno *tax = new mtg::annot::TaxonomyClsAnno();
    /*
     * Tree configuration:
     *      node 0 -> 1 2 3
     *      node 1 -> 4 5
     *      node 2 -> _
     *      node 3 -> 6
     *      node 4 -> 7 8
     */

    tax->rmq_data = {
            {0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0},
            {0, 1, 4, 4, 4, 4, 1, 1, 1, 0, 0, 0, 0, 3, 3, 0, 0},
            {0, 1, 4, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };
    tax->node_to_linearization_idx = {
            {0, 0},
            {1, 1},
            {2, 11},
            {3, 13},
            {4, 2},
            {5, 8},
            {6, 14},
            {7, 3},
            {8, 5},
    };
    tax->node_depth = {
            {0, 4},
            {1, 3},
            {2, 1},
            {3, 2},
            {4, 2},
            {5, 1},
            {6, 1},
            {7, 1},
            {8, 1},
    };

    struct query_lca {
        std::string test_id;
        uint32_t expected;
        std::vector<uint32_t> nodes;
    };

    std::vector<query_lca> queries = {
            {"test1", 0, {7, 6}},
            {"test2", 0, {1, 2}},
            {"test3", 0, {3, 4}},
            {"test4", 0, {1, 2, 5, 6}},
            {"test5", 2, {2}},
            {"test6", 3, {3, 6}},
            {"test6b", 3, {6, 3}},
            {"test7", 1, {7, 8, 5}},
            {"test8", 1, {4, 5}},
            {"test9", 4, {7, 8}},
            {"test10", 0, {0, 1, 2, 3, 4, 5, 6, 7, 8}},
    };

    for(const auto &it: queries) {
    EXPECT_EQ(make_pair(it.test_id, it.expected),
            make_pair(it.test_id, tax->find_lca(it.nodes)));
    }
}

TEST (TaxonomyTest, ClsAnno_ClassifierUpdateScoresAndLca) {
    mtg::annot::TaxonomyClsAnno tax_classifier;

    tax_classifier.root_node = 1;
    tax_classifier.node_parent = {         {1, 1},
                                    {2, 1},       {3, 1},
                                             {4, 3},    {5, 3},
                                        {6, 4}, {7, 4}
                                    };

    tsl::hopscotch_map<uint32_t, uint64_t> num_kmers_per_node = {
        {1, 20}, {2, 1}, {3, 15}, {4, 25}, {5, 6}, {6, 15}, {7, 3}  // leaves 2, 7 and 5 have a smaller number of kmers.
    };

    struct query_tax_map_update {
        std::string test_id;
        std::string description;
        uint64_t desired_number_kmers;
        vector<vector<uint32_t>> ordered_node_sets;
        tsl::hopscotch_map<uint32_t, uint64_t> expected_node_scores;
        tsl::hopscotch_set<uint32_t> expected_nodes_already_propagated;
        uint32_t expected_best_lca;
        uint32_t expected_best_lca_dist_to_root;
    };

    // All the lists in `ordered_node_sets` are covering the entire taxonomic tree.
    // Thus, the evaluation of `update_scores_and_lca` on any of those sets should return the same results.
    vector<vector<uint32_t>> ordered_node_sets = {
            {1, 2, 3, 4, 5, 6, 7},
            {7, 6, 5, 4, 3, 2, 1},
            {7, 4, 6, 3, 5, 1, 2},
            {4, 6, 7, 3, 5, 1, 2},
            {2, 5, 4, 6, 7, 3, 1},
            {2, 6, 7, 5},
            {6, 7, 5, 2},
            {6, 7, 5, 2, 1},
            {3, 5, 6, 7, 2}
    };

    std::vector<query_tax_map_update> tests = {
        {   "test1",
            "desired_number_kmers is equal to node_score[6]; expect LCA taxid = 6",
            75,
            ordered_node_sets,
            {{1, 85}, {2, 21}, {3, 84}, {4, 78}, {5, 41}, {6, 75}, {7, 63}},
            {1, 2, 3, 4, 5, 6, 7},
            6,
            4
        },
        {   "test2",
                "desired_number_kmers is equal to node_score[6]+1; expect LCA taxid = 4",
                76,
                ordered_node_sets,
                {{1, 85}, {2, 21}, {3, 84}, {4, 78}, {5, 41}, {6, 75}, {7, 63}},
                {1, 2, 3, 4, 5, 6, 7},
                4,
                3
        },
        {   "test3",
                "desired_number_kmers is equal to node_score[4]+1; expect LCA taxid = 3",
                79,
                ordered_node_sets,
                {{1, 85}, {2, 21}, {3, 84}, {4, 78}, {5, 41}, {6, 75}, {7, 63}},
                {1, 2, 3, 4, 5, 6, 7},
                3,
                2
        },
        {   "test4",
                "desired_number_kmers is equal to node_score[3]+1; expect LCA taxid = 1",
                85,
                ordered_node_sets,
                {{1, 85}, {2, 21}, {3, 84}, {4, 78}, {5, 41}, {6, 75}, {7, 63}},
                {1, 2, 3, 4, 5, 6, 7},
                1,
                1
        },
        {   "test5",
                "Check updated scores after processing only node 4",
                100,
                {{4}},
                {{4, 60}, {3, 60}, {1, 60}},
                {1, 3, 4},
                1,
                1
        },
        {   "test6",
                "Check updated scores after processing only the nodes 4 and 6",
                100,
                {{4, 6}, {6, 4}},
                {{6, 75}, {4, 75}, {3, 75}, {1, 75}},
                {1, 3, 4, 6},
                1,
                1
        },
        {   "test7",
                "Check updated scores after processing only the nodes 7 and 5",
                100,
                {{7, 5}, {5, 7}},
                {{7, 63}, {4, 63}, {5, 41}, {3, 69}, {1, 69}},
                {1, 3, 4, 5, 7},
                1,
                1
        },
        {   "test8",
                "Check updated scores after processing only the nodes 2, 6 and 7",
                100,
                {{2, 6, 7}, {2, 7, 6}, {6, 2, 7}, {6, 7, 2}, {7, 2, 6}, {7, 6, 2}},
                {{1, 79}, {2, 21}, {3, 78}, {4, 78}, {6, 75}, {7, 63}},
                {1, 2, 3, 4, 6, 7},
                1,
                1
        },
    };

    for (const auto &test: tests) {
        for (std::vector<uint32_t> nodes_set : test.ordered_node_sets) {
            tsl::hopscotch_set<uint32_t> nodes_already_propagated;
            tsl::hopscotch_map<uint32_t, uint64_t> node_scores;
            uint32_t best_lca = tax_classifier.root_node;
            uint32_t best_lca_dist_to_root = 1;

            for (uint64_t node: nodes_set) {
                tax_classifier.update_scores_and_lca(node, num_kmers_per_node, test.desired_number_kmers,
                                                     &node_scores, &nodes_already_propagated,
                                                     &best_lca, &best_lca_dist_to_root);
            }

            EXPECT_EQ(make_pair(test.test_id, test.expected_node_scores),
                    make_pair(test.test_id, node_scores));
            EXPECT_EQ(make_pair(test.test_id, test.expected_nodes_already_propagated),
                    make_pair(test.test_id, nodes_already_propagated));
            EXPECT_EQ(make_pair(test.test_id, test.expected_best_lca),
                    make_pair(test.test_id, best_lca));
            EXPECT_EQ(make_pair(test.test_id, test.expected_best_lca_dist_to_root),
                    make_pair(test.test_id, best_lca_dist_to_root));
        }
    }
}

}
