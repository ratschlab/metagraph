#include "gtest/gtest.h"

#define TESTING

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "annotation/taxonomy/tax_classifier.hpp"

namespace mtg {
namespace test {

TEST (TaxonomyTest, ClsAnno_DfsStatistics) {
    mtg::annot::TaxonomyClsAnno *tax = new mtg::annot::TaxonomyClsAnno();
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
    EXPECT_EQ(expected_node_depths, tax->node_depth);
    EXPECT_EQ(expected_node_to_linearization_idx, tax->node_to_linearization_idx);
}

TEST (TaxonomyTest, ClsAnno_RmqPreprocessing) {
    mtg::annot::TaxonomyClsAnno *tax = new mtg::annot::TaxonomyClsAnno();

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

    std::vector<uint32_t> linearization = {
            0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0
    };
    std::vector<std::vector<uint32_t> > expected_rmq = {
        {0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0},
        {0, 1, 4, 4, 4, 4, 1, 1, 1, 0, 0, 0, 0, 3, 3, 0, 0},
        {0, 1, 4, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };

    tax->rmq_preprocessing(linearization);
    EXPECT_EQ(expected_rmq, tax->rmq_data);
}

}
}
