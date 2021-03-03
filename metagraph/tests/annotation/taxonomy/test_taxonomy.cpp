#include "gtest/gtest.h"

#define private public

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#include <stdio.h>
#include <stdlib.h>

#include "annotation/taxonomy/taxonomic_db.hpp"
#include "../test_annotated_dbg_helpers.hpp"
#include "seq_io/sequence_io.hpp"


namespace mtg {
namespace test {

const std::string test_data_dir = "../tests/data/taxo_data/";
const std::string lookup_filepath = test_data_dir + "dumb.accession2taxid";
const std::string taxonomic_filepath = test_data_dir + "dumb_nodes.dmp";
const std::string input_filepath = test_data_dir + "taxo_input.fa";


tsl::hopscotch_set<std::string> get_all_labels_from_file(const std::string &filepath) {
    tsl::hopscotch_set<std::string> labels;
    seq_io::FastaParser fasta_parser(filepath);
    for (const seq_io::kseq_t &kseq : fasta_parser) {
        labels.insert(annot::TaxonomyDB::get_accession_version_from_label(kseq.name.s));
    }
    return labels;
}

void get_sequences_and_labels_from_file(const std::string &filepath,
                                        std::vector<std::string> &sequences,
                                        std::vector<std::string> &labels) {
    sequences.clear();
    labels.clear();
    seq_io::FastaParser fasta_parser(filepath);

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        labels.push_back(kseq.name.s);
        sequences.push_back(kseq.seq.s);
    }
}

TEST (TaxonomyTest, DfsStatistics) {
    tsl::hopscotch_set<std::string> labels = get_all_labels_from_file(input_filepath);
    annot::TaxonomyDB taxo(taxonomic_filepath, lookup_filepath, labels);
    taxo.node_depth.clear();
    taxo.node_depth.resize(9);
    taxo.node_to_linearization_idx.clear();
    taxo.node_to_linearization_idx.resize(9);

    std::vector<std::vector<uint64_t>> tree {
        {1, 2, 3},      // node 0 -> root
        {4, 5},         // node 1
        {},             // node 2
        {6},            // node 3
        {7, 8},         // node 4
        {},
        {},
        {},
        {},
    };

    std::vector<uint64_t> expected_linearization = {
        0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0
    };
    std::vector<uint64_t> expected_node_depths = {
        4, 3, 1, 2, 2, 1, 1, 1, 1
    };
    std::vector<uint64_t> expected_node_to_linearization_idx = {
        0, 1, 11, 13, 2, 8, 14, 3, 5
    };

    std::vector<uint64_t> tree_linearization;
    taxo.dfs_statistics(0, tree, tree_linearization);
    EXPECT_EQ(expected_linearization, tree_linearization);
    EXPECT_EQ(expected_node_depths, taxo.node_depth);
    EXPECT_EQ(expected_node_to_linearization_idx, taxo.node_to_linearization_idx);
}

TEST (TaxonomyTest, RmqPreprocessing) {
    tsl::hopscotch_set<std::string> labels = get_all_labels_from_file(input_filepath);
    annot::TaxonomyDB taxo(taxonomic_filepath, lookup_filepath, labels);

    for (uint64_t i = 0; i < taxo.rmq_data.size(); ++i) {
        taxo.rmq_data[i].clear();
    }
    taxo.rmq_data.clear();
    taxo.precalc_log2.clear();
    taxo.precalc_pow2.clear();
    taxo.node_depth = {4, 3, 1, 2, 2, 1, 1, 1, 1};

    std::vector<uint64_t> linearization = {
            0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0
    };
    std::vector<uint64_t> expected_pow2 = {1, 2, 4, 8, 16};
    std::vector<uint64_t> expected_log2 = {
        0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4
    };
    std::vector<std::vector<uint64_t> > expected_rmq = {
        {0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0},
        {0, 1, 4, 4, 4, 4, 1, 1, 1, 0, 0, 0, 0, 3, 3, 0, 0},
        {0, 1, 4, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };

    taxo.rmq_preprocessing(linearization);
    EXPECT_EQ(expected_pow2, taxo.precalc_pow2);
    EXPECT_EQ(expected_log2, taxo.precalc_log2);
    EXPECT_EQ(expected_rmq, taxo.rmq_data);
}

TEST (TaxonomyTest, FindLca) {
    tsl::hopscotch_set<std::string> labels = get_all_labels_from_file(input_filepath);
    annot::TaxonomyDB taxo(taxonomic_filepath, lookup_filepath, labels);

    /*
     * Tree configuration:
     *      node 0 -> 1 2 3
     *      node 1 -> 4 5
     *      node 2 -> _
     *      node 3 -> 6
     *      node 4 -> 7 8
     */

    taxo.rmq_data = {
        {0, 1, 4, 7, 4, 8, 4, 1, 5, 1, 0, 2, 0, 3, 6, 3, 0},
        {0, 1, 4, 4, 4, 4, 1, 1, 1, 0, 0, 0, 0, 3, 3, 0, 0},
        {0, 1, 4, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };
    taxo.node_to_linearization_idx = {0, 1, 11, 13, 2, 8, 14, 3, 5};
    taxo.precalc_pow2 = {1, 2, 4, 8, 16};
    taxo.precalc_log2 = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4};
    taxo.node_depth = {4, 3, 1, 2, 2, 1, 1, 1, 1};

    struct query_lca {
        std::string test_id;
        uint64_t expected;
        std::vector<uint64_t> nodes;
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
                  make_pair(it.test_id, taxo.find_lca(it.nodes)));
    }
}

TEST (TaxonomyTest, KmerToTaxidUpdate) {
    tsl::hopscotch_set<std::string> labels_set = get_all_labels_from_file(input_filepath);
    annot::TaxonomyDB taxo(taxonomic_filepath, lookup_filepath, labels_set);

    std::vector<std::string> all_sequences;
    std::vector<std::string> all_labels;

    /*
     * Tree configuration in '../tests/data/taxo_data/dumb_nodes.dmp'
     *
     *
     *                        0(normalized650)
     *                    5(n651)            \
     *               47(n568) 48(n503)       2(n310)
     *
     *
     *      For node  0 the normalized value is 650       ---> root
     *      For node  5 the normalized value is 651
     *      For node  2 the normalized value is 310       ---> SEQ1
     *      For node 47 the normalized value is 568       ---> SEQ2
     *      For node 48 the normalized value is 503       ---> SEQ3
     *
     */

    uint64_t k = 20;
    uint64_t SEQ1 = 0;
    uint64_t SEQ2 = 40;
    uint64_t SEQ3 = 41;

    get_sequences_and_labels_from_file(input_filepath, all_sequences, all_labels);
    std::mutex taxo_mutex;

    ASSERT_TRUE(taxo.taxonomic_map.size() == 0);

    // Iterating a hopscotch_map on linux and on darwin returns the objects in a different order.
    // Thus, the normalization ids are different and we need to compute all of them here.
    uint64_t normalized_taxid_seq1;
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(all_labels[SEQ1]),
        normalized_taxid_seq1));

    uint64_t normalized_taxid_seq2;
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(all_labels[SEQ2]),
        normalized_taxid_seq2));

    uint64_t normalized_taxid_seq3;
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(all_labels[SEQ3]),
        normalized_taxid_seq3));

    uint64_t normalized_lca_seq23 = taxo.find_lca({normalized_taxid_seq2,
                                                    normalized_taxid_seq3});
    uint64_t normalized_lca_seq123 = taxo.find_lca({normalized_taxid_seq1,
                                                    normalized_taxid_seq2,
                                                    normalized_taxid_seq3});

    struct query_taxo_map_update {
        std::string test_id;
        tsl::hopscotch_map<uint64_t, uint64_t> expected_freq_taxid;
        std::vector<uint64_t> seq_list;
    };

    std::vector<query_taxo_map_update> tests = {
        {"test1", {{normalized_taxid_seq1, 979}}, {SEQ1}},
        {"test2", {
                           {normalized_lca_seq123, 201},
                           {normalized_taxid_seq2, 778},
                           {normalized_taxid_seq1, 778}
                   }, {{SEQ1, SEQ2}}},
        {"test3", {
                           {normalized_lca_seq123, 217},
                           {normalized_lca_seq23, 650},
                           {normalized_taxid_seq1, 762},
                           {normalized_taxid_seq2, 128},
                           {normalized_taxid_seq3, 171}
                   }, {{SEQ1, SEQ2, SEQ3}}
        },
        {"test4", {
                           {normalized_lca_seq23, 792},
                           {normalized_taxid_seq2, 187},
                           {normalized_taxid_seq3, 187}
                   }, {{SEQ2, SEQ3}}},
    };

    for (const auto &test: tests) {
        // Clean taxonomic_map after the last test.
        for (uint64_t i = 0; i < taxo.taxonomic_map.size(); ++i) {
            taxo.taxonomic_map[i] = 0;
        }
        // Resize taxonomic_map to 0 for testing the piece of code which does taxonomic_map resizing.
        taxo.taxonomic_map.resize(0);

        std::vector<std::string> test_sequences;
        std::vector<std::string> test_labels;
        for (const auto &seq: test.seq_list) {
            test_sequences.push_back(all_sequences[seq]);
            test_labels.push_back(all_labels[seq]);
        }

        auto anno_graph_seq =
                build_anno_graph<DBGSuccinct>(k + 1, test_sequences, test_labels);

        taxo.kmer_to_taxid_map_update(anno_graph_seq->get_annotation(), taxo_mutex);
        tsl::hopscotch_map<uint64_t, uint64_t> frequencies_taxid;
        for (uint64_t i = 0; i < taxo.taxonomic_map.size(); ++i) {
            if (taxo.taxonomic_map[i] > 0) {
                ++frequencies_taxid[taxo.taxonomic_map[i]];
            }
        }

        EXPECT_EQ(make_pair(test.test_id, test.expected_freq_taxid),
                  make_pair(test.test_id, frequencies_taxid));
    }
}

}
}
