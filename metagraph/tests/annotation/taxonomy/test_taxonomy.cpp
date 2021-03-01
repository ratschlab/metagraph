#include "gtest/gtest.h"
#define private public

#include "annotation/taxonomy/taxonomic_db.hpp"
#include "../test_annotated_dbg_helpers.hpp"
#include "seq_io/sequence_io.hpp"

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>



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

    std::vector<std::vector<uint64_t>> tree;
    tree.push_back({1, 2, 3});  // node 0 -> root
    tree.push_back({4, 5});     // node 1
    tree.push_back({});         // node 2
    tree.push_back({6});        // node 3
    tree.push_back({7, 8});     // node 4

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

    std::vector<std::string> sequences;
    std::vector<std::string> labels;

    /*
     * Tree configuration in '../tests/data/taxo_data/dumb_nodes.dmp'
     *
     *
     *                        0(normalized891)
     *                   2(n906)            \
     *               5(n19)  26(n135)       63(n468)
     *
     *
     *      For node  0 the normalized value is 891      ---> root
     *      For node  2 the normalized value is 906
     *      For node  5 the normalized value is 19       ---> SEQ1
     *      For node 26 the normalized value is 135     ---> SEQ2
     *      For node 63 the normalized value is 468     ---> SEQ3
     *
     */

    int64_t k = 30;
    int64_t SEQ1 = 0;
    int64_t SEQ2 = 20;
    int64_t SEQ3 = 50;

    get_sequences_and_labels_from_file(input_filepath, sequences, labels);
    std::mutex taxo_mutex;

    ASSERT_TRUE(taxo.taxonomic_map.size() == 0);

    // Test the normalized taxids.
    uint64_t normalized_taxid;
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(labels[SEQ1]),
            normalized_taxid));
    EXPECT_EQ(19, normalized_taxid);
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(labels[SEQ2]),
        normalized_taxid));
    EXPECT_EQ(135, normalized_taxid);
    ASSERT_TRUE(taxo.get_normalized_taxid(taxo.get_accession_version_from_label(labels[SEQ3]),
            normalized_taxid));
    EXPECT_EQ(468, normalized_taxid);

//    Add sequence SEQ1 only to taxonomic_map.
    auto anno_graph_seq =
            build_anno_graph<DBGSuccinct>(k + 1,
                                          std::vector<std::string>{sequences[SEQ1]},
                                          std::vector<std::string>{labels[SEQ1]});
    taxo.kmer_to_taxid_map_update(anno_graph_seq->get_annotation(), taxo_mutex);
    tsl::hopscotch_map<uint64_t, uint64_t> frequencies_taxid;
    for (uint64_t i = 0; i < taxo.taxonomic_map.size(); ++i) {
        if (taxo.taxonomic_map[i] > 0) {
            ++frequencies_taxid[taxo.taxonomic_map[i]];
        }
    }
    tsl::hopscotch_map<uint64_t, uint64_t> expected_freq_taxid = {{19, 969}};
    EXPECT_EQ(expected_freq_taxid, frequencies_taxid);

//    Add sequences SEQ1 and SEQ2 to taxonomic_map.
    anno_graph_seq =
            build_anno_graph<DBGSuccinct>(k + 1,
                                          std::vector<std::string>{sequences[SEQ1], sequences[SEQ2]},
                                          std::vector<std::string>{labels[SEQ1], labels[SEQ2]});
    taxo.kmer_to_taxid_map_update(anno_graph_seq->get_annotation(), taxo_mutex);
    frequencies_taxid.clear();
    for (uint64_t i = 0; i < taxo.taxonomic_map.size(); ++i) {
        if (taxo.taxonomic_map[i] > 0) {
            ++frequencies_taxid[taxo.taxonomic_map[i]];
        }
    }
    expected_freq_taxid = {{906, 486}, {135, 484}, {19, 969}};
    EXPECT_EQ(expected_freq_taxid, frequencies_taxid);

//    Add sequences SEQ1, SEQ2 and SEQ3 to taxonomic_map.
    anno_graph_seq =
            build_anno_graph<DBGSuccinct>(k + 1,
                                          std::vector<std::string>{sequences[SEQ1], sequences[SEQ2], sequences[SEQ3]},
                                          std::vector<std::string>{labels[SEQ1], labels[SEQ2], labels[SEQ3]});
    taxo.kmer_to_taxid_map_update(anno_graph_seq->get_annotation(), taxo_mutex);
    frequencies_taxid.clear();
    for (uint64_t i = 0; i < taxo.taxonomic_map.size(); ++i) {
        if (taxo.taxonomic_map[i] > 0) {
            ++frequencies_taxid[taxo.taxonomic_map[i]];
        }
    }
    expected_freq_taxid = {{906, 790}, {891, 656}, {19, 648},
                            {468, 313}, {135, 502}};
    EXPECT_EQ(expected_freq_taxid, frequencies_taxid);
}

}
}
