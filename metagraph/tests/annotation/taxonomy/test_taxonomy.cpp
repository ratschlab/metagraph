#include "gtest/gtest.h"

#include "test_taxonomy_helpers.hpp"
#include "../../graph/all/test_dbg_helpers.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "annotation/taxonomy/taxonomic_db.hpp"
#include "annotation/taxonomy/taxo_classifier.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "../test_annotated_dbg_helpers.hpp"

namespace mtg {
namespace test {

const time_t FIXED_SEED = 2021;
const std::string dump_test_data_dir = "../tests/data/dump_test/";
const std::string taxo_db_filepath = dump_test_data_dir + "test.taxo";


// TODO add test for normal (not fast) taxo update
TEST (TaxonomyTest, ConstructTaxoDB_fast) {
    int64_t k = 30;
    int64_t available_ram = 4;
    double lca_coverage_threshold = 0.90;

    TaxonomyTestDataGenerator simulator(FIXED_SEED, dump_test_data_dir);
    ASSERT_TRUE(simulator.is_successful_simulation());

    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    simulator.get_all_extant_sequences(sequences, labels);


    auto anno_graph = build_anno_graph<DBGSuccinct>(k + 1, sequences, labels);
    EXPECT_EQ(sequences.size(), anno_graph->get_annotation().get_all_labels().size());

    annot::TaxonomyDB taxonomyDb(simulator.get_taxo_tree_filepath(),
                                 simulator.get_lookup_table_filepath(),
                                 simulator.get_fasta_headers_filepath());

    taxonomyDb.taxonomic_update_fast(*anno_graph, available_ram);
    taxonomyDb.export_to_file(taxo_db_filepath);

    std::vector<std::string> query_sequences;
    std::vector<uint64_t> expected_taxid;
    std::vector<uint64_t> num_children;
    simulator.get_query_sequences(query_sequences, expected_taxid, num_children);

    annot::TaxoClassifier classifier(taxo_db_filepath);

    uint64_t num_tip_correct = 0;
    uint64_t num_tip_total = 0;
    uint64_t num_internal_correct = 0;
    uint64_t num_internal_total = 0;
    for (uint64_t i = 0; i < query_sequences.size(); ++i) {
        if (num_children[i]) {
            ++num_internal_total;
        } else {
            ++num_tip_total;
        }
        uint64_t estimated = classifier.assign_class(anno_graph->get_graph(),
                                                     query_sequences[i],
                                                     lca_coverage_threshold);
        if (estimated == expected_taxid[i]) {
            if (num_children[i]) {
                ++num_internal_correct;
            } else {
                ++num_tip_correct;
            }
        }
    }
    ASSERT_TRUE(num_tip_total > 19000);
    ASSERT_TRUE(num_internal_total > 1800);

    double correct_tip_classification_rate = 1.0 * num_tip_correct / num_tip_total;
    double correct_internal_classification_rate = 1.0 * num_internal_correct / num_internal_total;

// Results for tip classification: rate=0.999947 -> 19013 of 19014
// Results for internal node classification: rate=0.393496 -> 726 of 1845

    ASSERT_TRUE(correct_tip_classification_rate > 0.99);
    ASSERT_TRUE(correct_internal_classification_rate > 0.39);
}

}
}
