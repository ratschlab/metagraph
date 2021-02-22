#ifndef __TEST_TAXONOMY_HELPERS_HPP__
#define __TEST_TAXONOMY_HELPERS_HPP__

#include "gtest/gtest.h"

#include <memory>
#include <vector>
#include <string>
#include <random>

#define protected public
#define private public
#include "graph/annotated_dbg.hpp"
#include <tsl/hopscotch_set.h>

namespace mtg {
namespace test {

/*
 * Use Birth-Death model for generating a TaxonomyDB.
 */
class TaxonomyTestDataGenerator {
  public:
    typedef std::uint64_t TaxId;

  private:
    struct Species {
        std::string genome;
        double time_sampling;
        uint64_t direct_ancestor;
        std::string accession_version;
        uint64_t num_children = 0;
    };

    std::vector<Species> species;
    tsl::hopscotch_set<std::string> genomes_already_existent;

    std::default_random_engine generator;
    std::exponential_distribution<> *birth_distribution;
    std::exponential_distribution<> *death_distribution;
    std::normal_distribution<> *mutation_distribution;
    std::uniform_real_distribution<> *unif_distribution;

    double no_deaths_in_the_first_x_time_percent = 0.1;
    double birth_rate = 0.5;
    double death_rate = 0.3;
    double mutation_mean_time = 0.03;
    double mutation_stddev_time = 0.01;
    double extinct_tip_sampling_prob = 0.05;
    double extant_tip_sampling_prob = 1;
    double query_seq_error_rate = 0.001;
    uint64_t simulation_time = 40;
    uint64_t seq_length = 1000;

    std::vector<char> valid_nucleotides = {'T', 'C', 'A', 'G'};

    bool successful_simulation = false;
    std::string dumb_fasta_fai_filepath = "dumb.fasta.fai";
    std::string dumb_taxo_tree_filepath = "dumb_nodes.dmp";
    std::string dumb_lookup_table_filepath = "dumb.accession2taxid";

    std::string get_random_sequence(const uint64_t seq_length);
    void single_nucleotide_mutation(std::string &seq);
    void simulate_evolution(std::string curr_seq, double curr_time,
                            TaxId last_sampled_ancestor);
    bool generate_fasta_fai(std::string filepath);
    bool generate_taxo_tree(std::string filepath);
    bool generate_lookup_table(std::string filepath);
    std::string get_acc_version(uint64_t x);

  public:
    TaxonomyTestDataGenerator(const time_t seed, const std::string test_dump_dirname);
    std::string get_taxo_tree_filepath();
    std::string get_lookup_table_filepath();
    std::string get_fasta_headers_filepath();
    void get_all_extant_sequences(std::vector<std::string> &sequences,
                                  std::vector<std::string> &labels);
    void get_query_sequences(std::vector<std::string> &query_sequences,
                             std::vector<TaxId> &required_taxid,
                             std::vector<uint64_t> &num_children);

    bool is_successful_simulation();
};

} // namespace test
} // namespace mtg

#endif // __TEST_TAXONOMY_HELPERS_HPP__
