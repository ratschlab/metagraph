#include "test_taxonomy_helpers.hpp"

#include <random>
#include "common/logger.hpp"
#include <fstream>
#include <string>
#include <filesystem>

namespace mtg {
namespace test {

using mtg::common::logger;

typedef TaxonomyTestDataGenerator::TaxId TaxId;

std::string TaxonomyTestDataGenerator::get_random_sequence(const uint64_t seq_length) {
    std::string seq;
    for (uint64_t i = 0; i < seq_length; ++i) {
        seq.push_back(valid_nucleotides[rand() % 4]);
    }
    return seq;
}

std::string TaxonomyTestDataGenerator::get_acc_version(uint64_t x) {
    return fmt::format("NC_{}.1", x);
};

void TaxonomyTestDataGenerator::single_nucleotide_mutation(std::string &seq) {
    seq[rand() % seq.size()] = valid_nucleotides[rand() % 4];
}

void TaxonomyTestDataGenerator::simulate_evolution(std::string curr_seq,
                                                   double curr_time,
                                                   TaxId last_sampled_ancestor) {
    if (curr_time > simulation_time) {
        if (genomes_already_existent.count(curr_seq) == 0 &&
            (*unif_distribution)(generator) <= extant_tip_sampling_prob) {
            species.push_back(Species{curr_seq, curr_time, last_sampled_ancestor,
                                      get_acc_version(species.size())});
            genomes_already_existent.insert(curr_seq);
            species[last_sampled_ancestor].num_children ++;
        }
        return;
    }

    double possible_birth = (*birth_distribution)(generator);
    double possible_death = (*death_distribution)(generator);
    if (possible_death < possible_birth &&
        curr_time + possible_death > simulation_time * no_deaths_in_the_first_x_time_percent) {
        return;
    }

    double birth_time = curr_time + possible_birth;
    while (true) {
        double next_mutation = (*mutation_distribution)(generator);
        if (next_mutation < 0) {
            continue;
        }
        if(curr_time + next_mutation < birth_time) {
            single_nucleotide_mutation(curr_seq);
        } else {
            break;
        }
        curr_time += next_mutation;
    }

    if (genomes_already_existent.count(curr_seq) == 0 &&
        (*unif_distribution)(generator) <= extinct_tip_sampling_prob) {
        species.push_back({curr_seq, curr_time, last_sampled_ancestor,
                           get_acc_version(species.size())});
        genomes_already_existent.insert(curr_seq);
        species[last_sampled_ancestor].num_children ++;
        last_sampled_ancestor = species.size() - 1;
    }

    simulate_evolution(curr_seq, birth_time, last_sampled_ancestor);
    simulate_evolution(curr_seq, birth_time, last_sampled_ancestor);
}

bool TaxonomyTestDataGenerator::generate_fasta_fai(std::string filepath) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        logger->error("Could not open file '{}'.", filepath);
        return false;
    }
    for (uint64_t i = 0; i < species.size(); ++i) {
        fout << fmt::format("gi|1234|ref|{}|\n", species[i].accession_version);
    }
    fout.close();
    return true;
}

bool TaxonomyTestDataGenerator::generate_taxo_tree(std::string filepath) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        logger->error("Could not open file '{}'.", filepath);
        return false;
    }
    for (uint64_t i = 0; i < species.size(); ++i) {
        fout << fmt::format("{}\t|\t{}\n", i, species[i].direct_ancestor);
    }
    fout.close();
    return true;
}

bool TaxonomyTestDataGenerator::generate_lookup_table(std::string filepath) {
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        logger->error("Could not open file '{}'.", filepath);
        return false;
    }
    fout << "accession\t accession.version\ttaxid\n";
    for (uint64_t i = 0; i < species.size(); ++i) {
        fout << fmt::format("none\t{}\t{}\n", species[i].accession_version, i);
    }
    fout.close();
    return true;
}

TaxonomyTestDataGenerator::TaxonomyTestDataGenerator(const time_t seed,
                                                     const std::string test_dump_dirname) {
    srand(seed);
    std::string seq_root = get_random_sequence(seq_length);

    birth_distribution = new std::exponential_distribution<> (birth_rate);
    death_distribution = new std::exponential_distribution<> (death_rate);
    mutation_distribution = new std::normal_distribution<double> (mutation_mean_time, mutation_stddev_time);
    unif_distribution = new std::uniform_real_distribution(0.0, 1.0);

    species.push_back(Species{seq_root, 0, 0, get_acc_version(0)});

    simulate_evolution(seq_root, 0, 0);

    if(std::filesystem::create_directories(test_dump_dirname)) {
        successful_simulation = false;
        logger->error("Failed to create directory '{}'", test_dump_dirname);
    }
    dumb_fasta_fai_filepath = test_dump_dirname + dumb_fasta_fai_filepath;
    if (! generate_fasta_fai(dumb_fasta_fai_filepath)) {
        successful_simulation = false;
        return;
    }

    dumb_taxo_tree_filepath = test_dump_dirname + dumb_taxo_tree_filepath;
    if (! generate_taxo_tree(dumb_taxo_tree_filepath)) {
        successful_simulation = false;
        return;
    }

    dumb_lookup_table_filepath = test_dump_dirname + dumb_lookup_table_filepath;
    if (! generate_lookup_table(dumb_lookup_table_filepath)) {
        successful_simulation = false;
        return;
    }

    successful_simulation = true;
}

void TaxonomyTestDataGenerator::get_query_sequences(std::vector<std::string> &query_sequences,
                                                    std::vector<TaxId> &required_taxid,
                                                    std::vector<uint64_t> &num_children) {
    for (uint64_t sp = 0; sp < species.size(); ++sp) {
        if (species[sp].num_children == 1) {
            // Very hard to be classified.
            continue;
        }
        std::string seq = species[sp].genome;
        for(uint64_t i = 0; i < seq.size(); ++i) {
            if ((*unif_distribution)(generator) <= query_seq_error_rate) {
                seq[i] = valid_nucleotides[rand() % 4];
            }
        }
        query_sequences.push_back(seq);
        required_taxid.push_back(sp);
        num_children.push_back(species[sp].num_children);
    }
}

void TaxonomyTestDataGenerator::get_all_extant_sequences(std::vector<std::string> &sequences,
                                                         std::vector<std::string> &labels) {
    for (auto x: species) {
        if (!x.num_children) {
            sequences.push_back(x.genome);
            labels.push_back(fmt::format("gi|1234|ref|{}|", x.accession_version));
        }
    }
}

bool TaxonomyTestDataGenerator::is_successful_simulation() {
    return successful_simulation;
}

std::string TaxonomyTestDataGenerator::get_taxo_tree_filepath() {
    return dumb_taxo_tree_filepath;
}

std::string TaxonomyTestDataGenerator::get_lookup_table_filepath() {
    return dumb_lookup_table_filepath;
}

std::string TaxonomyTestDataGenerator::get_fasta_headers_filepath() {
    return dumb_fasta_fai_filepath;
}

}
}
