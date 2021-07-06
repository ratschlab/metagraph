#include "tax_classifier.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/dac_vector.hpp>

#include "common/serialization.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;
using TaxId = TaxClassifier::TaxId;

void TaxClassifier::import_taxonomy(const std::string &taxdb_filepath) {
    Timer timer;
    logger->trace("Importing metagraph taxonomic data..");

    std::ifstream f(taxdb_filepath.c_str(), std::ios::in | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxdb file {}.", taxdb_filepath.c_str());
        std::exit(1);
    }

    if (!load_number_number_map(f, &this->node_parent)) {
        logger->error("Can't load serialized 'node_parent' from file {}.", taxdb_filepath.c_str());
        std::exit(1);
    }

    code_to_taxid.load(f);
    if (code_to_taxid.empty()) {
        logger->error("Can't load serialized 'code_to_taxid' from file {}.", taxdb_filepath.c_str());
        std::exit(1);
    }

    code.load(f);
    if (code.empty()) {
        logger->error("Can't load serialized 'code' from file {}.", taxdb_filepath.c_str());
        std::exit(1);
    }

    logger->trace("Finished taxdb import after {}s", timer.elapsed());
}

TaxClassifier::TaxClassifier(const std::string &taxdb_filepath,
                             const double lca_coverage_rate,
                             const double kmers_discovery_rate) {
    if (lca_coverage_rate <= 0 || lca_coverage_rate > 1) {
        logger->error("Error: current lca_coverage_rate is {}. Please modify its value to be a ratio: 0 < lca_coverage_rate <= 1.", lca_coverage_rate);
        exit(1);
    }
    assert(kmers_discovery_rate >= 0 && kmers_discovery_rate <= 1);
    Timer timer;
    logger->trace("Constructing Classifier object..");
    import_taxonomy(taxdb_filepath);

    // Find root_node.
    for (const pair<TaxId, TaxId> &it : this->node_parent) {
        if (it.first == it.second) {
            this->root_node = it.first;
            break;
        }
    }
    this->lca_coverage_rate = lca_coverage_rate;
    this->kmers_discovery_rate = kmers_discovery_rate;
    logger->trace("Finished TaxClassifier construction after {}s", timer.elapsed());
}

void TaxClassifier::update_scores_and_lca(const TaxId start_node,
                                          const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                                          const uint64_t desired_number_kmers,
                                          const TaxId root_node,
                                          const tsl::hopscotch_map<TaxId, TaxId> &node_parent,
                                          tsl::hopscotch_map<TaxId, uint64_t> *node_scores,
                                          tsl::hopscotch_set<TaxId> *nodes_already_propagated,
                                          TaxId *best_lca,
                                          uint64_t *best_lca_dist_to_root) {
    if (nodes_already_propagated->count(start_node)) {
        return;
    }
    uint64_t score_from_processed_parents = 0;
    uint64_t score_from_unprocessed_parents = num_kmers_per_node.at(start_node);

    // processed_parents represents the set of nodes on the path start_node->root that have already been processed in the previous iterations.
    std::vector<TaxId> processed_parents;
    std::vector<TaxId> unprocessed_parents;

    TaxId act_node = start_node;
    unprocessed_parents.push_back(act_node);

    while (act_node != root_node) {
        act_node = node_parent.at(act_node);
        if (!nodes_already_propagated->count(act_node)) {
            if (num_kmers_per_node.count(act_node)) {
                score_from_unprocessed_parents += num_kmers_per_node.at(act_node);
            }
            unprocessed_parents.push_back(act_node);
        } else {
            if (num_kmers_per_node.count(act_node)) {
                score_from_processed_parents += num_kmers_per_node.at(act_node);
            }
            processed_parents.push_back(act_node);
        }
    }
    // The score of all the nodes in 'processed_parents' will be updated with 'score_from_unprocessed_parents' only.
    // The nodes in 'unprocessed_parents' will be updated with the sum 'score_from_processed_parents + score_from_unprocessed_parents'.
    for (uint64_t i = 0; i < unprocessed_parents.size(); ++i) {
        TaxId &act_node = unprocessed_parents[i];
        (*node_scores)[act_node] =
                score_from_processed_parents + score_from_unprocessed_parents;
        nodes_already_propagated->insert(act_node);

        uint64_t act_dist_to_root =
                processed_parents.size() + unprocessed_parents.size() - i;

        // Test if the current node's score would be a better LCA result.
        if ((*node_scores)[act_node] >= desired_number_kmers
            && (act_dist_to_root > *best_lca_dist_to_root
            || (act_dist_to_root == *best_lca_dist_to_root && (*node_scores)[act_node] > (*node_scores)[*best_lca]))) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
    for (uint64_t i = 0; i < processed_parents.size(); ++i) {
        TaxId &act_node = processed_parents[i];
        (*node_scores)[act_node] += score_from_unprocessed_parents;

        uint64_t act_dist_to_root = processed_parents.size() - i;
        if ((*node_scores)[act_node] >= desired_number_kmers
            && (act_dist_to_root > *best_lca_dist_to_root
            || (act_dist_to_root == *best_lca_dist_to_root && (*node_scores)[act_node] > (*node_scores)[*best_lca]))) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
}

TaxId TaxClassifier::find_lca(const TaxId a,
                              const TaxId b,
                              const TaxId root_node,
                              const tsl::hopscotch_map<TaxId, TaxId> &node_parent) {
    std::vector<TaxId> ancestors_a;
    std::vector<TaxId> ancestors_b;

    TaxId curr_node = a;
    while (curr_node != root_node) {
        ancestors_a.push_back(curr_node);
        curr_node = node_parent.at(curr_node);
    }
    curr_node = b;
    while (curr_node != root_node) {
        ancestors_b.push_back(curr_node);
        curr_node = node_parent.at(curr_node);
    }

    TaxId lca = root_node;
    int idx_ancs_a = ancestors_a.size() - 1;
    int idx_ancs_b = ancestors_b.size() - 1;

    while (idx_ancs_a >= 0  && idx_ancs_b >= 0) {
        if (ancestors_a[idx_ancs_a] == ancestors_b[idx_ancs_b]) {
            lca = ancestors_a[idx_ancs_a];
        } else {
            break;
        }
        idx_ancs_a --;
        idx_ancs_b --;
    }
    return lca;
}

TaxId TaxClassifier::assign_class(const mtg::graph::DeBruijnGraph &graph,
                                  const std::string &sequence) const {
    tsl::hopscotch_map<TaxId, uint64_t> num_kmers_per_node;

    std::vector<uint64_t> forward_kmers;
    graph.map_to_nodes(sequence, [&](const uint64_t &i) {
        forward_kmers.push_back(i);
    });
    std::string reversed_sequence = sequence;
    reverse_complement(reversed_sequence.begin(), reversed_sequence.end());

    uint64_t total_kmers = forward_kmers.size();
    uint64_t backward_kmer_index = total_kmers;
    std::vector<uint64_t> backward_kmers(total_kmers);

    graph.map_to_nodes(reversed_sequence, [&](const uint64_t &i) {
        backward_kmers[--backward_kmer_index] = i;
    });

    uint64_t total_discovered_kmers = 0;

    // Find the LCA taxid for each kmer without any dependency on the orientation of the read.
    for (uint64_t i = 0; i < total_kmers; ++i) {
        if (forward_kmers[i] == 0 && backward_kmers[i] == 0) {
            continue;
        }

        TaxId curr_taxid;
        if (backward_kmers[i] == 0) {
            curr_taxid = this->code_to_taxid[this->code[forward_kmers[i] - 1]];
        } else if (forward_kmers[i] == 0) {
            curr_taxid = this->code_to_taxid[this->code[backward_kmers[i] - 1]];
        } else {
            // In case that both 'forward_taxid[i]' and 'backward_taxids[i]' are nonzero, compute the LCA.
            TaxId forward_taxid = this->code_to_taxid[this->code[forward_kmers[i] - 1]];
            TaxId backward_taxid = this->code_to_taxid[this->code[backward_kmers[i] - 1]];
            if (forward_taxid == 0) {
                curr_taxid = backward_taxid;
            } else if (backward_taxid == 0) {
                curr_taxid = forward_taxid;
            } else {
                curr_taxid = find_lca(forward_taxid, backward_taxid, this->root_node, this->node_parent);
            }
        }
        if (curr_taxid) {
            total_discovered_kmers += 1;
            num_kmers_per_node[curr_taxid]++;
        }
    }

    if (total_discovered_kmers <= this->kmers_discovery_rate * total_kmers) {
        return 0; // 0 is a wildcard for not enough discovered kmers.
    }

    tsl::hopscotch_set<TaxId> nodes_already_propagated;
    tsl::hopscotch_map<TaxId, uint64_t> node_scores;

    uint64_t desired_number_kmers = total_discovered_kmers * this->lca_coverage_rate;
    TaxId best_lca = this->root_node;
    uint64_t best_lca_dist_to_root = 1;

    // Update the nodes' score by iterating through all the nodes with nonzero kmers.
    for (const pair<TaxId, uint64_t> &node_pair : num_kmers_per_node) {
        TaxId start_node = node_pair.first;
        update_scores_and_lca(start_node, num_kmers_per_node, desired_number_kmers,
                              this->root_node, this->node_parent,
                              &node_scores, &nodes_already_propagated, &best_lca,
                              &best_lca_dist_to_root);
    }
    return best_lca;
}

} // namespace annot
} // namespace mtg
