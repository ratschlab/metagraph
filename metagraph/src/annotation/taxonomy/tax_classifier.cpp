#include "tax_classifier.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#include <vector>

#include "common/serialization.hpp"
#include "common/unix_tools.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;
using TaxId = TaxClassifier::TaxId;

void TaxClassifier::import_taxonomy(const std::string &filepath) {
    Timer timer;
    logger->trace("Importing metagraph taxonomic data..");

    std::ifstream f(filepath.c_str(), std::ios::in | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxdb file {}.", filepath.c_str());
        std::exit(1);
    }

    if (!load_number_number_map(f, &this->node_parent)) {
        logger->error("Can't load serialized 'node_parent' from file {}.", filepath.c_str());
        std::exit(1);
    }

    this->taxonomic_map.load(f);
    if (this->taxonomic_map.empty()) {
        logger->error("Can't load serialized 'taxonomic_map' from file {}.", filepath.c_str());
        std::exit(1);
    }
    logger->trace("Finished taxdb importing after {}s", timer.elapsed());
}

TaxClassifier::TaxClassifier(const std::string &filepath,
                             const double lca_coverage_rate,
                             const double kmers_discovery_rate) {
    if (lca_coverage_rate <= 0.5 || lca_coverage_rate > 1) {
        logger->error("Error: current lca_coverage_rate is {}. Please modify its value to be a fraction: 0.5 < lca_coverage_rate <= 1.", lca_coverage_rate);
        exit(1);
    }
    assert(kmers_discovery_rate >= 0 && kmers_discovery_rate <= 1);
    Timer timer;
    logger->trace("Constructing Classifier object..");
    import_taxonomy(filepath);

    // Find root_node.
    for (const pair<TaxId, TaxId> &it : this->node_parent) {
        if (it.first == it.second) {
            this->root_node = it.first;
            break;
        }
    }
    this->lca_coverage_rate = lca_coverage_rate;
    this->kmers_discovery_rate = kmers_discovery_rate;
    logger->trace("Finished TaxClassifier construction in {}s", timer.elapsed());
}

void TaxClassifier::update_scores_and_lca(const TaxId start_node,
                                          const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                                          const uint64_t desired_number_kmers,
                                          tsl::hopscotch_map<TaxId, uint64_t> *node_scores,
                                          tsl::hopscotch_set<TaxId> *nodes_already_propagated,
                                          TaxId *best_lca,
                                          uint64_t *best_lca_dist_to_root) const {
    if (nodes_already_propagated->count(start_node)) {
        return;
    }
    uint64_t score_from_processed_parents = 0;
    uint64_t score_from_unprocessed_parents = num_kmers_per_node.at(start_node);

    std::vector<TaxId> processed_parents;
    std::vector<TaxId> unprocessed_parents;

    TaxId act_node = start_node;
    unprocessed_parents.push_back(act_node);

    while (act_node != this->root_node) {
        act_node = this->node_parent.at(act_node);
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
    for (uint64_t i = 0; i < unprocessed_parents.size(); ++i) {
        TaxId &act_node = unprocessed_parents[i];
        (*node_scores)[act_node] =
                score_from_processed_parents + score_from_unprocessed_parents;
        nodes_already_propagated->insert(act_node);

        uint64_t act_dist_to_root =
                processed_parents.size() + unprocessed_parents.size() - i;
        if ((*node_scores)[act_node] >= desired_number_kmers &&
            act_dist_to_root > *best_lca_dist_to_root) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
    for (uint64_t i = 0; i < processed_parents.size(); ++i) {
        TaxId &act_node = processed_parents[i];
        (*node_scores)[act_node] += score_from_unprocessed_parents;

        uint64_t act_dist_to_root = processed_parents.size() - i;
        if ((*node_scores)[act_node] >= desired_number_kmers &&
            act_dist_to_root > *best_lca_dist_to_root) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
}

TaxId TaxClassifier::assign_class(const mtg::graph::DeBruijnGraph &graph,
                                  const std::string &sequence) const {
    tsl::hopscotch_map<TaxId, uint64_t> num_kmers_per_node;
    uint64_t total_discovered_kmers = 0;
    uint64_t total_nonzero_kmers = 0;

    graph.map_to_nodes(sequence, [&](const uint64_t &i) {
        if (i <= 0) {
            return;
        }
        total_nonzero_kmers ++;
        if (this->taxonomic_map[i - 1]) {
            num_kmers_per_node[this->taxonomic_map[i - 1]]++;
            total_discovered_kmers ++;
        }
    });

    if (total_discovered_kmers < this->kmers_discovery_rate * total_nonzero_kmers) {
        return 0; // 0 is a wildcard for not enough discovered kmers.
    }

    tsl::hopscotch_set<TaxId> nodes_already_propagated;
    tsl::hopscotch_map<TaxId, uint64_t> node_scores;

    uint64_t desired_number_kmers = total_discovered_kmers * this->lca_coverage_rate;
    TaxId best_lca = this->root_node;
    uint64_t best_lca_dist_to_root = 1;
    for (const pair<TaxId, uint64_t> &node_pair : num_kmers_per_node) {
        TaxId start_node = node_pair.first;
        update_scores_and_lca(start_node, num_kmers_per_node, desired_number_kmers,
                              &node_scores, &nodes_already_propagated, &best_lca,
                              &best_lca_dist_to_root);
    }
    return best_lca;
}

} // namespace annot
} // namespace mtg
