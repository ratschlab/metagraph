#include "taxo_classifier.hpp"

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

typedef TaxoClassifier::TaxId TaxId;

void TaxoClassifier::import_taxonomy(const std::string &filepath) {
    Timer timer;
    logger->trace("Importing metagraph taxonomic data..");

    std::ifstream f(filepath.c_str(), std::ios::in | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }

    if (!load_number_number_map(f, &node_parent)) {
        logger->error("Can't load serialized 'node_parent' from file '{}'.", filepath.c_str());
        std::exit(1);
    }

    taxonomic_map.load(f);
    if (taxonomic_map.empty()) {
        logger->error("Can't load serialized 'taxonomic_map' from file '{}'.", filepath.c_str());
        std::exit(1);
    }
    f.close();
    logger->trace("Finished with importing metagraph taxonomicDB after '{}' sec", timer.elapsed());
}

TaxoClassifier::TaxoClassifier(const std::string &filepath) {
    Timer timer;
    logger->trace("Constructing Classifier object..");
    import_taxonomy(filepath);

    for (const auto &it: node_parent) {
        if (it.first == it.second) {
            root_node = it.first;
            break;
        }
    }
    logger->trace("Finished the Taxonomic Classifier's constructor in '{}' sec", timer.elapsed());
}

void TaxoClassifier::update_scores_and_lca(const TaxId start_node,
                                           const tsl::hopscotch_map<TaxId, uint64_t> &num_kmers_per_node,
                                           const uint64_t &desired_number_kmers,
                                           tsl::hopscotch_map<TaxId, uint64_t> &node_scores,
                                           tsl::hopscotch_set<TaxId> &nodes_already_propagated,
                                           TaxId &best_lca,
                                           uint64_t &best_lca_dist_to_root) {
    if (nodes_already_propagated.count(start_node)) {
        return;
    }
    uint64_t score_from_processed_parents = 0;
    uint64_t score_from_unprocessed_parents = num_kmers_per_node.at(start_node);

    std::vector<TaxId> processed_parents;
    std::vector<TaxId> unprocessed_parents;

    TaxId act_node = start_node;
    unprocessed_parents.push_back(act_node);

    while (act_node != root_node) {
        act_node = node_parent[act_node];
        if (!nodes_already_propagated.count(act_node)) {
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
        node_scores[act_node] =
                score_from_processed_parents + score_from_unprocessed_parents;
        nodes_already_propagated.insert(act_node);

        uint64_t act_dist_to_root =
                processed_parents.size() + unprocessed_parents.size() - i;
        if (node_scores[act_node] >= desired_number_kmers &&
            act_dist_to_root > best_lca_dist_to_root) {
            best_lca = act_node;
            best_lca_dist_to_root = act_dist_to_root;
        }
    }
    for (uint64_t i = 0; i < processed_parents.size(); ++i) {
        TaxId &act_node = processed_parents[i];
        node_scores[act_node] += score_from_unprocessed_parents;

        uint64_t act_dist_to_root = processed_parents.size() - i;
        if (node_scores[act_node] >= desired_number_kmers &&
            act_dist_to_root > best_lca_dist_to_root) {
            best_lca = act_node;
            best_lca_dist_to_root = act_dist_to_root;
        }
    }
}

TaxId TaxoClassifier::assign_class(const mtg::graph::DeBruijnGraph &graph,
                                   const std::string &sequence,
                                   const double &lca_coverage_threshold) {
    if (lca_coverage_threshold <= 0.5) {
        logger->error("Error: received lca coverage threshold {} <= 0.5. Please modify its value to be > 0.5 for having a unique taxid lca solution.", lca_coverage_threshold);
        exit(1);
    }
    tsl::hopscotch_map<TaxId, uint64_t> num_kmers_per_node;
    uint64_t total_kmers = 0;

    graph.map_to_nodes(sequence, [&](const auto &i) {
      if (i > 0) {
          // We need this i-1, because of the way how annotation cmd is implemented.
          num_kmers_per_node[taxonomic_map[i - 1]]++;
          total_kmers++;
      }
    });

    tsl::hopscotch_set<TaxId> nodes_already_propagated;
    tsl::hopscotch_map<TaxId, uint64_t> node_scores;

    uint64_t desired_number_kmers = total_kmers * lca_coverage_threshold;
    TaxId best_lca = root_node;
    uint64_t best_lca_dist_to_root = 1;
    for (const auto &node_pair: num_kmers_per_node) {
        TaxId start_node = node_pair.first;
        update_scores_and_lca(start_node, num_kmers_per_node, desired_number_kmers,
                              node_scores, nodes_already_propagated, best_lca,
                              best_lca_dist_to_root);
    }
    return best_lca;
}

}
}
