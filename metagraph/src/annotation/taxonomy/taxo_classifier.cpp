#include "taxo_classifier.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/serialization.hpp"
#include "common/unix_tools.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace annot {

using mtg::common::logger;

typedef TaxoClassifier::NormalizedTaxId NormalizedTaxId;

void TaxoClassifier::import_taxonomy(const std::string &filepath) {
    Timer timer;
    logger->trace("Importing metagraph taxonomic data..");

    std::ifstream f(filepath.c_str(), std::ios::in | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }

    if (!load_number_number_map(f, &taxonomic_map)) {
        logger->error("Can't load serialized 'taxonomic_map' from file '{}'.", filepath.c_str());
        std::exit(1);
    }

    if (!load_string_vector(f, &node_to_acc_version)) {
        logger->error("Can't load serialized 'node_to_acc_version' from file '{}'.", filepath.c_str());
        std::exit(1);
    }
    if (!load_number_vector(f, &node_parent)) {
        logger->error("Can't load serialized 'node_parent' from file '{}'.", filepath.c_str());
        std::exit(1);
    }

    f.close();
    logger->trace("Finished with importing metagraph taxonomic data after '{}' sec", timer.elapsed());
}

TaxoClassifier::TaxoClassifier(const std::string &filepath) {
    Timer timer;
    logger->trace("Constructing Classifier object..");
    import_taxonomy(filepath);

    for (uint64_t i = 0; i < node_parent.size(); ++i) {
        if (i == node_parent[i]) {
            root_node = i;
            break;
        }
    }
    logger->trace("Finished the Taxonomic Classifier's constructor in '{}' sec", timer.elapsed());
}

std::string TaxoClassifier::assign_class(const mtg::graph::DeBruijnGraph &graph,
                                         const std::string &sequence,
                                         const double &lca_coverage_threshold) {
    tsl::hopscotch_map<NormalizedTaxId, uint64_t> raw_node_matches;
    uint64_t total_kmers = 0;

    graph.map_to_nodes(sequence, [&](const auto &i) {
        if (i > 0) {
            // We need this i-1, because of the way how annotation cmd is implemented.
            raw_node_matches[taxonomic_map[i - 1]]++;
            total_kmers++;
        }
    });
    tsl::hopscotch_set<NormalizedTaxId> nodes_already_propagated;

    uint64_t desired_number_kmers = total_kmers * lca_coverage_threshold;
    NormalizedTaxId best_lca = root_node;
    uint64_t best_lca_dist_to_root = 1;
    for (const auto &node_pair: raw_node_matches) {
        uint64_t start_node = node_pair.first;
        if (nodes_already_propagated.count(start_node)) {
            continue;
        }
        uint64_t score_from_processed_parents = 0;
        uint64_t score_from_unprocessed_parents = raw_node_matches[start_node];

        std::vector<NormalizedTaxId> processed_parents;
        std::vector<NormalizedTaxId> unprocessed_parents;

        NormalizedTaxId act_node = start_node;
        unprocessed_parents.push_back(act_node);

        while (act_node != root_node) {
            act_node = node_parent[act_node];
            if (!nodes_already_propagated.count(act_node)) {
                if (raw_node_matches.count(act_node)) {
                    score_from_unprocessed_parents += raw_node_matches[act_node];
                }
                unprocessed_parents.push_back(act_node);
            } else {
                if (raw_node_matches.count(act_node)) {
                    score_from_processed_parents += raw_node_matches[act_node];
                }
                processed_parents.push_back(act_node);
            }
        }

        tsl::hopscotch_map<NormalizedTaxId, uint64_t> node_scores;
        for (uint64_t i = 0; i < unprocessed_parents.size(); ++i) {
            NormalizedTaxId &act_node = unprocessed_parents[i];
            node_scores[act_node] = score_from_processed_parents +
                                      score_from_unprocessed_parents;
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
            NormalizedTaxId &act_node = processed_parents[i];
            node_scores[act_node] += score_from_unprocessed_parents;

            uint64_t act_dist_to_root = processed_parents.size() - i;
            if (node_scores[act_node] >= desired_number_kmers &&
                act_dist_to_root > best_lca_dist_to_root) {
                best_lca = act_node;
                best_lca_dist_to_root = act_dist_to_root;
            }
        }
    }
    return node_to_acc_version[best_lca];
}

}
}
