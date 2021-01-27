#include "classifier.hpp"

#include "common/logger.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream> // TODO delete
#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

typedef Classifier::TaxoLabel TaxoLabel;

void Classifier::import_taxonomy(
        const std::string &filepath,
        std::vector<TaxoLabel> &linearization
    ) {
    std::ifstream f(filepath.c_str(), std::ios::in);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }
    uint64_t size_map;
    uint64_t key;
    uint64_t value;
    f >> size_map;
    for (uint64_t i = 0; i < size_map; ++i) {
        f >> key >> value;
        taxonomic_map[key] = value;
    }

    f >> num_nodes;
    index_to_label.resize(num_nodes);
    for (uint64_t i = 0; i < num_nodes; ++i) {
        f >> index_to_label[i];
    }
    linearization.resize(2 * num_nodes - 1);
    for (uint64_t i = 0; i < 2 * num_nodes - 1; ++i) {
        f >> linearization[i];
    }
    f.close();
}

Classifier::Classifier(const std::string &filepath) {
    std::cout << "inside classifier constructor \n" << std::endl;
    std::vector<TaxoLabel> linearization;
    import_taxonomy(filepath, linearization);
    root_node = linearization[0];

    std::vector<bool> visited_node(num_nodes);
    node_depth.resize(num_nodes);
    node_parent.resize(num_nodes);
    for (uint64_t i = 1; i < linearization.size(); ++i) {
        uint64_t act = linearization[i];
        uint64_t prv = linearization[i - 1];
        if (!visited_node[act]) {
//            Previous node was the parent.
            node_parent[act] = prv;
            visited_node[act] = true;
        } else {
//            previous node was a child;
            if (node_depth[prv] + 1 > node_depth[act]) {
                node_depth[act] = node_depth[prv] + 1;
            }
        }
    }

//    for(uint i = 0; i < num_nodes; ++i) {
//        std::cout << "d[" << i << "]=" << node_depth[i] << "\n";
//    }
}

std::string Classifier::assign_class(const mtg::graph::DeBruijnGraph &graph,
                               const std::string &sequence,
                               const double &lca_coverage_threshold) {
    std::cout << "inside assign_class" << std::endl;
    tsl::hopscotch_map<TaxoLabel, uint64_t> raw_node_score;
    uint64_t total_kmers = 0;

    std::cout << "sequence=" << sequence << std::endl;
//    like in src/graph/annotated_dbg.cpp line 53
    graph.map_to_nodes(sequence, [&](const auto &i) {
        if (i > 0) {
//            std::cout << "map to node " << i << std::endl;
            raw_node_score[taxonomic_map[i]]++;
            total_kmers++;
        }
    });

    std::cout << "before\n";
    for (auto &it: raw_node_score) {
        std::cout << "taxon_before=" << it.first << " -> " << it.second << "\n";
    }
    tsl::hopscotch_set<TaxoLabel> nodes_already_propagated;
    tsl::hopscotch_map<TaxoLabel, uint64_t> prc_node_score;

    uint64_t desired_number_kmers = total_kmers * lca_coverage_threshold;
    TaxoLabel best_lca = root_node;
    for (const auto &node_pair: raw_node_score) {
        uint64_t node = node_pair.first;
        if (nodes_already_propagated.count(node)) {
            continue;
        }
        uint64_t first_processed_parent = node;
        uint64_t amount_from_processed_parents = 0;
        uint64_t amount_from_unprocessed_parents = raw_node_score[node];

        TaxoLabel it_node = node;
        do {
            it_node = node_parent[it_node];
            if (!nodes_already_propagated.count(it_node)) {
                amount_from_unprocessed_parents += raw_node_score[it_node];
            } else {
                first_processed_parent = it_node;
                amount_from_processed_parents = raw_node_score[it_node];
                break;
            }
        } while (it_node != root_node);
        while (it_node != root_node) {
            it_node = node_parent[it_node];
            amount_from_processed_parents += raw_node_score[it_node];
        }
        it_node = node;
        while (it_node != first_processed_parent) {
            prc_node_score[it_node] = amount_from_processed_parents +
                                        amount_from_unprocessed_parents;
            if (!nodes_already_propagated.count(it_node)) {
                nodes_already_propagated.insert(it_node);
            }
            if (prc_node_score[it_node] >= desired_number_kmers &&
                node_depth[it_node] < node_depth[best_lca]) {
                best_lca = it_node;
            }
            it_node = node_parent[it_node];
        }
        prc_node_score[it_node] += amount_from_unprocessed_parents;
        if (prc_node_score[it_node] >= desired_number_kmers &&
            node_depth[it_node] < node_depth[best_lca]) {
            best_lca = it_node;
        }
        while (it_node != root_node) {
            it_node = node_parent[it_node];
            prc_node_score[it_node] += amount_from_unprocessed_parents;
            if (prc_node_score[it_node] >= desired_number_kmers &&
                node_depth[it_node] < node_depth[best_lca]) {
                best_lca = it_node;
            }
        }
    }

    std::cout << "after\n";
    for (auto &it: prc_node_score) {
        std::cout << "taxon_after=" << it.first << " -> " << it.second << "\n";
    }
    std::cout << "\n\n\n";
    return index_to_label[best_lca];
}

}
}