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
#include "common/serialization.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

typedef Classifier::TaxoLabel TaxoLabel;

void Classifier::import_taxonomy(
        const std::string &filepath,
        std::vector<TaxoLabel> &linearization
    ) {
    std::cout << "in inport" << std::endl;
    std::ifstream f(filepath.c_str(), std::ios::in | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }
    std::cout << "in inport 2" << std::endl;
//    uint64_t size_map;
//    uint64_t key;
//    uint64_t value;
//    f >> size_map;
//    for (uint64_t i = 0; i < size_map; ++i) {
//        f >> key >> value;
//        taxonomic_map[key] = value;
//    }
    if (!load_number_number_map(f, &taxonomic_map)) {
        logger->error("Can't load serialized 'taxonomic_map' from file '{}'.", filepath.c_str());
        std::exit(1);
    }

    std::vector<uint> distrib(30);
    for (auto &it: taxonomic_map) {
        distrib[it.second]++;
//        std::cout << it.first << " " << it.second << "\n";
    }
    for (int i = 0; i < 30; ++i) {
        std::cout << distrib[i] << " ";
    }

    std::cout << "in inport 3" << std::endl;
//    f >> num_nodes;
//    index_to_label.resize(num_nodes);
//    for (uint64_t i = 0; i < num_nodes; ++i) {
//        f >> index_to_label[i];
//    }
    if (!load_string_vector(f, &index_to_label)) {
        logger->error("Can't load serialized 'index_to_label' from file '{}'.", filepath.c_str());
        std::exit(1);
    }
    num_nodes = index_to_label.size();

    for (auto &it: index_to_label) {
        std::cout << it << " ";
    }
    std::cout << "\n";


    std::cout << "in inport 4" << std::endl;
//    linearization.resize(2 * num_nodes - 1);
//    for (uint64_t i = 0; i < 2 * num_nodes - 1; ++i) {
//        f >> linearization[i];
//    }
//    F_ = load_number_vector_raw<edge_index>(instream);

    linearization = load_number_vector_raw<TaxoLabel>(f);

    for (auto &it : linearization) {
        std::cout << it << " ";
    }
    std::cout << "\n";

    std::cout << "in inport fin" << std::endl;
//    if (!linearization = load_number_vector_raw<TaxoLabel>(f)) {
//        logger->error("Can't load serialized 'linearization' from file '{}'.", filepath.c_str());
//        std::exit(1);
//    }
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
        std::cout << "in this for i=" << i << std::endl;
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
            std::cout << i << " ";
            raw_node_score[taxonomic_map[i - 1]]++;
            total_kmers++;
        }
    });
    std::cout << "\n\n";

    std::cout << "before\n";
    for (const auto &it: raw_node_score) {
        std::cout << "taxon_before=" << it.first << " -> " << it.second << "\n";
    }
    tsl::hopscotch_set<TaxoLabel> nodes_already_propagated;
    tsl::hopscotch_map<TaxoLabel, uint64_t> prc_node_score;

    uint64_t desired_number_kmers = total_kmers * lca_coverage_threshold;
    TaxoLabel best_lca = root_node;
    for (const auto &node_pair: raw_node_score) {
        uint64_t node = node_pair.first;
        std::cout << "node in for = " << node << "\n";
        if (nodes_already_propagated.count(node)) {
            std::cout << "skip\n";
            continue;
        }
        uint64_t amount_from_processed_parents = 0;
        uint64_t amount_from_unprocessed_parents = raw_node_score[node];

        std::vector<TaxoLabel> processed_parents;
        std::vector<TaxoLabel> unprocessed_parents;

        TaxoLabel it_node = node;
        unprocessed_parents.push_back(it_node);
        do {
            it_node = node_parent[it_node];
            if (!nodes_already_propagated.count(it_node)) {
                if (raw_node_score.count(it_node)) {
                    amount_from_unprocessed_parents += raw_node_score[it_node];
                }
                unprocessed_parents.push_back(it_node);
            } else {
                if (raw_node_score.count(it_node)) {
                    amount_from_processed_parents += raw_node_score[it_node];
                }
                processed_parents.push_back(it_node);
            }
        } while (it_node != root_node);

        std::cout << " amount_from_processed_parents=" << amount_from_processed_parents <<
                " amount_from_unprocessed_parents=" << amount_from_unprocessed_parents <<
                "\n";

        for (const auto &it_node: unprocessed_parents) {
            prc_node_score[it_node] = amount_from_processed_parents +
                                      amount_from_unprocessed_parents;
            nodes_already_propagated.insert(it_node);
            if (prc_node_score[it_node] >= desired_number_kmers &&
                node_depth[it_node] < node_depth[best_lca]) {
                best_lca = it_node;
            }
        }
        for (const auto &it_node: processed_parents) {
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