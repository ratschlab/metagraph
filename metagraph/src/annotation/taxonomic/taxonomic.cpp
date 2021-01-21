#include "taxonomic.hpp"

#include "common/logger.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream> // TODO delete


namespace mtg {
namespace annot {

using mtg::common::logger;

typedef Taxonomy::TaxoLabel TaxoLabel;

void Taxonomy::calculate_node_depth(const TaxoLabel &node, const ChildrenList &tree) {
    uint64_t depth = 0;
    for(auto child: tree[node]) {
        calculate_node_depth(child, tree);
        if (node_depth[child] > depth) {
            depth = node_depth[child];
        }
    }
    node_depth[node] = depth + 1;
}

void Taxonomy::dfs_linearization(
        const TaxoLabel &node,
        const ChildrenList &tree,
        std::vector<TaxoLabel> &tree_linearization
) {
    linerization_idx[node] = tree_linearization.size();
    tree_linearization.push_back(node);
    for(auto child: tree[node]) {
        dfs_linearization(child, tree, tree_linearization);
        tree_linearization.push_back(node);
    }
}

void Taxonomy::read_tree(const std::string &tree_filepath, Taxonomy::ChildrenList &tree) {
    std::ifstream f(tree_filepath);
    std::string parent;
    std::string child;

    TaxoLabel next_available_label = 0;
    while (f >> parent >> child) {
        if (!label_to_index.count(parent)) {
            label_to_index[parent] = next_available_label++;
            tree.push_back(std::vector<TaxoLabel>());
        }
        if (!label_to_index.count(child)) {
            label_to_index[child] = next_available_label++;
            tree.push_back(std::vector<TaxoLabel>());
        }
        tree[label_to_index[parent]].push_back(label_to_index[child]);
    }
    f.close();

    std::vector<uint64_t> node_indegree(tree.size());
    for (uint node = 0; node < tree.size(); ++node) {
        for (uint i = 0; i < tree[node].size(); ++i) {
            node_indegree[tree[node][i]]++;
        }
    }

    for (TaxoLabel i = 0; i < tree.size(); ++i) {
        if (node_indegree[i] == 0) {
            root_node = i;
            break;
        }
    }
}

void Taxonomy::rmq_preprocessing(const ChildrenList &tree) {
    std::vector<TaxoLabel> tree_linearization;
    linerization_idx.resize(tree.size());
    dfs_linearization(root_node, tree, tree_linearization);

    uint num_rmq_rows = log (tree_linearization.size()) + 2;
    rmq_data = new std::unique_ptr<TaxoLabel[]>[num_rmq_rows];
    for (uint i = 0; i < num_rmq_rows; ++i) {
        rmq_data[i] = std::make_unique<TaxoLabel[]>(tree_linearization.size());
    }

//    Copy tree_linearization on rmq line 0.
    for (uint i = 0; i < tree_linearization.size(); ++i) {
        rmq_data[0][i] = tree_linearization[i];
    }

    uint delta = 1;
    for (uint row = 1; row < num_rmq_rows; ++row) {
        for (uint i = 0; i + delta < tree_linearization.size(); ++i) {
            if (node_depth[rmq_data[row-1][i]] > node_depth[rmq_data[row-1][i + delta]]) {
                rmq_data[row][i] = rmq_data[row-1][i];
            } else {
                rmq_data[row][i] = rmq_data[row-1][i + delta];
            }
        }
        delta *= 2;
    }

    for (uint row = 0; row < num_rmq_rows; ++row) {
        for (uint i = 0; i < tree_linearization.size(); ++i) {
            std::cout << rmq_data[row][i] << " ";
        }
        std::cout << "\n";
    }

}

Taxonomy::Taxonomy(const std::string &tree_filepath) {

    if (!std::filesystem::exists(tree_filepath)) {
//        TODO check extension
        logger->error("Can't open taxonomic tree file '{}'.", tree_filepath);
        std::exit(1);
    }
    Taxonomy::ChildrenList tree;
    read_tree(tree_filepath, tree);
    node_depth.resize(tree.size());
    calculate_node_depth(root_node, tree);
    rmq_preprocessing(tree);

    std::cout << "\n\n\n";

    for (uint node = 0; node < tree.size(); ++node) {
        std::cout << node << " -> ";
        for (uint i = 0; i < tree[node].size(); ++i) {
            std::cout << tree[node][i] << " ";
        }
        std::cout << std::endl;
    }

}

}
}
