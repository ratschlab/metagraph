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
    linearization_idx[node] = tree_linearization.size();
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
            index_to_label.push_back(parent);
            tree.push_back(std::vector<TaxoLabel>());
        }
        if (!label_to_index.count(child)) {
            label_to_index[child] = next_available_label++;
            index_to_label.push_back(child);
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
    linearization_idx.resize(tree.size());
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

    precalc_log.resize(tree_linearization.size());
    precalc_pow2.push_back(1);
    for (uint i = 2; i < tree_linearization.size(); ++i) {
        precalc_log[i] = 1 + precalc_log[i/2];
        if(precalc_log[i] > precalc_log[i-1]) {
            precalc_pow2.push_back(i);
        }
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

TaxoLabel Taxonomy::find_lca(const std::vector<TaxoLabel> &labels) {
//    std::cout << "\nin find_lca with: ";
//    for (auto it: labels) {
//        std::cout << it << " ";
//    }
//    std::cout << "\n";
    uint64_t left_idx = linearization_idx[labels[0]];
    uint64_t right_idx = linearization_idx[labels[0]];
    for (const auto &label: labels) {
        if (linearization_idx[label] < left_idx) {
            left_idx = linearization_idx[label];
        }
        if (linearization_idx[label] > right_idx) {
            right_idx = linearization_idx[label];
        }
    }
//    std::cout << "left_idx=" << left_idx << " right_idx=" << right_idx << "\n";
    uint64_t log_dist = precalc_log[right_idx - left_idx];
    uint64_t left_lca = rmq_data[log_dist][left_idx];
    uint64_t right_lca = rmq_data[log_dist][right_idx - precalc_pow2[log_dist] + 1];
//    std::cout << "left_lca=" << left_lca << " right_lca=" << right_lca << "\n";

    if (node_depth[left_lca] > node_depth[right_lca]) {
        return left_lca;
    }
    return right_lca;

}

TaxoLabel Taxonomy::find_lca(const std::vector<Label> &raw_labels) {
    std::vector<TaxoLabel> labels;
    for (const auto &label: raw_labels) {
        if (! label_to_index.count(label)) {
            std::cerr << "\ncannot find label: " << label << " in label_to_index\n";
            exit(1);
        }
        labels.push_back(label_to_index[label]);
    }
    return find_lca(labels);
}

TaxoLabel Taxonomy::find_lca(const TaxoLabel &label1, const TaxoLabel &label2) {
    return find_lca(std::vector<TaxoLabel> {label1, label2});
}

void Taxonomy::update_row_indices(const std::vector<row_index> &indices, const std::vector<Label> &labels) {
    std:: cout << "here update_row_indices :\n";
    for (auto i: indices) {
        std::cout << i << " ";
    }
    std::cout << "\n";
    for (auto i: labels) {
        std::cout << i << " ";
    }

    TaxoLabel lca = find_lca(labels);
    std::cout << "\n lca =" << lca << "\n";

    for (const auto &kmer: indices) {
        if (taxonomic_map.count(kmer)) {
            taxonomic_map[kmer] = find_lca(taxonomic_map[kmer], lca);
//            std::cout << " upd (" << kmer << ") -> " << taxonomic_map[kmer] << "\n";
        } else {
            taxonomic_map[kmer] = lca;
//            std::cout << " fst (" << kmer << ") -> " << taxonomic_map[kmer] << "\n";
        }
    }
}

void Taxonomy::export_to_file(const std::string &filepath) {
    std::ofstream f(filepath.c_str(), std::ios::out);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }

    std::vector<uint> distrib(30);
    f << taxonomic_map.size() << "\n";
    for (auto &it: taxonomic_map) { // Certainly can be serialized.
        f << it.first << " " << it.second << "\n";
        distrib[it.second] ++;
    }
    for (int i = 0; i < 30; ++i) {
        std::cout << "distrib["<<i<<"]=" << distrib[i] << "\n";
    }

    f << index_to_label.size() << "\n";
    for (auto &it: index_to_label) {
        f << it << " ";
    }
    f << "\n";

    uint64_t num_nodes = index_to_label.size();
    for (uint64_t i = 0; i < 2 * num_nodes - 1; ++i) {
        f << rmq_data[0][i] << " ";
    }
    f << "\n";
    f.close();
}

}
}
