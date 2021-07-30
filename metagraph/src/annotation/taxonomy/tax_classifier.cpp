#include "tax_classifier.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "common/unix_tools.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

TaxonomyClsAnno::TaxonomyClsAnno(const graph::AnnotatedDBG &anno,
                                 const double lca_coverage_rate,
                                 const double kmers_discovery_rate,
                                 const std::string &tax_tree_filepath,
                                 const std::string &label_taxid_map_filepath) : _anno_matrix(&anno) {
    _lca_coverage_rate = lca_coverage_rate;
    _kmers_discovery_rate = kmers_discovery_rate;

    if (!std::filesystem::exists(tax_tree_filepath)) {
        logger->error("Can't open taxonomic tree file {}.", tax_tree_filepath);
        std::exit(1);
    }

    bool require_accversion_to_taxid_map = false;
    assign_label_type(_anno_matrix->get_annotation().get_all_labels()[0], &require_accversion_to_taxid_map);

    Timer timer;
    if (require_accversion_to_taxid_map) {
        logger->trace("Parsing label_taxid_map file..");
        read_accversion_to_taxid_map(label_taxid_map_filepath, _anno_matrix);
        logger->trace("Finished label_taxid_map file in {}s", timer.elapsed());
    }

    timer.reset();
    logger->trace("Parsing taxonomic tree..");
    ChildrenList tree;
    read_tree(tax_tree_filepath, &tree);
    logger->trace("Finished taxonomic tree read in {}s.", timer.elapsed());

    timer.reset();
    logger->trace("Calculating tree statistics..");
    std::vector<TaxId> tree_linearization;
    dfs_statistics(root_node, tree, &tree_linearization);
    logger->trace("Finished tree statistics calculation in {}s.", timer.elapsed());

    timer.reset();
    logger->trace("Starting rmq preprocessing..");
    rmq_preprocessing(tree_linearization);
    logger->trace("Finished rmq preprocessing in {}s.", timer.elapsed());
}

void TaxonomyClsAnno::read_tree(const std::string &tax_tree_filepath,
                                ChildrenList *tree) {
    std::ifstream f(tax_tree_filepath);
    if (!f.good()) {
        logger->error("Failed to open Taxonomic Tree file {}.", tax_tree_filepath);
        exit(1);
    }

    std::string line;
    tsl::hopscotch_map<TaxId, TaxId> full_parents_list;
    while (getline(f, line)) {
        if (line == "") {
            logger->error("The Taxonomic Tree file contains empty lines. Please make sure that this file was not manually modified: {}.",
                          tax_tree_filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("The Taxonomic tree filepath contains incomplete lines. Please make sure that this file was not manually modified: {}.",
                          tax_tree_filepath);
            exit(1);
        }
        uint32_t act = static_cast<uint32_t>(std::stoull(parts[0]));
        uint32_t parent = static_cast<uint32_t>(std::stoull(parts[2]));
        full_parents_list[act] = parent;
        this->node_parent[act] = parent;
    }

    std::vector<TaxId> relevant_taxids;
    // 'considered_relevant_taxids' is used to make sure that there are no duplications in 'relevant_taxids'.
    tsl::hopscotch_set<TaxId> considered_relevant_taxids;

    if (this->accversion_to_taxid_map.size()) {
        // Store only the taxonomic nodes that exists in the annotation matrix.
        for (const pair<std::string, TaxId>  &it : this->accversion_to_taxid_map) {
            relevant_taxids.push_back(it.second);
            considered_relevant_taxids.insert(it.second);
        }
    } else {
        // If 'this->accversion_to_taxid_map' is empty, store the entire taxonomic tree.
        for (auto it : full_parents_list) {
            relevant_taxids.push_back(it.first);
            considered_relevant_taxids.insert(it.first);
        }
    }
    assert(relevant_taxids.size());

    uint64_t num_taxid_failed = 0; // num_taxid_failed is used for logging only.
    for (uint32_t i = 0; i < relevant_taxids.size(); ++i) {
        const TaxId taxid = relevant_taxids[i];
        if (!full_parents_list.count(taxid)) {
            num_taxid_failed += 1;
            continue;
        }

        if (considered_relevant_taxids.find(full_parents_list[taxid]) == considered_relevant_taxids.end()) {
            relevant_taxids.push_back(full_parents_list[taxid]);
            considered_relevant_taxids.insert(full_parents_list[taxid]);
        }

        // Check if the current taxid is the root.
        if (taxid == full_parents_list[taxid]) {
            this->root_node = taxid;
        }
    }
    if (num_taxid_failed) {
        logger->warn("During the tax_tree_filepath {} parsing, {} taxids were not found out of {} evaluations.",
                     tax_tree_filepath, num_taxid_failed, relevant_taxids.size());
    }

    // Construct the output tree.
    for (const TaxId &taxid : relevant_taxids) {
        if (taxid == this->root_node) {
            continue;
        }
        (*tree)[full_parents_list[taxid]].push_back(taxid);
    }
}

void TaxonomyClsAnno::dfs_statistics(const TaxId node,
                                     const ChildrenList &tree,
                                     std::vector<TaxId> *tree_linearization) {
    this->node_to_linearization_idx[node] = tree_linearization->size();
    tree_linearization->push_back(node);
    uint32_t depth = 0;
    for (const TaxId &child : tree.at(node)) {
        dfs_statistics(child, tree, tree_linearization);
        tree_linearization->push_back(node);
        if (this->node_depth[child] > depth) {
            depth = this->node_depth[child];
        }
    }
    this->node_depth[node] = depth + 1;
}

void TaxonomyClsAnno::rmq_preprocessing(const std::vector<TaxId> &tree_linearization) {
    uint32_t num_rmq_rows = log2(tree_linearization.size()) + 1;

    this->rmq_data.resize(num_rmq_rows);
    for (uint32_t i = 0; i < num_rmq_rows; ++i) {
        this->rmq_data[i].resize(tree_linearization.size());
    }

    // Copy tree_linearization to rmq[0].
    for (uint32_t i = 0; i < tree_linearization.size(); ++i) {
        this->rmq_data[0][i] = tree_linearization[i];
    }

    // Delta represents the size of the RMQ's sliding window (always a power of 2).
    uint32_t delta = 1;
    for (uint32_t row = 1; row < num_rmq_rows; ++row) {
        for (uint32_t i = 0; i + delta < tree_linearization.size(); ++i) {
            // rmq_data[row][i] covers an interval of size delta=2^row and returns the node with the maximal depth among positions [i, i+2^row-1] in the linearization.
            // According to 'this->dfs_statistics()': node_depth[leaf] = 1 and node_depth[root] = maximum distance to a leaf.
            if (this->node_depth[this->rmq_data[row - 1][i]] >
                this->node_depth[this->rmq_data[row - 1][i + delta]]) {
                this->rmq_data[row][i] = this->rmq_data[row - 1][i];
            } else {
                this->rmq_data[row][i] = this->rmq_data[row - 1][i + delta];
            }
        }
        delta *= 2;
    }

    // Compute fast tables for log2 and pow2.
    this->fast_log2.resize(tree_linearization.size());
    this->fast_pow2.push_back(1);
    for (uint32_t i = 2; i < tree_linearization.size(); ++i) {
        this->fast_log2[i] = 1 + this->fast_log2[i/2];
        if (this->fast_log2[i] > this->fast_log2[i-1]) {
            this->fast_pow2.push_back(i);
        }
    }
}

TaxId TaxonomyClsAnno::assign_class(const std::string &sequence) const {
    std::cerr << "assign class not implemented " << sequence << "\n\n";
    return 0;
}

} // namespace annot
} // namespace mtg
