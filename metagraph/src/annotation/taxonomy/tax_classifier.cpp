#include "tax_classifier.hpp"

#include <string>
#include <vector>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

bool TaxonomyBase::get_taxid_from_label(const std::string &label, TaxId *taxid) const {
    if (label_type_ == TAXID) {
        *taxid = std::stoul(utils::split_string(label, "|")[1]);
        return true;
    } else if (label_type_ == GEN_BANK) {
        std::string acc_version = get_accession_version_from_label(label);
        auto it = accversion_to_taxid_map_.find(acc_version);
        if (it == accversion_to_taxid_map_.end()) {
            return false;
        }
        *taxid = it->second;
        return true;
    }

    logger->error("Error: Could not get the taxid for label {}", label);
    exit(1);
}

std::string TaxonomyBase::get_accession_version_from_label(const std::string &label) const {
    if (label_type_ == TAXID) {
        return utils::split_string(utils::split_string(label, "|")[2], " ")[0];
    } else if (label_type_ == GEN_BANK) {
        return utils::split_string(label, "|")[3];;
    }

    logger->error("Error: Could not get the accession version for label {}", label);
    exit(1);
}

// TODO improve this by parsing the compressed ".gz" version (or use https://github.com/pmenzel/taxonomy-tools)
void TaxonomyBase::read_accversion_to_taxid_map(const std::string &filepath,
                                                const graph::AnnotatedDBG *anno_matrix = NULL) {
    std::ifstream f(filepath);
    if (!f.good()) {
        logger->error("Error: Failed to open accession to taxid map table {}. \n"
                      "In the cases when the taxid is not specified in the label string, "
                      "the acc_version to taxid lookup table filepath must be given as a flag.", filepath);
        exit(1);
    }

    std::string line;
    getline(f, line);
    if (!utils::starts_with(line, "accession\taccession.version\ttaxid\t")) {
        logger->error("Error: The accession to taxid map table is not in the standard (*.accession2taxid) format {}",
                      filepath);
        exit(1);
    }

    tsl::hopscotch_set<std::string> input_accessions;
    if (anno_matrix != NULL) {
        for (const std::string &label : anno_matrix->get_annotation().get_all_labels()) {
            input_accessions.insert(get_accession_version_from_label(label));
        }
    }

    while (getline(f, line)) {
        if (line == "") {
            logger->error("Error: The accession to taxid map table contains empty lines. "
                          "Please make sure that this file was not manually modified {}", filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("Error: The accession to taxid map table contains incomplete lines. "
                          "Please make sure that this file was not manually modified {}", filepath);
            exit(1);
        }
        if (input_accessions.size() == 0 || input_accessions.count(parts[1])) {
            accversion_to_taxid_map_[parts[1]] = std::stoul(parts[2]);
        }
    }
}

TaxonomyClsAnno::TaxonomyClsAnno(const graph::AnnotatedDBG &anno,
                                 const double lca_coverage_rate,
                                 const double kmers_discovery_rate,
                                 const std::string &tax_tree_filepath,
                                 const std::string &label_taxid_map_filepath)
             : TaxonomyBase(lca_coverage_rate, kmers_discovery_rate),
               anno_matrix_(&anno) {
    if (!std::filesystem::exists(tax_tree_filepath)) {
        logger->error("Error: Can't open taxonomic tree file {}", tax_tree_filepath);
        exit(1);
    }

    // Take one sample label and find the label type.
    std::string sample_label = anno_matrix_->get_annotation().get_all_labels()[0];

    // If true, require_accversion_to_taxid_map means that the taxid is not mentioned as part of the label. Thus, the program
    // needs to parse an additional accession_version to taxid lookup file (optional argument: label_taxid_map_filepath).
    bool require_accversion_to_taxid_map;
    if (utils::starts_with(sample_label, "gi|")) {
        // e.g.   >gi|1070643132|ref|NC_031224.1| Arthrobacter phage Mudcat, complete genome
        label_type_ = GEN_BANK;
        require_accversion_to_taxid_map = true;
    } else if (utils::starts_with(utils::split_string(sample_label, ":")[1], "taxid|")) {
        // e.g.   >kraken:taxid|2016032|NC_047834.1 Alteromonas virus vB_AspP-H4/4, complete genome
        label_type_ = TAXID;
        require_accversion_to_taxid_map = false;
    } else {
        logger->error("Error: Can't determine the type of the given label {}. "
                      "Make sure the labels are in a recognized format.", sample_label);
        exit(1);
    }

    Timer timer;
    if (require_accversion_to_taxid_map) {
        logger->trace("Parsing label_taxid_map file...");
        read_accversion_to_taxid_map(label_taxid_map_filepath, anno_matrix_);
        logger->trace("Finished label_taxid_map file in {} sec", timer.elapsed());
    }

    timer.reset();
    logger->trace("Parsing taxonomic tree...");
    ChildrenList tree;
    read_tree(tax_tree_filepath, &tree);
    logger->trace("Finished taxonomic tree read in {} sec.", timer.elapsed());

    timer.reset();
    logger->trace("Calculating tree statistics...");
    std::vector<TaxId> tree_linearization;
    dfs_statistics(root_node_, tree, &tree_linearization);
    logger->trace("Finished tree statistics calculation in {} sec.", timer.elapsed());

    timer.reset();
    logger->trace("Starting rmq preprocessing...");
    rmq_preprocessing(tree_linearization);
    logger->trace("Finished rmq preprocessing in {} sec.", timer.elapsed());
}

void TaxonomyClsAnno::read_tree(const std::string &tax_tree_filepath, ChildrenList *tree) {
    std::ifstream f(tax_tree_filepath);
    if (!f.good()) {
        logger->error("Error: Failed to open Taxonomic Tree file {}", tax_tree_filepath);
        exit(1);
    }

    std::string line;
    tsl::hopscotch_map<TaxId, TaxId> full_parents_list;
    while (getline(f, line)) {
        if (line == "") {
            logger->error("Error: The Taxonomic Tree file contains empty lines. "
                          "Please make sure that this file was not manually modified: {}",
                          tax_tree_filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("Error: The Taxonomic tree filepath contains incomplete lines. "
                          "Please make sure that this file was not manually modified: {}",
                          tax_tree_filepath);
            exit(1);
        }
        uint32_t act = std::stoul(parts[0]);
        uint32_t parent = std::stoul(parts[2]);
        full_parents_list[act] = parent;
        node_parent_[act] = parent;
    }

    std::vector<TaxId> relevant_taxids;
    // 'considered_relevant_taxids' is used to make sure that there are no duplications in 'relevant_taxids'.
    tsl::hopscotch_set<TaxId> considered_relevant_taxids;

    if (accversion_to_taxid_map_.size()) {
        // Store only the taxonomic nodes that exists in the annotation matrix.
        for (const std::pair<std::string, TaxId> &pair : accversion_to_taxid_map_) {
            relevant_taxids.push_back(pair.second);
            considered_relevant_taxids.insert(pair.second);
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

        if (not considered_relevant_taxids.count(full_parents_list[taxid])) {
            relevant_taxids.push_back(full_parents_list[taxid]);
            considered_relevant_taxids.insert(full_parents_list[taxid]);
        }

        // Check if the current taxid is the root.
        if (taxid == full_parents_list[taxid]) {
            root_node_ = taxid;
        }
    }
    if (num_taxid_failed) {
        logger->warn("During the tax_tree_filepath {} parsing, {} taxids were not found out of {} total evaluations.",
                     tax_tree_filepath, num_taxid_failed, relevant_taxids.size());
    }

    // Construct the output tree.
    for (const TaxId &taxid : relevant_taxids) {
        if (taxid == root_node_) {
            continue;
        }
        (*tree)[full_parents_list[taxid]].push_back(taxid);
    }
}

void TaxonomyClsAnno::dfs_statistics(const TaxId node,
                                     const ChildrenList &tree,
                                     std::vector<TaxId> *tree_linearization) {
    node_to_linearization_idx_[node] = tree_linearization->size();
    tree_linearization->push_back(node);
    uint32_t depth = 0;

    auto it = tree.find(node);
    if (it != tree.end()) {
        for (const TaxId &child : it->second) {
            dfs_statistics(child, tree, tree_linearization);
            tree_linearization->push_back(node);
            if (node_depth_[child] > depth) {
                depth = node_depth_[child];
            }
        }
    }
    node_depth_[node] = depth + 1;
}

void TaxonomyClsAnno::rmq_preprocessing(const std::vector<TaxId> &tree_linearization) {
    uint32_t num_rmq_rows = sdsl::bits::hi(tree_linearization.size()) + 1;

    rmq_data_.resize(num_rmq_rows);
    for (uint32_t i = 0; i < num_rmq_rows; ++i) {
        rmq_data_[i].resize(tree_linearization.size());
    }

    // Copy tree_linearization to rmq[0].
    for (uint32_t i = 0; i < tree_linearization.size(); ++i) {
        rmq_data_[0][i] = tree_linearization[i];
    }

    // Delta represents the size of the RMQ's sliding window (always a power of 2).
    uint32_t delta = 1;
    for (uint32_t row = 1; row < num_rmq_rows; ++row) {
        for (uint32_t i = 0; i + delta < tree_linearization.size(); ++i) {
            // rmq_data[row][i] covers an interval of size delta=2^row and returns the node with the maximal depth
            // among positions [i, i+2^row-1] in the linearization.
            // According to 'this->dfs_statistics()':
            //     node_depth[leaf] = 1 and node_depth[root] = maximum distance to a leaf.

            if (node_depth_[rmq_data_[row - 1][i]] >
                node_depth_[rmq_data_[row - 1][i + delta]]) {
                rmq_data_[row][i] = rmq_data_[row - 1][i];
            } else {
                rmq_data_[row][i] = rmq_data_[row - 1][i + delta];
            }
        }
        delta *= 2;
    }
}

std::vector<TaxId> TaxonomyClsAnno::get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const {
    throw std::runtime_error("get_lca_taxids_for_seq TaxonomyClsAnno not implemented. Received seq size"
                             + to_string(sequence.size()) + to_string(reversed));
}

std::vector<TaxId> TaxonomyClsImportDB::get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const {
    throw std::runtime_error("get_lca_taxids_for_seq TaxonomyClsImportDB not implemented. Received seq size"
                             + to_string(sequence.size()) + to_string(reversed));
}

TaxId TaxonomyClsAnno::find_lca(const std::vector<TaxId> &taxids) const {
    throw std::runtime_error("find_lca TaxonomyClsAnno not implemented. Received taxids size"
                             + to_string(taxids.size()));
}

TaxId TaxonomyClsImportDB::find_lca(const std::vector<TaxId> &taxids) const {
    throw std::runtime_error("find_lca TaxonomyClsImportDB not implemented. Received taxids size"
                             + to_string(taxids.size()));
}

} // namespace annot
} // namespace mtg
