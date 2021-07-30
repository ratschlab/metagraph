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

std::string TaxonomyBase::get_accession_version_from_label(const std::string &label) const {
    switch (label_type_) {
        case TAXID:
            return utils::split_string(utils::split_string(label, "|")[2], " ")[0];
        case GEN_BANK:
            return utils::split_string(label, "|")[3];
    }

    logger->error("Error: Could not get the accession version for label {}", label);
    exit(1);
}

// TODO improve this by parsing the compressed ".gz" version (or use https://github.com/pmenzel/taxonomy-tools)
void TaxonomyBase::read_accversion_to_taxid_map(const std::string &filepath,
                                                const graph::AnnotatedDBG *anno_matrix) {
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
            // e.g. of nucl.accession2taxid file:
            //
            // A00001	A00001.1	10641	58418
            //
            // Thus, parts[1] represents the accession version and parts[2] the corresponding taxid.
            accversion_to_taxid_map_[parts[1]] = std::stoul(parts[2]);
        }
    }
}

TaxonomyClsAnno::TaxonomyClsAnno(const graph::AnnotatedDBG &anno,
                                 const std::string &tax_tree_filepath,
                                 double lca_coverage_rate,
                                 double kmers_discovery_rate,
                                 const std::string &label_taxid_map_filepath)
             : TaxonomyBase(lca_coverage_rate, kmers_discovery_rate),
               anno_matrix_(&anno) {
    if (!std::filesystem::exists(tax_tree_filepath)) {
        logger->error("Error: Can't open taxonomic tree file {}", tax_tree_filepath);
        exit(1);
    }

    // Take one sample label and find the label type.
    std::string sample_label = anno_matrix_->get_annotation().get_all_labels()[0];

    if (utils::starts_with(sample_label, "gi|")) {
        // e.g.   >gi|1070643132|ref|NC_031224.1| Arthrobacter phage Mudcat, complete genome
        label_type_ = GEN_BANK;
    } else if (utils::starts_with(utils::split_string(sample_label, ":")[1], "taxid|")) {
        // e.g.   >kraken:taxid|2016032|NC_047834.1 Alteromonas virus vB_AspP-H4/4, complete genome
        label_type_ = TAXID;
    } else {
        logger->error("Error: Can't determine the type of the given label {}. "
                      "Make sure the labels are in a recognized format.", sample_label);
        exit(1);
    }

    Timer timer;
    if (label_type_ == GEN_BANK) {
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
        // e.g. of nodes.dmp file:
        //
        // 2	|	131567	|	superkingdom	|		|	0	|	0
        //
        // Thus, parts[0] represents the child taxid and parts[2] the parent taxid.
        node_parent_[std::stoul(parts[0])] = std::stoul(parts[2]);
    }

    std::vector<TaxId> relevant_taxids;
    // 'considered_relevant_taxids' is used to make sure that there are no duplications in 'relevant_taxids'.
    tsl::hopscotch_set<TaxId> considered_relevant_taxids;

    if (accversion_to_taxid_map_.size()) {
        // Store only the taxonomic nodes that exists in the annotation matrix.
        for (const auto &[_, taxid] : accversion_to_taxid_map_) {
            relevant_taxids.push_back(taxid);
            considered_relevant_taxids.insert(taxid);
        }
    } else {
        // If 'this->accversion_to_taxid_map' is empty, store the entire taxonomic tree.
        for (const auto &[child, _] : node_parent_) {
            relevant_taxids.push_back(child);
            considered_relevant_taxids.insert(child);
        }
    }
    assert(relevant_taxids.size());

    uint64_t num_taxid_failed = 0; // num_taxid_failed is used for logging only.
    for (uint32_t i = 0; i < relevant_taxids.size(); ++i) {
        TaxId taxid = relevant_taxids[i];
        auto it_taxid_parent = node_parent_.find(taxid);
        if (it_taxid_parent == node_parent_.end()) {
            num_taxid_failed += 1;
            continue;
        }
        TaxId taxid_parent = it_taxid_parent->second;

        if (not considered_relevant_taxids.count(taxid_parent)) {
            relevant_taxids.push_back(taxid_parent);
            considered_relevant_taxids.insert(taxid_parent);
        }

        // Check if the current taxid is the root.
        if (taxid == taxid_parent) {
            root_node_ = taxid;
        }
    }
    if (num_taxid_failed) {
        logger->warn("During the tax_tree_filepath {} parsing, {} taxids were not found out of {} total evaluations",
                     tax_tree_filepath, num_taxid_failed, relevant_taxids.size());
    }

    // Construct the output tree.
    for (const TaxId &taxid : relevant_taxids) {
        if (taxid == root_node_)
            continue;
        auto it_taxid_parent = node_parent_.find(taxid);
        if (it_taxid_parent != node_parent_.end()) {
            (*tree)[it_taxid_parent->second].push_back(taxid);
        }
    }
}

void TaxonomyClsAnno::dfs_statistics(TaxId node,
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

} // namespace annot
} // namespace mtg
