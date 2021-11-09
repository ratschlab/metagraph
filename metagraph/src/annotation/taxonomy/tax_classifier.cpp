#include "tax_classifier.hpp"

#include <string>
#include <vector>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "common/unix_tools.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

bool TaxonomyBase::get_taxid_from_label(const std::string &label, TaxId *taxid) const {
    if (label_type_ == TAXID) {
        *taxid = std::stoul(utils::split_string(label, "|")[1]);
        return true;
    } else if (label_type_ == GEN_BANK) {
        auto it_acc_version_taxid = accversion_to_taxid_map_.find(get_accession_version_from_label(label));
        if (it_acc_version_taxid == accversion_to_taxid_map_.end()) {
            return false;
        }
        *taxid = it_acc_version_taxid->second;
        return true;
    }

    logger->error("Error: Could not get the taxid for label {}", label);
    exit(1);
}

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

    uint32_t num_taxid_failed = 0; // num_taxid_failed is used for logging only.
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

std::vector<TaxId> TaxonomyClsAnno::get_lca_taxids_for_seq(const std::string_view &sequence, bool reversed) const {
    // num_kmers represents the total number of kmers.
    uint32_t num_kmers = 0;

    // 'kmer_idx' and 'kmer_val' are storing the indexes and values of all the nonzero kmers in the given read.
    // The list of kmers (kmer_val) will be further sent to "matrix.getrows()" method;
    // The list of indexes (kmer_idx) will be used to link each row from "matrix.getrows()" to the corresponding kmer index.
    std::vector<uint32_t> kmer_idx;
    std::vector<node_index> kmer_val;

    if (sequence.size() >= std::numeric_limits<uint32_t>::max()) {
        logger->error("Error: The given sequence contains more than 2^32 bp.");
        exit(1);
    }

    std::shared_ptr<const graph::SequenceGraph> anno_graph = anno_matrix_->get_graph_ptr();
    anno_graph->map_to_nodes(sequence, [&](node_index i) {
        num_kmers++;
        if (!i || i >= anno_graph->max_index())
            return;
        kmer_val.push_back(i - 1);
        kmer_idx.push_back(num_kmers - 1);
    });

    // Compute the LCA taxid for each nonzero kmer in the given read.
    const auto unique_matrix_rows = anno_matrix_->get_annotation().get_matrix().get_rows(kmer_val);
    //TODO make sure that this function works even if we have duplications in 'rows'. Then, delete this error catch.
    if (kmer_val.size() != unique_matrix_rows.size()) {
        throw std::runtime_error("Error: The current implementation doesn't work in case of multiple occurrences"
                                 " of the same kmer in one read.");
    }

    if (unique_matrix_rows.size() >= std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("Error: There must be less than 2^32 unique rows in one anno matrix query. "
                                 "Please reduce the query batch size.");
    }
    const auto &label_encoder = anno_matrix_->get_annotation().get_label_encoder();

    TaxId taxid;
    uint32_t curr_kmer_identifier = 0;
    std::vector<TaxId> curr_kmer_taxids;
    std::vector<TaxId> lca_taxids(num_kmers);

    for (auto row : unique_matrix_rows) {
        for (auto cell : row) {
            if (get_taxid_from_label(label_encoder.decode(cell), &taxid)) {
                curr_kmer_taxids.push_back(taxid);
            }
        }
        if (curr_kmer_taxids.size() != 0) {
            if (not reversed) {
                lca_taxids[kmer_idx[curr_kmer_identifier]] = find_lca(curr_kmer_taxids);
            } else {
                lca_taxids[num_kmers - 1 - kmer_idx[curr_kmer_identifier]] = find_lca(curr_kmer_taxids);
            }
        }
        curr_kmer_identifier++;
        curr_kmer_taxids.clear();
    }

    return lca_taxids;
}

TaxId TaxonomyBase::assign_class(const std::string &sequence) const {
    std::vector<TaxId> forward_taxids = get_lca_taxids_for_seq(sequence, false);

    std::string reversed_sequence(sequence);
    reverse_complement(reversed_sequence.begin(), reversed_sequence.end());
    std::vector<TaxId> backward_taxids = get_lca_taxids_for_seq(reversed_sequence, true);

    tsl::hopscotch_map<TaxId, uint32_t> num_kmers_per_taxid;

    // num_discovered_kmers represents the number of nonzero kmers according to the forward and/or reversed read.
    uint32_t num_discovered_kmers = 0;
    // num_total_kmers is equal to the total number of zero/nonzero kmers in the read: size_read - k + 1
    uint32_t num_total_kmers = forward_taxids.size();

    // Find the LCA taxid for each kmer without considering the orientation of the read.
    // In case that both forward and reversed read have a nonzero kmer, then assign the LCA of those 2 kmers.
    for (uint32_t i = 0; i < num_total_kmers; ++i) {
        if (forward_taxids[i] == 0 && backward_taxids[i] == 0) {
            continue;
        }
        TaxId curr_taxid;
        if (backward_taxids[i] == 0) {
            curr_taxid = forward_taxids[i];
        } else if (forward_taxids[i] == 0) {
            curr_taxid = backward_taxids[i];
        } else {
            // In case that both 'forward_taxid[i]' and 'backward_taxids[i]' are nonzero, compute the LCA.
            TaxId forward_taxid = forward_taxids[i];
            TaxId backward_taxid = backward_taxids[i];
            if (forward_taxid == 0) {
                curr_taxid = backward_taxid;
            } else if (backward_taxid == 0) {
                curr_taxid = forward_taxid;
            } else {
                curr_taxid = find_lca({forward_taxid, backward_taxid});
            }
        }
        if (curr_taxid) {
            num_discovered_kmers ++;
            num_kmers_per_taxid[curr_taxid]++;
        }
    }

    if (num_discovered_kmers <= kmers_discovery_rate_ * num_total_kmers) {
        return 0; // 0 is a wildcard for not enough discovered kmers.
    }

    tsl::hopscotch_set<TaxId> taxid_already_propagated;
    tsl::hopscotch_map<TaxId, uint32_t> taxid_scores;

    uint32_t min_required_kmers = num_discovered_kmers * lca_coverage_rate_;
    TaxId best_lca = root_node_;
    uint32_t best_lca_dist_to_root = 1;

    // Update the nodes' score by iterating through all the nodes with nonzero kmers.
    for (const auto &[taxid, _] : num_kmers_per_taxid) {
        this->update_scores_and_lca(taxid, num_kmers_per_taxid, min_required_kmers, &taxid_scores,
                                    &taxid_already_propagated, &best_lca, &best_lca_dist_to_root);
    }
    return best_lca;
}

void TaxonomyBase::update_scores_and_lca(TaxId start_taxid,
                                         const tsl::hopscotch_map<TaxId, uint32_t> &num_kmers_per_taxid,
                                         uint32_t min_required_kmers,
                                         tsl::hopscotch_map<TaxId, uint32_t> *taxid_scores,
                                         tsl::hopscotch_set<TaxId> *taxid_already_propagated,
                                         TaxId *best_lca,
                                         uint32_t *best_lca_dist_to_root) const {
    if (taxid_already_propagated->count(start_taxid)) {
        return;
    }
    uint32_t score_from_processed_parents = 0;
    uint32_t score_from_unprocessed_parents = num_kmers_per_taxid.at(start_taxid);

    // processed_parents represents the set of nodes on the path start_taxid->root that have already been processed in the previous iterations.
    std::vector<TaxId> processed_parents;
    std::vector<TaxId> unprocessed_parents;

    TaxId act_node = start_taxid;
    unprocessed_parents.push_back(act_node);

    while (act_node != root_node_) {
        act_node = node_parent_.at(act_node);
        auto num_kmers_act_node_it = num_kmers_per_taxid.find(act_node);
        if (!taxid_already_propagated->count(act_node)) {
            if (num_kmers_act_node_it != num_kmers_per_taxid.end()) {
                score_from_unprocessed_parents += num_kmers_act_node_it->second;
            }
            unprocessed_parents.push_back(act_node);
        } else {
            if (num_kmers_act_node_it != num_kmers_per_taxid.end()) {
                score_from_processed_parents += num_kmers_act_node_it->second;
            }
            processed_parents.push_back(act_node);
        }
    }
    // The score of all the nodes in 'processed_parents' will be updated with 'score_from_unprocessed_parents' only.
    // The nodes in 'unprocessed_parents' will be updated with the sum 'score_from_processed_parents + score_from_unprocessed_parents'.
    for (uint32_t i = 0; i < unprocessed_parents.size(); ++i) {
        TaxId &act_node = unprocessed_parents[i];
        (*taxid_scores)[act_node] =
                score_from_processed_parents + score_from_unprocessed_parents;
        taxid_already_propagated->insert(act_node);

        uint32_t act_dist_to_root =
                processed_parents.size() + unprocessed_parents.size() - i;

        // Test if the current node's score would be a better LCA result.
        if ((*taxid_scores)[act_node] >= min_required_kmers
            && (act_dist_to_root > *best_lca_dist_to_root
                || (act_dist_to_root == *best_lca_dist_to_root
                    && (*taxid_scores)[act_node] > (*taxid_scores)[*best_lca])
            )
                ) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
    for (uint32_t i = 0; i < processed_parents.size(); ++i) {
        TaxId &act_node = processed_parents[i];
        (*taxid_scores)[act_node] += score_from_unprocessed_parents;

        uint32_t act_dist_to_root = processed_parents.size() - i;
        if ((*taxid_scores)[act_node] >= min_required_kmers
            && (act_dist_to_root > *best_lca_dist_to_root
                || (act_dist_to_root == *best_lca_dist_to_root
                    && (*taxid_scores)[act_node] > (*taxid_scores)[*best_lca])
            )
                ) {
            *best_lca = act_node;
            *best_lca_dist_to_root = act_dist_to_root;
        }
    }
}

TaxId TaxonomyClsAnno::find_lca(const std::vector<TaxId> &taxids) const {
    if (taxids.empty()) {
        logger->error("Error: Can't find LCA for an empty set of normalized taxids.");
        exit(1);
    }
    uint32_t left_idx = node_to_linearization_idx_.at(taxids[0]);
    uint32_t right_idx = node_to_linearization_idx_.at(taxids[0]);

    for (const TaxId &taxid : taxids) {
        uint32_t curr_idx = node_to_linearization_idx_.at(taxid);
        if (curr_idx < left_idx) {
            left_idx = curr_idx;
        }
        if (curr_idx > right_idx) {
            right_idx = curr_idx;
        }
    }
    // The node with maximum node_depth in 'linearization[left_idx : right_idx+1]' is the LCA of the given set.

    // Find the maximum node_depth between the 2 overlapping intervals of size 2^log_dist.
    uint32_t log_dist = sdsl::bits::hi(right_idx - left_idx);
    if (rmq_data_.size() <= log_dist) {
        logger->error("Error: the RMQ was not precomputed before the LCA queries.");
        exit(1);
    }

    uint32_t left_lca = rmq_data_[log_dist][left_idx];
    uint32_t right_lca = rmq_data_[log_dist][right_idx - (1 << log_dist) + 1];

    return node_depth_.at(left_lca) > node_depth_.at(right_lca) ? left_lca : right_lca;
}

TaxId TaxonomyClsAnno::assign_class_toplabels(const std::string &sequence, double label_fraction) const {
    // Get all the labels with a frequency higher than 'label_fraction' among the kmers in the forward read.
    std::vector<std::string> labels_discovered = anno_matrix_->get_labels(sequence, label_fraction);

    std::string reversed_sequence = sequence;
    reverse_complement(reversed_sequence.begin(), reversed_sequence.end());
    // Get all the labels with a frequency higher than 'label_fraction' among the kmers in the reversed read.
    std::vector<std::string> labels_discovered_rev = anno_matrix_->get_labels(reversed_sequence, label_fraction);

    // Usually, only one of the two sets ('labels_discovered', 'labels_discovered_rev') will be nonempty.

    std::vector<TaxId> curr_taxids;
    for (uint32_t i = 0; i < labels_discovered.size(); ++i) {
        TaxId act;
        if(get_taxid_from_label(labels_discovered[i], &act)) {
            curr_taxids.push_back(act);
        }
    }
    for (uint32_t i = 0; i < labels_discovered_rev.size(); ++i) {
        TaxId act;
        if(get_taxid_from_label(labels_discovered_rev[i], &act)) {
            curr_taxids.push_back(act);
        }
    }
    // Wildcard 0 for not being able to assign a taxid.
    return curr_taxids.size() == 0 ? 0 : find_lca(curr_taxids);
}

} // namespace annot
} // namespace mtg
