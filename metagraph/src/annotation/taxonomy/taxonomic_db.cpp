#include "taxonomic_db.hpp"
#include "tax_classifier.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <sdsl/int_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <vector>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "common/serialization.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/threads/threading.hpp"

#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

using NormalizedTaxId = TaxonomyDB::NormalizedTaxId;
using TaxId = TaxonomyDB::TaxId;

uint64_t TaxonomyDB::num_get_taxid_calls = 0;
uint64_t TaxonomyDB::num_get_taxid_calls_failed = 0;

void TaxonomyDB::dfs_statistics(const NormalizedTaxId node,
                                const ChildrenList &tree,
                                std::vector<NormalizedTaxId> *tree_linearization) {
    this->node_to_linearization_idx[node] = tree_linearization->size();
    tree_linearization->push_back(node);
    uint64_t depth = 0;
    for (const NormalizedTaxId &child : tree[node]) {
        dfs_statistics(child, tree, tree_linearization);
        tree_linearization->push_back(node);
        if (this->node_depth[child] > depth) {
            depth = this->node_depth[child];
        }
    }
    this->node_depth[node] = depth + 1;
}

void TaxonomyDB::read_tree(const std::string &tax_tree_filepath,
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
            logger->error("The Taxonomic Tree file contains empty lines. Please make sure that this file was not manually modified: {}",
                          tax_tree_filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("The Taxonomic tree filepath contains incomplete lines. Please make sure that this file was not manually modified: {}",
                          tax_tree_filepath);
            exit(1);
        }
        uint64_t act = static_cast<uint64_t>(std::stoull(parts[0]));
        uint64_t parent = static_cast<uint64_t>(std::stoull(parts[2]));
        full_parents_list[act] = parent;
        this->node_parent[act] = parent;
    }

    std::queue<TaxId> relevant_taxids;
    for (const pair<AccessionVersion, TaxId>  &it : this->label_taxid_map) {
        relevant_taxids.push(it.second);
    }

    // If this->label_taxid_map (the set of labels from where the reads were generated) is empty, add and process all the taxids from the tree.
    if (this->label_taxid_map.size() == 0) {
        for (auto it : full_parents_list) {
            relevant_taxids.push(it.first);
        }
    }
    assert(relevant_taxids.size());

    // We want to treat the case when "taxid==0" as a wildcard for not initialised value.
    // Thus, 'this->denormalized_taxid' and 'this->normalized_taxid' must start their indexing from 1.
    this->denormalized_taxid.push_back(0);

    uint64_t num_nodes = 0;
    uint64_t num_taxid_failed = 0; // num_taxid_failed used for logging only.
    while (relevant_taxids.size()) {
        const TaxId taxid = relevant_taxids.front();
        relevant_taxids.pop();
        if (this->normalized_taxid.count(taxid)) {
            continue;
        }
        if (!full_parents_list.count(taxid)) {
            num_taxid_failed += 1;
            continue;
        }

        this->normalized_taxid[taxid] = ++num_nodes;
        this->denormalized_taxid.push_back(taxid);

        // Check if the current taxid is the root.
        if (taxid == full_parents_list[taxid]) {
            this->root_node = this->normalized_taxid[taxid];
            continue;
        }
        relevant_taxids.push(full_parents_list[taxid]);
    }
    if (num_taxid_failed) {
        logger->warn("During the tax_tree_filepath {} parsing, {} taxids were not found out of {} evaluations.",
                     tax_tree_filepath, num_taxid_failed, num_nodes);
    }

    // Construct the output tree.
    (*tree).resize(num_nodes + 1);
    for (const pair<TaxId, NormalizedTaxId> &it : this->normalized_taxid) {
        TaxId taxid = it.first;
        if (this->normalized_taxid[taxid] == this->root_node) {
            continue;
        }
        (*tree)[this->normalized_taxid[full_parents_list[taxid]]].push_back(
                this->normalized_taxid[taxid]);
    }
    assert(this->normalized_taxid.size() + 1 == this->denormalized_taxid.size()); // denormalized_taxid contains the '0' wildcard.
}

void TaxonomyDB::rmq_preprocessing(const std::vector<NormalizedTaxId> &tree_linearization) {
    uint num_rmq_rows = log2(tree_linearization.size()) + 1;

    this->rmq_data.resize(num_rmq_rows);
    for (uint64_t i = 0; i < num_rmq_rows; ++i) {
        this->rmq_data[i].resize(tree_linearization.size());
    }

    // Copy tree_linearization to rmq[0].
    for (uint64_t i = 0; i < tree_linearization.size(); ++i) {
        this->rmq_data[0][i] = tree_linearization[i];
    }

    // Delta represents the size of the rmq's sliding window (it is a power of 2).
    uint64_t delta = 1;
    for (uint row = 1; row < num_rmq_rows; ++row) {
        for (uint64_t i = 0; i + delta < tree_linearization.size(); ++i) {
            // rmq_data[row][i] has size 2^row and is computed as the max between the 2 halves of size 2^{row-1}.
            // According to 'this->dfs_statistics()': node_depth[leaf] = 1 and node_depth[root] = MAX_dist_to_tip.
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
    for (uint64_t i = 2; i < tree_linearization.size(); ++i) {
        this->fast_log2[i] = 1 + this->fast_log2[i/2];
        if (this->fast_log2[i] > this->fast_log2[i-1]) {
            this->fast_pow2.push_back(i);
        }
    }
}

// TODO - decide if we want to keep only 'get_taxid_from_label()' function and delete 'get_accession_version_from_label()' and 'get_normalized_taxid()'.
std::string TaxonomyDB::get_accession_version_from_label(const std::string &label) {
    std::vector<std::string> label_parts = utils::split_string(label, "|");
    if (label_parts.size() <= 3 || label_parts[3] == "") {
        logger->error("Failed to get accession from fasta label. Please make sure that the labels in the annotation matrix are in the standard ncbi format.");
        exit(1);
    }
    return label_parts[3];
}

bool TaxonomyDB::get_normalized_taxid(const std::string accession_version,
                                      NormalizedTaxId *taxid) const {
    num_get_taxid_calls += 1;
    if (!this->label_taxid_map.count(accession_version)) {
        // accession_version not in 'this->label_taxid_map'
        num_get_taxid_calls_failed += 1;
        return false;
    }
    if (!this->normalized_taxid.count(this->label_taxid_map.at(accession_version))) {
        // taxid (corresponding to the current accession_version) not in the taxonomic tree.
        num_get_taxid_calls_failed += 1;
        return false;
    }

    *taxid = this->normalized_taxid.at(this->label_taxid_map.at(accession_version));
    return true;
}

bool TaxonomyDB::get_normalized_taxid_from_label(const std::string &label, NormalizedTaxId *taxid) const {
    TaxId act = static_cast<uint64_t>(std::stoull(utils::split_string(label, "|")[1]));
    if (!this->normalized_taxid.count(act)) {
        return false;
    }
    *taxid = this->normalized_taxid.at(act);
    return true;
}

// TODO improve this by parsing the compressed ".gz" version (or use https://github.com/pmenzel/taxonomy-tools)
void TaxonomyDB::read_label_taxid_map(const std::string &label_taxid_map_filepath,
                                      const tsl::hopscotch_set<AccessionVersion> &input_accessions) {
    std::ifstream f(label_taxid_map_filepath);
    if (!f.good()) {
        logger->error("Failed to open accession to taxid map table {}",
                      label_taxid_map_filepath);
        exit(1);
    }

    std::string line;
    getline(f, line);
    if (!utils::starts_with(line, "accession\taccession.version\ttaxid\t")) {
        logger->error("The accession to taxid map table is not in the standard (*.accession2taxid) format {}.", label_taxid_map_filepath);
        exit(1);
    }

    while (getline(f, line)) {
        if (line == "") {
            logger->error("The accession to taxid map table contains empty lines. Please make sure that this file was not manually modified {}", label_taxid_map_filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("The accession to taxid map table contains incomplete lines. Please make sure that this file was not manually modified {}", label_taxid_map_filepath);
            exit(1);
        }
        if (input_accessions.count(parts[1])) {
            this->label_taxid_map[parts[1]] = static_cast<TaxId>(std::stoull(parts[2]));
        }
    }
}

TaxonomyDB::TaxonomyDB(const std::string &tax_tree_filepath,
                       const std::string &label_taxid_map_filepath,
                       const tsl::hopscotch_set<AccessionVersion> &input_accessions) {
    if (!std::filesystem::exists(tax_tree_filepath)) {
        logger->error("Can't open taxonomic tree file {}.", tax_tree_filepath);
        std::exit(1);
    }
    Timer timer;
    if (label_taxid_map_filepath != "") {
        logger->trace("Parsing label_taxid_map file..");
        read_label_taxid_map(label_taxid_map_filepath, input_accessions);
        logger->trace("Finished label_taxid_map file in {}s", timer.elapsed());
    }

    timer.reset();
    logger->trace("Parsing taxonomic tree..");
    TaxonomyDB::ChildrenList tree;
    read_tree(tax_tree_filepath, &tree);
    logger->trace("Finished taxonomic tree read in {}s", timer.elapsed());

    timer.reset();
    logger->trace("Calculating tree statistics..");
    std::vector<NormalizedTaxId> tree_linearization;
    this->node_to_linearization_idx.resize(tree.size());
    this->node_depth.resize(tree.size());
    dfs_statistics(root_node, tree, &tree_linearization);
    logger->trace("Finished tree statistics calculation in {}s", timer.elapsed());

    timer.reset();
    logger->trace("Starting rmq preprocessing..");
    rmq_preprocessing(tree_linearization);
    logger->trace("Finished rmq preprocessing in {}s", timer.elapsed());
}

NormalizedTaxId TaxonomyDB::find_normalized_lca(const std::vector<NormalizedTaxId> &taxids) const {
    if (taxids.empty()) {
        logger->error("Internal error: Can't find LCA for an empty set of normalized taxids.");
        std::exit(1);
    }
    uint64_t left_idx = this->node_to_linearization_idx[taxids[0]];
    uint64_t right_idx = this->node_to_linearization_idx[taxids[0]];

    // The node with maximum node_depth in 'linearization[left_idx : right_idx+1]' is the LCA of the given set.

    for (const NormalizedTaxId &taxid : taxids) {
        if (this->node_to_linearization_idx[taxid] < left_idx) {
            left_idx = this->node_to_linearization_idx[taxid];
        }
        if (this->node_to_linearization_idx[taxid] > right_idx) {
            right_idx = this->node_to_linearization_idx[taxid];
        }
    }

    // Find the maximum node_depth between the 2 overlapping intervals of size 2^log_dist.
    uint64_t log_dist = this->fast_log2[right_idx - left_idx];
    uint64_t left_lca = this->rmq_data[log_dist][left_idx];
    uint64_t right_lca =
            this->rmq_data[log_dist][right_idx - this->fast_pow2[log_dist] + 1];

    if (this->node_depth[left_lca] > this->node_depth[right_lca]) {
        return left_lca;
    }
    return right_lca;
}

void TaxonomyDB::kmer_to_taxid_map_update(const annot::MultiLabelEncoded<std::string> &annot) {
    std::vector<std::string> all_labels = annot.get_all_labels();
    for (const std::string &label : all_labels) {
        AccessionVersion accession_version = get_accession_version_from_label(label);
        NormalizedTaxId taxid;
        if (!get_normalized_taxid(accession_version, &taxid)) {
            continue;
        }
        annot.call_objects(label, [&](const KmerId &index) {
            if (index <= 0) {
                return;
            }
            if (index >= this->taxonomic_map.size()) {
                // Double the size of this->taxonomic_map until the current 'index' fits inside.
                uint64_t new_size = this->taxonomic_map.size();
                while (index >= new_size) {
                    new_size = new_size * 2 + 1;
                }
                // Create a aux_taxonomic_map because taxonomic_map.resize() does not preserve the data.
                sdsl::int_vector<> aux_taxonomic_map = this->taxonomic_map;
                this->taxonomic_map.resize(new_size);
                for (uint64_t i = 0; i < this->taxonomic_map.size(); ++i) {
                    this->taxonomic_map[i] = 0;
                }
                for (uint64_t i = 0; i < aux_taxonomic_map.size(); ++i) {
                    this->taxonomic_map[i] = aux_taxonomic_map[i];
                }
            }

            if (this->taxonomic_map[index] == 0) {
                this->taxonomic_map[index] = taxid;
            } else {
                this->taxonomic_map[index] = find_normalized_lca(
                        std::vector<NormalizedTaxId>{this->taxonomic_map[index], taxid});
          }
        });
    }
}

void TaxonomyDB::export_to_file(const std::string &filepath) {
    if (num_get_taxid_calls_failed) {
        logger->warn("During the annotation matrix parsing, {} taxids were not found out of a total of {} evaluations.",
                     num_get_taxid_calls_failed, num_get_taxid_calls);
    }

    Timer timer;
    logger->trace("Exporting metagraph taxonomic data..");

    std::ofstream f(filepath.c_str(), std::ios::out | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file {}.", filepath.c_str());
        std::exit(1);
    }
    serialize_number_number_map(f, this->node_parent);
    std::vector<std::pair<uint64_t, uint64_t>> taxid_frequencies(1);

    // Compute taxid_frequencies.
    for (size_t i = 0; i < this->taxonomic_map.size(); ++i) {
        if (this->taxonomic_map[i] == 0) {
            continue;
        }
        uint64_t taxid = this->taxonomic_map[i];
        while (taxid >= taxid_frequencies.size()) {
            taxid_frequencies.resize(2 * taxid_frequencies.size());
        }
        ++taxid_frequencies[taxid].first;
        taxid_frequencies[taxid].second = taxid;
    }
    // Sort the taxids in descending frequencies order for better compression.
    std::sort(taxid_frequencies.begin(), taxid_frequencies.end(),
        [](const std::pair<uint64_t, uint64_t> &a, const std::pair<uint64_t, uint64_t> &b){
            if (a.first != b.first) {
                return a.first > b.first;
            }
            return a.second > b.second;
        });

    // Construct 'code_to_taxid' and 'code' such that 'taxonomic_map[i]=code_to_taxid[code[i]]'.

    sdsl::int_vector<> code_to_taxid(taxid_frequencies.size());
    for (uint64_t i = 0; i < taxid_frequencies.size(); ++i) {
        if (taxid_frequencies[i].first == 0) {
            break;
        }
        code_to_taxid[i] = taxid_frequencies[i].second;
    }

    tsl::hopscotch_map<TaxId, TaxId> taxid_to_code;
    for (uint64_t i = 0; i < code_to_taxid.size(); ++i) {
        taxid_to_code[code_to_taxid[i]] = i;
    }

    sdsl::int_vector<> code_vector(this->taxonomic_map.size());
    uint64_t code_vector_cnt = 0;
    for (TaxId taxid : this->taxonomic_map) {
        code_vector[code_vector_cnt++] = taxid_to_code[taxid];
    }
    sdsl::dac_vector_dp<sdsl::rrr_vector<>> code(code_vector);

    // Denormalize code_to_taxid.
    for (uint64_t i = 0; i < code_to_taxid.size(); ++i) {
        code_to_taxid[i] = this->denormalized_taxid[code_to_taxid[i]];
    }
    code_to_taxid.serialize(f);

    code.serialize(f);
    serialize_number_number_map(f, node_parent);

    logger->trace("Finished exporting metagraph taxonomic data after {}s",
                  timer.elapsed());
}

TaxId TaxonomyDB::assign_class_getrows(const graph::AnnotatedDBG &anno,
                                       const std::string &sequence,
                                       const double lca_coverage_rate,
                                       const double kmers_discovery_rate) const {
    std::vector<NormalizedTaxId> curr_taxids;
    auto callback_cell = [&](const std::string &label) {
        NormalizedTaxId taxid;
        if(this->get_normalized_taxid_from_label(label, &taxid)) {
            curr_taxids.push_back(taxid);
        }
    };


    // num_kmers represents the total number of kmers in one read (read.len - graph.k + 1).
    uint64_t num_kmers = 0;

    // 'poz_forward' and 'forward_kmers' are storing the indexes and values of all the nonzero kmers in the given read.
    // 'forward_kmers' will be further sent to "matrix.getrows()" method;
    // 'poz_forward' will be used to associate one row from "matrix.getrows()" with the corresponding kmer index.
    vector<uint64_t> poz_forward;
    std::vector<node_index> forward_kmers;
    anno.get_graph_ptr()->map_to_nodes(sequence, [&](node_index i) {
        num_kmers++;
        if (i <= 0 || i >= anno.get_graph_ptr()->max_index()) {
            return;
        }
        forward_kmers.push_back(i - 1);
        poz_forward.push_back(num_kmers - 1);
    });

    // Compute the LCA normalized taxid for each nonzero kmer in the given read.
    std::vector<uint64_t> forward_taxids(num_kmers);

    uint64_t cnt = 0;
    anno.call_annotated_rows(forward_kmers, callback_cell,[&]() {
        if (curr_taxids.size() != 0) {
            forward_taxids[poz_forward[cnt++]] = this->denormalized_taxid[find_normalized_lca(curr_taxids)];
        }
        curr_taxids.clear();
    });

    std::string reversed_sequence = sequence;
    reverse_complement(reversed_sequence.begin(), reversed_sequence.end());

    // 'poz_backward' and 'backward_kmers' are storing the indexes and values of all the nonzero kmers in the given read reversed.
    // 'backward_kmers' will be further sent to "matrix.getrows()" method;
    // 'poz_backward' will be used to associate one row from "matrix.getrows()" with the corresponding kmer index.
    std::vector<node_index> backward_kmers;
    vector<uint64_t> poz_backward;
    cnt = 0;
    anno.get_graph_ptr()->map_to_nodes(reversed_sequence, [&](node_index i) {
        cnt++;
        if (i <= 0 || i >= anno.get_graph_ptr()->max_index()) {
            return;
        }
        backward_kmers.push_back(i - 1);
        poz_backward.push_back(cnt - 1);
    });

    // Compute the LCA normalized taxid for each nonzero kmer in the given read reversed.
    std::vector<uint64_t> backward_taxids(num_kmers);
    cnt = 0;

    anno.call_annotated_rows(backward_kmers, callback_cell, [&]() {
        if (curr_taxids.size() != 0) {
            // We will want to compare backward_taxids[i] with forward_taxids[i],
            // then we need to compute 'num_kmers - 1 - poz_backward[cnt++]'.
            backward_taxids[num_kmers - 1 - poz_backward[cnt++]] =
                    this->denormalized_taxid[find_normalized_lca(curr_taxids)];
        }
        curr_taxids.clear();
    });

    tsl::hopscotch_map<TaxId, uint64_t> num_kmers_per_node;

    // total_discovered_kmers represents the number of nonzero kmers according to both forward and reversed read.
    uint64_t total_discovered_kmers = 0;

    // Find the LCA taxid for each kmer without any dependency on the orientation of the read.
    for (uint64_t i = 0; i < backward_taxids.size(); ++i) {
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
                curr_taxid = TaxClassifier::find_lca(forward_taxid, backward_taxid, denormalized_taxid[this->root_node], this->node_parent);
                // todo - decide if it is worth to compute the normalization for all taxids. If not, we can delete one of  TaxClassifier::find_lca() or this.find_normalized_lca().
            }
        }
        if (curr_taxid) {
            total_discovered_kmers += 1;
            num_kmers_per_node[curr_taxid]++;
        }
    }

    if (total_discovered_kmers <= kmers_discovery_rate * num_kmers) {
        return 0; // 0 is a wildcard for not enough discovered kmers.
    }

    tsl::hopscotch_set<TaxId> nodes_already_propagated;
    tsl::hopscotch_map<TaxId, uint64_t> node_scores;

    uint64_t desired_number_kmers = total_discovered_kmers * lca_coverage_rate;
    TaxId best_lca = denormalized_taxid[this->root_node];
    uint64_t best_lca_dist_to_root = 1;

    // Update the nodes' score by iterating through all the nodes with nonzero kmers.
    for (const pair<TaxId, uint64_t> &node_pair : num_kmers_per_node) {
        TaxId start_node = node_pair.first;
        TaxClassifier::update_scores_and_lca(start_node, num_kmers_per_node, desired_number_kmers,
                              denormalized_taxid[this->root_node], this->node_parent,
                              &node_scores, &nodes_already_propagated, &best_lca,
                              &best_lca_dist_to_root);
    }
    return best_lca;
}


TaxId TaxonomyDB::assign_class_toplabels(const graph::AnnotatedDBG &anno,
                                         const std::string &sequence,
                                         const double label_fraction
                                         ) const {
    // Get all the labels with a frequency better than 'label_fraction' among the kmers in the forward read.
    std::vector<std::string> labels_discovered = anno.get_labels(sequence, label_fraction);

    std::string reversed_sequence = sequence;
    reverse_complement(reversed_sequence.begin(), reversed_sequence.end());
    // Get all the labels with a frequency better than 'label_fraction' among the kmers in the reversed read.
    std::vector<std::string> labels_discovered_rev = anno.get_labels(reversed_sequence, label_fraction);

    // Usually, only one of the two sets ('labels_discovered', 'labels_discovered_rev') will be nonempty.

    std::vector<NormalizedTaxId> curr_taxids;
    for (uint64_t i = 0; i < labels_discovered.size(); ++i) {
        TaxId act;
        if(this->get_normalized_taxid_from_label(labels_discovered[i], &act)) {
            curr_taxids.push_back(act);
        }
    }
    for (uint64_t i = 0; i < labels_discovered_rev.size(); ++i) {
        TaxId act;
        if(this->get_normalized_taxid_from_label(labels_discovered_rev[i], &act)) {
            curr_taxids.push_back(act);
        }
    }
    if (curr_taxids.size() == 0) {
        return 0; // Wildcard for not being able to assign a taxid.
    }
    return this->denormalized_taxid[find_normalized_lca(curr_taxids)];
}

} // namespace annot
} // namespace mtg
