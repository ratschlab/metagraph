#include "taxonomic_db.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <sdsl/int_vector.hpp>
#include <vector>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "common/serialization.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace annot {

using mtg::common::logger;

typedef TaxonomyDB::NormalizedTaxId NormalizedTaxId;

uint64_t TaxonomyDB::num_get_taxid_calls = 0;
uint64_t TaxonomyDB::num_get_taxid_calls_failed = 0;

void TaxonomyDB::dfs_statistics(const NormalizedTaxId &node, const ChildrenList &tree,
                                std::vector<NormalizedTaxId> &tree_linearization) {
    node_to_linearization_idx[node] = tree_linearization.size();
    tree_linearization.push_back(node);
    uint64_t depth = 0;
    for(const auto &child: tree[node]) {
        dfs_statistics(child, tree, tree_linearization);
        tree_linearization.push_back(node);
        if (node_depth[child] > depth) {
            depth = node_depth[child];
        }
    }
    node_depth[node] = depth + 1;
}

void TaxonomyDB::read_tree(const std::string &taxo_tree_filepath,
                           ChildrenList &tree, NormalizedTaxId &root_node) {
    std::ifstream f(taxo_tree_filepath);
    std::string line;

    tsl::hopscotch_map<TaxId, TaxId> full_parents_list;
    while (getline(f, line) ) {
        std::vector<std::string> parts = utils::split_string(line, "\t");
        uint64_t act = static_cast<uint64_t>(std::stoull(parts[0]));
        uint64_t parent = static_cast<uint64_t>(std::stoull(parts[2]));
        full_parents_list[act] = parent;
    }
    f.close();

    std::queue<TaxId> relevant_taxids;
    for (const auto &it: lookup_table) {
        relevant_taxids.push(it.second);
    }

    uint64_t num_nodes = 0;
    // num_taxid_failed used for logging only.
    uint64_t num_taxid_failed = 0;
    while (relevant_taxids.size()) {
        const TaxId taxid = relevant_taxids.front();
        relevant_taxids.pop();
        if (normalized_taxid.count(taxid)) {
            continue;
        }
        if (! full_parents_list.count(taxid)) {
            logger->warn("Taxid {} cannot be found in the taxonomic tree.", taxid);
            num_taxid_failed += 1;
            continue;
        }
        normalized_taxid[taxid] = num_nodes++;
        denormalized_taxid.push_back(taxid);
        if (taxid == full_parents_list[taxid]) {
            root_node = normalized_taxid[taxid];
            continue;
        }
        relevant_taxids.push(full_parents_list[taxid]);
    }
    if (num_taxid_failed) {
        logger->warn("Could not find {} taxid out of the total number of {} evaluated taxids.", num_taxid_failed, num_nodes);
    }

    tree.resize(num_nodes);
    for (const auto &it: normalized_taxid) {
        TaxId taxid = it.first;
        if (normalized_taxid[taxid] == root_node) {
            continue;
        }
        tree[normalized_taxid[full_parents_list[taxid]]].push_back(
                normalized_taxid[taxid]);
    }
}

void TaxonomyDB::rmq_preprocessing(const std::vector<NormalizedTaxId> &tree_linearization) {
    uint num_rmq_rows = log2(tree_linearization.size()) + 1;
    rmq_data.resize(num_rmq_rows);
    for (uint64_t i = 0; i < num_rmq_rows; ++i) {
        rmq_data[i].resize(tree_linearization.size());
    }

    // Copy tree_linearization on rmq[0].
    for (uint64_t i = 0; i < tree_linearization.size(); ++i) {
        rmq_data[0][i] = tree_linearization[i];
    }

    uint delta = 1;
    for (uint row = 1; row < num_rmq_rows; ++row) {
        for (uint64_t i = 0; i + delta < tree_linearization.size(); ++i) {
            if (node_depth[rmq_data[row - 1][i]] > node_depth[rmq_data[row - 1][i + delta]]) {
                rmq_data[row][i] = rmq_data[row - 1][i];
            } else {
                rmq_data[row][i] = rmq_data[row - 1][i + delta];
            }
        }
        delta *= 2;
    }

    precalc_log2.resize(tree_linearization.size());
    precalc_pow2.push_back(1);
    for (uint64_t i = 2; i < tree_linearization.size(); ++i) {
        precalc_log2[i] = 1 + precalc_log2[i/2];
        if(precalc_log2[i] > precalc_log2[i-1]) {
            precalc_pow2.push_back(i);
        }
    }
}

std::string TaxonomyDB::get_accession_version_from_label(const std::string &label) {
    return  utils::split_string(label, "|")[3];
}

// TODO improve this by parsing the compressed ".gz" version (or use https://github.com/pmenzel/taxonomy-tools)
void TaxonomyDB::read_lookup_table(const std::string &lookup_table_filepath,
                                   const tsl::hopscotch_set<AccessionVersion> &input_accessions) {
    //   TODO Print an error when the file format (e.g. '\t') is broken.
    std::ifstream f(lookup_table_filepath);
    std::string line;
    while (getline(f, line) ) {
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if(input_accessions.count(parts[1])) {
            TaxId act = static_cast<TaxId>(std::stoull(parts[2]));
            lookup_table[parts[1]] = act;
        }
    }
    f.close();
}

TaxonomyDB::TaxonomyDB(const std::string &taxo_tree_filepath,
                       const std::string &lookup_table_filepath,
                       const tsl::hopscotch_set<AccessionVersion> &input_accessions) {

    if (!std::filesystem::exists(taxo_tree_filepath)) {
        logger->error("Can't open taxonomic tree file '{}'.", taxo_tree_filepath);
        std::exit(1);
    }
    if (!std::filesystem::exists(lookup_table_filepath)) {
        logger->error("Can't open taxonomic lookup table file '{}'.", lookup_table_filepath);
        std::exit(1);
    }

    Timer timer;
    logger->trace("Parsing lookup table..");
    read_lookup_table(lookup_table_filepath, input_accessions);
    logger->trace("Finished parsing tookup table in '{}' sec", timer.elapsed());

    timer.reset();
    logger->trace("Parsing taxonomic tree..");
    TaxonomyDB::ChildrenList tree;
    NormalizedTaxId root_node;
    read_tree(taxo_tree_filepath, tree, root_node);
    logger->trace("Finished parsing taxonomic tree in '{}' sec", timer.elapsed());

    timer.reset();
    logger->trace("Calculating tree statistics..");
    std::vector<NormalizedTaxId> tree_linearization;
    node_to_linearization_idx.resize(tree.size());
    node_depth.resize(tree.size());
    dfs_statistics(root_node, tree, tree_linearization);
    logger->trace("Finished calculating tree statistics in '{}' sec", timer.elapsed());

    timer.reset();
    logger->trace("Starting rmq preprocessing..");
    rmq_preprocessing(tree_linearization);
    logger->trace("Finished rmq preprocessing in '{}' sec", timer.elapsed());
}

NormalizedTaxId TaxonomyDB::find_lca(const std::vector<NormalizedTaxId> &taxids) const {
    if (taxids.empty()) {
        logger->error("Can't find LCA for an empty set of taxids.");
        std::exit(1);
    }
    uint64_t left_idx = node_to_linearization_idx[taxids[0]];
    uint64_t right_idx = node_to_linearization_idx[taxids[0]];
    for (const auto &taxid: taxids) {
        if (node_to_linearization_idx[taxid] < left_idx) {
            left_idx = node_to_linearization_idx[taxid];
        }
        if (node_to_linearization_idx[taxid] > right_idx) {
            right_idx = node_to_linearization_idx[taxid];
        }
    }

    uint64_t log_dist = precalc_log2[right_idx - left_idx];
    uint64_t left_lca = rmq_data[log_dist][left_idx];
    uint64_t right_lca = rmq_data[log_dist][right_idx - precalc_pow2[log_dist] + 1];

    if (node_depth[left_lca] > node_depth[right_lca]) {
        return left_lca;
    }
    return right_lca;
}

bool TaxonomyDB::get_normalized_taxid(const std::string accession_version, NormalizedTaxId &taxid) const {
    num_get_taxid_calls += 1;
    if (! lookup_table.count(accession_version)) {
        // accession_version not in the lookup_table
        num_get_taxid_calls_failed += 1;
        return false;
    }
    if (! normalized_taxid.count(lookup_table.at(accession_version))) {
        // taxid (corresponding to the current accession_version) not in the taxonomic tree.
        return false;
    }

    taxid = normalized_taxid.at(lookup_table.at(accession_version));
    return true;
}

void TaxonomyDB::kmer_to_taxid_map_update(const annot::MultiLabelEncoded<std::string> &annot,
                                          std::mutex &taxo_mutex) {
    for (const auto &label: annot.get_all_labels()) {
        AccessionVersion accession_version = get_accession_version_from_label(label);
        uint64_t taxid;
        if(!get_normalized_taxid(accession_version, taxid)) {
            continue;
        }
        annot.call_objects(label, [&](const auto &index) {
          if (index <= 0) {
              return;
          }
          if (index >= taxonomic_map.size()) {
              // Double the size of taxonomic_map until the current 'index' fits inside.
              taxo_mutex.lock();
              uint64_t new_size = taxonomic_map.size();
              while (index >= new_size) {
                  new_size = new_size * 2 + 1;
              }
              sdsl::int_vector<> aux_taxonomic_map = taxonomic_map;
              taxonomic_map.resize(new_size);
              for(uint64_t i = 0; i < taxonomic_map.size(); ++i) {
                  taxonomic_map[i] = 0;
              }
              for(uint64_t i = 0; i < aux_taxonomic_map.size(); ++i) {
                  taxonomic_map[i] = aux_taxonomic_map[i];
              }
              taxo_mutex.unlock();
          }

          if (taxonomic_map[index] == 0) {
              taxo_mutex.lock();
              taxonomic_map[index] = taxid;
              taxo_mutex.unlock();
          } else {
              auto lca = find_lca(std::vector<uint64_t>{taxonomic_map[index], taxid});
              taxo_mutex.lock();
              taxonomic_map[index] = lca;
              taxo_mutex.unlock();
          }
        });
    }
}

void TaxonomyDB::export_to_file(const std::string &filepath) {
    if (num_get_taxid_calls_failed) {
        logger->warn("Failed to get taxid {} times out of the total number of {} calls.",
                     num_get_taxid_calls_failed, num_get_taxid_calls);
    }

    Timer timer;
    logger->trace("Exporting metagraph taxonomic data..");

    std::ofstream f(filepath.c_str(), std::ios::out | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }

    const std::vector<NormalizedTaxId> &linearization = rmq_data[0];
    tsl::hopscotch_map<TaxId, TaxId> node_parent;

    node_parent[denormalized_taxid[linearization[0]]] = denormalized_taxid[linearization[0]];

    // Start the iteration from 1, because we are computing linearization[i - 1].
    for (uint64_t i = 1; i < linearization.size(); ++i) {
        uint64_t act = linearization[i];
        uint64_t prv = linearization[i - 1];
        if (node_parent.count(denormalized_taxid[act])) {
            // Prv is act's child.
            continue;
        }
        // Prv is act's parent.
        node_parent[denormalized_taxid[act]] = denormalized_taxid[prv];
    }
    serialize_number_number_map(f, node_parent);

    // Denormalize taxids in taxonomic map.
    for (size_t i = 0; i < taxonomic_map.size(); ++i) {
        if (taxonomic_map[i]) {
            taxonomic_map[i] = denormalized_taxid[taxonomic_map[i]];
        }
    }
    taxonomic_map.serialize(f);

    f.close();
    logger->trace("Finished exporting metagraph taxonomic data after '{}' sec", timer.elapsed());
}

}
}
