#include "taxonomic_db.hpp"

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
                           ChildrenList *tree,
                           NormalizedTaxId *root_node) {
    std::ifstream f(tax_tree_filepath);
    if (!f.good()) {
        logger->error("Failed to open Taxonomic Tree file {}", tax_tree_filepath);
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
    }

    std::queue<TaxId> relevant_taxids;
    for (const pair<AccessionVersion, TaxId>  &it : this->label_taxid_map) {
        relevant_taxids.push(it.second);
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
        if (taxid == full_parents_list[taxid]) {
            (*root_node) = this->normalized_taxid[taxid];
            continue;
        }
        relevant_taxids.push(full_parents_list[taxid]);
    }
    if (num_taxid_failed) {
        logger->warn("During the tax_tree_filepath {} parsing, there could not be found {} taxids associated to the given accession versions (total number of evaluations: {}).",
                     tax_tree_filepath, num_taxid_failed, num_nodes);
    }

    (*tree).resize(num_nodes + 1);
    for (const pair<TaxId, NormalizedTaxId> &it : this->normalized_taxid) {
        TaxId taxid = it.first;
        if (this->normalized_taxid[taxid] == *root_node) {
            continue;
        }
        (*tree)[this->normalized_taxid[full_parents_list[taxid]]].push_back(
                this->normalized_taxid[taxid]);
    }
    assert(this->normalized_taxid.size() + 1 == this->denormalized_taxid.size()); // denormalized_taxid contains a wildcard.
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

    uint64_t delta = 1;
    for (uint row = 1; row < num_rmq_rows; ++row) {
        for (uint64_t i = 0; i + delta < tree_linearization.size(); ++i) {
            if (this->node_depth[this->rmq_data[row - 1][i]] >
                this->node_depth[this->rmq_data[row - 1][i + delta]]) {
                this->rmq_data[row][i] = this->rmq_data[row - 1][i];
            } else {
                this->rmq_data[row][i] = this->rmq_data[row - 1][i + delta];
            }
        }
        delta *= 2;
    }

    this->fast_log2.resize(tree_linearization.size());
    this->fast_pow2.push_back(1);
    for (uint64_t i = 2; i < tree_linearization.size(); ++i) {
        this->fast_log2[i] = 1 + this->fast_log2[i/2];
        if (this->fast_log2[i] > this->fast_log2[i-1]) {
            this->fast_pow2.push_back(i);
        }
    }
}

// std::string TaxonomyDB::get_accession_version_from_label(const std::string &label) {
//     std::vector<std::string> label_parts = utils::split_string(label, "|");
//     if (label_parts.size() <= 3 || label_parts[3] == "") {
//         logger->error("Failed to get accession from fasta label. Please make sure that the labels in the annotation matrix are in the standard ncbi format.");
//         exit(1);
//     }
//     return label_parts[3];
// }

std::string TaxonomyDB::get_accession_version_from_label(const std::string &label) {
    auto aux = utils::split_string(label, "|")[2];
    return utils::split_string(aux, " ")[0];
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

    std::cerr << "input_accessions.size=" << input_accessions.size() << "\n\n\n";

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
        if (input_accessions.size() == 0 || input_accessions.count(parts[1])) {
            this->label_taxid_map[parts[1]] = static_cast<TaxId>(std::stoull(parts[2]));
        }
    }
}

TaxonomyDB::TaxonomyDB(const std::string &tax_tree_filepath,
                       const std::string &label_taxid_map_filepath,
                       const tsl::hopscotch_set<AccessionVersion> &input_accessions) {
//    if (!input_accessions.size()) {
//        logger->error("Can't construct TaxonomyDB for an empty set of accession versions.");
//        std::exit(1);
//    }
    if (!std::filesystem::exists(tax_tree_filepath)) {
        logger->error("Can't open taxonomic tree file {}.", tax_tree_filepath);
        std::exit(1);
    }
    // if (!std::filesystem::exists(label_taxid_map_filepath)) {
    //     logger->error("Can't open label_taxid_map file {}.", label_taxid_map_filepath);
    //     std::exit(1);
    // }

    Timer timer;
    if (label_taxid_map_filepath != "") {
        logger->trace("Parsing label_taxid_map file..");
        read_label_taxid_map(label_taxid_map_filepath, input_accessions);
        logger->trace("Finished label_taxid_map file in {}s", timer.elapsed());
    }

    timer.reset();
    logger->trace("Parsing taxonomic tree..");
    TaxonomyDB::ChildrenList tree;
    NormalizedTaxId root_node;
    read_tree(tax_tree_filepath, &tree, &root_node);
    logger->trace("Finished parsing taxonomic tree in {}s", timer.elapsed());

    timer.reset();
    logger->trace("Calculating tree statistics..");
    std::vector<NormalizedTaxId> tree_linearization;
    this->node_to_linearization_idx.resize(tree.size());
    this->node_depth.resize(tree.size());
    dfs_statistics(root_node, tree, &tree_linearization);
    logger->trace("Finished calculating tree statistics in {}s", timer.elapsed());

    timer.reset();
    logger->trace("Starting rmq preprocessing..");
    rmq_preprocessing(tree_linearization);
    logger->trace("Finished rmq preprocessing in {}s", timer.elapsed());
}

NormalizedTaxId TaxonomyDB::find_lca(const std::vector<NormalizedTaxId> &taxids) const {
    if (taxids.empty()) {
        logger->error("Internal error: Can't find LCA for an empty set of taxids.");
        std::exit(1);
    }
    uint64_t left_idx = this->node_to_linearization_idx[taxids[0]];
    uint64_t right_idx = this->node_to_linearization_idx[taxids[0]];
    for (const NormalizedTaxId &taxid : taxids) {
        if (this->node_to_linearization_idx[taxid] < left_idx) {
            left_idx = this->node_to_linearization_idx[taxid];
        }
        if (this->node_to_linearization_idx[taxid] > right_idx) {
            right_idx = this->node_to_linearization_idx[taxid];
        }
    }

    uint64_t log_dist = this->fast_log2[right_idx - left_idx];
    uint64_t left_lca = this->rmq_data[log_dist][left_idx];
    uint64_t right_lca =
            this->rmq_data[log_dist][right_idx - this->fast_pow2[log_dist] + 1];

    if (this->node_depth[left_lca] > this->node_depth[right_lca]) {
        return left_lca;
    }
    return right_lca;
}

bool TaxonomyDB::get_normalized_taxid(const std::string accession_version,
                                      NormalizedTaxId *taxid) const {
    num_get_taxid_calls += 1;
    if (!this->label_taxid_map.count(accession_version)) {
        // accession_version not in the this->label_taxid_map
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
                NormalizedTaxId lca =
                        find_lca(std::vector<NormalizedTaxId>{this->taxonomic_map[index],
                                                         taxid});
                this->taxonomic_map[index] = lca;
          }
        });
    }
}

void TaxonomyDB::export_to_file(const std::string &filepath) {
    if (num_get_taxid_calls_failed) {
        logger->warn("During the annotation matrix parsing, there could not be found {} taxids associated to the labels in the anno matrix (total evaluations:  {}).",
                     num_get_taxid_calls_failed, num_get_taxid_calls);
    }

    Timer timer;
    logger->trace("Exporting metagraph taxonomic data..");

    std::ofstream f(filepath.c_str(), std::ios::out | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file {}.", filepath.c_str());
        std::exit(1);
    }

    const std::vector<NormalizedTaxId> &linearization = this->rmq_data[0];
    tsl::hopscotch_map<TaxId, TaxId> node_parent;

    node_parent[this->denormalized_taxid[linearization[0]]] =
            this->denormalized_taxid[linearization[0]];

    // Start the iteration from 1, because we are computing linearization[i - 1].
    for (uint64_t i = 1; i < linearization.size(); ++i) {
        uint64_t act = linearization[i];
        uint64_t prv = linearization[i - 1];
        if (node_parent.count(this->denormalized_taxid[act])) {
            // Prv is act's child.
            continue;
        }
        // Prv is act's parent.
        node_parent[this->denormalized_taxid[act]] = this->denormalized_taxid[prv];
    }
    serialize_number_number_map(f, node_parent);

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
    std::sort(taxid_frequencies.begin(), taxid_frequencies.end(),
        [](const std::pair<uint64_t, uint64_t> &a, const std::pair<uint64_t, uint64_t> &b){
            if (a.first != b.first) {
                return a.first > b.first;
            }
            return a.second > b.second;
        });

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

TaxId TaxonomyDB::assign_class(const graph::AnnotatedDBG &anno, const std::string &sequence) const {
	std::vector<std::string> labels_discovered = anno.get_labels(sequence, 0.3);

    std::vector<NormalizedTaxId> curr_taxids;

    std::cerr << "this->normalized_taxid.size=" << this->normalized_taxid.size() << " ";

	for (uint64_t i = 0; i < labels_discovered.size(); ++i) {
        std::string act_str = utils::split_string(labels_discovered[i], "|")[1];
        TaxId act = static_cast<uint64_t>(std::stoull(act_str));
        std::cerr << " act_str=" << act_str << " act=" << act;
        if (this->normalized_taxid.count(act)) {
            curr_taxids.push_back(this->normalized_taxid.at(act));
        }
	}
    if (curr_taxids.size() == 0) {
        std::cerr << "could not find valid labels\n";
        return 0;
    }
    std::cerr << " " << this->denormalized_taxid[find_lca(curr_taxids)] << "\n";
	return this->denormalized_taxid[find_lca(curr_taxids)];
}

} // namespace annot
} // namespace mtg
