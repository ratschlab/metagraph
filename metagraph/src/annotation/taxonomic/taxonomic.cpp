#include "taxonomic.hpp"

#include "common/logger.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream> // TODO delete
#include "common/serialization.hpp"
#include "seq_io/sequence_io.hpp"
#include "common/utils/string_utils.hpp"
#include "common/unix_tools.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

typedef Taxonomy::TaxNormalizedId TaxNormalizedId;


void Taxonomy::calculate_node_depth(const TaxNormalizedId &node, const ChildrenList &tree) {
//    std::cout << "calculate_node_depth for node=" << node << std::endl;
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
        const TaxNormalizedId &node,
        const ChildrenList &tree,
        std::vector<TaxNormalizedId> &tree_linearization
) {
    linearization_idx[node] = tree_linearization.size();
    tree_linearization.push_back(node);
    for(auto child: tree[node]) {
        dfs_linearization(child, tree, tree_linearization);
        tree_linearization.push_back(node);
    }
}

void Taxonomy::read_tree(const std::string &tree_filepath,
                         ChildrenList &tree) {
    std::ifstream f(tree_filepath);
    std::string line;

    tsl::hopscotch_map<TaxId, TaxId> full_parents_list;
    uint cnt = 0;
    while (getline(f, line) ) {
        std::vector<std::string> parts = utils::split_string(line, "\t");
        uint64_t act = static_cast<uint64_t>(std::stoll(parts[0]));
        uint64_t parent = static_cast<uint64_t>(std::stoll(parts[2]));

        if (full_parents_list.count(act)) {
            std::cout << "something already existent in full_parents_list!!" << "\n";
        }
        full_parents_list[act] = parent;
        if (cnt < 20) {
            std::cout << "tree -> " << parts[0] << " " << parts[2] << "\n";
            std::cout << "tree -0 " << act << " " << parent << "\n";
            std::cout << " ---> " << full_parents_list[act] << "\n\n";
        }
        if (parts[0] == parts[2]) {
            std::cout << "tree equal -> " << parts[0] << "\n";
        }
        cnt += 1;
    }
    f.close();

    std::queue<TaxId> relevant_taxids;
    for (const auto &it: reversed_lookup_table) {
        relevant_taxids.push(it.first);
    }
    std::cout << "init size relevant_taxids = " << relevant_taxids.size() << "\n";

    uint64_t num_nodes = 0;
    while (relevant_taxids.size()) {
        TaxId taxid = relevant_taxids.front();

        if (taxid == 1) {
            std::cout << " check 1 -> " << full_parents_list[taxid] << "\n";
        }

        relevant_taxids.pop();
        if (normalized_tax_it.count(taxid)) {
            continue;
        }
        if (! full_parents_list.count(taxid)) {
            std::cout << "WARNING taxid " << taxid << " cannot be found in the taxonomic tree\n";
//            todo add logs
            continue;
        }
        normalized_tax_it[taxid] = num_nodes++;
        if (taxid == full_parents_list[taxid]) {
            root_node = normalized_tax_it[taxid];
            std::cout << "root_node = " << root_node << "\n";
            continue;
        }
        relevant_taxids.push(full_parents_list[taxid]);
        if (reversed_lookup_table.count(taxid)) {
            lookup_table[reversed_lookup_table[taxid]] = normalized_tax_it[taxid];
        }
    }
    std::cout << "num_nodes = " << num_nodes << "\n";

    tree.resize(num_nodes);
    for (const auto &it: normalized_tax_it) {
        TaxId taxid = it.first;
        if (normalized_tax_it[taxid] == root_node) {
            continue;
        }
        tree[normalized_tax_it[full_parents_list[taxid]]].push_back(
                normalized_tax_it[taxid]);
    }
}

void Taxonomy::rmq_preprocessing(const ChildrenList &tree) {
    std::vector<TaxNormalizedId> tree_linearization;
    linearization_idx.resize(tree.size());
    dfs_linearization(root_node, tree, tree_linearization);

    uint num_rmq_rows = log (tree_linearization.size()) + 2;
//    rmq_data = std::make_unique<std::vector<TaxNormalizedId>>(num_rmq_rows);
    rmq_data.resize(num_rmq_rows); // = new std::unique_ptr<std::vector<TaxNormalizedId>>[num_rmq_rows];
    for (uint64_t i = 0; i < num_rmq_rows; ++i) {
//        rmq_data[i] = std::make_unique<std::vector<TaxNormalizedId>>();
        rmq_data[i].resize(tree_linearization.size());
    }

//    Copy tree_linearization on rmq line 0.
    for (uint64_t i = 0; i < tree_linearization.size(); ++i) {
//        rmq_data[0]->push_back(std::move(tree_linearization[i]));
        rmq_data[0][i] = tree_linearization[i];
    }

    uint delta = 1;
    for (uint64_t row = 1; row < num_rmq_rows; ++row) {
        for (uint64_t i = 0; i + delta < tree_linearization.size(); ++i) {
            if (node_depth[rmq_data[row-1][i]] > node_depth[rmq_data[row-1][i + delta]]) {
//                rmq[row]->push_back(std::move(rmq_data[row-1][i]));
                rmq_data[row][i] = rmq_data[row-1][i];
            } else {
//                rmq[row]->push_back(std::move(rmq_data[row-1][i + delta]));
                rmq_data[row][i] = rmq_data[row-1][i + delta];
            }
        }
        delta *= 2;
    }

    precalc_log.resize(tree_linearization.size());
    precalc_pow2.push_back(1);
    for (uint64_t i = 2; i < tree_linearization.size(); ++i) {
        precalc_log[i] = 1 + precalc_log[i/2];
        if(precalc_log[i] > precalc_log[i-1]) {
            precalc_pow2.push_back(i);
        }
    }

    for (uint64_t row = 0; row < num_rmq_rows; ++row) {
        for (uint i = 0; i < tree_linearization.size(); ++i) {
//            std::cout << rmq_data[row][i] << " ";
            if (rmq_data[row][i] > normalized_tax_it.size()) {
                std::cout << "rmq_data["<<row<<"]["<<i<<"] = " <<rmq_data[row][i] << std::endl;
                exit(1);
            }
        }
//        std::cout << "\n";
    }

    std::cout << "num_rmq_rows=" << num_rmq_rows << "\n";
//    exit(0);
}

void Taxonomy::get_input_accessions(const std::string &fasta_fai,
                                    tsl::hopscotch_set<AccessionVersion> &input_accessions) {
    std::cout << "in get_input_accessions\n";
    std::ifstream f(fasta_fai);
    std::string line;

    uint cnt = 0;
    while (getline(f, line) ) {
        std::string fasta_header = utils::split_string(line, "\t")[0];
        std::string accession_version = utils::split_string(fasta_header, "|")[3];
        input_accessions.insert(accession_version);
//        if (cnt < 10) {
//            std::cout << "found accession = " << accession_version << "\n";
//        }
        cnt += 1;
    }
    f.close();

    std::cout << "finished get_input_accessions (cnt=" << cnt << ")\n";
}

void Taxonomy::parse_lookup_table(const std::string &lookup_table_filepath,
                                  const tsl::hopscotch_set<AccessionVersion> &input_accessions) {
    std::ifstream f(lookup_table_filepath);
//    Better use a small bash script or https://github.com/pmenzel/taxonomy-tools for parsing ".gz" directly.
    std::string line;
//    uint64_t cnt = 0;
    while (getline(f, line) ) {
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if(input_accessions.count(parts[1])) {
            reversed_lookup_table[std::stoull(parts[2])] = parts[1];

//            todo assert if too few parts
//            std::cout << "found " << parts[1] << "\n";
//            if (cnt < 10) {
//                std::cout << parts[1] <<" " << parts[2] << "\n";
//                cnt += 1;
//            }
        }
    }
    std::cout << " input_accessions.size=" << input_accessions.size() << "\n" <<
            " reversed_lookup_table.size=" << reversed_lookup_table.size() << "\n";

//  todo  index_to_label.resize()
    f.close();
}

Taxonomy::Taxonomy(const std::string &tree_filepath,
                   const std::string &lookup_table_filepath,
                   const std::string &fasta_fai) {

    if (!std::filesystem::exists(tree_filepath)) {
//        TODO check extension
        logger->error("Can't open taxonomic tree file '{}'.", tree_filepath);
        std::exit(1);
    }
    if (!std::filesystem::exists(lookup_table_filepath)) {
//        TODO check extension
        logger->error("Can't open taxonomic lookup table file '{}'.", lookup_table_filepath);
        std::exit(1);
    }
    if (!std::filesystem::exists(fasta_fai)) {
//        TODO check extension
        logger->error("Can't open '.fasta.fai' file '{}'.", fasta_fai);
        std::exit(1);
    }
    Timer timer;
    logger->trace("Parsing input accessions");
    tsl::hopscotch_set<AccessionVersion> input_accessions;
    get_input_accessions(fasta_fai, input_accessions);
    logger->trace("Finished parsing input accessions after '{}' sec", timer.elapsed());
    timer.reset();
    logger->trace("Parsing lookup table");
    parse_lookup_table(lookup_table_filepath, input_accessions);
    logger->trace("Finished parsing tookup table after '{}' sec", timer.elapsed());
    timer.reset();
    logger->trace("Parsing taxonomic tree");
    Taxonomy::ChildrenList tree;
    read_tree(tree_filepath, tree);
    logger->trace("Finished parsing taxonomic tree after '{}' sec", timer.elapsed());
    node_depth.resize(tree.size());
    calculate_node_depth(root_node, tree);
    logger->trace("Starting rmq preprocessing");
    rmq_preprocessing(tree);
    logger->trace("Finished rmq preprocessing after '{}' sec", timer.elapsed());

    std::cout << "\n\n\n";

//    for (uint node = 0; node < tree.size(); ++node) {
//        std::cout << node << " -> ";
//        for (uint i = 0; i < tree[node].size(); ++i) {
//            std::cout << tree[node][i] << " ";
//        }
//        std::cout << std::endl;
//    }
}

TaxNormalizedId Taxonomy::find_lca(const std::vector<TaxNormalizedId> &labels) {
//    std::cout << "\nin find_lca with: ";
//    for (auto it: labels) {
//        std::cout << it << " ";
//    }
//    std::cout << "\n";
    if (labels[0] == 6437 && labels[1] == 3511) {
        std::cout << "on in" << std::endl;
    }
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
    assert(left_idx >= 0);
    assert(right_idx >= 0);
    if (labels[0] == 6437 && labels[1] == 3511) {
        std::cout << "## rmq.size=" << rmq_data.size() << "\n";
        for (uint64_t row = 0; row < rmq_data.size(); ++row) {
            std::cout << "\t## rmq["<<row<<"].size=" << rmq_data[row].size() << "\n";
            for (uint i = 0; i < rmq_data[row].size(); ++i) {
                if (rmq_data[row][i] > normalized_tax_it.size()) {
                    std::cout << "rmq_data[" << row << "][" << i
                              << "] = " << rmq_data[row][i] << std::endl;
                    exit(1);
                }
            }
        }
        std::cout << "rmq_data.size=" << rmq_data.size();
        std::cout << "rmq_data[12].size=" << rmq_data[12].size();
        std::cout << "rmq val = " << rmq_data[12][14487] << "\n";
        std::cout << "left_idx=" << left_idx << " right_idx=" << right_idx <<  std::endl;
    }
//    std::cout << "left_idx=" << left_idx << " right_idx=" << right_idx << "\n";
    uint64_t log_dist = precalc_log[right_idx - left_idx];
    assert(log_dist < rmq_data.size());
    if (labels[0] == 6437 && labels[1] == 3511) {
        std::cout << "log_dist=" << log_dist << std::endl;
    }
    assert(left_idx < rmq_data[log_dist].size());
    uint64_t left_lca = rmq_data[log_dist][left_idx];

    if (labels[0] == 6437 && labels[1] == 3511) {
        std::cout << "right_idx - precalc_pow2[log_dist] + 1=" << right_idx - precalc_pow2[log_dist] + 1 << std::endl;
    }
    assert(right_idx - precalc_pow2[log_dist] + 1 < rmq_data[log_dist].size());
    uint64_t right_lca = rmq_data[log_dist][right_idx - precalc_pow2[log_dist] + 1];
//    std::cout << "left_lca=" << left_lca << " right_lca=" << right_lca << "\n";


    if (labels[0] == 6437 && labels[1] == 3511) {
        std::cout << "off in left_lca=" << left_lca << " right_lca=" << right_lca <<  std::endl;
        std::cout << "node_depth[left_lca]=" << node_depth[left_lca] << std::endl;
        std::cout << "node_depth[right_lca]=" << node_depth[right_lca] << std::endl;
    }

    if (node_depth[left_lca] > node_depth[right_lca]) {
        return left_lca;
    }
    return right_lca;

}

bool Taxonomy::find_lca(const std::vector<AccessionVersion> &raw_labels,
                                   TaxNormalizedId &lca) {
    std::vector<TaxNormalizedId> labels;
    for (const auto &label: raw_labels) {
        std::string accession_version = utils::split_string(label, "|")[3];
        std::cout << "accession_version=" << accession_version << " label =" << label << "\n";
        if (! lookup_table.count(accession_version)) {
            logger->warn("Accession version {} cannot be found in the lookup_table.", accession_version);
            return false;
        }
        labels.push_back(lookup_table[accession_version]);
    }
    lca = find_lca(labels);
    return true;
}

TaxNormalizedId Taxonomy::find_lca(const TaxNormalizedId &label1, const TaxNormalizedId &label2) {
    return find_lca(std::vector<TaxNormalizedId> {label1, label2});
}

void Taxonomy::update_kmers_lca(const std::vector<KmerId> &indices,
                                const TaxNormalizedId &lca) {
    std:: cout << "here update_row_indices :\n";
//    for (auto i: indices) {
//        std::cout << i << " ";
//    }
//    std::cout << "\n";
//    for (auto i: labels) {
//        std::cout << i << " ";
//    }
//
//    TaxNormalizedId lca = find_lca(labels);
    std::cout << "\n lca =" << lca << "\n";

    if(lca == 3511) {
        std::cout << "on" << std::endl;
    }

    for (const auto &kmer: indices) {
        if (taxonomic_map.count(kmer)) {
//            if(lca == 3511) {
//                std::cout << "kmer=" << kmer << std::endl;
//                std::cout << " old taxonomic_map[kmer]=" << taxonomic_map[kmer] << std::endl;
//            }
            taxonomic_map[kmer] = find_lca(taxonomic_map[kmer], lca);
//            std::cout << " upd (" << kmer << ") -> " << taxonomic_map[kmer] << "\n";
        } else {
//            if(lca == 3511) {
//                std::cout << "kmer=" << kmer << " no old " << std::endl;
//            }
            taxonomic_map[kmer] = lca;
//            std::cout << " fst (" << kmer << ") -> " << taxonomic_map[kmer] << "\n";
        }
    }
    if(lca == 3511) {
        std::cout << "off" << std::endl;
    }

}

void Taxonomy::export_to_file(const std::string &filepath) {
    Timer timer;
    logger->trace("Exporting metagraph taxonomic data");

    std::ofstream f(filepath.c_str(), std::ios::out | std::ios::binary);
    if (!f.is_open()) {
        logger->error("Can't open taxonomic file '{}'.", filepath.c_str());
        std::exit(1);
    }

//    std::vector<uint> distrib(30);
    serialize_number_number_map(f, taxonomic_map);
//    f << taxonomic_map.size() << "\n";
//    for (auto &it: taxonomic_map) { // Certainly can be serialized.
//        f << it.first << " " << it.second << "\n";
//        distrib[it.second] ++;
//    }
//    for (int i = 0; i < 30; ++i) {
//        std::cout << "distrib["<<i<<"]=" << distrib[i] << "\n";
//    }

    serialize_string_vector(f, index_to_label);
//    f << index_to_label.size() << "\n";
//    for (auto &it: index_to_label) {
//        f << it << " ";
//    }
//    f << "\n";

    serialize_number_vector_raw(f, rmq_data[0]);
//    uint64_t num_nodes = index_to_label.size();
//    for (uint64_t i = 0; i < 2 * num_nodes - 1; ++i) {
//        f << rmq_data[0][i] << " ";
//    }
//    f << "\n";
    f.close();
    logger->trace("Finished exporting metagraph taxonomic data after '{}' sec", timer.elapsed());
}

}
}
