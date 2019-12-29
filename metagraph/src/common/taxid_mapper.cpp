#include "taxid_mapper.hpp"

#include <iostream>
#include <sstream>

#include <zlib.h>

#include "common/utils/string_utils.hpp"
#include "common/serialization.hpp"


std::string TaxIDMapper::parse_label(const std::string &gb) {
    auto split = utils::split_string(gb, "|");
    if (split.size() > 1) {
        auto version_split = utils::split_string(split[1], ".");
        if (version_split.size() > 1) {
            version_split.pop_back();
            return utils::join_strings(version_split, ".");
        } else {
            return split[1];
        }
    } else {
        auto version_split = utils::split_string(gb, ".");
        if (version_split.size() > 1) {
            version_split.pop_back();
            return utils::join_strings(version_split, ".");
        } else {
            return gb;
        }
    }
}

bool TaxIDMapper::parse_accession2taxid(const std::string &accession2taxid) {
    gzFile input_p = gzopen(accession2taxid.c_str(), "rb");
    if (input_p == Z_NULL)
        return false;

    char buf[1024];
    char *line;

    std::string accession, ignore;
    taxid_t taxid;

    // check if reading header works
    if ((line = gzgets(input_p, buf, sizeof(buf))) == NULL)
        return false;

    while ((line = gzgets(input_p, buf, sizeof(buf))) != NULL) {
        std::istringstream sin(line);

        if (!(sin >> accession >> ignore >> taxid))
            return false;

        auto insert = gb_to_taxid_.emplace(accession, taxid);
        if (!insert.second)
            std::cerr << "Warning: duplicate accession ID " << accession << std::endl;
    }

    gzclose(input_p);
    return true;
}

bool TaxIDMapper::parse_nodes(const std::string &nodes) {
    gzFile input_p = gzopen(nodes.c_str(), "rb");
    if (input_p == Z_NULL)
        return false;

    char buf[1024];
    char *line;
    std::string block1, block2;
    uint64_t taxid;
    uint64_t parentid;

    while ((line = gzgets(input_p, buf, sizeof(buf))) != NULL) {
        std::istringstream sin(line);

        if (!(sin >> taxid >> block1 >> parentid >> block1 >> block1 >> block2))
            return false;

        if (block2 != "|")
            block1 += " " + block2;

        auto insert_tree = taxid_to_parent_.emplace(taxid, parentid);
        if (!insert_tree.second)
            std::cerr << "Warning: duplicate taxid in taxonomy tree "
                      << taxid << " " << insert_tree.first->second
                      << std::endl;

        auto insert_rank = taxid_to_rank_label_.emplace(taxid, block1);
        if (!insert_rank.second)
            std::cerr << "Warning: duplicate taxid in rank map "
                      << taxid << " " << block1
                      << std::endl;
    }

    gzclose(input_p);
    return true;
}

TaxIDMapper::taxid_t TaxIDMapper::gb_to_taxid(const std::string &gb) const {
    auto find = gb_to_taxid_.find(parse_label(gb));
    if (find == gb_to_taxid_.end())
        return 0;

    return find->second;
}

TaxIDMapper::taxid_t TaxIDMapper::get_parent(const taxid_t &id) const {
    if (id == 1)
        return 1;

    auto find = taxid_to_parent_.find(id);
    if (find == taxid_to_parent_.end())
        return 0;

    return find->second;
}

std::string TaxIDMapper::get_rank_label(const taxid_t &id) const {
    auto find = taxid_to_rank_label_.find(id);
    if (find == taxid_to_rank_label_.end())
        return "";

    return find->second;
}

TaxIDMapper::taxid_t TaxIDMapper
::get_ancestor_with_rank_label(taxid_t id,
                               const std::string &rank_label) const {
    assert(rank_label.length());

    do {
        auto cur_rank_label = get_rank_label(id);
        if (!cur_rank_label.length())
            return 0;

        if (cur_rank_label == rank_label)
            return id;

        id = get_parent(id);
    } while (id != 1);

    return rank_label == "no rank" ? id : 0;
}

std::pair<TaxIDMapper::taxid_t, std::string> TaxIDMapper
::get_ancestor_labeled(taxid_t id) const {
    do {
        auto cur_rank_label = get_rank_label(id);
        if (!cur_rank_label.length())
            return { 0, "" };

        if (cur_rank_label != "no rank")
            return { id, cur_rank_label };

        id = get_parent(id);
    } while (id != 1);

    return { 0, "" };
}

bool TaxIDMapper::load(std::ifstream &in) {
    if (!in.good())
        return false;

    if (!load_string_number_map(in, &gb_to_taxid_)) {
        std::cerr << "Failed to load Accession ID to Taxonomy ID map" << std::endl;
        return false;
    }

    if (!load_number_number_map(in, &taxid_to_parent_)) {
        std::cerr << "Failed to load taxonomy tree" << std::endl;
        return false;
    }

    if (!load_number_string_map(in, &taxid_to_rank_label_)) {
        std::cerr << "Failed to load Taxonomy ID to rank label map" << std::endl;
        return false;
    }

    return true;
}

void TaxIDMapper::serialize(std::ofstream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump mapper to a bad stream");

    serialize_string_number_map(out, gb_to_taxid_);
    serialize_number_number_map(out, taxid_to_parent_);
    serialize_number_string_map(out, taxid_to_rank_label_);
}
