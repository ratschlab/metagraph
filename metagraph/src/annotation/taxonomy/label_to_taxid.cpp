#include "tax_classifier.hpp"

#include "common/utils/string_utils.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace annot {

using mtg::common::logger;

void TaxonomyBase::assign_label_type(const std::string &label, bool *require_accversion_to_taxid_map) {
    if (!utils::starts_with(label, ">gi|")) {
        // e.g.   >gi|1070643132|ref|NC_031224.1| Arthrobacter phage Mudcat, complete genome
        this->label_type = GEN_BANK;
        *require_accversion_to_taxid_map = true;
    } else if (!utils::starts_with(label, ">kraken:")) {
        // e.g.   >kraken:taxid|2016032|NC_047834.1 Alteromonas virus vB_AspP-H4/4, complete genome
        this->label_type = KRAKEN;
        *require_accversion_to_taxid_map = false;
    } else {
        logger->error("Can't determine the type of the given label {}. Please make sure that the labels are in a recognized format.", label);
        std::exit(1);
    }
}

bool TaxonomyBase::get_taxid_from_label(const std::string &label, TaxId *taxid) const {
    if (this->label_type == KRAKEN) {
        *taxid = static_cast<uint32_t>(std::stoull(utils::split_string(label, "|")[1]));
        return true;
    } else if (TaxonomyBase::label_type == GEN_BANK) {
        std::string acc_version = this->get_accession_version_from_label(label);
        if (not this->accversion_to_taxid_map.count(acc_version)) {
            return false;
        }
        *taxid = this->accversion_to_taxid_map.at(acc_version);
        return true;
    }

    logger->error("Run get_taxid_from_label() for unknown label {}.", label);
    std::exit(1);
}

std::string TaxonomyBase::get_accession_version_from_label(const std::string &label) const {
    if (this->label_type == KRAKEN) {
        return utils::split_string(utils::split_string(label, "|")[2], " ")[0];
    } else if (this->label_type == GEN_BANK) {
        return utils::split_string(label, "|")[3];;
    }

    logger->error("Run get_accession_version_from_label() for unknown label {}.", label);
    std::exit(1);
}

// TODO improve this by parsing the compressed ".gz" version (or use https://github.com/pmenzel/taxonomy-tools)
void TaxonomyBase::read_accversion_to_taxid_map(const std::string &filepath,
                                                const graph::AnnotatedDBG *anno_matrix = NULL) {
    std::ifstream f(filepath);
    if (!f.good()) {
        logger->error("Failed to open accession to taxid map table {}", filepath);
        exit(1);
    }

    std::string line;
    getline(f, line);
    if (!utils::starts_with(line, "accession\taccession.version\ttaxid\t")) {
        logger->error("The accession to taxid map table is not in the standard (*.accession2taxid) format {}.", filepath);
        exit(1);
    }

    tsl::hopscotch_set<std::string> input_accessions;
    if (anno_matrix != NULL) {
        for (const std::string &accversion : anno_matrix->get_annotation().get_all_labels()) {
            input_accessions.insert(accversion);
        }
    }

    while (getline(f, line)) {
        if (line == "") {
            logger->error("The accession to taxid map table contains empty lines. Please make sure that this file was not manually modified {}.", filepath);
            exit(1);
        }
        std::vector<std::string> parts = utils::split_string(line, "\t");
        if (parts.size() <= 2) {
            logger->error("The accession to taxid map table contains incomplete lines. Please make sure that this file was not manually modified {}.", filepath);
            exit(1);
        }
        if (input_accessions.size() == 0 || input_accessions.count(parts[1])) {
            this->accversion_to_taxid_map[parts[1]] = static_cast<TaxId>(std::stoull(parts[2]));
        }
    }
}

} // namespace annot
} // namespace mtg
