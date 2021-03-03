#include "transform_anno_tax.hpp"

#include <filesystem>
#include <tsl/hopscotch_set.h>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/taxonomy/taxonomic_db.hpp"
#include "common/threads/threading.hpp"
#include "config/config.hpp"
#include "cli/load/load_annotation.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

std::vector<std::string> get_anno_filenames(const std::vector<std::string> &paths) {
    std::vector<std::string> filenames;

    for(const auto &path: paths) {
        if (! filesystem::is_directory(path)) {
            if (!utils::ends_with(path, ".annodbg")) {
                logger->warn("File '{}' is not an '.annodbg' file.", path);
                continue;
            }
            filenames.push_back(path);
            continue;
        }
        logger->trace("Looking for anno files in directory: '{}'", path);
        for (auto &rec_file: filesystem::recursive_directory_iterator(path)) {
            if (utils::ends_with(rec_file.path().string(), ".annodbg")) {
                logger->trace("Found anno file '{}'", rec_file.path().string());
                filenames.push_back(rec_file.path().string());
            }
        }
    }
    return filenames;
}

void call_taxonomy_map_updates(annot::TaxonomyDB &taxonomy,
                               const std::vector<std::string> &filenames,
                               cli::Config *config) {
    std::mutex taxo_mutex;
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (uint64_t i = 0; i < filenames.size(); ++i) {
        const auto &file = filenames[i];
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annot =
                cli::initialize_annotation(file, *config);
        if (!annot->load(file)) {
            logger->error("Cannot load annotations from file '{}'", file);
            exit(1);
        }
        taxonomy.kmer_to_taxid_map_update(*annot, taxo_mutex);
    }
}

int transform_anno_taxo(Config *config) {
    assert(config);
    const auto &filenames = get_anno_filenames(config->fnames);
    uint64_t mem_bytes = config->memory_available * 1e9;

    tsl::hopscotch_set<std::string> all_labels;

    // Get all labels from all the annotation files.
    for (const auto &file: filenames) {
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annot =
                cli::initialize_annotation(file, *config);
        if (!annot->load(file)) {
            logger->error("Cannot load annotations from file '{}'", file);
            exit(1);
        }
        auto file_labels = annot->get_all_labels();
        for (const auto &it: file_labels) {
            all_labels.insert(annot::TaxonomyDB::get_accession_version_from_label(it));
        }
    }

    annot::TaxonomyDB taxonomy(config->taxonomic_tree, config->lookup_table, all_labels);

    // Load as many annotation matrix files as we can fit in memory.
    for (uint32_t i = 0; i < filenames.size(); ) {
        logger->trace("Loading a new batch of files...");
        size_t mem_bytes_left = mem_bytes;
        std::vector<std::string> file_batch;
        for ( ; i < filenames.size(); ++i) {
            uint64_t file_size = std::filesystem::file_size(filenames[i]);

            // If there is only one file in tha batch and its size exceeds 'mem_bytes', still construct a batch of size 1.
            if (file_size > mem_bytes_left && file_batch.size() > 0) {
                break;
            }

            if (file_size > mem_bytes) {
                logger->warn(
                        "File {} requires {} MB, more memory than the available {} MB. Can't optimize the processing of this file.",
                        filenames[i], file_size / 1e6, mem_bytes / 1e6);
            }

            mem_bytes_left -= file_size;
            file_batch.push_back(filenames[i]);
        }
        call_taxonomy_map_updates(taxonomy, file_batch, config);
    }

    taxonomy.export_to_file(config->outfbase + ".taxo");
    return 0;
}

}
}