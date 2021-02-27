#include "transform_anno_tax.hpp"

#include "common/logger.hpp"
#include "config/config.hpp"

#include <tsl/hopscotch_set.h>
#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/taxonomy/taxonomic_db.hpp"
#include "cli/load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

int transform_anno_taxo(Config *config) {
    assert(config);
    const auto &filenames = config->fnames;
    uint64_t mem_bytes = config->memory_available * 1e9;

    tsl::hopscotch_set<std::string> all_labels;

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

    // load as many annotation matrix files as we can fit in memory.
    for (uint32_t i = 0; i < filenames.size(); ) {
        logger->trace("Loading columns for batch-conversion...");
        size_t mem_bytes_left = mem_bytes;
        std::vector<std::string> file_batch;
        for ( ; i < filenames.size(); ++i) {
            uint64_t file_size = std::filesystem::file_size(filenames[i]);
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
        taxonomy.kmer_to_taxid_map_update(file_batch, config);
    }

    taxonomy.export_to_file(config->outfbase + ".taxo");
    return 0;
}

}
}