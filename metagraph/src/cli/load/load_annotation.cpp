#include "load_annotation.hpp"

#include <string>

#include "common/logger.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;


Config::AnnotationType parse_annotation_type(const std::string &filename) {
    if (utils::ends_with(filename, anno::ColumnCompressed<>::kExtension)) {
        return Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, anno::RowCompressed<>::kExtension)) {
        return Config::AnnotationType::RowCompressed;

    } else if (utils::ends_with(filename, anno::MultiBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::BRWT;

    } else if (utils::ends_with(filename, anno::BinRelWT_sdslAnnotator::kExtension)) {
        return Config::AnnotationType::BinRelWT_sdsl;

    } else if (utils::ends_with(filename, anno::BinRelWTAnnotator::kExtension)) {
        return Config::AnnotationType::BinRelWT;

    } else if (utils::ends_with(filename, anno::RowFlatAnnotator::kExtension)) {
        return Config::AnnotationType::RowFlat;

    } else if (utils::ends_with(filename, anno::RainbowfishAnnotator::kExtension)) {
        return Config::AnnotationType::RBFish;

    } else if (utils::ends_with(filename, annotate::RbBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::RbBRWT;

    } else {
        logger->error("Unknown annotation format in '{}'", filename);
        exit(1);
    }
}

std::unique_ptr<anno::MultiLabelEncoded<std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      size_t column_compressed_num_columns_cached,
                      bool row_compressed_sparse,
                      uint64_t num_rows) {
    std::unique_ptr<anno::MultiLabelEncoded<std::string>> annotation;

    switch (anno_type) {
        case Config::ColumnCompressed: {
            annotation.reset(
                new anno::ColumnCompressed<>(num_rows, column_compressed_num_columns_cached)
            );
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new anno::RowCompressed<>(num_rows, row_compressed_sparse));
            break;
        }
        case Config::BRWT: {
            annotation.reset(new anno::MultiBRWTAnnotator());
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new anno::BinRelWT_sdslAnnotator());
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new anno::BinRelWTAnnotator());
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new anno::RowFlatAnnotator());
            break;
        }
        case Config::RBFish: {
            annotation.reset(new anno::RainbowfishAnnotator());
            break;
        }
        case Config::RbBRWT: {
            annotation.reset(new annotate::RbBRWTAnnotator());
            break;
        }
    }

    return annotation;
}

} // namespace cli
} // namespace mtg
