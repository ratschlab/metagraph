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
    if (utils::ends_with(filename, annot::ColumnCompressed<>::kExtension)) {
        return Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, annot::RowDiffAnnotator::kExtension)) {
        return Config::AnnotationType::RowDiff;

    } else if (utils::ends_with(filename, annot::RowCompressed<>::kExtension)) {
        return Config::AnnotationType::RowCompressed;

    } else if (utils::ends_with(filename, annot::MultiBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::BRWT;

    } else if (utils::ends_with(filename, annot::RowDiffBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::RowDiffBRWT;

    } else if (utils::ends_with(filename, annot::BinRelWT_sdslAnnotator::kExtension)) {
        return Config::AnnotationType::BinRelWT_sdsl;

    } else if (utils::ends_with(filename, annot::BinRelWTAnnotator::kExtension)) {
        return Config::AnnotationType::BinRelWT;

    } else if (utils::ends_with(filename, annot::RowFlatAnnotator::kExtension)) {
        return Config::AnnotationType::RowFlat;

    } else if (utils::ends_with(filename, annot::RainbowfishAnnotator::kExtension)) {
        return Config::AnnotationType::RBFish;

    } else if (utils::ends_with(filename, annot::RbBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::RbBRWT;

    } else {
        logger->error("Unknown annotation format in '{}'", filename);
        exit(1);
    }
}

std::unique_ptr<annot::MultiLabelEncoded<std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      size_t column_compressed_num_columns_cached,
                      bool row_compressed_sparse,
                      uint64_t num_rows) {
    std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation;

    switch (anno_type) {
        case Config::ColumnCompressed: {
            annotation.reset(
                new annot::ColumnCompressed<>(num_rows, column_compressed_num_columns_cached)
            );
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new annot::RowCompressed<>(num_rows, row_compressed_sparse));
            break;
        }
        case Config::RowDiff: {
            annotation.reset(new annot::RowDiffAnnotator());
            break;
        }
        case Config::BRWT: {
            annotation.reset(new annot::MultiBRWTAnnotator());
            break;
        }
        case Config::RowDiffBRWT: {
            annotation.reset(new annot::RowDiffBRWTAnnotator());
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new annot::BinRelWT_sdslAnnotator());
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new annot::BinRelWTAnnotator());
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new annot::RowFlatAnnotator());
            break;
        }
        case Config::RBFish: {
            annotation.reset(new annot::RainbowfishAnnotator());
            break;
        }
        case Config::RbBRWT: {
            annotation.reset(new annot::RbBRWTAnnotator());
            break;
        }
    }

    return annotation;
}

} // namespace cli
} // namespace mtg
