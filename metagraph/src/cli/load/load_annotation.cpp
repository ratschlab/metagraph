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

const uint64_t kBytesInGigabyte = 1'000'000'000;


Config::AnnotationType parse_annotation_type(const std::string &filename) {
    if (utils::ends_with(filename, annot::ColumnCompressed<>::kExtension)) {
        return Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, annot::ColumnCoordAnnotator::kExtension)) {
        return Config::AnnotationType::ColumnCoord;

    } else if (utils::ends_with(filename, annot::RowDiffCoordAnnotator::kExtension)) {
        return Config::AnnotationType::RowDiffCoord;

    } else if (utils::ends_with(filename, annot::RowDiffColumnAnnotator::kExtension)) {
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

    } else if (utils::ends_with(filename, annot::RowSparseAnnotator::kExtension)) {
        return Config::AnnotationType::RowSparse;

    } else if (utils::ends_with(filename, annot::RowDiffRowSparseAnnotator ::kExtension)) {
        return Config::AnnotationType::RowDiffRowSparse;

    } else if (utils::ends_with(filename, annot::RainbowfishAnnotator::kExtension)) {
        return Config::AnnotationType::RBFish;

    } else if (utils::ends_with(filename, annot::RbBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::RbBRWT;

    } else if (utils::ends_with(filename, annot::IntMultiBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::IntBRWT;

    } else if (utils::ends_with(filename, annot::IntRowDiffBRWTAnnotator::kExtension)) {
        return Config::AnnotationType::IntRowDiffBRWT;

    } else {
        logger->error("Unknown annotation format in '{}'", filename);
        exit(1);
    }
}

std::unique_ptr<annot::MultiLabelEncoded<std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      size_t column_compressed_num_columns_cached,
                      bool row_compressed_sparse,
                      uint64_t num_rows,
                      const std::string &swap_dir,
                      double memory_available_gb,
                      uint8_t count_width) {
    std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation;

    switch (anno_type) {
        case Config::ColumnCompressed: {
            annotation.reset(
                new annot::ColumnCompressed<>(num_rows, column_compressed_num_columns_cached,
                                              swap_dir, memory_available_gb * kBytesInGigabyte,
                                              count_width)
            );
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new annot::RowCompressed<>(num_rows, row_compressed_sparse));
            break;
        }
        case Config::RowDiff: {
            annotation.reset(new annot::RowDiffColumnAnnotator());
            break;
        }
        case Config::RowSparse: {
            annotation.reset(new annot::RowSparseAnnotator());
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
        case Config::RowDiffRowSparse: {
            annotation.reset(new annot::RowDiffRowSparseAnnotator());
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
        case Config::IntBRWT: {
            annotation.reset(new annot::IntMultiBRWTAnnotator());
            break;
        }
        case Config::IntRowDiffBRWT: {
            annotation.reset(new annot::IntRowDiffBRWTAnnotator());
            break;
        }
        case Config::ColumnCoord: {
            annotation.reset(new annot::ColumnCoordAnnotator());
            break;
        }
        case Config::RowDiffCoord: {
            annotation.reset(new annot::RowDiffCoordAnnotator());
            break;
        }
    }

    return annotation;
}

} // namespace cli
} // namespace mtg
