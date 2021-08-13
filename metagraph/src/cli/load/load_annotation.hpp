#ifndef __LOAD_ANNOTATION_HPP__
#define __LOAD_ANNOTATION_HPP__

#include <string>

#include "annotation/representation/base/annotation.hpp"
#include "cli/config/config.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"


namespace mtg {
namespace cli {

Config::AnnotationType parse_annotation_type(const std::string &filename);

std::unique_ptr<annot::MultiLabelEncoded<std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      size_t column_compressed_num_columns_cached = 1,
                      bool row_compressed_sparse = false,
                      uint64_t num_rows = 0,
                      const std::string &swap_dir = "",
                      double memory_available_gb = 1,
                      uint8_t count_width = 8);

std::unique_ptr<annot::RowCompressed<std::string>>
initialize_annotation_dd(Config::AnnotationType anno_type,
                      bool row_compressed_sparse = false,
                      uint64_t num_rows = 0);

inline std::unique_ptr<annot::MultiLabelEncoded<std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      const Config &config,
                      uint64_t num_rows = 0) {
    return initialize_annotation(anno_type, config.num_columns_cached, config.sparse,
                                 num_rows, config.tmp_dir, config.memory_available,
                                 config.count_width);
}

inline std::unique_ptr<annot::RowCompressed<std::string>>
initialize_annotation_dd(Config::AnnotationType anno_type,
                         const uint64_t max_index) {
    return initialize_annotation_dd(anno_type, true, max_index);
}

template <typename... Args>
inline std::unique_ptr<annot::MultiLabelEncoded<std::string>>
initialize_annotation(const std::string &filename, const Args &... args) {
    return initialize_annotation(parse_annotation_type(filename), args...);
}

template <typename... Args>
inline std::unique_ptr<annot::RowCompressed<std::string>>
initialize_annotation_dd(const std::string &filename, const uint64_t max_index22) {
    return initialize_annotation_dd(parse_annotation_type(filename), max_index22);
}

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATION_HPP__
