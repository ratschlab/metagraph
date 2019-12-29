#ifndef __LOAD_ANNOTATION_HPP__
#define __LOAD_ANNOTATION_HPP__

#include <string>

#include "annotation/annotate.hpp"
#include "cli/config.hpp"


Config::AnnotationType parse_annotation_type(const std::string &filename);

std::unique_ptr<annotate::MultiLabelEncoded<uint64_t, std::string>>
initialize_annotation(Config::AnnotationType anno_type,
                      const Config &config,
                      uint64_t num_rows);

inline auto initialize_annotation(const std::string &filename,
                                  const Config &config) {
    return initialize_annotation(parse_annotation_type(filename), config, 0);
}

#endif // __LOAD_ANNOTATION_HPP__
