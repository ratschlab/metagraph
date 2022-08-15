#pragma once

#include <fstream>
#include <string>

namespace mtg {
namespace common {

class IfstreamWithNameAndOffset : public std::ifstream {
    std::string fname;
    uint64_t offset {};

  public:
    explicit IfstreamWithNameAndOffset(const std::string &fname, std::ios_base::openmode mode)
        : std::ifstream(fname, mode), fname(fname) {}

    const std::string &get_name() const { return fname; }

    void set_offset(uint64_t offset) { this->offset = offset; }

    uint64_t get_offset() const { return offset; }
};

} // namespace common
} // namespace mtg
