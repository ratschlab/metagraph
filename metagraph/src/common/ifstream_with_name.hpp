#pragma once

#include <fstream>
#include <string>

namespace mtg {
namespace common {

class IfstreamWithName : public std::ifstream {
  public:
    IfstreamWithName(const std::string &fname, std::ios_base::openmode mode)
        : std::ifstream(fname, mode), fname_(fname) {}

    const std::string& get_name() const { return fname_; }

    private:
      std::string fname_;
};

} // namespace common
} // namespace mtg
