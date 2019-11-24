#pragma once

#include "common/chunked_wait_queue.hpp"

#include <string>
#include <vector>

namespace mg {
namespace common {
template <typename T>
class FileMerger {
  public:
    FileMerger(const std::vector<std::string> &sources, const std::string &out_file);
    void merge(ChunkedWaitQueue<T>* result);

  private:
    std::vector<std::string> sources_;
    std::string out_file_;
};

template <typename T>
void merge_files(const std::vector<std::string> &sources, const std::string &out_file,
        ChunkedWaitQueue<T>* result);

} // namespace common
} // namespace mg
