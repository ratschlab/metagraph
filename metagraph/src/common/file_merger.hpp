#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "common/logger.hpp"
#include "utils/template_utils.hpp"

namespace mg {
namespace common {
/**
 * Given a list of n source files, containing ordered elements of type T, merge the n
 * sources into a single (ordered) list.
 * @param sources the files containing sorted lists of type T
 * @param on_new_item callback to invoke when a new element was merged
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T>
void merge_files(const std::vector<std::string> sources,
                 std::function<void(const T &)> on_new_item) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::fstream> chunk_files(sources.size());

    std::priority_queue<std::pair<T, uint32_t>, std::vector<std::pair<T, uint32_t>>, utils::GreaterFirst>
            merge_heap;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        T data_item;
        if (chunk_files[i].good()) {
            chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(data_item));
            merge_heap.push({ data_item, i });
        } else {
            logger->error("Error: Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
    }
    uint64_t totalSize = 0;

    // initialized to suppress maybe-uninitialized warnings in GCC
    T last_written = {};

    bool has_written = false;
    while (!merge_heap.empty()) {
        std::pair<T, uint32_t> smallest = merge_heap.top();
        merge_heap.pop();
        if (!has_written || smallest.first != last_written) {
            has_written = true;
            on_new_item(smallest.first);
            last_written = smallest.first;
            totalSize++;
        }
        if (chunk_files[smallest.second].good()) {
            T data_item;
            if (chunk_files[smallest.second].read(reinterpret_cast<char *>(&data_item),
                                                  sizeof(data_item))) {
                merge_heap.push({ data_item, smallest.second });
            }
        }
    }
    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });
}

} // namespace common
} // namespace mg
