#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "utils/template_utils.hpp"
#include "common/logger.hpp"

namespace mg {
namespace common {
/**
 * Given a list of n source files, containing ordered elements of type T, merge the n
 * sources into a single (ordered) list.
 * @param sources the files containing sorted lists of type T
 * @param on_new_item callback to invoke when a new element was merged
 *
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T>
uint64_t merge_files(const std::vector<std::string> sources,
                     std::function<void(const T &)> on_new_item) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::ifstream> chunk_files(sources.size());
    uint64_t num_elements = 0;

    std::priority_queue<std::pair<T, uint32_t>,
                        std::vector<std::pair<T, uint32_t>>,
                        utils::GreaterFirst>
                    merge_heap;
    T data_item;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        // TODO: use a bigger buffer in stream to reduce the number
        // of random reads or insert into the heap a batch of elements at a time
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Error: Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        if (chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(T))) {
            merge_heap.emplace(data_item, i);
            num_elements++;
        }
    }

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
        }

        if (chunk_files[smallest.second]
             && chunk_files[smallest.second].read(reinterpret_cast<char *>(&data_item),
                                                  sizeof(T))) {
            merge_heap.emplace(data_item, smallest.second);
            num_elements++;
        }
    }

    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });

    return num_elements;
}

} // namespace common
} // namespace mg
