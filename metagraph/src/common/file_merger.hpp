#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"

namespace mg {
namespace common {

/**
 * Heap implemented as a sorted vector.
 */
// Note: profiling shows that using a sorted vector instead of a std::priority queue is
// faster up to heaps with 1000 elements.  Using an unsorted vector (faster insert,
// slower pop()) is ~40% slower. Preventing duplicates in the heap so that we don't
// need to test for dupes at pop time is ~60% slower.
template <typename T, class Compare = std::less<T>>
class VectorHeap {
    /** The heap stores triplets of the form <Element, Count, SourceIndex> */
    using value_type = std::pair<T, uint32_t>;

  public:
    void emplace(T el, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), value_type(el, idx),
                                   [this](const value_type &a, const value_type &b) {
                                       return compare_(a.first, b.first);
                                   });
        els.emplace(it, el, idx);
    }

    value_type pop() {
        value_type result = els.back();
        els.pop_back();
        return result;
    }

    bool empty() { return els.empty(); }

  private:
    // elements stored in decreasing order of the first tuple member
    std::vector<value_type> els;
    Compare compare_ = Compare();
};


/**
 * Given a list of n source files, containing ordered elements of type T, merge the n
 * sources into a single (ordered) list and delete the original files.
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
    uint64_t num_elements_read = 0;

    VectorHeap<T, std::greater<T>> merge_heap;
    T data_item;
    // profiling indicates that setting a larger buffer slightly increases performance
    char *buffer = new char[sources.size() * 1024 * 1024];
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].rdbuf()->pubsetbuf((buffer + i * 1024 * 1024), 1024 * 1024);
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Error: Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        if (chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(T))) {
            merge_heap.emplace(data_item, i);
            num_elements_read++;
        }
    }

    std::optional<T> last_written;
    while (!merge_heap.empty()) {
        std::pair<T, uint32_t> smallest = merge_heap.pop();

        if (!last_written.has_value() || smallest.first != last_written.value()) {
            on_new_item(smallest.first);
            last_written = smallest.first;
        }

        if (chunk_files[smallest.second]
            && chunk_files[smallest.second].read(reinterpret_cast<char *>(&data_item),
                                                 sizeof(T))) {
            merge_heap.emplace(data_item, smallest.second);
            num_elements_read++;
        }
    }

    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });

    delete[] buffer;
    return num_elements_read;
}


/**
 * Given a list of n source files, containing ordered pairs of  <element, count>,
 * merge the n sources (and the corresponding counts) into a single list, ordered by el
 * and delete the original files.
 * If two pairs have the same first element, the counts are added together.
 * @param sources the files containing sorted lists of pairs of type <T, C>
 * @param on_new_item callback to invoke when a new element was merged
 *
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
// Note: Profiling indicates that merging elements within the heap is ~30% slower than
// inserting duplicates and merging them when popping from the heap, as we do now
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> sources,
                     std::function<void(const std::pair<T, C> &)> on_new_item) {
    // Convenience type for storing a pair and it's source file index
    using El = std::pair<T, C>;
    using CountedEl = std::pair<El, uint32_t>;
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::ifstream> chunk_files(sources.size());
    uint64_t num_elements = 0;

    VectorHeap<El, utils::GreaterFirst> merge_heap;
    El data_item;
    std::unique_ptr<char[]> buffer(new char[sources.size() * 1024 * 1024]);
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].rdbuf()->pubsetbuf((buffer.get() + i * 1024 * 1024), 1024 * 1024);
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        if (chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(El))) {
            merge_heap.emplace(data_item, i);
            num_elements++;
        }
    }

    std::optional<CountedEl> current;
    while (!merge_heap.empty()) {
        CountedEl smallest = merge_heap.pop();

        if (current.has_value() && smallest.first.first != current.value().first.first) {
            on_new_item(current.value().first);
            current = smallest;
        } else {
            if (current.has_value()) {
                current.value().first.second += smallest.first.second;
            } else {
                current = smallest;
            }
        }

        uint32_t chunk_index = smallest.second;
        if (chunk_files[chunk_index]
            && chunk_files[chunk_index].read(reinterpret_cast<char *>(&data_item),
                                             sizeof(El))) {
            merge_heap.emplace(data_item, chunk_index);
            num_elements++;
        }
    }
    if (current.has_value()) {
        on_new_item(current.value().first);
    }

    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });

    return num_elements;
}

} // namespace common
} // namespace mg
