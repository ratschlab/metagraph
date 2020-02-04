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
 * Profiling shows that using a sorted vector instead of a std::priority queue is
 * faster if the queue has less than ~45 elements, which is the case for us, as each
 * element represents a 1GB chunk, and SRAs typically expand to ~15 chunks.
 */
template <typename T>
class VectorHeap {
    /** The heap stores triplets of the form <Element, Count, SourceIndex> */
    using value_type = std::pair<T, uint32_t>;

  public:
    VectorHeap(size_t size) { els.resize(size); }
    void emplace(T el, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), value_type(el, idx),
                                   [](const value_type &a, const value_type &b) {
                                       return a.first > b.first;
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
    uint64_t num_elements = 0;

    VectorHeap<T> merge_heap(sources.size());
    T data_item;
    // profiling iidicates that setting a larger buffer slightly increases performance
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
            num_elements++;
        }
    }

    // initialized to suppress maybe-uninitialized warnings in GCC
    T last_written = {};

    bool has_written = false;

    while (!merge_heap.empty()) {
        std::pair<T, uint32_t> smallest = merge_heap.pop();

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

    delete[] buffer;
    return num_elements;
}

/**
 * A heap that merges elements that are equal by adding their counts.
 * The heap uses a vector as the underlying data structure, thus making it efficient
 * only for a small (<100) number of elements.
 * //TODO(ddanciu): profile is using an actual binary heap is faster in practice
 * @tparam T the actual heap element type; this is the type tested for equality
 * @tparam C the count type for elements in the heap (typically some unsigned integer)
 */
 //TODO(ddanciu) - instead of using MergingHeap, we can use a simple
 // std::priority_queue which stores duplicate elements and merge elements when popping
 // Profile which one is faster. Also profile if using an  std::vector as in MergeQueue
 // is faster than using an std::queue.
template <typename T, typename C>
class MergingHeap {
    /** The heap stores triplets of the form <Element, Count, SourceIndex> */
    using value_type = std::tuple<T, C, uint32_t> ;

  public:
    bool emplace(T el, C count, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), value_type(el, count, idx),
                                   [](const value_type &a, const value_type &b) {
                                       return std::get<0>(a) > std::get<0>(b);
                                   });
        if (it != els.end() && el == std::get<0>(*it)) {
            C &value = std::get<1>(*it);
            if (value < std::numeric_limits<C>::max() - count) {
                value += count;
            } else {
                value = std::numeric_limits<C>::max();
            }
            return true;
        }
        els.emplace(it, el, count, idx);
        return false;
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
};

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
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> sources,
                     std::function<void(const std::pair<T, C> &)> on_new_item) {
    // Convenience type for storing a pair and it's source file index
    using CountedEl = std::tuple<T, C, uint32_t>;
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::ifstream> chunk_files(sources.size());
    uint64_t num_elements = 0;

    MergingHeap<T, C> merge_heap;
    std::pair<T, C> data_item;
    char *buffer = new char[sources.size() * 1024 * 1024];
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].rdbuf()->pubsetbuf((buffer + i * 1024 * 1024), 1024 * 1024);
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        bool found = true;
        while (found) {
            if (chunk_files[i].read(reinterpret_cast<char *>(&data_item),
                                    sizeof(std::pair<T, C>))) {
                found = merge_heap.emplace(data_item.first, data_item.second, i);
                num_elements++;
            } else {
                found = false;
            }
        }
    }

    while (!merge_heap.empty()) {
        CountedEl smallest = merge_heap.pop();
        on_new_item({ std::get<0>(smallest), std::get<1>(smallest) });

        bool found = true;
        uint32_t chunk_index = std::get<2>(smallest);
        while (found && chunk_files[chunk_index]
               && chunk_files[chunk_index].read(reinterpret_cast<char *>(&data_item),
                                                sizeof(std::pair<T, C>))) {
            found = merge_heap.emplace(data_item.first, data_item.second, chunk_index);
            num_elements++;
        }
    }

    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });
    delete[] buffer;

    return num_elements;
}

} // namespace common
} // namespace mg
