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
 *
 * The heap uses a vector as the underlying data structure, thus making it efficient
 * only for a small (<1000) number of elements.
 * @tparam T the actual heap element type; this is the type tested for equality
 */
// Note: profiling shows that using a sorted vector instead of a std::priority queue is
// faster up to ~1000 elements.  Using an unsorted vector (faster insert,
// slower pop()) is ~40% slower. Preventing duplicates in the heap so that we don't
// need to test for dupes at pop time is ~60% slower.
// Note: Profiling indicates that merging elements within the heap is ~30% slower than
// inserting duplicates and merging them when popping from the heap, as we do now
template <typename T, class Compare = std::greater<T>>
class MergeHeap {
    /** The heap stores pairs <Element, SourceIndex> */
    using value_type = std::pair<T, uint32_t>;

  public:
    void emplace(T el, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), el,
                                   [this](const value_type &p, const T &v) {
                                       return compare_(p.first, v);
                                   });
        els.emplace(it, el, idx);
    }

    const value_type& top() const { return els.back(); }

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
uint64_t merge_files(const std::vector<std::string> &sources,
                     std::function<void(const T &)> on_new_item,
                     bool cleanup = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::ifstream> chunk_files(sources.size());
    uint64_t num_elements_read = 0;

    MergeHeap<T> merge_heap;
    T data_item;
    // profiling indicates that setting a larger buffer slightly increases performance
    std::unique_ptr<char[]> buffer(new char[sources.size() * 1024 * 1024]);
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].rdbuf()->pubsetbuf((buffer.get() + i * 1024 * 1024), 1024 * 1024);
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        if (chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(T))) {
            merge_heap.emplace(data_item, i);
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
        if (cleanup) {
            std::for_each(sources.begin(), sources.end(),
                          [](const std::string &s) { std::filesystem::remove(s); });
        }
        return num_elements_read;
    }

    T last_written = merge_heap.top().first;
    on_new_item(last_written);

    while (!merge_heap.empty()) {
        auto [smallest, chunk_index] = merge_heap.pop();

        if (smallest != last_written) {
            on_new_item(smallest);
            last_written = smallest;
        }

        if (chunk_files[chunk_index]
            && chunk_files[chunk_index].read(reinterpret_cast<char *>(&data_item),
                                             sizeof(T))) {
            merge_heap.emplace(data_item, chunk_index);
            num_elements_read++;
        }
    }

    if (cleanup) {
        std::for_each(sources.begin(), sources.end(),
                      [](const std::string &s) { std::filesystem::remove(s); });
    }

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
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> &sources,
                     std::function<void(const std::pair<T, C> &)> on_new_item,
                     bool cleanup = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::ifstream> chunk_files(sources.size());
    uint64_t num_elements = 0;

    MergeHeap<std::pair<T, C>, utils::GreaterFirst> merge_heap;
    std::pair<T, C> data_item;
    std::unique_ptr<char[]> buffer(new char[sources.size() * 1024 * 1024]);
    for (uint32_t i = 0; i < sources.size(); ++i) {
        chunk_files[i].rdbuf()->pubsetbuf((buffer.get() + i * 1024 * 1024), 1024 * 1024);
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        if (!chunk_files[i].good()) {
            logger->error("Unable to open chunk file '{}'", sources[i]);
            std::exit(EXIT_FAILURE);
        }
        if (chunk_files[i].read(reinterpret_cast<char *>(&data_item),
                                sizeof(std::pair<T, C>))) {
            merge_heap.emplace(data_item, i);
            num_elements++;
        }
    }

    if (merge_heap.empty()) {
        if (cleanup) {
            std::for_each(sources.begin(), sources.end(),
                          [](const std::string &s) { std::filesystem::remove(s); });
        }
        return num_elements;
    }

    // initialize the smallest element
    std::pair<T, C> current = { merge_heap.top().first.first, 0 };

    while (!merge_heap.empty()) {
        auto [smallest, chunk_index] = merge_heap.pop();

        if (smallest.first != current.first) {
            on_new_item(current);
            current = smallest;
        } else {
            if (current.second < std::numeric_limits<C>::max() - smallest.second) {
                current.second += smallest.second;
            } else {
                current.second = std::numeric_limits<C>::max();
            }
        }

        if (chunk_files[chunk_index]
            && chunk_files[chunk_index].read(reinterpret_cast<char *>(&data_item),
                                             sizeof(std::pair<T, C>))) {
            merge_heap.emplace(data_item, chunk_index);
            num_elements++;
        }
    }
    on_new_item(current);

    if (cleanup) {
        std::for_each(sources.begin(), sources.end(),
                      [](const std::string &s) { std::filesystem::remove(s); });
    }

    return num_elements;
}

} // namespace common
} // namespace mg
