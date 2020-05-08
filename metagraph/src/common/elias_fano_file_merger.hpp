#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "common/elias_fano.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"

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
 * sources into a single (ordered) list of type T and delete the original files.
 * @tparam T the type of the  elements to be merged (typically a 64/128 or 256-bit k-mer)
 * @param sources the files containing sorted lists of type T
 * @param on_new_item callback to invoke when a new element was merged
 * @param remove_sources if true, remove source files after merging
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 0;

    MergeHeap<T> merge_heap;
    std::optional<T> data_item;

    std::vector<EliasFanoDecoder<T>> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.emplace_back(sources[i], remove_sources);
        data_item = decoders.back().next();
        if (data_item.has_value()) {
            merge_heap.emplace(data_item.value(), i);
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
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

        if ((data_item = decoders[chunk_index].next()).has_value()) {
            merge_heap.emplace(data_item.value(), chunk_index);
            num_elements_read++;
        }
    }

    return num_elements_read;
}

// TODO: these two `merge_files` are almost identical. Merge them into one.
//       Implement the merging mechanism  (remove duplicates, increment
//       counters) in the caller?
/**
 * Given a list of n source files, containing ordered pairs of  <element, count>,
 * merge the n sources (and the corresponding counts) into a single list, ordered by el
 * and delete the original files.
 * If two pairs have the same first element, the counts are added together.
 * @param sources the files containing sorted lists of pairs of type <T, C>
 * @param on_new_item callback to invoke when a new element was merged
 * @param remove_sources if true, remove source files after merging
 *
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 0;

    MergeHeap<std::pair<T, C>, utils::GreaterFirst> merge_heap;
    std::optional<std::pair<T, C>> data_item;
    std::vector<EliasFanoDecoder<std::pair<T, C>>> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.emplace_back(sources[i], remove_sources);
        data_item = decoders.back().next();
        if (data_item.has_value()) {
            merge_heap.emplace(data_item.value(), i);
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
        return num_elements_read;
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

        if ((data_item = decoders[chunk_index].next()).has_value()) {
            merge_heap.emplace(data_item.value(), chunk_index);
            num_elements_read++;
        }
    }
    on_new_item(current);

    return num_elements_read;
}

/**
 * Merges the Ts in #source with the Ts in #source_zero_count. This is no different than
 * calling merge() for all files.
 * @param source name of a source file containing Elias-Fano encoded INTs
 * @param source_zero_count name of a soruce file containing Elias-Fano encoded INTs corresponding to dummy k-mers
 * @param on_new_item callback to invoke for each merged item
 * @param remove_sources if true, the #source and #source_zero_count files will be removed
 */
template <typename T>
uint64_t merge_dummy(const std::string &source,
                     std::vector<std::string> source_zero_count,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    source_zero_count.push_back(source);
    return merge_files(source_zero_count, on_new_item, remove_sources);
}

/**
 * Merges the <T, C> pairs in #source with the Ts in #source_zero_count. The INTs in
 * source_zero_count will be assigned a count of 0.
 */
template <typename T, typename C>
uint64_t merge_dummy(const std::string &source,
                     const std::vector<std::string> &source_zero_count,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    uint64_t num_elements_read = 0;

    EliasFanoDecoder<std::pair<T,C>> decoder(source, remove_sources);
    std::optional<std::pair<T, C>> data_item = decoder.next();

    merge_files(source_zero_count,
                [&data_item, &on_new_item, &num_elements_read, &decoder](const T &v) {
                    num_elements_read++;
                    while (data_item.has_value() && data_item.value().first < v) {
                        num_elements_read++;
                        on_new_item(data_item.value());
                        data_item = decoder.next();
                    }
                    assert(!data_item.has_value() || data_item.value().first != v);
                    on_new_item({ v, 0 });
                });

    // add the leftover contents from #source
    while (data_item.has_value()) {
        on_new_item(data_item.value());
        data_item = decoder.next();
    }

    return num_elements_read;
}


} // namespace common
} // namespace mg
