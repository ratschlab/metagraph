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
 * Decoder that reads data from several sorted files and merges it into a single sorted
 * stream.
 * @tparam T the type of data being stored
 */
template <typename T>
class MergeDecoder {
  public:
    MergeDecoder(const std::vector<std::string> &source_names, bool remove_sources) {
        sources_.reserve(source_names.size());
        for (uint32_t i = 0; i < source_names.size(); ++i) {
            sources_.emplace_back(source_names[i], remove_sources);
            std::optional<T> data_item = sources_.back().next();
            if (data_item.has_value()) {
                heap_.emplace(data_item.value(), i);
            }
        }
    }

    std::optional<T> top() {
        if (heap_.empty()) {
            return std::nullopt;
        }
        return heap_.top().first;
    }

    std::optional<T> next() {
        if (heap_.empty()) {
            return std::nullopt;
        }
        auto [result, chunk_index] = heap_.pop();
        std::optional<T> data_item = sources_[chunk_index].next();
        if (data_item.has_value()) {
            heap_.emplace(data_item.value(), chunk_index);
        }
        return result;
    }

  private:
    std::vector<EliasFanoDecoder<T>> sources_;
    common::MergeHeap<T> heap_;
};

/**
 * Merges Elias-Fano sorted compressed files into a single stream.
 */
template <typename T>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    MergeDecoder<T> decoder = MergeDecoder<T>(sources, remove_sources);
    size_t num_elements_read = 0;
    while (decoder.top().has_value()) {
        on_new_item(decoder.next().value());
        num_elements_read++;
    }
    return num_elements_read;
}

/**
 * Given a list of n source files, containing ordered pairs of  <element, count>,
 * merge the n sources (and the corresponding counts) into a single list, ordered by el.
 * If two pairs have the same first element, the counts are added together.
 * @param sources the files containing sorted lists of pairs of type <T, C>
 * @param on_new_item callback to invoke when a new element was merged
 * @param remove_sources if true, remove source files after merging
 *
 * @return the total number of (merged) elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    MergeDecoder<std::pair<T, C>> decoder
            = MergeDecoder<std::pair<T, C>>(sources, remove_sources);
    if (!decoder.top().has_value()) {
        return 0;
    }
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 1;
    std::pair<T, C> current = decoder.next().value();
    while (decoder.top().has_value()) {
        std::pair<T,C> last = decoder.top().value();
        if (current != last) {
            on_new_item(current);
            num_elements_read++;
            current = last;
        } else {
            if (current.second < std::numeric_limits<C>::max() - last.second) {
                current.second += last.second;
            } else {
                current.second = std::numeric_limits<C>::max();
            }
        }
    }
    on_new_item(current);

    return num_elements_read;
}

/**
 * Merges the <T, C> pairs in #sources with the Ts in #source_no_count. The INTs in
 * source_no_count will be assigned a count of 0.
 */
template <typename T>
uint64_t merge_dummy(const std::vector<std::string> &source,
                     std::vector<std::string> source_no_count,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    source_no_count.insert(source_no_count.end(), source.begin(), source.end());
    return merge_files(source_no_count, on_new_item, remove_sources);
}

/**
 * Merges the <T, C> pairs in #sources with the Ts in #source_no_count. The INTs in
 * source_no_count will be assigned a count of 0.
 */
template <typename T, typename C>
uint64_t merge_dummy(const std::vector<std::string> &sources,
                     const std::vector<std::string> &sources_no_count,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    uint64_t num_elements_read = 0;
    MergeDecoder<std::pair<T, C>> decoder
            = MergeDecoder<std::pair<T, C>>(sources, remove_sources);
    MergeDecoder<T> decoder_no_count = MergeDecoder<T>(sources_no_count, remove_sources);
    while (decoder.top().has_value()) {
        std::pair<T,C> next = decoder.next().value();
        while (decoder_no_count.top().has_value() && decoder_no_count.top().value() < next.first) {
            on_new_item({decoder_no_count.next().value(), 0U});
            num_elements_read++;
        }
        on_new_item(next);
        num_elements_read++;
    }
    while (decoder_no_count.top().has_value()) {
        on_new_item({ decoder_no_count.next().value(), 0U });
        num_elements_read++;
    }

    return num_elements_read;
}

} // namespace common
} // namespace mg
