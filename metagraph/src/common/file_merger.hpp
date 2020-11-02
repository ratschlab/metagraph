#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "common/logger.hpp"

namespace mtg {
namespace common {

/**
 * Heap implemented as a sorted vector.
 *
 * The heap uses a vector as the underlying data structure, thus making it efficient
 * only for a small (<1000) number of elements.
 * @tparam T the actual heap element type; this is the type tested for equality
 */
// Note: profiling shows that using a sorted vector instead of a std::priority_queue is
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
    inline void emplace(T el, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), el,
                                   [this](const value_type &p, const T &v) {
                                     return compare_(p.first, v);
                                   });
        els.emplace(it, el, idx);
    }

    inline const value_type& top() const { return els.back(); }

    inline value_type pop() {
        value_type result = els.back();
        els.pop_back();
        return result;
    }

    inline bool empty() const { return els.empty(); }

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
    typedef T value_type;

    MergeDecoder(const std::vector<std::string> &source_names, bool remove_sources) {
        sources_.reserve(source_names.size());
        for (uint32_t i = 0; i < source_names.size(); ++i) {
            sources_.emplace_back(source_names[i], std::ios::binary);
            T v;
            if (sources_.back().read(&v, sizeof(T))) {
                heap_.emplace(v, i);
            }
        }
    }

    inline bool empty() const { return heap_.empty(); }

    inline const T& top() const {
#ifndef NDEBUG
        if (heap_.empty())
            throw std::runtime_error("Popping an empty MergeDecoder");
#endif
        return heap_.top().first;
    }

    inline T pop() {
#ifndef NDEBUG
        if (heap_.empty())
            throw std::runtime_error("Popping an empty MergeDecoder");
#endif
        auto [result, source_index] = heap_.pop();
        std::optional<T> data_item = sources_[source_index].next();
        if (data_item.has_value()) {
            heap_.emplace(data_item.value(), source_index);
        }
        return result;
    }

  private:
    std::vector<std::ifstream> sources_;
    common::MergeHeap<T> heap_;
};

/**
 * Merges Elias-Fano sorted compressed files into a single stream.
 */
template <typename T>
void merge_files(const std::vector<std::string> &sources,
                 const std::function<void(const T &)> &on_new_item,
                 bool remove_sources = true) {
    MergeDecoder<T> decoder(sources, remove_sources);
    if (decoder.empty())
        return;

    T last = decoder.pop();
    while (!decoder.empty()) {
        T curr = decoder.pop();
        if (curr != last) {
            on_new_item(last);
            last = curr;
        }
    }
    on_new_item(last);
}

} // namespace common
} // namespace mtg

