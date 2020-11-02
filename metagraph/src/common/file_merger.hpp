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
        els.push_back({ el, idx });
        std::push_heap(els.begin(), els.end(), cmp_);
    }

    inline const value_type& top() const { return els.front(); }

    inline value_type pop() {
        value_type result = els.front();
        std::pop_heap(els.begin(), els.end(), cmp_);
        els.pop_back();
        return result;
    }

    inline bool empty() const { return els.empty(); }

  private:
    // elements stored in decreasing order of the first tuple member
    std::vector<value_type> els;
    Compare compare_ = Compare();
    std::function<int(const value_type &a, const value_type &b)> cmp_ =
            [&](const value_type &a, const value_type &b) {
                return compare_(a.first, b.first);
            };
};


/**
 * Decoder that reads data from several sorted files and merges it into a single sorted
 * stream.
 * @tparam T the type of data being stored
 */
template <typename T>
class FileMerger {
  public:
    typedef T value_type;

    FileMerger(const std::vector<std::string> &source_names) {
        sources_.reserve(source_names.size());
        for (uint32_t i = 0; i < source_names.size(); ++i) {
            sources_.emplace_back(source_names[i], std::ios::binary);
            T v;
            if (sources_.back().read(reinterpret_cast<char *>(&v), sizeof(T))) {
                heap_.emplace(v, i);
            }
        }
    }

    inline bool empty() const { return heap_.empty(); }

    inline const T& top() const {
#ifndef NDEBUG
        if (heap_.empty())
            throw std::runtime_error("Popping an empty FileMerger");
#endif
        return heap_.top().first;
    }

    inline T pop() {
#ifndef NDEBUG
        if (heap_.empty())
            throw std::runtime_error("Popping an empty FileMerger");
#endif
        auto [result, source_index] = heap_.pop();
        T v;
        if (sources_[source_index].read(reinterpret_cast<char *>(&v), sizeof(T))) {
            heap_.emplace(v, source_index);
        }
        return result;
    }

  private:
    std::vector<std::ifstream> sources_;
    common::MergeHeap<T> heap_;
};

/**
 * Merges sorted file into a single stream.
 */
template <typename T>
void merge_files(const std::vector<std::string> &sources,
                 const std::function<void(const T &)> &on_new_item) {
    FileMerger<T> decoder(sources);
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

