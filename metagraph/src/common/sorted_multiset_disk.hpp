#pragma once

#include <mutex>
#include <shared_mutex>
#include <iostream>
#include <vector>
#include <cassert>

#include <ips4o.hpp>
#include "common/sorted_set_disk_base.hh"
#include "common/threads/chunked_wait_queue.hpp"


namespace mg {
namespace common {

/**
 * Thread safe data storage that is able to sort and count elements from the underlying
 * Container using external storage to limit the amount of required memory.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is sorted, merged, and written to disk,
 * each disk write into a different file. The #data() method returns the
 * globally sorted data.
 *
 * @tparam T the type of the elements that are being stored, sorted and counted,
 * typically #KMerBOSS instances
 * @param C the type used to count the number of appearences of each element of type T
 */
template <typename T, typename C = uint8_t>
class SortedMultisetDisk : public SortedDiskBase<std::pair<T, C>> {
  public:
    typedef std::pair<T, C> value_type;
    typedef Vector<value_type> storage_type;
    typedef ChunkedWaitQueue<value_type> result_type;
    typedef typename storage_type::iterator Iterator;

    /**
     * Constructs a SortedMultisetDisk instance and initializes its buffers sizes to the
     * value specified in #reserved_num_elements.
     * @param cleanup function to run each time a chunk is written to disk; typically
     * performs cleanup operations, such as removing redundant dummy source k-mers
     * @param num_threads the number of threads to use by the sorting algorithm
     * @param chunk_file_prefix the prefix of the temporary files where chunks are
     * written before being merged
     * @param container_size the size of the in-memory container that is written
     * to disk when full
     */
    SortedMultisetDisk(
            std::function<void(storage_type *)> cleanup = [](storage_type *) {},
            size_t num_threads = 1,
            size_t reserved_num_elements = 1e6,
            const std::string &chunk_file_prefix = "/tmp/chunk_",
            std::function<void(const value_type &)> on_item_pushed
            = [](const value_type &) {},
            size_t num_last_elements_cached = 100)
        : SortedDiskBase<std::pair<T, C>>(cleanup,
                                          num_threads,
                                          reserved_num_elements,
                                          chunk_file_prefix,
                                          on_item_pushed,
                                          num_last_elements_cached) {}

    static constexpr uint64_t max_count() { return std::numeric_limits<C>::max(); }

    /**
     * Insert the data between #begin and #end into the buffer. If the buffer is
     * full, the data is sorted, counted and written to disk, after which the
     * buffer is cleared.
     */
    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        std::optional<size_t> offset = this->prepare_insert(begin, end);

        // different threads will insert to different chunks of memory, so it's okay
        // (and desirable) to allow concurrent inserts
        std::shared_lock<std::shared_timed_mutex> multi_insert_lock(this->multi_insert_mutex_);
        if (offset) {
            if constexpr (std::is_same<T, std::remove_cv_t<std::remove_reference_t<decltype(*begin)>>>::value) {
                std::transform(begin, end, this->data_.begin() + offset.value(),
                               [](const T &value) { return std::make_pair(value, C(1)); });
            } else {
                std::copy(begin, end, this->data_.begin() + offset.value());
            }
        }
    }

  protected:
    virtual void sort_and_merge_duplicates(storage_type *vector, size_t num_threads) const {
        assert(vector);
        ips4o::parallel::sort(
                vector->begin(), vector->end(),
                [](const value_type &first, const value_type &second) {
                    return first.first < second.first;
                },
                this->num_threads_);

        auto first = vector->begin();
        auto last = vector->end();

        auto dest = first;

        while (++first != last) {
            if (first->first == dest->first) {
                if (first->second < max_count() - dest->second) {
                    dest->second += first->second;
                } else {
                    dest->second = max_count();
                }
            } else {
                *++dest = std::move(*first);;
            }
        }

        vector->erase(++dest, this->data_.end());

        this->cleanup_(vector);
    }
};

} // namespace common
} // namespace mg
