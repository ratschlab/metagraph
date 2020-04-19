#pragma once

#include <cassert>
#include <functional>
#include <optional>
#include <shared_mutex>
#include <string>

#include <ips4o.hpp>

#include "common/sorted_set_disk_base.hpp"

namespace mg {
namespace common {

/**
 * Specialization of SortedSetDiskBase that is able to both sort and count elements.
 *
 * @tparam T the type of the elements that are being stored, sorted and counted,
 * typically k-mers
 * @param INT the integer representation of the element being sorted (this representation
 * is used for writing to disk, using Elias-Fano compression)
 * @param C the type used to count the multiplicity of each value in the multi-set
 */
template <typename T, typename C = uint8_t>
class SortedMultisetDisk : public SortedSetDiskBase<std::pair<T, C>> {
  public:
    typedef T key_type;
    typedef C count_type;
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
     * @param tmp_dir the prefix of the temporary files where chunks are
     * written before being merged
     * @param container_size the size of the in-memory container that is written
     * to disk when full
     */
    SortedMultisetDisk(
            std::function<void(storage_type *)> cleanup = [](storage_type *) {},
            size_t num_threads = 1,
            size_t reserved_num_elements = 1e6,
            const std::filesystem::path &tmp_dir = "/tmp/",
            size_t max_disk_space_bytes = 1e9,
            size_t num_last_elements_cached = 100)
        : SortedSetDiskBase<value_type>(cleanup,
                                        num_threads,
                                        reserved_num_elements,
                                        tmp_dir,
                                        max_disk_space_bytes,
                                        num_last_elements_cached) {}

    static constexpr uint64_t max_count() { return std::numeric_limits<C>::max(); }

    /**
     * Insert the data between #begin and #end into the buffer.
     */
    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        Iterator original_end = end;
        end = this->safe_advance(begin, original_end, this->buffer_size());
        while (begin != end) {
            // acquire the mutex to restrict the number of writing threads
            std::unique_lock<std::mutex> exclusive_lock(this->mutex_);

            size_t offset = this->prepare_insert(begin, end);

            std::shared_lock<std::shared_timed_mutex> multi_insert_lock(
                    this->multi_insert_mutex_);
            // different threads will insert to different chunks of memory, so it's okay
            // (and desirable) to allow concurrent inserts
            exclusive_lock.unlock();
            if constexpr (std::is_same<T, std::decay_t<decltype(*begin)>>::value) {
                std::transform(begin, end, this->data_.begin() + offset, [](const T &value) {
                    return std::make_pair(value, static_cast<C>(1));
                });
            } else {
                std::copy(begin, end, this->data_.begin() + offset);
            }
            begin = end;
            end = this->safe_advance(begin, original_end, this->buffer_size());
        }
    }

  private:
    virtual void sort_and_remove_duplicates(storage_type *vector,
                                            size_t num_threads) const override {
        assert(vector);
        ips4o::parallel::sort(
                vector->begin(), vector->end(),
                [](const value_type &first, const value_type &second) {
                    return first.first < second.first;
                },
                num_threads);

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
                *++dest = std::move(*first);
            }
        }

        vector->erase(++dest, this->data_.end());

        this->cleanup_(vector);
    }
};

} // namespace common
} // namespace mg
