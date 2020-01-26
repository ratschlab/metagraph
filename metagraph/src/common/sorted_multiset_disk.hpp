#pragma once

#include <mutex>
#include <shared_mutex>
#include <iostream>
#include <vector>
#include <cassert>

#include <ips4o.hpp>

#include "common/threads/chunked_wait_queue.hpp"
#include "common/file_merger.hpp"
#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/vector.hpp"

namespace mg {
namespace common {

// Thread safe data storage for counting
/**
 * Thread safe data storage that is able to sort and count elements from the underlying
 * Container using external storage to limit the amount of required memory.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is sorted and written to disk,
 * each disk write into a different file. The #data() method returns the
 * globally sorted data.
 *
 * @tparam T the type of the elements that are being stored, sorted and counted,
 * typically #KMerBOSS instances
 * @param C the type used to count the number of appearences of each element of type T
 */
template <typename T, typename C = uint8_t>
class SortedMultisetDisk {
  public:
    typedef T key_type;
    typedef C count_type;
    typedef std::pair<T, C> value_type;
    typedef Vector<value_type> storage_type;
    typedef ChunkedWaitQueue<value_type> result_type;

    /**
     * Constructs a SortedMultisetDisk instance and initializes its buffers to the value
     * specified in #reserved_num_elements.
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
        : num_threads_(num_threads),
          chunk_file_prefix_(chunk_file_prefix),
          merge_queue_(reserved_num_elements, num_last_elements_cached, on_item_pushed),
          cleanup_(cleanup) {
        if (reserved_num_elements == 0) {
            logger->error("SortedSetDisk buffer cannot have size 0");
            std::exit(EXIT_FAILURE);
        }
        try_reserve(reserved_num_elements);
    }

    ~SortedMultisetDisk() {
        merge_queue_.shutdown(); // make sure the data was processed
    }

    static constexpr uint64_t max_count() { return std::numeric_limits<C>::max(); }

    /**
     * Insert the data between #begin and #end into the buffer. If the buffer is
     * full, the data is sorted, counted and written to disk, after which the
     * buffer is cleared.
     */
    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        assert(begin <= end);

        uint64_t batch_size = end - begin;

        if (!batch_size)
            return;

        // acquire the mutex to restrict the number of writing threads
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        assert(batch_size <= data_.capacity()
               && "Batch size exceeded the buffer's capacity.");

        if (data_.size() + batch_size > data_.capacity()) { // time to write to disk
            std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
            shrink_data();
            dump_to_file_async();
        }

        size_t offset = data_.size();
        // resize to the future size after insertion (no reallocation will happen)
        data_.resize(data_.size() + batch_size);
        std::shared_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        // different threads will insert to different chunks of memory, so it's okay
        // (and desirable) to allow concurrent inserts
        exclusive_lock.unlock();

        if constexpr(std::is_same<T, std::remove_cv_t<
                std::remove_reference_t<
                        decltype(*begin)>>>::value) {
            std::transform(begin, end, data_.begin() + offset,
                           [](const T &value) { return std::make_pair(value, C(1)); });
        } else {
            std::copy(begin, end, data_.begin() + offset);
        }
    }

    size_t buffer_size() const { return data_.capacity(); }

    /**
     * Returns the globally sorted and counted data. Typically called once all the data
     * was inserted via insert().
     */
    ChunkedWaitQueue<value_type> &data() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);

        if (!is_merging_) {
            is_merging_ = true;
            // write any residual data left
            if (!data_.empty()) {
                sort_and_merge_duplicates(&data_, num_threads_);
                dump_to_file_async();
            }
            async_worker_.join(); // make sure all pending data was written
            // TODO(ddanciu): instead of writing to file, pass buffer to merge_func to
            //  avoid writing/reading to disk for small graphs
            start_merging();
        }
        return merge_queue_;
    }

    void start_merging() {
        async_worker_.enqueue([this]() {
            std::vector<std::string> file_names(chunk_count_);
            for (size_t i = 0; i < chunk_count_; ++i) {
                file_names[i] = chunk_file_prefix_ + std::to_string(i);
            }
            merge_files<T, C>(file_names,
                              [this](const value_type &v) { merge_queue_.push(v); });
            merge_queue_.shutdown();
        });
    }

    void clear(std::function<void(const T &)> on_item_pushed = [](const T &) {}) {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        is_merging_ = false;
        chunk_count_ = 0;
        data_.resize(0); // this makes sure the buffer is not reallocated
        data_dump_.resize(0);
        merge_queue_.reset(on_item_pushed);
    }

  private:
    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erasing duplicate values...");

        size_t old_size = data_.size();
        sort_and_merge_duplicates(&data_, num_threads_);
        sorted_end_ = data_.size();

        logger->trace("...done. Size reduced from {} to {}, {}MiB", old_size,
                      data_.size(), (data_.size() * sizeof(value_type) >> 20));
    }

    template <class Array>
    void sort_and_merge_duplicates(Array *vector, size_t num_threads) const {
        assert(vector);
        ips4o::parallel::sort(
                vector->begin(), vector->end(),
                [](const value_type &first, const value_type &second) {
                    return first.first < second.first;
                },
                num_threads_);

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

        vector->erase(++dest, data_.end());

        cleanup_(vector);
    }

    /**
     * Dumps the given data to a file, synchronously.
     * @tparam storage_type the container type for the data being dumped (typically
     * std::vector of folly::Vector)
     * @param[in] chunk_count the current chunk number
     * @param[in,out] data the data to be dumped. At the end of the call, the data will be
     * empty, but its allocated memory will remain unaffected
     */
    template <class storage_type>
    static void dump_to_file(const std::string &file_name, storage_type *data) {
        std::fstream binary_file = std::fstream(file_name, std::ios::out | std::ios::binary);
        if (!binary_file) {
            logger->error("Error: Creating chunk file '{}' failed", file_name);
            std::exit(EXIT_FAILURE);
        }
        if (!binary_file.write((char *)&((*data)[0]), sizeof((*data)[0]) * data->size())) {
            logger->error("Error: Writing to '{}' failed", file_name);
            std::exit(EXIT_FAILURE);
        }
        binary_file.close();
        data->resize(0);
    }

    /**
     * Write the current chunk to disk, asynchronously.
     * While the current chunk is being written to disk, the class may accept new
     * insertions.
     */
    void dump_to_file_async() {
        async_worker_.join(); // wait for other thread to finish writing
        assert(!data_.empty());
        std::string file_name = chunk_file_prefix_ + std::to_string(chunk_count_);

        data_dump_.swap(data_);
        async_worker_.enqueue(dump_to_file<storage_type>, file_name, &data_dump_);
        chunk_count_++;
    }

    void try_reserve(size_t size, size_t min_size = 0) {
        size = std::max(size, min_size);
        size_t original_size = size;
        while (size > min_size) {
            try {
                data_.reserve(size);
                data_dump_.reserve(size);
                if (size != original_size) {
                    logger->warn("SortedSetDisk: Requested {}MiB, but only found {}MiB",
                                 (original_size * sizeof(T)) >> 20,
                                 (size * sizeof(T)) >> 20);
                }
                return;
            } catch (const std::bad_alloc &exception) {
                size = min_size + (size - min_size) * 2 / 3;
            }
        }
        data_.reserve(min_size);
        data_dump_.reserve(min_size);
    }

    /**
     * Current number of chunks written to disk. We expect this to be at most in
     * the order of thousands, so a 32 bit integer should suffice for storage.
     */
    uint32_t chunk_count_ = 0;
    /**
     * Hold the data filled in via #insert.
     */
    storage_type data_;
    /**
     * Buffer containing the data that is currently being dumped to disk (while #data_
     * is being filled with new information.
     */
    storage_type data_dump_;

    size_t num_threads_;

    std::string chunk_file_prefix_;
    /**
     * True if the data merging thread was started, and data started flowing into the #merge_queue_.
     */
    bool is_merging_ = false;

    // indicate the end of the preprocessed distinct and sorted values
    uint64_t sorted_end_ = 0;

    /**
     * Ensures mutually exclusive access (and thus thread-safety) to #data.
     */
    mutable std::mutex mutex_;
    /**
     * Mutex that can be acquired by multiple threads that are appending to
     * non-overlapping areas of #data_
     */
    mutable std::shared_timed_mutex multi_insert_mutex_;

    /**
     * Thread pool for writing to disk and  for merging data from disk. Since writing
     * to disk happens before the merging, a single thread is needed.
     */
    ThreadPool async_worker_ = ThreadPool(1, 1);

    ChunkedWaitQueue<value_type> merge_queue_;

    std::function<void(storage_type *)> cleanup_;
};

} // namespace common
} // namespace mg
