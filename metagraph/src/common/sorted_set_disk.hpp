#ifndef __SORTED_SET_DISK_HPP__
#define __SORTED_SET_DISK_HPP__

#include <cassert>
#include <fstream> //TODO(ddanciu) - try boost mmapped instead
#include <mutex>
#include <queue>
#include <shared_mutex>

#include <ips4o.hpp>

#include "common/chunked_wait_queue.hpp"
#include "common/file_merger.hpp"
#include "common/logger.hpp"
#include "common/threading.hpp"
#include "utils/vectors.hpp"

namespace mg {
namespace common {

/**
 * Thread safe data storage that is able to sort and extract distinct elements
 * from the underlying Container using external storage to avoid memory overlow.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is sorted and written to disk,
 * each disk write into a different file. The #data() method returns the
 * globally sorted data.
 *
 * @tparam T the type of the elements that are being stored and sorted,
 * typically #KMerBOSS instances
 * */
template <typename T>
class SortedSetDisk {
  public:
    typedef T key_type;
    typedef T value_type;
    typedef Vector<T> storage_type;
    typedef ChunkedWaitQueue<T> result_type;

    /**
     * Constructs a SortedSetDisk instance and initializes its buffer to the value
     * specified in CONTAINER_SIZE_BYTES.
     * @param cleanup function to run each time a chunk is written to disk; typically
     * performs cleanup operations, such as removing redundant dummy source k-mers
     * @param num_threads the number of threads to use by the sorting algorithm
     * @param chunk_file_prefix the prefix of the temporary files where chunks are
     * written before being merged
     * @param container_size the size of the in-memory container that is written
     * to disk when full
     */
    SortedSetDisk(
            std::function<void(storage_type *)> cleanup = [](storage_type *) {},
            size_t num_threads = 1,
            size_t reserved_num_elements = 1e6,
            const std::string &chunk_file_prefix = "/tmp/chunk_",
            std::function<void(const T &)> on_item_pushed = [](const T &) {},
            size_t num_last_elements_cached = 100)
        : num_threads_(num_threads),
          chunk_file_prefix_(chunk_file_prefix),
          merge_queue_(reserved_num_elements ,
                       num_last_elements_cached,
                       on_item_pushed),
          cleanup_(cleanup) {
        if (reserved_num_elements == 0) {
            logger->error("SortedSetDisk buffer cannot have size 0");
            std::exit(EXIT_FAILURE);
        }
        try_reserve(reserved_num_elements);
    }

    ~SortedSetDisk() {
        merge_queue_.shutdown(); // make sure the data was processed
    }

    /**
     * Insert the data between #begin and #end into the buffer. If the buffer is
     * full, the data is sorted, de-duped and written to disk, after which the
     * buffer is cleared.
     */
    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        assert(begin <= end);

        uint64_t batch_size = end - begin;

        // acquire the mutex to restrict the number of writing threads
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        if (data_.size() + batch_size > data_.capacity()) { // time to write to disk
            std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
            // TODO(ddanciu) - test if it's worth it to keep adding in memory if the
            // shrinking was significant
            shrink_data();

            dump_to_file_async();
        }

        size_t offset = data_.size();
        // resize to the future size after insertion (no reallocation will happen)
        data_.resize(data_.size() + batch_size);
        std::shared_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        exclusive_lock.unlock();
        // different threads will insert to different chunks of memory, so it's okay
        // (and desirable) to allow concurrent inserts
        std::copy(begin, end, data_.begin() + offset);
    }

    /**
     * Returns the globally sorted and de-duped data. Typically called once all
     * the data was inserted via insert().
     */
    ChunkedWaitQueue<T> &data() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);

        if (!is_merging_) {
            is_merging_ = true;
            // write any residual data left
            if (!data_.empty()) {
                sort_and_remove_duplicates(&data_, num_threads_);
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
            merge_files<T>(file_names, [this](const T &v) { merge_queue_.push(v); });
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

    template <class Array>
    void sort_and_remove_duplicates(Array *vector, size_t num_threads) const {
        assert(vector);

        ips4o::parallel::sort(vector->begin(), vector->end(),
                              std::less<typename Array::value_type>(), num_threads);
        // remove duplicates
        auto unique_end = std::unique(vector->begin(), vector->end());
        vector->erase(unique_end, vector->end());

        cleanup_(vector); // typically removes source dummy k-mers
    }

    size_t buffer_size() const { return data_.capacity(); }

  private:
    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erasing duplicate values...");

        size_t old_size = data_.size();
        sort_and_remove_duplicates(&data_, num_threads_);

        logger->trace("...done. Size reduced from {} to {}, {}MiB", old_size,
                      data_.size(), (data_.size() * sizeof(T) >> 20));
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
     * insertions which will be written into the other #data_ buffer.
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

    ChunkedWaitQueue<T> merge_queue_;

    std::function<void(storage_type *)> cleanup_;
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_DISK_HPP__
