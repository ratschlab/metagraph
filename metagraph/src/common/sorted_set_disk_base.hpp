#ifndef __SORTED_SET_DISK_BASE_HPP__
#define __SORTED_SET_DISK_BASE_HPP__

#include <cassert>
#include <filesystem>
#include <functional>
#include <mutex>
#include <iostream>
#include <optional>
#include <shared_mutex>
#include <vector>

#include <ips4o.hpp>

#include "common/file_merger.hpp"
#include "common/threads/chunked_wait_queue.hpp"
#include "common/vector.hpp"

namespace mg {
namespace common {

/**
 * Abstract thread safe data storage that is able to sort  elements from the
 * underlying Container using external storage to limit the amount of required memory.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is processed using
 * #sort_and_remove_duplicates and written to disk, each disk write into a different file.
 * The #data() method returns the globally sorted data.
 *
 * @tparam T the type of the elements that are being stored and sorted,
 * typically #KMerBOSS instances, or <#KmerBOSS, count> pairs.
 */
template <typename T>
class SortedSetDiskBase {
    typedef T value_type;
    typedef Vector<value_type> storage_type;
    typedef ChunkedWaitQueue<value_type> result_type;
    typedef typename storage_type::iterator Iterator;

    static constexpr uint32_t MERGE_L1_COUNT = 4;

    /**
     * The number of elements in the merge queue. Should be large enough to reduce lock
     * contention, but small enough to not lose time with memory allocation.
     */
    static constexpr size_t QUEUE_EL_COUNT = 30'000;

  public:
    SortedSetDiskBase(std::function<void(storage_type *)> cleanup,
                      size_t num_threads,
                      size_t reserved_num_elements,
                      const std::filesystem::path &tmp_dir,
                      std::function<void(const T &)> on_item_pushed,
                      size_t num_last_elements_cached)
        : num_threads_(num_threads),
          chunk_file_prefix_(tmp_dir / "chunk_"),
          merge_queue_(std::min(reserved_num_elements, QUEUE_EL_COUNT),
                       num_last_elements_cached,
                       on_item_pushed),
          cleanup_(cleanup) {
        std::filesystem::create_directory(tmp_dir);
        if (reserved_num_elements == 0) {
            logger->error("SortedSetDisk buffer cannot have size 0");
            std::exit(EXIT_FAILURE);
        }
        try_reserve(reserved_num_elements);
    }

    virtual ~SortedSetDiskBase() {
        merge_queue_.shutdown(); // make sure the data was processed
    }

    size_t buffer_size() const { return data_.capacity(); }

    /**
     * Returns the globally sorted and counted data. Typically called once all the data
     * was inserted via insert().
     */
    ChunkedWaitQueue<T> &data() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);

        if (!is_merging_) {
            is_merging_ = true;
            // write any residual data left
            if (!data_.empty()) {
                sort_and_remove_duplicates(&data_, num_threads_);
                dump_to_file_async(true /* is_done */);
            }
            start_merging();
        }
        return merge_queue_;
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

  protected:
    virtual void sort_and_remove_duplicates(storage_type *vector,
                                            size_t num_threads) const = 0;

    void start_merging() {
        async_worker_.join(); // make sure all pending data was written
        async_merge_l1_.join(); // make sure all L1 merges are done

        async_worker_.enqueue([this]() {
            const uint32_t l1_chunk_count
                    = (chunk_count_ + MERGE_L1_COUNT - 1) / MERGE_L1_COUNT;
            std::vector<std::string> file_names(l1_chunk_count);
            for (size_t i = 0; i < l1_chunk_count; ++i) {
                file_names[i] = merged_l1_name(chunk_file_prefix_, i);
            }
            std::function<void(const value_type &)> on_new_item
                    = [this](const value_type &v) { merge_queue_.push(v); };
            merge_files(file_names, on_new_item);
            merge_queue_.shutdown();
        });
    }

    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erasing duplicate values...");

        size_t old_size = data_.size();
        sort_and_remove_duplicates(&data_, num_threads_);

        logger->trace("...done. Size reduced from {} to {}, {}MiB", old_size,
                      data_.size(), (data_.size() * sizeof(T) >> 20));
    }

    static inline std::string merged_l1_name(const std::string &prefix, uint32_t count) {
        return prefix + "m" + std::to_string(count);
    }

    template <class value_type>
    static void merge_l1(const std::string &chunk_file_prefix, uint32_t chunk_count) {
        const std::string &merged_l1_file_name
                = merged_l1_name(chunk_file_prefix, chunk_count / MERGE_L1_COUNT);
        std::fstream merged_file(merged_l1_file_name, std::ios::binary | std::ios::out);
        const uint32_t to_merge_count = chunk_count % MERGE_L1_COUNT + 1;
        std::vector<std::string> to_merge(to_merge_count);
        for (uint32_t i = 0; i < to_merge_count; ++i) {
            to_merge[i] = chunk_file_prefix + std::to_string(chunk_count - i);
        }
        logger->trace("Starting merging last {} chunks into {}", to_merge_count,
                      merged_l1_file_name);
        std::function<void(const value_type &)> on_new_item
                = [&merged_file, &merged_l1_file_name](const value_type &v) {
                      if (!merged_file.write(reinterpret_cast<const char *>(&v), sizeof(T))) {
                          std::cerr << "Error: Writing of merged data to "
                                    << merged_l1_file_name << " failed." << std::endl;
                          std::exit(EXIT_FAILURE);
                      }
                  };
        merge_files(to_merge, on_new_item, true /* clean up */);
        logger->trace("Merging last {} chunks into {} done", to_merge_count,
                      merged_l1_file_name);
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
    static void dump_to_file(ThreadPool *threadPool,
                             const std::string &chunk_file_prefix,
                             storage_type *data,
                             uint32_t chunk_count,
                             bool is_done) {
        std::string file_name = chunk_file_prefix + std::to_string(chunk_count);
        std::fstream binary_file(file_name, std::ios::out | std::ios::binary);
        if (!binary_file) {
            logger->error("Error: Creating chunk file '{}' failed", file_name);
            std::exit(EXIT_FAILURE);
        }
        if (!binary_file.write((char *)&((*data)[0]), sizeof((*data)[0]) * data->size())) {
            logger->error("Error: Writing to '{}' failed", file_name);
            std::exit(EXIT_FAILURE);
        }
        binary_file.close();
        if (is_done || (chunk_count + 1) % MERGE_L1_COUNT == 0) {
            threadPool->enqueue(merge_l1<typename storage_type::value_type>,
                                chunk_file_prefix, chunk_count);
        }
        data->resize(0);
    }

    /**
     * Write the current chunk to disk, asynchronously.
     * While the current chunk is being written to disk, the class may accept new
     * insertions.
     * @param is_done true if no more data is added (this is the last call to this
     * function)
     */
    void dump_to_file_async(bool is_done) {
        async_worker_.join(); // wait for other thread to finish writing
        assert(!data_.empty());

        data_dump_.swap(data_);
        async_worker_.enqueue(dump_to_file<storage_type>, &async_merge_l1_,
                              chunk_file_prefix_, &data_dump_, chunk_count_, is_done);
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

    /** Advances #it by step or points to #end, whichever comes first. */
    template <typename Iterator>
    Iterator safe_advance(const Iterator &it, const Iterator &end, size_t step) {
        assert(it <= end);
        return static_cast<size_t>(end - it) < this->buffer_size() ? end : it + step;
    }

    /**
     * Prepare inserting the data between #begin and #end into the buffer. If the
     * buffer is full, the data is processed using #shink_data(), then written to disk,
     * after which the buffer is cleared.
     */
    template <class Iterator>
    size_t prepare_insert(Iterator begin, Iterator end) {
        uint64_t batch_size = end - begin;

        assert(begin <= end && batch_size <= data_.capacity());

        if (data_.size() + batch_size > data_.capacity()) { // time to write to disk
            std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
            shrink_data();

            dump_to_file_async(false /* is_done */);
        }

        size_t offset = data_.size();
        // resize to the future size after insertion (no reallocation will happen)
        data_.resize(data_.size() + batch_size);

        return offset;
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
    const std::string chunk_file_prefix_;
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

    /**
     * Thread pool for doing the "level 1" merging, i.e. merging #MERGE_L1_COUNT chunk
     * files at a time, while new files are still being added to the data structure.
     */
    // TODO: this is suboptimal, as ThreadPool imposes a max of 5*threads tasks
    ThreadPool async_merge_l1_ = ThreadPool(1, 100);

    ChunkedWaitQueue<T> merge_queue_;

    std::function<void(storage_type *)> cleanup_;
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_DISK_BASE_HPP__
