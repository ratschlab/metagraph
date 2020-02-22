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
          reserved_num_elements_(reserved_num_elements),
          chunk_file_prefix_(tmp_dir / "chunk_"),
          cleanup_(cleanup),
          on_item_pushed_(on_item_pushed),
          num_last_elements_cached_(num_last_elements_cached) {
        std::filesystem::create_directory(tmp_dir);
        if (reserved_num_elements == 0) {
            logger->error("SortedSetDisk buffer cannot have size 0");
            std::exit(EXIT_FAILURE);
        }
        try_reserve(reserved_num_elements);
    }

    size_t buffer_size() const { return data_.capacity(); }

    /**
     * Returns the globally sorted and counted data. Typically called once all the data
     * was inserted via insert().
     */
    ChunkedWaitQueue<T> &data() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        ChunkedWaitQueue<T> result(std::min(data_.size(), QUEUE_EL_COUNT),
                                   num_last_elements_cached_, on_item_pushed_);
        // write any residual data left
        if (!data_.empty()) {
            sort_and_remove_duplicates(&data_, num_threads_);
            dump_to_file_async(true /* is_done */);
        }
        data_.reserve(0);
        data_dump_.reserve(0);
        start_merging(&result);
        return result;
    }

    void clear(std::function<void(const T &)> on_item_pushed = [](const T &) {},
            const std::filesystem::path& tmp_path = "/tmp/") {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        try_reserve(reserved_num_elements_);
        chunk_count_ = 0;
        on_item_pushed_ = on_item_pushed;
        chunk_file_prefix_ = tmp_path / "chunk_";
        std::filesystem::create_directory(tmp_path);
    }

  protected:
    virtual void sort_and_remove_duplicates(storage_type *vector,
                                            size_t num_threads) const = 0;

    void start_merging(ChunkedWaitQueue<T>* queue) {
        async_worker_.join(); // make sure all pending data was written
        async_merge_l1_.join(); // make sure all L1 merges are done

        async_worker_.enqueue([this, queue]() {
            std::vector<std::string> file_names;
            for (size_t i = 0; i < merged_index_; ++i) {
                file_names.push_back(merged_l1_name(chunk_file_prefix_, i));
            }
            for (size_t i = MERGE_L1_COUNT * merged_index_; i < chunk_count_; ++i) {
                file_names.push_back(chunk_file_prefix_ + std::to_string(i));
            }
            std::function<void(const value_type &)> on_new_item
                    = [queue](const value_type &v) { queue->push(v); };
            merge_files(file_names, on_new_item);
            queue->shutdown();
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
    static void merge_l1(const std::string &chunk_file_prefix,
                         uint32_t chunk_count,
                         uint32_t *merged_index) {
        const std::string &merged_l1_file_name
                = merged_l1_name(chunk_file_prefix, chunk_count / MERGE_L1_COUNT);
        std::fstream merged_file(merged_l1_file_name, std::ios::binary | std::ios::out);
        const uint32_t to_merge_count = chunk_count % MERGE_L1_COUNT + 1;
        std::vector<std::string> to_merge(to_merge_count);
        for (uint32_t i = 0; i < to_merge_count; ++i) {
            to_merge[i] = chunk_file_prefix + std::to_string(chunk_count - i);
        }
        if (to_merge_count == 1) { // small optimization if only 1 file to merge
            std::filesystem::rename(to_merge[0], merged_l1_file_name);
            return;
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
        (*merged_index)++;
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
                             bool is_done,
                             uint32_t *merged_index) {
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
        if (is_done) {
            threadPool->clear();
        } else if ((chunk_count + 1) % MERGE_L1_COUNT == 0) {
            threadPool->enqueue(merge_l1<typename storage_type::value_type>,
                                chunk_file_prefix, chunk_count, merged_index);
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
                              chunk_file_prefix_, &data_dump_, chunk_count_, is_done,
                              &merged_index_);
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
     * The number of L1 merges that were successfully performed.
     */
    __uint32_t merged_index_ = 0;

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

    size_t reserved_num_elements_;

    const std::string chunk_file_prefix_;

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

    std::function<void(storage_type *)> cleanup_;

    /**
     * Function to call when an item is pushed into the final sorted set.
     */
    std::function<void(const T &)> on_item_pushed_;

    /**
     * Number of elements cached by the sorted set for backwards iteration.
     */
    size_t num_last_elements_cached_;
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_DISK_BASE_HPP__
