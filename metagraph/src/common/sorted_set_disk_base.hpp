#ifndef __SORTED_SET_DISK_BASE_HPP__
#define __SORTED_SET_DISK_BASE_HPP__

#include <atomic>
#include <cassert>
#include <filesystem>
#include <functional>
#include <mutex>
#include <iostream>
#include <optional>
#include <shared_mutex>
#include <vector>

#include <ips4o.hpp>

#include "common/elias_fano.hpp"
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
template <typename T, typename INT>
class SortedSetDiskBase {
    typedef T value_type;
    typedef Vector<value_type> storage_type;
    typedef ChunkedWaitQueue<value_type> result_type;
    typedef typename storage_type::iterator Iterator;

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
                      size_t max_disk_space_bytes,
                      std::function<void(const T &)> on_item_pushed,
                      size_t num_last_elements_cached)
        : num_threads_(num_threads),
          reserved_num_elements_(reserved_num_elements),
          max_disk_space_bytes_(max_disk_space_bytes),
          on_item_pushed_(on_item_pushed),
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
                dump_to_file(true /* is_done */);
            }
            Vector<T>().swap(data_); // free up the (usually very large) buffer
            start_merging();
        }
        return merge_queue_;
    }

    /**
     * Returns the files to be merged - useful if the caller prefers to do the merging.
     */
    std::vector<std::string> files_to_merge() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        // write any residual data left
        if (!data_.empty()) {
            sort_and_remove_duplicates(&data_, num_threads_);
            dump_to_file(true /* is_done */);
        }
        return get_file_names();
    }

    void clear(const std::filesystem::path &tmp_path = "/tmp/") {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        is_merging_ = false;
        chunk_count_ = 0;
        l1_chunk_count_ = 0;
        total_chunk_size_bytes_ = 0;
        try_reserve(reserved_num_elements_);
        data_.resize(0); // this makes sure the buffer is not reallocated
        chunk_file_prefix_ = tmp_path / "chunk_";
        std::filesystem::create_directory(tmp_path);
        merge_queue_.reset(on_item_pushed_);
    }

  protected: // TODO: move most of these methods to private before submitting
    virtual void sort_and_remove_duplicates(storage_type *vector,
                                            size_t num_threads) const = 0;

    virtual void start_merging() = 0;

    std::vector<std::string> get_file_names() {
        async_merge_l1_.join(); // make sure all L1 merges are done
        std::vector<std::string> file_names;
        if (merged_all_) {
            file_names.push_back(merged_all_name(chunk_file_prefix_));
        }

        for (size_t i = 0; i < l1_chunk_count_; ++i) {
            file_names.push_back(merged_l1_name(chunk_file_prefix_, i));
        }
        for (size_t i = MERGE_L1_COUNT * l1_chunk_count_; i < chunk_count_; ++i) {
            file_names.push_back(chunk_file_prefix_ + std::to_string(i));
        }
        return file_names;
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

    static inline std::string merged_all_name(const std::string &prefix) {
        return prefix + "_all";
    }

    /**
     * Dumps the given data to a file, synchronously. If the maximum allowd disk size
     * is reached, all chunks will be merged into a single chunk in an effort to reduce
     * disk space.
     * @param is_done if this is the last chunk being dumped
     */
    virtual void dump_to_file(bool is_done) = 0;

    void try_reserve(size_t size, size_t min_size = 0) {
        size = std::max(size, min_size);
        size_t original_size = size;
        while (size > min_size) {
            try {
                data_.reserve(size);
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

            dump_to_file(false /* is_done */);
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
    std::atomic<uint32_t > l1_chunk_count_ = 0;

    /**
     * Hold the data filled in via #insert.
     */
    storage_type data_;

    size_t num_threads_;

    size_t reserved_num_elements_;

    size_t max_disk_space_bytes_;

    std::function<void(const T &)> on_item_pushed_;

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
     * Thread for merging data from disk.
     */
    ThreadPool async_worker_ = ThreadPool(1, 1);

    ChunkedWaitQueue<T> merge_queue_;

    /**
     * Thread pool for doing the "level 1" merging, i.e. merging #MERGE_L1_COUNT chunk
     * files at a time, while new files are still being added to the data structure.
     */
    ThreadPool async_merge_l1_ = ThreadPool(1, 100);

    std::function<void(storage_type *)> cleanup_;

    std::atomic<size_t> total_chunk_size_bytes_ = 0;

    bool merged_all_ = false;
  private:
    /** Number of chunks for "level 1" intermediary merging. */
    static constexpr uint32_t MERGE_L1_COUNT = 4;
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_DISK_BASE_HPP__
