#ifndef __SORTED_SET_DISK_BASE_HPP__
#define __SORTED_SET_DISK_BASE_HPP__

#include <atomic>
#include <cassert>
#include <filesystem>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <vector>

#include "common/threads/chunked_wait_queue.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace common {

/**
 * Abstract thread safe data storage that is able to sort  elements from the
 * underlying Container using external storage to limit the amount of required memory.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is processed using
 * #sort_and_dedupe and written to disk, each disk write into a different file.
 * The #data() method returns the globally sorted data.
 *
 * @tparam T the type of the elements that are being stored and sorted,
 * typically #KMerBOSS instances, or <#KmerBOSS, count> pairs.
 * @tparam INT the corresponding integer representation of T, used for compressed storage
 * on disk
 */
template <typename T>
class SortedSetDiskBase {
    typedef T value_type;
    typedef Vector<T> storage_type;
    typedef ChunkedWaitQueue<T> result_type;
    typedef typename storage_type::iterator Iterator;

    /**
     * The number of elements in the merge queue. Should be large enough to reduce lock
     * contention, but small enough to not lose time with memory allocation.
     */
    static constexpr size_t QUEUE_EL_COUNT = 30'000;

  public:
    SortedSetDiskBase(size_t num_threads,
                      size_t reserved_num_elements,
                      const std::filesystem::path &tmp_dir,
                      size_t max_disk_space_bytes);

    virtual ~SortedSetDiskBase() {
        // remove the files that have not been requested to merge
        for (const auto &chunk_file : get_file_names()) {
            std::filesystem::remove(chunk_file);
        }
        async_worker_.join(); // make sure the data was processed
    }

    size_t buffer_size() const { return data_.capacity(); }

    /**
     * Returns the globally sorted and counted data. Typically called once all the data
     * was inserted via insert().
     */
    ChunkedWaitQueue<T>& data(bool free_buffer = true);

    /**
     * Returns the files to be merged - useful if the caller prefers to do the merging.
     */
    std::vector<std::string> files_to_merge();

    /**
     * Clears the set, preparing it to be re-used for another merge. Creating a new
     * sorted set may be expensive when #data_ is large. In these cases, prefer calling
     * #clear and re-using the buffer.
     */
    void clear(const std::filesystem::path &tmp_path = "/tmp/");

    /**
     * Insert already sorted data into the set. This data is written directly to a
     * designated chunk without being sorted and de-duped.
     * Use this method along with #insert when some of the data that is inserted into the
     * set is known to be sorted.
     * Note: if calling #insert_sorted multiple times, #data must be globally sorted
     */
    void insert_sorted(const std::vector<T> &data);

  protected:
    /** Advances #it by step or points to #end, whichever comes first. */
    template <typename Iterator>
    Iterator safe_advance(const Iterator &it, const Iterator &end, size_t step) {
        assert(it <= end);
        return static_cast<size_t>(end - it) < buffer_size() ? end : it + step;
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
     * Hold the data filled in via #insert.
     */
    storage_type data_;

    /**
     * Ensures mutually exclusive access (and thus thread-safety) to #data.
     */
    mutable std::mutex mutex_;
    /**
     * Mutex that can be acquired by multiple threads that are appending to
     * non-overlapping areas of #data_
     */
    mutable std::shared_timed_mutex multi_insert_mutex_;

    size_t num_threads_;

  private:
    virtual void sort_and_dedupe() = 0;

    void start_merging_async();

    void shrink_data();

    /**
     * Dumps the given data to a file, synchronously. If the maximum allowed disk size
     * is reached, all chunks will be merged into a single chunk in an effort to reduce
     * disk space.
     * @param is_done if this is the last chunk being dumped
     */
    void dump_to_file(bool is_done);

    void try_reserve(size_t size, size_t min_size = 0);

    /**
     * Current number of chunks written to disk. We expect this to be at most in
     * the order of thousands, so a 32 bit integer should suffice for storage.
     */
    uint32_t chunk_count_ = 0;

    /**
     * The number of L1 merges that were successfully performed.
     */
    std::atomic<uint32_t> l1_chunk_count_ = 0;

    size_t reserved_num_elements_;

    size_t max_disk_space_bytes_;

    std::string chunk_file_prefix_;

    /**
     * True if the data merging thread was started, and data started flowing into the #merge_queue_.
     */
    bool is_merging_ = false;

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

    std::atomic<size_t> total_chunk_size_bytes_ = 0;

    uint32_t merged_all_count_ = 0;

    /** Number of chunks for "level 1" intermediary merging. */
    static constexpr uint32_t MERGE_L1_COUNT = 4;

    static std::string merged_l1_name(const std::string &prefix, uint32_t count) {
        return prefix + "m" + std::to_string(count);
    }

    static std::string merged_all_name(const std::string &prefix, uint32_t count) {
        return prefix + "all_" + std::to_string(count);
    }

    static void merge_l1(const std::string &chunk_file_prefix,
                         uint32_t chunk_count,
                         std::atomic<uint32_t> *l1_chunk_count,
                         std::atomic<size_t> *total_size);

    static void merge_all(const std::string &out_file,
                          const std::vector<std::string> &to_merge);

    std::vector<std::string> get_file_names();
};

} // namespace common
} // namespace mtg

#endif // __SORTED_SET_DISK_BASE_HPP__
