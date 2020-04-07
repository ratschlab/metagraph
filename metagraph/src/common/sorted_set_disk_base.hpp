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
 * @tparam INT the corresponding integer representation of T, used for compressed storage
 * on disk
 */
template <typename T, typename INT>
class SortedSetDiskBase {
    typedef T value_type;
    typedef INT int_type;
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
                      size_t num_last_elements_cached,
                      std::function<int_type(const value_type &v)> to_int)
        : num_threads_(num_threads),
          reserved_num_elements_(reserved_num_elements),
          max_disk_space_bytes_(max_disk_space_bytes),
          on_item_pushed_(on_item_pushed),
          chunk_file_prefix_(tmp_dir / "chunk_"),
          merge_queue_(std::min(reserved_num_elements, QUEUE_EL_COUNT),
                       num_last_elements_cached,
                       on_item_pushed),
          cleanup_(cleanup),
          to_int_(to_int) {
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

    /**
     * Clears the set, preparing it to be re-used for another merge. Creating a new
     * sorted set may be expensive when #data_ is large. In these cases, prefere calling
     * #clear and re-using the sorted.
     */
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

    void start_merging() {
        const std::vector<std::string> file_names = this->get_file_names();
        this->async_worker_.enqueue([file_names, this]() {
          std::function<void(const value_type &)> on_new_item
                  = [this](const value_type &v) { this->merge_queue_.push(v); };
          if constexpr (utils::is_pair<T> {}) {
              merge_files<typename T::first_type, typename T::second_type,
                          typename INT::first_type>(file_names, on_new_item);
          } else {
              merge_files<T, INT>(file_names, on_new_item);
          }
          this->merge_queue_.shutdown();
        });
    }

    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erasing duplicate values...");

        size_t old_size = data_.size();
        sort_and_remove_duplicates(&data_, num_threads_);

        logger->trace("...done. Size reduced from {} to {}, {}MiB", old_size,
                      data_.size(), (data_.size() * sizeof(T) >> 20));
    }

    /**
     * Dumps the given data to a file, synchronously. If the maximum allowed disk size
     * is reached, all chunks will be merged into a single chunk in an effort to reduce
     * disk space.
     * @param is_done if this is the last chunk being dumped
     */
    void dump_to_file(bool is_done) {
        assert(!this->data_.empty());

        std::string file_name
                = this->chunk_file_prefix_ + std::to_string(this->chunk_count_);

        EliasFanoEncoder<INT> encoder(this->data_.size(),
                                      utils::get_first(to_int_(this->data_.front())),
                                      utils::get_first(to_int_(this->data_.back())),
                                      file_name);
        for (const auto &v : this->data_) {
            encoder.add(to_int_(v));
        }
        this->total_chunk_size_bytes_ += encoder.finish();
        this->data_.resize(0);
        if (is_done) {
            this->async_merge_l1_.clear();
        } else if (this->total_chunk_size_bytes_ > this->max_disk_space_bytes_) {
            this->async_merge_l1_.clear();
            this->async_merge_l1_.join();
            std::string all_merged_file = this->chunk_file_prefix_ + "_all.tmp";
            this->chunk_count_++; // needs to be incremented, so that get_file_names()
            // returns correct values
            this->merge_all(all_merged_file, this->get_file_names(), to_int_);
            this->merged_all_ = true;
            this->total_chunk_size_bytes_ = std::filesystem::file_size(all_merged_file);
            if (this->total_chunk_size_bytes_ > this->max_disk_space_bytes_ * 0.8) {
                logger->critical("Disk space reduced by < 20%. Giving up.");
                std::exit(EXIT_FAILURE);
            }
            std::filesystem::rename(all_merged_file,
                                    this->merged_all_name(this->chunk_file_prefix_));
            this->chunk_count_ = 0;
            this->l1_chunk_count_ = 0;
            return;
        } else if ((this->chunk_count_ + 1) % MERGE_L1_COUNT == 0) {
            this->async_merge_l1_.enqueue(this->merge_l1, this->chunk_file_prefix_,
                                          this->chunk_count_, &this->l1_chunk_count_,
                                          &this->total_chunk_size_bytes_, to_int_);
        }
        this->chunk_count_++;
    }

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
    std::atomic<uint32_t> l1_chunk_count_ = 0;

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

    std::function<int_type(const value_type &v)> to_int_;

  private:
    /** Number of chunks for "level 1" intermediary merging. */
    static constexpr uint32_t MERGE_L1_COUNT = 4;

    static inline std::string merged_l1_name(const std::string &prefix, uint32_t count) {
        return prefix + "m" + std::to_string(count);
    }

    static inline std::string merged_all_name(const std::string &prefix) {
        return prefix + "_all";
    }

    static void merge_l1(const std::string &chunk_file_prefix,
                         uint32_t chunk_count,
                         std::atomic<uint32_t> *l1_chunk_count,
                         std::atomic<size_t> *total_size,
                         std::function<int_type(const value_type &v)> to_int) {
        const std::string &merged_l1_file_name
                = SortedSetDiskBase<value_type, INT>::merged_l1_name(
                        chunk_file_prefix, chunk_count / MERGE_L1_COUNT);
        std::vector<std::string> to_merge(MERGE_L1_COUNT);
        for (uint32_t i = 0; i < MERGE_L1_COUNT; ++i) {
            to_merge[i] = chunk_file_prefix + std::to_string(chunk_count - i);
            *total_size -= static_cast<int64_t>(std::filesystem::file_size(to_merge[i]));
        }
        logger->trace("Starting merging last {} chunks into {}", MERGE_L1_COUNT,
                      merged_l1_file_name);
        EliasFanoEncoderBuffered<int_type> encoder(merged_l1_file_name, 1000);
        std::function<void(const value_type &v)> on_new_item
                = [&encoder, to_int](const value_type &v) { encoder.add(to_int(v)); };
        if constexpr (utils::is_pair<T> {}) {
            merge_files<typename T::first_type, typename T::second_type, typename INT::first_type>(
                    to_merge, on_new_item);
        } else {
            merge_files<T, INT>(to_merge, on_new_item);
        }
        encoder.finish();

        (*l1_chunk_count)++;
        logger->trace("Merging last {} chunks into {} done", MERGE_L1_COUNT,
                      merged_l1_file_name);
        *total_size += std::filesystem::file_size(merged_l1_file_name);
    }

    static void merge_all(const std::string &out_file,
                          const std::vector<std::string> &to_merge,
                          std::function<int_type(const value_type &v)> to_int) {
        std::ofstream merged_count_file(out_file + ".count",
                                        std::ios::binary | std::ios::out);
        logger->trace(
                "Max allocated disk capacity exceeded. Starting merging all {} chunks "
                "into {}",
                to_merge.size(), out_file);
        EliasFanoEncoderBuffered<int_type> encoder(out_file, 1000);
        std::function<void(const value_type &v)> on_new_item
                = [&encoder, to_int](const value_type &v) { encoder.add(to_int(v)); };
        if constexpr (utils::is_pair<T> {}) {
            merge_files<typename T::first_type, typename T::second_type, typename INT::first_type>(
                    to_merge, on_new_item);
        } else {
            merge_files<T, INT>(to_merge, on_new_item);
        }
        encoder.finish();
        logger->trace("Merging all {} chunks into {} of size {:.0f}MiB done",
                      to_merge.size(), out_file, std::filesystem::file_size(out_file) / 1e6);
    }

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
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_DISK_BASE_HPP__
