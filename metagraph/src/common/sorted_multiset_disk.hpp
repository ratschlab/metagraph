#pragma once

#include <cassert>
#include <functional>
#include <optional>
#include <shared_mutex>
#include <string>

#include <ips4o.hpp>

#include "common/sorted_set_disk_base.hpp"
#include "common/threads/chunked_wait_queue.hpp"


namespace mg {
namespace common {

/**
 * Specialization of SortedSetDiskBase that is able to both sort and count elements.
 *
 * @tparam T the type of the elements that are being stored, sorted and counted,
 * typically k-mers
 * @param C the type used to count the multiplicity of each value in the multi-set
 */
template <typename T, typename INT = T, typename C = uint8_t>
class SortedMultisetDisk : public SortedSetDiskBase<std::pair<T, C>, INT> {
  public:
    typedef T key_type;
    typedef C count_type;
    typedef std::pair<T, C> value_type;
    typedef std::pair<INT, C> int_pair;
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
            std::function<void(const value_type &)> on_item_pushed
            = [](const value_type &) {},
            size_t num_last_elements_cached = 100,
            std::function<int_pair(const value_type &v)> to_int
            = [](const value_type &v) { return std::make_pair(INT(v.first), v.second); })
        : SortedSetDiskBase<value_type, INT>(cleanup,
                                             num_threads,
                                             reserved_num_elements,
                                             tmp_dir,
                                             max_disk_space_bytes,
                                             on_item_pushed,
                                             num_last_elements_cached),
          to_int_(to_int) {}

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
    virtual void start_merging() override {
        const std::vector<std::string> file_names = this->get_file_names();
        this->async_worker_.enqueue([file_names, this]() {
            std::function<void(const value_type &)> on_new_item
                    = [this](const value_type &v) { this->merge_queue_.push(v); };
            merge_files<T, C, INT>(file_names, on_new_item);
            this->merge_queue_.shutdown();
        });
    }

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

    /**
     * Dumps the given data to a file, synchronously. If the maximum allowed disk size
     * is reached, all chunks will be merged into a single chunk in an effort to reduce
     * disk space.
     * @param is_done if this is the last chunk being dumped
     */
    virtual void dump_to_file(bool is_done) override {
        assert(!this->data_.empty());

        std::string file_name
                = this->chunk_file_prefix_ + std::to_string(this->chunk_count_);

        EliasFanoEncoder<int_pair> encoder(this->data_.size(), to_int_(this->data_.front()).first,
                                           to_int_(this->data_.back()).first, file_name);
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
            this->async_merge_l1_.enqueue(merge_l1, this->chunk_file_prefix_,
                                          this->chunk_count_, &this->l1_chunk_count_,
                                          &this->total_chunk_size_bytes_, to_int_);
        }
        this->chunk_count_++;
    }

    static void merge_l1(const std::string &chunk_file_prefix,
                         uint32_t chunk_count,
                         std::atomic<uint32_t> *l1_chunk_count,
                         std::atomic<size_t> *total_size,
                         std::function<int_pair(const value_type &v)> to_int) {
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
        EliasFanoEncoderBuffered<int_pair> encoder(merged_l1_file_name, 1000);
        merge_files<T, C, INT>(to_merge, [&encoder, to_int](const std::pair<T, C> &v) {
            encoder.add(to_int(v));
        });
        encoder.finish();

        (*l1_chunk_count)++;
        logger->trace("Merging last {} chunks into {} done", MERGE_L1_COUNT,
                      merged_l1_file_name);
        *total_size += std::filesystem::file_size(merged_l1_file_name);
    }

    static void merge_all(const std::string &out_file,
                          const std::vector<std::string> &to_merge,
                          std::function<int_pair(const value_type &v)> to_int) {
        std::ofstream merged_count_file(out_file + ".count",
                                        std::ios::binary | std::ios::out);
        logger->trace(
                "Max allocated disk capacity exceeded. Starting merging all {} chunks "
                "into {}",
                to_merge.size(), out_file);
        EliasFanoEncoderBuffered<int_pair> encoder(out_file, 1000);
        merge_files<T, C, INT>(to_merge, [&encoder, to_int](const std::pair<T, C> &v) {
            encoder.add(to_int(v));
        });
        encoder.finish();
        logger->trace("Merging all {} chunks into {} of size {:.0f}MiB done",
                      to_merge.size(), out_file, std::filesystem::file_size(out_file) / 1e6);
    }

  private:
    /** Number of chunks for "level 1" intermediary merging. */
    static constexpr uint32_t MERGE_L1_COUNT = 4000;

    std::function<int_pair(const value_type &v)> to_int_;
};

} // namespace common
} // namespace mg
