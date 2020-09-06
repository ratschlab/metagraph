#include "sorted_set_disk_base.hpp"

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/elias_fano.hpp"
#include "common/elias_fano_file_merger.hpp"


namespace mtg {
namespace common {

template <typename T>
SortedSetDiskBase<T>::SortedSetDiskBase(size_t num_threads,
                                        size_t reserved_num_elements,
                                        const std::filesystem::path &tmp_dir,
                                        size_t max_disk_space_bytes)
    : num_threads_(num_threads),
      reserved_num_elements_(reserved_num_elements),
      max_disk_space_bytes_(max_disk_space_bytes),
      chunk_file_prefix_(tmp_dir/"chunk_"),
      merge_queue_(std::min(reserved_num_elements, QUEUE_EL_COUNT)) {
    if (reserved_num_elements == 0) {
        logger->error("SortedSetDisk buffer cannot have size 0");
        std::exit(EXIT_FAILURE);
    }
    try_reserve(reserved_num_elements);
}

template <typename T>
ChunkedWaitQueue<T>& SortedSetDiskBase<T>::data(bool free_buffer) {
    std::unique_lock<std::mutex> exclusive_lock(mutex_);
    std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);

    if (!is_merging_) {
        is_merging_ = true;
        // write any residual data left
        if (!data_.empty()) {
            sort_and_dedupe();
            dump_to_file(true /* is_done */);
        }
        if (free_buffer) {
            Vector<T>().swap(data_); // free up the (usually very large) buffer
        }
        assert(data_.empty());
        start_merging_async();
        chunk_count_ = 0;
        l1_chunk_count_ = 0;
        total_chunk_size_bytes_ = 0;
    }
    return merge_queue_;
}

template <typename T>
std::vector<std::string> SortedSetDiskBase<T>::files_to_merge() {
    std::unique_lock<std::mutex> exclusive_lock(mutex_);
    std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
    // write any residual data left
    if (!data_.empty()) {
        sort_and_dedupe();
        dump_to_file(true /* is_done */);
    }
    return get_file_names();
}

template <typename T>
void SortedSetDiskBase<T>::clear(const std::filesystem::path &tmp_path, bool remove_files) {
    std::unique_lock<std::mutex> exclusive_lock(mutex_);
    std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
    is_merging_ = false;
    if (remove_files) {
        // remove the files that have not been requested to merge
        for (const auto &chunk_file : get_file_names()) {
            std::filesystem::remove(chunk_file);
        }
    }
    chunk_count_ = 0;
    l1_chunk_count_ = 0;
    total_chunk_size_bytes_ = 0;
    try_reserve(reserved_num_elements_);
    data_.resize(0); // this makes sure the buffer is not reallocated
    chunk_file_prefix_ = tmp_path/"chunk_";
    std::filesystem::create_directory(tmp_path);
}

template <typename T>
void SortedSetDiskBase<T>::start_merging_async() {
    // wait until the previous merging job is done
    async_worker_.join();
    // now the queue can be reinitialized and used in the next merge
    merge_queue_.reset();
    const std::vector<std::string> file_names = get_file_names();
    async_worker_.enqueue([file_names, this]() {
        std::function<void(const T &)> on_new_item
                = [this](const T &v) { merge_queue_.push(v); };
        merge_files(file_names, on_new_item, false);
        merge_queue_.shutdown();
    });
}

template <typename T>
void SortedSetDiskBase<T>::shrink_data() {
    logger->trace("Allocated capacity exceeded, erasing duplicate values...");

    size_t old_size = data_.size();
    sort_and_dedupe();

    logger->trace("...done. Size reduced from {} to {}, {} MiB", old_size,
                  data_.size(), (data_.size() * sizeof(T) >> 20));
}

template <typename T>
void SortedSetDiskBase<T>::dump_to_file(bool is_done) {
    assert(!data_.empty());

    std::string file_name = chunk_file_prefix_ + std::to_string(chunk_count_);

    EliasFanoEncoder<T> encoder(data_.size(), utils::get_first(data_.front()),
                                utils::get_first(data_.back()), file_name);
    for (const auto &v : data_) {
        encoder.add(v);
    }
    total_chunk_size_bytes_ += encoder.finish();
    data_.resize(0);
    if (is_done) {
        async_merge_l1_.remove_waiting_tasks();
    } else if (total_chunk_size_bytes_ > max_disk_space_bytes_) {
        async_merge_l1_.remove_waiting_tasks();
        async_merge_l1_.join();
        std::string all_merged_file = merged_all_name(chunk_file_prefix_, merged_all_count_);
        // increment chunk_count, so that get_file_names() returns correct values
        chunk_count_++;
        merge_all(all_merged_file, get_file_names());

        total_chunk_size_bytes_ = std::filesystem::file_size(all_merged_file);
        if (total_chunk_size_bytes_ > max_disk_space_bytes_ * 0.8) {
            logger->critical("Disk space reduced by < 20%. Giving up.");
            std::exit(EXIT_FAILURE);
        }
        merged_all_count_++;
        chunk_count_ = 0;
        l1_chunk_count_ = 0;
        return;
    } else if ((chunk_count_ + 1) % MERGE_L1_COUNT == 0) {
        async_merge_l1_.enqueue(merge_l1, chunk_file_prefix_, chunk_count_,
                                &l1_chunk_count_, &total_chunk_size_bytes_);
    }
    chunk_count_++;
}

template <typename T>
void SortedSetDiskBase<T>::try_reserve(size_t size, size_t min_size) {
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

template <typename T>
void SortedSetDiskBase<T>::merge_l1(const std::string &chunk_file_prefix,
                                    uint32_t chunk_count,
                                    std::atomic<uint32_t> *l1_chunk_count,
                                    std::atomic<size_t> *total_size) {
    const std::string &merged_l1_file_name
            = SortedSetDiskBase<T>::merged_l1_name(chunk_file_prefix,
                                                   chunk_count / MERGE_L1_COUNT);
    std::vector<std::string> to_merge(MERGE_L1_COUNT);
    for (uint32_t i = 0; i < MERGE_L1_COUNT; ++i) {
        to_merge[i] = chunk_file_prefix + std::to_string(chunk_count - i);
        *total_size -= static_cast<int64_t>(std::filesystem::file_size(to_merge[i]));
    }
    logger->trace("Starting merging last {} chunks into {}", MERGE_L1_COUNT,
                  merged_l1_file_name);
    EliasFanoEncoderBuffered<T> encoder(merged_l1_file_name, 1000);
    std::function<void(const T &v)> on_new_item
            = [&encoder](const T &v) { encoder.add(v); };
    merge_files(to_merge, on_new_item);
    encoder.finish();

    (*l1_chunk_count)++;
    logger->trace("Merging last {} chunks into {} done", MERGE_L1_COUNT,
                  merged_l1_file_name);
    *total_size += std::filesystem::file_size(merged_l1_file_name);
}

template <typename T>
void SortedSetDiskBase<T>::merge_all(const std::string &out_file,
                                     const std::vector<std::string> &to_merge) {
    logger->trace(
            "Max allocated disk capacity exceeded. Starting merging all {} chunks "
            "into {}",
            to_merge.size(), out_file);
    EliasFanoEncoderBuffered<T> encoder(out_file, 1000);
    std::function<void(const T &v)> on_new_item
            = [&encoder](const T &v) { encoder.add(v); };
    merge_files(to_merge, on_new_item);
    encoder.finish();
    logger->trace("Merging all {} chunks into {} of size {:.0f}MiB done",
                  to_merge.size(), out_file, std::filesystem::file_size(out_file) / 1e6);
}

template <typename T>
std::vector<std::string> SortedSetDiskBase<T>::get_file_names() {
    async_merge_l1_.join(); // make sure all L1 merges are done
    std::vector<std::string> file_names;
    if (merged_all_count_ > 0) {
        file_names.push_back(merged_all_name(chunk_file_prefix_, merged_all_count_ - 1));
    }

    for (size_t i = 0; i < l1_chunk_count_; ++i) {
        file_names.push_back(merged_l1_name(chunk_file_prefix_, i));
    }
    for (size_t i = MERGE_L1_COUNT * l1_chunk_count_; i < chunk_count_; ++i) {
        file_names.push_back(chunk_file_prefix_ + std::to_string(i));
    }
    return file_names;
}

template class SortedSetDiskBase<uint64_t>;
template class SortedSetDiskBase<sdsl::uint128_t>;
template class SortedSetDiskBase<sdsl::uint256_t>;
template class SortedSetDiskBase<std::pair<uint64_t, uint8_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint128_t, uint8_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint256_t, uint8_t>>;
template class SortedSetDiskBase<std::pair<uint64_t, uint16_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint128_t, uint16_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint256_t, uint16_t>>;
template class SortedSetDiskBase<std::pair<uint64_t, uint32_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint128_t, uint32_t>>;
template class SortedSetDiskBase<std::pair<sdsl::uint256_t, uint32_t>>;

} // namespace common
} // namespace mtg
