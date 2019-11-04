#pragma once

#include <mutex>
#include <shared_mutex>

#include <ips4o.hpp>

#include "utils.hpp"

const std::string output_dir = "/tmp/"; // TODO - use a flag instead

// Thread safe data storage to extract distinct elements
template <typename T, class Container = Vector<T>,
          class Cleaner = utils::NoCleanup>
class SortedSetExternal {
public:
  /**
   * Size of the underlying data structure storing the kmers. The class is
   * guaranteed to flush to disk when the newly added data would exceed this
   * size.
   */
  static constexpr size_t CONTAINER_SIZE_BYTES = 1e9; // 1 GB
  static_assert(std::is_same_v<T, typename Container::value_type>);

  typedef T key_type;
  typedef T value_type;
  typedef Container storage_type;

  SortedSetExternal(size_t num_threads = 1, bool verbose = false)
      : container_size_(CONTAINER_SIZE_BYTES), num_threads_(num_threads),
        verbose_(verbose) {
    try {
      try_reserve(CONTAINER_SIZE_BYTES);
    } catch (const std::bad_alloc &exception) {
      std::cerr << "ERROR: Not enough memory for SortedSetExternal. Requested"
                << CONTAINER_SIZE_BYTES << " bytes" << std::endl;
      exit(1);
    }
  }

  template <class Iterator> void insert(Iterator begin, Iterator end) {
    assert(begin <= end);

    uint64_t batch_size = end - begin;

    // acquire the mutex to restrict the number of writing threads
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);

    if (data_.size() + batch_size > data_.capacity()) { // time to write to disk
      shrink_data();

      std::fstream binary_file = std::fstream(
          output_dir + "file" + std::to_string(chunk_count_) + ".binary",
          std::ios::out | std::ios::binary);
      binary_file.write((char *)&data_[0], data_.size());
      binary_file.close();
      data_.clear();
      chunk_count_++;
    }

    std::copy(begin, end, data_.end());
  }

  void reserve(size_t size) {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    try_reserve(size);
  }

  storage_type &data() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    if (sorted_end_ != data_.size()) {
      sort_and_remove_duplicates(&data_, num_threads_);
      sorted_end_ = data_.size();
    }

    return data_;
  }

  void clear() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    data_ = storage_type();
    sorted_end_ = 0;
  }

  template <class Array>
  void sort_and_remove_duplicates(Array *vector, size_t num_threads) const {
    assert(vector);

    ips4o::parallel::sort(vector->begin(), vector->end(),
                          std::less<typename Array::value_type>(), num_threads);
    // remove duplicates
    auto unique_end = std::unique(vector->begin(), vector->end());
    vector->erase(unique_end, vector->end());

    Cleaner::cleanup(vector);
  }

private:
  void shrink_data() {
    if (verbose_) {
      std::cout << "Allocated capacity exceeded, erase duplicate values..."
                << std::flush;
    }

    size_t old_size = data_.size();
    sort_and_remove_duplicates(&data_, num_threads_);
    sorted_end_ = data_.size();

    if (verbose_) {
      std::cout << " done. Size reduced from " << old_size << " to "
                << data_.size() << ", " << (data_.size() * sizeof(T) >> 20)
                << "Mb" << std::endl;
    }
  }

  void try_reserve(size_t size, size_t min_size = 0) {
    if constexpr (std::is_same_v<utils::DequeStorage<T>, storage_type>) {
      data_.try_reserve(size, min_size);

    } else {
      size = std::max(size, min_size);

      while (size > min_size) {
        try {
          data_.reserve(size);
          return;
        } catch (const std::bad_alloc &exception) {
          size = min_size + (size - min_size) * 2 / 3;
        }
      }
      data_.reserve(min_size);
    }
  }

  uint64_t chunk_count_ = 0;
  size_t container_size_;
  storage_type data_;
  size_t num_threads_;
  bool verbose_;

  // indicate the end of the preprocessed distinct and sorted values
  uint64_t sorted_end_ = 0;

  mutable std::mutex mutex_resize_;
  mutable std::shared_timed_mutex mutex_copy_;
};
