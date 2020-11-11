#ifndef __SORTED_SET_DISK_HPP__
#define __SORTED_SET_DISK_HPP__

#include "sorted_set_disk_base.hpp"


namespace mtg {
namespace common {

/**
 * Thread safe data storage that is able to sort and extract distinct elements
 * from the underlying Container using external storage to avoid memory overflow.
 * Data is pushed into the data structure in chunks using the #insert() method.
 * Once the internal buffer is full, the data is sorted and written to disk,
 * each disk write into a different file. The #data() method returns the
 * globally sorted data.
 *
 * @tparam T the type of the elements that are being stored and sorted,
 * typically #KMerBOSS instances
 */
template <typename T>
class SortedSetDisk : public SortedSetDiskBase<T> {
  public:
    typedef T key_type;
    typedef T value_type;
    typedef Vector<T> storage_type;
    typedef ChunkedWaitQueue<T> result_type;

    /**
     * Constructs a SortedSetDisk instance and initializes its buffers to the value
     * specified in #reserved_num_elements.
     * @param num_threads the number of threads to use by the sorting algorithm
     * @param tmp_dir the prefix of the temporary files where chunks are
     * written before being merged
     * @param container_size the size of the in-memory container that is written
     * to disk when full
     */
    SortedSetDisk(size_t num_threads = 1,
                  size_t reserved_num_elements = 1e6,
                  const std::filesystem::path &tmp_dir = "/tmp/",
                  size_t max_disk_space_bytes = 1e9,
                  size_t merge_count = 4)
        : SortedSetDiskBase<T>(num_threads,
                               reserved_num_elements,
                               tmp_dir,
                               max_disk_space_bytes,
                               merge_count) {}

    /**
     * Insert the data between #begin and #end into the buffer. If the buffer is
     * full, the data is sorted, de-duped and written to disk, after which the
     * buffer is cleared.
     */
    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        while (begin != end) {
            Iterator batch_end = this->safe_advance(begin, end, this->buffer_size());
            // acquire the mutex to restrict the number of writing threads
            std::unique_lock<std::mutex> exclusive_lock(this->mutex_);
            size_t offset = this->prepare_insert(begin, batch_end);

            std::shared_lock<std::shared_timed_mutex> multi_insert_lock(
                    this->multi_insert_mutex_);
            // different threads will insert to different chunks of memory, so it's okay
            // (and desirable) to allow concurrent inserts
            exclusive_lock.unlock();
            std::copy(begin, batch_end, this->data_.begin() + offset);
            begin = batch_end;
        }
    }

  private:
    virtual void sort_and_dedupe() override;
};

} // namespace common
} // namespace mtg

#endif // __SORTED_SET_DISK_HPP__
