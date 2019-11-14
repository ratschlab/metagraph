#pragma once

#include "common/chunked_wait_queue.hpp"
#include "common/deque_vector.hpp"
#include "common/threading.hpp"
#include "common/vectors.hpp"

#include <ips4o.hpp>

#include <cassert>
#include <fstream> //TODO(ddanciu) - try boost mmapped instead
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <shared_mutex>

const std::string output_dir = "/tmp/"; // TODO(ddanciu) - use a flag instead

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
    /**
     * Size of the underlying data structure storing the kmers. The class is
     * guaranteed to flush to disk when the newly added data would exceed this
     * size.
     */
    static constexpr size_t CONTAINER_SIZE_BYTES = 1e3; // 1 GB
    //TODO:set back to 1e9 before submitting, also below

    /** Size of the merge queue's underlying circular buffer */
    static constexpr size_t MERGE_QUEUE_SIZE = 1e3; // 1GB
    /** Number of elements that can be iterated backwards in the merge queue */
    static constexpr size_t MERGE_QUEUE_BACKWARDS_COUNT = 100;

    typedef T key_type;
    typedef T value_type;
    typedef Vector<T> storage_type;

    /**
     * Constructs a SortedSetDisk instance and initializes its buffer to the value
     * specified in CONTAINER_SIZE_BYTES.
     * @param num_threads the number of threads to use by the sorting algorithm
     * @param verbose true if verbose logging is desired
     * @param container_size the size of the in-memory container that is written
     * to disk when full
     */
    SortedSetDisk(
            std::function<void(storage_type *)> cleanup = [](storage_type *) {},
            size_t num_threads = 1,
            bool verbose = false,
            size_t container_size = CONTAINER_SIZE_BYTES)
        : num_threads_(num_threads),
          verbose_(verbose),
          cleanup_(cleanup),
          merge_queue_(MERGE_QUEUE_SIZE, MERGE_QUEUE_SIZE / 10, MERGE_QUEUE_BACKWARDS_COUNT) {
        try {
            try_reserve(container_size);
        } catch (const std::bad_alloc &exception) {
            std::cerr << "ERROR: Not enough memory for SortedSetDisk. Requested"
                      << CONTAINER_SIZE_BYTES << " bytes" << std::endl;
            exit(1);
        }
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
        if (data_->size() + batch_size > data_->capacity()) { // time to write to disk
            std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
            // TODO(ddanciu) - test if it's worth it to keep adding in memory if the
            // shrinking was significant
            shrink_data();

            dump_to_file_async();
        }

        size_t offset = data_->size();
        // resize to the future size after insertion (no reallocation will happen)
        data_->resize(data_->size() + batch_size);
        std::shared_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        exclusive_lock.unlock();
        // different threads will insert to different chunks of memory, so it's okay
        // (and desirable) to allow concurrent inserts
        std::copy(begin, end, data_->begin() + offset);
    }

    storage_type &data() {
        // TODO(ddanciu) - implement an adaptor from ChunkedWaitQueue to the expected
        //  data structures in KmerStorage
        throw std::runtime_error("Function not yet implemented");
    }

    /**
     * Returns the globally sorted and de-duped data. Typically called once all
     * the data was inserted via insert().
     */
    threads::ChunkedWaitQueue<T> &dataStream() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);

        if (!is_merging_) {
            is_merging_ = true;
            // write any residual data left
            if (!data_->empty()) {
                sort_and_remove_duplicates(data_, num_threads_);
                dump_to_file_async();
            }
            if (write_to_disk_future_.valid()) {
                write_to_disk_future_.get(); // make sure all pending data was written
            }

            start_merging();
        }
        return merge_queue_;
    }

    static void merge_data(uint32_t chunk_count, threads::ChunkedWaitQueue<T> *merge_queue) {
        // start merging disk chunks by using a heap to store the current element
        // from each chunk
        std::vector<std::fstream> chunk_files(chunk_count);
        std::fstream sorted_file(output_dir + "sorted.bin", std::ios::binary | std::ios::out);

        auto comp = [](const auto &a, const auto &b) { return a.first > b.first; };

        std::priority_queue<std::pair<T, uint32_t>, std::vector<std::pair<T, uint32_t>>, decltype(comp)>
                merge_heap(comp);

        for (uint32_t i = 0; i < chunk_count; ++i) {
            const std::string file_name
                    = output_dir + "chunk_" + std::to_string(i) + ".bin";
            chunk_files[i].open(file_name, std::ios::in | std::ios::binary);
            T data_item;
            if (chunk_files[i].good()) {
                chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(data_item));
                merge_heap.push({ data_item, i });
            } else {
                throw std::runtime_error("Unable to open chunk file " + file_name);
            }
        }
        uint64_t totalSize = 0;
        // init with any value that is not the top
        T last_written;
        bool has_written = false;
        while (!merge_heap.empty()) {
            std::pair<T, uint32_t> smallest = merge_heap.top();
            merge_heap.pop();
            if (!has_written || smallest.first != last_written) {
                has_written = true;
                merge_queue->push_front(smallest.first);

                last_written = smallest.first;
                totalSize++;
            }
            if (chunk_files[smallest.second].good()) {
                T data_item;
                if (chunk_files[smallest.second].read(reinterpret_cast<char *>(&data_item),
                                                      sizeof(data_item))) {
                    // TODO(ddanciu): apply logic from cleanup_ to remove redundant
                    // dummy BOSS kmers
                    merge_heap.push({ data_item, smallest.second });
                }
            }
        }
        merge_queue->shutdown();
    }

    std::future<void> start_merging() {
        return thread_pool_.enqueue(merge_data, chunk_count_, &merge_queue_);
    }

    void clear() {
        std::unique_lock<std::mutex> exclusive_lock(mutex_);
        std::unique_lock<std::shared_timed_mutex> multi_insert_lock(multi_insert_mutex_);
        data_->resize(0); // this makes sure the buffer is not reallocated
    }

    template <class Array>
    void sort_and_remove_duplicates(Array *vector, size_t num_threads) const {
        assert(vector);

        ips4o::parallel::sort(vector->begin(), vector->end(),
                              std::less<typename Array::value_type>(), num_threads);
        // remove duplicates
        auto unique_end = std::unique(vector->begin(), vector->end());
        vector->erase(unique_end, vector->end());

        cleanup_(vector);
    }

    void reserve(size_t size) {
        std::cerr << "SortedSetDisk: Ignoring reserving size " << size << std::endl;
    }

  private:
    void shrink_data() {
        if (verbose_) {
            std::cout << "Allocated capacity exceeded, erasing duplicate values..."
                      << std::flush;
        }

        size_t old_size = data_->size();
        sort_and_remove_duplicates(data_, num_threads_);

        if (verbose_) {
            std::cout << " done. Size reduced from " << old_size << " to " << data_->size()
                      << ", " << (data_->size() * sizeof(T) >> 20) << "Mb" << std::endl;
        }
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
    static void dump_to_file_sync(uint32_t chunk_count, storage_type *data) {
        std::fstream binary_file
                = std::fstream(output_dir + "chunk_" + std::to_string(chunk_count)
                                       + ".bin",
                               std::ios::out | std::ios::binary);
        if (!binary_file.write((char *)&((*data)[0]), sizeof((*data)[0]) * data->size())) {
            throw std::runtime_error("Writing of chunk no " + std::to_string(chunk_count)
                                     + " failed.");
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
        if (write_to_disk_future_.valid()) {
            write_to_disk_future_.wait(); // wait for other thread to finish writing
        }
        if (data_->data() == data_first.data()) {
            write_to_disk_future_ = thread_pool_.enqueue(dump_to_file_sync<storage_type>,
                                                         chunk_count_, &data_first);
            data_ = &data_second;
        } else {
            write_to_disk_future_ = thread_pool_.enqueue(dump_to_file_sync<storage_type>,
                                                         chunk_count_, &data_second);
            data_ = &data_first;
        }
        chunk_count_++;
    }

    void try_reserve(size_t size, size_t min_size = 0) {
        if constexpr (std::is_same_v<DequeStorage<T>, storage_type>) {
            data_first.try_reserve(size, min_size);
            data_second.try_reserve(size, min_size);

        } else {
            size = std::max(size, min_size);

            while (size > min_size) {
                try {
                    data_first.reserve(size);
                    data_second.reserve(size);
                    return;
                } catch (const std::bad_alloc &exception) {
                    size = min_size + (size - min_size) * 2 / 3;
                }
            }
            data_first.reserve(min_size);
            data_second.reserve(min_size);
        }
    }

    /**
     * Current number of chunks written to disk. We expect this to be at most in
     * the order of thousands, so a 32 bit integer should suffice for storage.
     */
    uint32_t chunk_count_ = 0;
    /**
     * Alternating buffers into which data is inserted using #insert(). While
     * one buffer is being written to disk using #dump_to_file() the other one
     * is inserted into using #insert(). Once the insert buffer is full, we
     * wait for the disk write operation to finish (if needed) and then swap
     * buffers.
     */
    storage_type data_first, data_second;
    /**
     * Reference to the buffer data is currently written into (either #data1_
     * or #data2_)
     */
    storage_type *data_ = &data_first;
    size_t num_threads_;
    /**
     * True if the data merging thread was started, and data started flowing into the #merge_queue_.
     */
    bool is_merging_ = false;
    bool verbose_;

    /**
     * Removes redundant elements from container while extracting unique ones.
     */
    std::function<void(storage_type *)> cleanup_;

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
     * Thread pool with two threads: one for writing to disk and one for merging data from disk.
     * Each thread has a single task (write to disk and merge from disk, respectively).
     */
    ThreadPool thread_pool_ = ThreadPool(2, 1);
    /**
     * Future that signals whether the current writing to disk job was done.
     */
    std::future<void> write_to_disk_future_;

    threads::ChunkedWaitQueue<T> merge_queue_;
};
