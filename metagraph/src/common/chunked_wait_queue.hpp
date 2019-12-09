#pragma once

#include <cassert>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <optional>
#include <stack>
#include <thread>
#include <vector>

#include "common/threading.hpp"

namespace mg {
namespace common {

/**
 *  A ChunkedWaitQueue is a data structure that allows safe zero-copy data transfer from
 * one thread to another. There are a few fundamental differences between a normal
 * WaitQueue and this implementation:
 *  1. Users have no direct control over how elements are popped from the queue.
 *  Elements are popped from the queue in chunks whenever the single iterator of the queue
 * is far enough from the oldest element.
 *  2. The queue exposes a single iterator instance.  The Iterator instance allows
 * iterating the elements of the queue forward. The iterator will remove elements from the
 * queue as they become old. Limited iteration backwards is permitted with a maximum of
 * #fence_size elements from the most recently read element.
 *  4. Since the queue exposes a single iterator, it expects a single reader thread.
 * Multiple writers may push information into the queue.
 *
 *  The ChunkedWaitQueue uses a pre-allocated circular array to store the elements in
 * order to avoid resize and heap allocations at runtime.
 *
 * Writers can push_front a value to the queue, the reader can access elements in order
 * via a single iterator exposed by the class. The reader/iterator will block if there
 * are no available elements, the writers will block if the queue is full i.e. if it
 * contains buffer_size number of elements.
 *
 *  When the queue is shutdown, the reader will unblock. At destruction time, the queue is
 * shut down (if not already shut down) and the destructor will wait until all elements
 * were read (i.e. the iterator reaches past the last available element).
 *
 * The class will optionally write all enqueued elements to a file specified via
 * #set_out_file. The #push_front operation will block until an output file is specified.
 *
 * The class is not copyable or copy constructible.
 */
template <typename T, typename Alloc = std::allocator<T>>
class ChunkedWaitQueue {
  public:
    class Iterator;
    // typedefs for STL compatibility
    typedef size_t size_type;
    typedef T value_type;
    typedef Iterator iterator;

    ChunkedWaitQueue(const ChunkedWaitQueue &other) = delete;
    ChunkedWaitQueue &operator=(const ChunkedWaitQueue &) = delete;

    /**
     * Constructs a WaitQueue with the given size parameters.
     * @param buffer_size the size of the buffer (number of elements) used internally by
     * the queue. The actual memory footprint is buffer_size*sizeof(T)
     * @param fence_size the number of elements that can be iterated backwards. The queue
     * will always keep at least fence_size elements behind the furthest iterator for this
     * purpose. Must be smaller than #buffer_size, and in practice it's orders of
     * magnitude smaller.
     */
    explicit ChunkedWaitQueue(size_type buffer_size, size_type fence_size)
        : chunk_size_(std::min(std::max(1UL, buffer_size / 3), buffer_size - fence_size)),
          fence_size_(fence_size),
          buffer_size_(buffer_size),
          queue_(buffer_size),
          end_iterator_(Iterator(this, buffer_size)),
          is_shutdown_(false) {
        assert(fence_size < buffer_size);
        output_file_.status = OutputFile::Status::NotSet;
    }

    /**
     * Destroys the wait queue. Notifies all readers and waits until all data was read.
     */
    ~ChunkedWaitQueue() {
        shutdown();
        std::unique_lock<std::mutex> lock(mutex_);
        empty_.wait(lock, [this] { return empty() || iterator_ == end(); });
    }

    /**
     * Resets the queue to an empty state.
     * Undefined behavior if the queue is reset while iterating over it.
     */
    void reset() {
        shutdown();
        first_ = 0;
        last_ = buffer_size_;
        is_shutdown_ = false;
        iterator_.idx_ = 0;
        output_file_.status = OutputFile::Status::NotSet;
    }

    /**
     * Returns true if the queue is empty. This can happen only before an element is
     * added. Once an element is added, the queue will never be empty again (because
     * the queue will always keep at least #fence_size_ elements for backwards iteration).
     */
    bool empty() const { return last_ == buffer_size_; }

    /**
     * Returns true if the buffer of the queue is full. #push_front operations will block
     * until #iterator() advances far enough in the queue to allow garbage collecting
     * the older elements via #pop_chunk().
     */
    bool full() const {
        return last_ != buffer_size_ && first_ == (last_ + 1) % queue_.size();
    }

    /**
     * Returns the capacity of the queue's internal buffer.
     */
    size_type capacity() const { return queue_.capacity(); }

    /**
     * Returns the number of elements in the buffer. Sort of irrelevant, as the number of
     * elements is not under user control.
     */
    size_type size() {
        if (last_ == buffer_size_) {
            return 0;
        }
        return first_ <= last_ ? last_ - first_ + 1 : queue_.size() + last_ - first_ + 1;
    }

    /**
     * Close the queue and notify any blocked readers.
     */
    void shutdown() {
        std::unique_lock<std::mutex> lock(mutex_);
        if (is_shutdown_) {
            return;
        }
        is_shutdown_ = true;
        not_empty_.notify_all();
        file_write_pool_.join();
        if (output_file_.status == OutputFile::Status::Set && output_file_.stream.is_open()) {
            output_file_.stream.close();
        }
    }

    /**
     * Enqueues x by *moving* it into the queue, blocks when full or when no output
     * file has yet been set.
     * Note that this function receives its parameter by value, so make sure you
     * std::move it into the queue if the copy construction is expensive.
     */
    void push(value_type x) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [this] {
            return !full() && output_file_.status != OutputFile::Status::NotSet;
        });
        bool was_all_read = iterator_.no_more_elements();
        last_ = (last_ == buffer_size_) ? 0 : (last_ + 1) % queue_.size();
        queue_[last_] = std::move(x);
        if (output_file_.status == OutputFile::Status::Set) {
            auto write_to_file = [this](const T v) {
                if (!output_file_.stream.write(reinterpret_cast<const char *>(&v),
                                               sizeof(value_type))) {
                    std::cerr << "Error: Writing of merged data failed." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            };
            file_write_pool_.enqueue(write_to_file, queue_[last_]);
        }
        if (was_all_read) { // queue was empty or all items were read
            not_empty_.notify_one();
        }
    }

    /**
     * Added for compatibility with Vector. Note that this method returns the class'
     * unique iterator, which will point to the first element only before starting
     * iterating.
     */
    Iterator &begin() { return iter(); }

    /**
     * Special iterator indicating the end of the queue - the end is reached when the
     * queue was shut down *and* the iterator past the last element in the queue.
     */
    Iterator &end() { return end_iterator_; }

    /**
     * Sets the name of the output file.
     * The queue's #push_front will block until a (possibly empty) output file was
     * specified.
     */
    void set_out_file(const std::string &output_name) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (output_file_.status == OutputFile::Status::Set && output_file_.stream.is_open()) {
            file_write_pool_.join();
            output_file_.stream.close();
        }
        if (output_name == "") {
            output_file_.status = OutputFile::Status::SetNoOutput;
        } else {
            output_file_.stream
                    = std::fstream(output_name, std::ios::binary | std::ios::out);
            if (!output_file_.stream) {
                std::cerr << "Error: Could not open output file " + output_name << std::endl;
                std::exit(EXIT_FAILURE);
            }
            output_file_.status = OutputFile::Status::Set;
        }
        // we can now start pushing elements into the queue
        not_full_.notify_all();
    }

  private:
    const size_type chunk_size_;
    const size_type fence_size_;
    const size_type buffer_size_;

    std::vector<T, Alloc> queue_;

    /**
     * mutex_ used for synchronizing access to the queue. Both #mutex_ and #not_empty_
     * need to be mutable in order to simulate const iterators (see #iterator_)
     */
    std::mutex mutex_;
    /**
     * Signals that the queue is not empty (an element was added) or that it
     * was shut down.
     */
    std::condition_variable not_empty_;
    std::condition_variable empty_;

    /**
     * Signals that the queue is ready to accept new elements, i.e. it is not full and
     * an output file was set vai #set_out_file.
     */
    std::condition_variable not_full_;

    size_type first_ = 0;
    size_type last_ = buffer_size_;

    /** The queue may optionally write its input to a file */
    struct OutputFile {
        enum class Status { NotSet, SetNoOutput, Set } status;
        std::fstream stream;
    } output_file_;

    // TODO(ddanciu): consider using a WaitQueue instead, to reduce per-task overhead
    ThreadPool file_write_pool_ = ThreadPool(1, 100);

    Iterator iterator_ = Iterator(this);
    Iterator end_iterator_;

    bool is_shutdown_;

  private:
    void pop_chunk() {
        if (size() < chunk_size_) { // nothing to pop
            return;
        }

        const bool was_full = full();

        first_ = (first_ + chunk_size_) % queue_.size();

        if (was_full) { // notify waiting writer that it can start writing again
            not_full_.notify_one();
        }
    }

    /**
     * Returns the iterator of the queue. Note that a queue only has one iterator, so
     * multiple calls to this method will return the same iterator instance.
     * Should not be called from within multiple threads.
     */
    Iterator &iter() {
        std::unique_lock<std::mutex> l(mutex_);
        // if the queue is empty, we don't know if elements will be added at a later time
        // so we need to wait until either an element is added or shutdown() is called
        if (empty()) {
            not_empty_.wait(l, [this]() { return is_shutdown_ || !empty(); });
        }
        if (empty() && is_shutdown_) {
            return end_iterator_;
        }

        return iterator_;
    }
};

/**
 * Defines the (unique) const iterator of a ChunkedWaitQueue. The iterator traverses the
 * queue forward and has support for a limited window of backward operations.
 * As the iterator is moving forward, it will remove older chunks from the queue in
 * order to make room for new elements. An iterator with a nullptr parent denotes the
 * end of the queue.
 */
template <typename T, typename Alloc>
class ChunkedWaitQueue<T, Alloc>::Iterator {
    friend ChunkedWaitQueue;
  public:
    Iterator(const Iterator &other) = delete; // non construction-copyable
    Iterator &operator=(const Iterator &) = delete; // non copyable

    /**
     * Creates an iterator for the given queue.
     * @param parent the ChunkedWaitQueue this class iterates
     */
    Iterator(ChunkedWaitQueue *parent) : parent_(parent) {}

    /**
     * Creates an iterator for the given queue.
     * @param parent the ChunkedWaitQueue this class iterates
     */
    Iterator(ChunkedWaitQueue *parent, size_t idx) : parent_(parent), idx_(idx) {}

    /**
     * Returns the element currently pointed at by the iterator.
     * Undefined behavior if the iterator is pointing at the past-the-end element.
     */
    T operator*() const {
#ifdef DEBUG
        if (idx_ == parent_->buffer_size_) {
            std::cerr << "Attempting to dereference past-the-end iterator." << std::endl;
            std::exit(EXIT_FAILURE);
        }
#endif
        return parent_->queue_[idx_];
    }

    /**
     * Moves the iterator to the next element in the queue.
     * Blocks if no elements are available. If the iterator moved past
     * parent_->chunk_size_ + parent_->fence_size_ from the oldest element, the oldest
     * parent_->chunk_size_ elements are cleaned up from the queue to make room for
     * new elements.
     * @return a pointer to the next element or the special "end iterator" if there
     * are no elements left and the queue was shut down.
     */
    Iterator &operator++() {
        std::unique_lock<std::mutex> l(parent_->mutex_);
        // make some room, if possible
        if (index_dist() + 1 >= parent_->chunk_size_ + parent_->fence_size_) {
            parent_->pop_chunk();
        }
        if (no_more_elements()) {
            parent_->not_empty_.wait(l, [this]() {
                return parent_->is_shutdown_ || !no_more_elements();
            });
            if (parent_->is_shutdown_ && no_more_elements()) { // reached the end
                // notify waiting destructor that object is ready to be destroyed
                parent_->empty_.notify_all();
                idx_ = parent_->buffer_size_;
                return *this;
            }
            // the queue may have filled up while we were sleeping; make some room
            if (index_dist() + 1 >= parent_->chunk_size_ + parent_->fence_size_) {
                parent_->pop_chunk();
            }
        }
        idx_ = (idx_ + 1) % parent_->queue_.size();
        return *this;
    }

    /**
     * Moves the iterator to the previous element in the queue.
     * @throw std::runtime_error if attempting to move before the first element
     */
    Iterator &operator--() {
        std::unique_lock<std::mutex> l(parent_->mutex_);
        if (idx_ == parent_->first_) { // underflow
            throw std::runtime_error("Attempting to move before the first element.");
        }
        if (idx_ == parent_->buffer_size_) {
            idx_ = parent_->last_;
        } else {
            (idx_ > 0) ? idx_-- : idx_ = parent_->queue_.size() - 1;
        }
        return *this;
    }

    /**
     * Returns true if the iterator is pointing at the first available element. Note
     * that this is not necessarily the same as pointing to the first element in the
     * queue - as this element my no longer be available.
     * Useful when iterating backwards - this method returns false if the iterator
     * reached the oldest available element.
     * */
    bool at_begin() {
        return (idx_ == parent_->first_);
    }

    void push_pos() {
        saved_indexes_.push(idx_);
    }

    void pop_pos() {
        idx_ = saved_indexes_.top();
        saved_indexes_.pop();
    }

    bool operator==(const Iterator &other) {
        return parent_ == other.parent_ && idx_ == other.idx_;
    }
    bool operator!=(const Iterator &other) { return !(*this == other); }

    void reset() {
        idx_ = 0;
        saved_indexes_.swap(std::stack<T>());
    }

  private:
    ChunkedWaitQueue *parent_;
    size_type idx_ = 0;
    std::stack<size_type> saved_indexes_;

  private:
    size_type index_dist() {
        return idx_ >= parent_->first_ ? idx_ - parent_->first_
                                       : parent_->size() + idx_ - parent_->first_;
    }

    bool no_more_elements() { return parent_->empty() || idx_ == parent_->last_; }
};

} // namespace common
} // namespace mg
