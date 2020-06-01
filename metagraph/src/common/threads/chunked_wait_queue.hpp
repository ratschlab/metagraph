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

#include "common/logger.hpp"
#include "threading.hpp"

namespace mtg {
namespace common {

/**
 *  A ChunkedWaitQueue is a data structure that allows safe data transfer from a single
 *  writer thread to a single reader thread. There are a few fundamental differences
 *  between a normal WaitQueue and this implementation:
 *  1. Users have no direct control over how elements are popped from the queue.
 *  Elements are popped from the queue in chunks whenever the single iterator of the queue
 * is far enough from the oldest element.
 *  2. The queue exposes a single iterator instance.  The Iterator instance allows
 * iterating the elements of the queue forward. The iterator will remove elements from the
 * queue as they become old.
 *  3. To reduce the amount of locking, the queue expects one reader and one writer
 *  thread.  Elements are "committed" (made visible to the reader thread) by the #push()
 *  method in batches. The writer thread will only lock once every
 * #write_buf_.capacity() pushes to the queue, so the larger this value the less
 * locking, but also the later the reader thread will see the committed values.
 *
 * The ChunkedWaitQueue uses a pre-allocated circular array to store the elements in
 * order to avoid heap allocations at runtime.
 *
 * Writers can #push a value to the queue, the reader can access elements in order
 * via a single iterator exposed by the class. The reader/iterator will block if there
 * are no available elements, the writers will block if the queue is full i.e. if it
 * contains buffer_size number of elements.
 *
 *  When the queue is shutdown, the reader will unblock. At destruction time, the queue is
 * shut down (if not already shut down) and the destructor will wait until all elements
 * were read (i.e. the iterator reaches past the last available element).
 *
 * The class is not copyable or copy-constructible.
 */
template <typename T, typename Alloc = std::allocator<T>>
class ChunkedWaitQueue {
  public:
    class Iterator;
    // typedefs for STL compatibility
    typedef size_t size_type;
    typedef T value_type;
    typedef Iterator& iterator;

    static constexpr size_t WRITE_BUF_SIZE = 10000;

    ChunkedWaitQueue(const ChunkedWaitQueue &other) = delete;
    ChunkedWaitQueue &operator=(const ChunkedWaitQueue &) = delete;

    /**
     * Constructs a WaitQueue with the given size parameters.
     * @param buffer_size the size of the buffer (number of elements) used internally by
     * the queue. The actual memory footprint is buffer_size*sizeof(T)
     */
    explicit ChunkedWaitQueue(size_type buffer_size)
        : chunk_size_(std::max(1UL, buffer_size / 3)),
          buffer_size_(buffer_size),
          buffer_(buffer_size),
          end_iterator_(Iterator(this, std::min(WRITE_BUF_SIZE, chunk_size_), buffer_size)),
          is_shutdown_(false) {
        write_buf_.reserve(std::min(WRITE_BUF_SIZE, chunk_size_));
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
        iterator_.reset();
    }

    /**
     * Returns the capacity of the queue's internal buffer.
     */
    size_type buffer_size() const { return buffer_.size(); }

    /**
     * Close the queue and notify any blocked readers.
     */
    void shutdown() {
        std::unique_lock<std::mutex> lock(mutex_);
        if (is_shutdown_) {
            return;
        }
        can_flush_.wait(lock, [this] { return can_flush(); });
        flush();
        is_shutdown_ = true;
        not_empty_.notify_all();
    }

    /**
     * Enqueues x by *moving* it into the queue, blocks when full.
     */
    void push(value_type&& x) {
        write_buf_.push_back(std::move(x));
        if (write_buf_.size() == write_buf_.capacity()) {
            std::unique_lock<std::mutex> lock(mutex_);
            can_flush_.wait(lock, [this] { return can_flush(); });
            flush();
        }
    }

    /**
     * Pushes x into the queue, blocks when full.
     */
    void push(const value_type &x) {
        write_buf_.push_back(x);
        if (write_buf_.size() == write_buf_.capacity()) {
            std::unique_lock<std::mutex> lock(mutex_);
            can_flush_.wait(lock, [this] { return can_flush(); });
            flush();
        }
    }

    /**
     * Added for compatibility with Vector. Note that this method returns the class'
     * unique iterator, which will point to the first element only before starting
     * iterating.
     */
    // TODO: construct iterator and return it instead of returning a reference
    Iterator& begin() const { return const_cast<ChunkedWaitQueue *>(this)->iter(); }

    /**
     * Special iterator indicating the end of the queue - the end is reached when the
     * queue was shut down *and* the iterator past the last element in the queue.
     */
    // TODO: construct iterator and return it instead of returning a reference
    Iterator& end() const { return const_cast<ChunkedWaitQueue*>(this)->end_iterator_; }

  private:
    const size_type chunk_size_;
    const size_type buffer_size_;

    std::vector<T, Alloc> buffer_;

    std::vector<T> write_buf_;

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
     * Signals that the queue is ready to accept new elements.
     */
    std::condition_variable can_flush_;

    size_type first_ = 0;
    size_type last_ = buffer_size_;

    Iterator iterator_ = Iterator(this, std::min(WRITE_BUF_SIZE, chunk_size_));
    Iterator end_iterator_;

    bool is_shutdown_;

  private:
    /**
     * Returns true if the queue is empty.
     */
    bool empty() const { return last_ == buffer_size_; }

    /** Returns the number of elements in the buffer. */
    size_type size() {
        if (last_ == buffer_size_) {
            return 0;
        }
        return first_ <= last_ ? last_ - first_ + 1 : buffer_size_ + last_ - first_ + 1;
    }

    /**
     * Returns true if there is enough space in #queue_ to flush #write_buf_ into it.
     * If #can_flush() is false, #push() operations will block until #iterator() advances
     * far enough in the queue to allow garbage collecting the older elements via
     * #pop_chunk().
     */
    bool can_flush() const {
        return empty() || (last_ < first_ && first_ - last_ > write_buf_.size())
                || (last_ >= first_
                    && buffer_size_ - last_ + first_ > write_buf_.size());
    }

    void pop_chunk() {
        const bool could_flush = can_flush();
        first_ += chunk_size_;
        if (first_ >= buffer_size_) {
            first_ -= buffer_size_;
        }

        if (!could_flush && can_flush()) {
            // notify waiting writer that it can start writing  again
            can_flush_.notify_one();
        }
    }

    /** Write the contents of buf to the queue - the caller must hold the mutex. */
    void flush() {
        bool was_all_read = !iterator_.can_increment();
        for (auto &v : write_buf_) {
            if (++last_ >= buffer_size_) {
                last_ = 0;
            }
            buffer_[last_] = std::move(v);
        }
        if (was_all_read) { // queue was empty or all items were read
            not_empty_.notify_all();
        }
        write_buf_.resize(0);
    }


    /**
     * Returns the iterator of the queue. Note that a queue only has one iterator, so
     * multiple calls to this method will return the same iterator instance.
     * Should not be called from within multiple threads.
     */
    Iterator &iter() {
        std::unique_lock<std::mutex> lock(mutex_);
        // if the queue is empty, we don't know if elements will be added at a later time
        // so we need to wait until either an element is added or shutdown() is called
        if (empty()) {
            not_empty_.wait(lock, [this]() { return is_shutdown_ || !empty(); });
        }
        if (empty() && is_shutdown_) {
            return end_iterator_;
        }
        iterator_.init();
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
    Iterator(ChunkedWaitQueue *parent, size_t read_buf_size)
        : Iterator(parent, read_buf_size, 0) {}

    /**
     * Creates an iterator for the given queue at a specific position.
     * @param parent the ChunkedWaitQueue this class iterates
     */
    Iterator(ChunkedWaitQueue *parent, size_t read_buf_size, size_t idx)
        : queue_(parent), idx_(idx), read_buf_size_(read_buf_size) {
    }

    /**
     * Returns the element currently pointed at by the iterator.
     * Undefined behavior if the iterator is pointing at the past-the-end element.
     */
    T operator*() const { return read_buf_[read_buf_idx_]; }

    /**
     * Moves the iterator to the next element in the queue.
     * Blocks if no elements are available. If the iterator moved past
     * parent_->chunk_size_ from the oldest element, the oldest parent_->chunk_size_
     * elements are cleaned up from the queue to make room for new elements.
     * @return a pointer to the next element or the special "end iterator" if there
     * are no elements left and the queue was shut down.
     */
    Iterator &operator++() {
        read_buf_idx_++;
        if (read_buf_idx_ < read_buf_.size()) {
            return *this;
        }
        read_buf_idx_ = 0; // moved here to reduce size of critical section
        std::unique_lock<std::mutex> l(queue_->mutex_);
        // make some room, if possible
        if (elements_read() > queue_->chunk_size_) {
            queue_->pop_chunk();
        }
        if (!can_read_from_queue()) {
            queue_->not_empty_.wait(l, [this]() {
                return queue_->is_shutdown_ || can_read_from_queue();
            });
            if (queue_->is_shutdown_ && !can_increment()) { // reached the end
                idx_ = queue_->buffer_size_;
                // notify waiting destructor that object is ready to be destroyed
                queue_->empty_.notify_all();
                return *this;
            }
            // the queue may have filled up while we were sleeping; make some room
            if (elements_read() > queue_->chunk_size_) {
                queue_->pop_chunk();
            }
        }
        read_buf_.resize(read_buf_size_);
        size_t i;
        for (i = 0; i < read_buf_.size() && idx_ != queue_->last_ ; ++i) {
            if (++idx_ == queue_->buffer_size_) {
                idx_ = 0;
            }
            read_buf_[i] = queue_->buffer_[idx_];
        }
        if (i < read_buf_.size()) { // only happens if queue was shut down
            assert(queue_->is_shutdown_ && !can_increment());
            read_buf_.resize(i);
        }
        return *this;
    }

    bool operator==(const Iterator &other) const {
        return queue_ == other.queue_ && idx_ == other.idx_;
    }
    bool operator!=(const Iterator &other) const { return !(*this == other); }

    void reset() {
        idx_ = 0;
        read_buf_.resize(0);
        read_buf_idx_ = 0;
    }

  private:
    ChunkedWaitQueue *queue_;
    size_type idx_ = 0;
    std::vector<T> read_buf_;
    size_t read_buf_idx_ = 0;
    size_t read_buf_size_;

  private:
    /** Returns the number of elements between the current element and the oldest */
    size_type elements_read() const {
        // idx_ points to the current element, hence the '+1' in the result
        return idx_ >= queue_->first_ ? idx_ - queue_->first_ + 1
                                      : queue_->buffer_size_ + idx_ - queue_->first_ + 1;
    }

    /** Returns true if the iterator can be incremented without blocking */
    bool can_increment() const {
        return !queue_->empty() && idx_ != queue_->last_ && idx_ != queue_->buffer_size_;
    }

    /** Returns true if there are enough unread elements to fill #read_buf_ */
    bool can_read_from_queue() const {
        if (idx_ == queue_->buffer_size_) {
            return false; // this is the end iterator
        }
        uint32_t unread_element_count = idx_ <= queue_->last_
                ? queue_->last_ - idx_
                : queue_->buffer_size_ + queue_->last_ - idx_;

        return unread_element_count >= read_buf_size_;
    }

    /**
     * Initializes the read iterator's buffer. Must be called before iterator is used for
     * the first time, but after the parent has enough elements to initialize the buffer.
     */
    void init() {
        uint32_t el_count = std::min(read_buf_size_, queue_->last_ + 1);
        assert(el_count == read_buf_size_ || queue_->is_shutdown_);
        read_buf_.resize(el_count);
        std::copy_n(queue_->buffer_.begin(), el_count, read_buf_.begin());
        idx_ = el_count - 1;
    }
};

} // namespace common
} // namespace mtg
