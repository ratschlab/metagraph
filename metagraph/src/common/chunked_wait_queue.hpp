#pragma once

#include <cassert>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

namespace threads {

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
 */
template <typename T, typename Alloc = std::allocator<T>>
class ChunkedWaitQueue {
  public:
    // typedefs for STL compatibility
    typedef size_t size_type;
    typedef T value_type;

    /**
     * Constructs a WaitQueue with the given size parameters.
     * @param buffer_size the size of the buffer used internally by the queue
     * @param chunk_size the size of the chunks that are popped back when the iterator goes past
     * chunk_size_ + fence_size_ elements past the first element in the buffer
     * @param fence_size the number of elements that can be iterated backwards. The queue will
     * always keep at least fence_size elements behind the furthest iterator for this purpose.
     */
    explicit ChunkedWaitQueue(size_type buffer_size, size_type chunk_size, size_type fence_size)
        : chunk_size_(chunk_size), fence_size_(fence_size), is_shutdown_(false) {
        queue_ = std::vector<T>(buffer_size);
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
     * Returns true if the queue is empty. This can happen only before an element is
     * added. Once an element is added, the queue will never be empty again (because
     * the queue will always keep at least #fence_size_ elements for backwards iteration).
     */
    bool empty() const { return last_ == INVALID_IDX; }

    /**
     * Returns true if the buffer of the queue is full. #push_front operations will block
     * until #iterator() advances far enough in the queue to allow garbage collecting
     * the older elements via #pop_chunk().
     */
    bool full() const {
        return last_ != INVALID_IDX && first_ == (last_ + 1) % queue_.size();
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
        if (last_ == INVALID_IDX) {
            return 0;
        }
        return first_ <= last_ ? last_ - first_ + 1 : queue_.size() + last_ - first_ + 1;
    }

    /**
     * Close the queue and notify any blocked readers.
     */
    void shutdown() {
        std::unique_lock<std::mutex> lock(mutex_);
        is_shutdown_ = true;
        not_empty_.notify_all();
    }

    /**
     * Enqueues x by *moving* it into the queue, blocks when full.
     * Note that this function receives its parameter by value, so make sure you
     * std::move it into the queue if the copy construction is expensive.
     */
    void push_front(value_type x) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [this] { return !full(); });
        bool was_all_read = empty();
        last_ = (last_ + 1) % queue_.size();
        queue_[last_] = std::move(x);
        if (was_all_read) { // queue was empty or all items were read
            not_empty_.notify_one();
        }
    }

    class Iterator;

    /**
     * Special iterator indicating the end of the queue - the end is reached when the
     * queue was shut down *and* the iterator past the last element in the queue.
     */
    Iterator &end() { return end_iterator; }

    /**
     * Returns the iterator of the queue. Note that a queue only has one iterator, so
     * multiple calls to this method will return the same iterator instance.
     * Should not be called from within multiple threads.
     */
    Iterator &iterator() {
        std::unique_lock<std::mutex> l(mutex_);
        // if the queue is empty, we don't know if elements will be added at a later time
        // so we need to wait until either an element is added or shutdown() is called
        if (empty()) {
            not_empty_.wait(l, [this]() { return is_shutdown_ || !empty(); });
        }
        if (empty() && is_shutdown_) {
            return end_iterator;
        }

        return iterator_;
    }

  private:
    static Iterator end_iterator;
    static constexpr size_type INVALID_IDX = std::numeric_limits<size_type>::max();
    size_type chunk_size_;
    size_type fence_size_;
    std::vector<T, Alloc> queue_;
    std::mutex mutex_;
    std::condition_variable empty_;
    /**
     * Signals that the queue is not empty (an element was added) or that it
     * was shut down.
     */
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    size_type first_ = 0;
    size_type last_ = INVALID_IDX;

    Iterator iterator_ = Iterator(this);

    bool is_shutdown_;

  private:
    ChunkedWaitQueue(const ChunkedWaitQueue &other) = delete; // non construction-copyable
    ChunkedWaitQueue &operator=(const ChunkedWaitQueue &) = delete; // non copyable

    void pop_chunk() {
        // this should only be called if we *know* a chunk can be popped safely
        assert(size() >= chunk_size_);

        const bool was_full = full();

        first_ = (first_ + chunk_size_) % queue_.size();

        if (was_full) { // notify waiting writer that it can start writing again
            not_full_.notify_one();
        }
    }
};


template <typename T, typename Alloc>
typename ChunkedWaitQueue<T, Alloc>::Iterator ChunkedWaitQueue<T, Alloc>::end_iterator
        = ChunkedWaitQueue<T, Alloc>::Iterator(nullptr);

/**
 * Defines the (unique) iterator of a ChunkedWaitQueue. The iterator traverses the
 * queue forward and has support for a limited window of backward operations.
 * As the iterator is moving forward, it will remove older chunks from the queue in
 * order to make room for new elements. An iterator with a nullptr parent denotes the
 * end of the queue.
 */
template <typename T, typename Alloc>
class ChunkedWaitQueue<T, Alloc>::Iterator {
  public:
    Iterator(const Iterator &other) = delete; // non construction-copyable
    Iterator &operator=(const Iterator &) = delete; // non copyable

    /**
     * Creates an iterator for the given queue.
     * @param parent the ChunkedWaitQueue this class iterates
     */
    Iterator(ChunkedWaitQueue *parent) : parent_(parent) {}

    /**
     * Returns the element currently pointed at by the iterator.
     * Undefined behavior if the iterator is pointing at the past-the-end element.
     */
    T operator*() {
#ifdef DEBUG
        if (parent_ == nullptr) {
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
        if (no_more_elements()) {
            parent_->not_empty_.wait(l, [this]() {
                return parent_->is_shutdown_ || !no_more_elements();
            });
            if (parent_->is_shutdown_ && no_more_elements()) { // reached the end
                // notify waiting destructor that object is ready to be destroyed
                parent_->empty_.notify_all();
                parent_ = nullptr;
                idx_ = 0;
                return *this;
            }
        }
        idx_ = (idx_ + 1) % parent_->queue_.size();
        if (index_dist() >= parent_->chunk_size_ + parent_->fence_size_) {
            parent_->pop_chunk();
        }
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
        (idx_ > 0) ? idx_-- : idx_ = parent_->queue_.size() - 1;
        return *this;
    }

    bool operator==(const Iterator &other) {
        return parent_ == other.parent_ && idx_ == other.idx_;
    }
    bool operator!=(const Iterator &other) { return !(*this == other); }

  private:
    size_type idx_ = 0;
    ChunkedWaitQueue *parent_;

  private:
    size_type index_dist() {
        return idx_ >= parent_->first_ ? idx_ - parent_->first_
                                       : parent_->size() + idx_ - parent_->first_;
    }

    bool no_more_elements() { return parent_->empty() || idx_ == parent_->last_; }
};

} // namespace threads
