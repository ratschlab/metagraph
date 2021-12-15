#ifndef __BATCH_WORKER_HPP__
#define __BATCH_WORKER_HPP__

#include <cassert>
#include <functional>
#include <limits>
#include <vector>


// A buffer for accumulating data and processing them in batches.
// The data processing is triggered once the buffer capacity is exceeded
// or if the total cost of all data accumulated to that moment exceeds
// a pre-defined threshold.
template <typename T,
          typename CostType = long long int,
          class Buffer = std::vector<T>>
class BatchAccumulator {
  public:
    typedef T value_type;

    BatchAccumulator() {}

    BatchAccumulator(const std::function<void(Buffer&&)> &process_batch,
                     size_t batch_size_limit,
                     CostType batch_cost_limit = std::numeric_limits<CostType>::max(),
                     size_t reserved_size = 0)
          : process_batch_(process_batch),
            batch_size_limit_(batch_size_limit),
            batch_cost_limit_(batch_cost_limit),
            reserved_size_(reserved_size) {
        static_assert(std::is_same<typename Buffer::value_type, T>::value);
        assert(batch_size_limit_);
        assert(batch_cost_limit_ >= data_cost_);

        data_.reserve(reserved_size_);
    }

    ~BatchAccumulator() {
        if (data_.size())
            process_batch_(std::move(data_));
    }

    // push a data point with zero cost
    inline void push(const T &data_point) {
        data_.push_back(data_point);

        if (data_.size() == batch_size_limit_)
            process_buffered_without_check();
    }

    // push a data point with zero cost
    template <typename... Args>
    inline void push(Args&&... args) {
        data_.emplace_back(std::forward<Args>(args)...);

        if (data_.size() == batch_size_limit_)
            process_buffered_without_check();
    }

    // Push a data point into the buffer and increment the batch cost.
    // If the total batch will exceed the limit, process the batch of
    // data points accumulated to the moment first and then insert
    // the new data point.
    template <typename... Args>
    inline void push_and_pay(CostType cost, Args&&... args) {
        if (data_cost_ + cost > batch_cost_limit_ && data_.size())
            process_buffered_without_check();

        data_cost_ += cost;

        push(std::forward<Args>(args)...);
    }

    // process data accumulated to the moment in batch
    inline void process_all_buffered() {
        if (data_.size())
            process_buffered_without_check();
    }

  private:
    void process_buffered_without_check() {
        assert(data_.size());

        Buffer data_batch;
        data_batch.swap(data_);

        process_batch_(std::move(data_batch));

        data_.reserve(reserved_size_);
        data_cost_ = 0;

        assert(data_.empty());
    }

    std::function<void(Buffer&&)> process_batch_;
    size_t batch_size_limit_;
    CostType batch_cost_limit_;
    size_t reserved_size_;

    Buffer data_;
    CostType data_cost_ = 0;
};

#endif // __BATCH_WORKER_HPP__
