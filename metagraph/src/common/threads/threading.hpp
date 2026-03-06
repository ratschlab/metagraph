#ifndef __THREADING_HPP__
#define __THREADING_HPP__

#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <atomic>


void set_num_threads(unsigned int num_threads);
unsigned int get_num_threads();

inline size_t get_chunk_size(size_t total_size,
                             size_t max_chunk_size,
                             size_t num_threads,
                             bool one_chunk_if_single_thread = true) {
    num_threads = std::max<size_t>(1, num_threads);
    return (one_chunk_if_single_thread && num_threads < 2)
            ? std::max<size_t>(1, total_size)
            : std::max<size_t>(1, std::min<size_t>(max_chunk_size,
                                                    (total_size + num_threads - 1) / num_threads));
}


/**
 * A Thread Pool for parallel execution of tasks with arbitrary parameters
 *
 * The implementation is based on:
 * https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
 */
class ThreadPool {
  public:
    ThreadPool(size_t num_workers, size_t max_num_tasks);
    explicit ThreadPool(size_t num_workers)
        : ThreadPool(num_workers, num_workers * 5) {}

    template <class F, typename... Args>
    auto enqueue(F&& f, Args&&... args) {
        return emplace(false, false, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <class F, typename... Args>
    auto force_enqueue(F&& f, Args&&... args) {
        return emplace(false, true, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <class F, typename... Args>
    auto force_enqueue_front(F&& f, Args&&... args) {
        return emplace(true, true, std::forward<F>(f), std::forward<Args>(args)...);
    }

    void join();

    void remove_waiting_tasks();

    ~ThreadPool();

  private:
    void initialize(size_t num_threads);

    std::vector<std::thread> workers;
    size_t num_waiting_;
    std::deque<std::function<void()>> tasks;
    size_t max_num_tasks_;

    std::mutex queue_mutex;
    std::condition_variable empty_condition;
    std::condition_variable full_condition;
    std::condition_variable all_waiting;

    bool joining_;
    bool stop_;

    template <class F, typename... Args>
    auto emplace(bool front, bool force, F&& f, Args&&... args) {
        using return_type = decltype(f(args...));
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::shared_future<return_type> future(task->get_future());

        auto wrapped_task = [task,future]() {
            (*task)();
            future.get(); // re-thrown exceptions (if any) from packaged_task
        };

        if (!workers.size()) {
            wrapped_task();
            return future;
        } else {
            std::unique_lock<std::mutex> lock(queue_mutex);
            full_condition.wait(lock, [this,force]() {
                return this->tasks.size() < this->max_num_tasks_ || force;
            });
            if (front) {
                tasks.emplace_front(std::move(wrapped_task));
            } else {
                tasks.emplace_back(std::move(wrapped_task));
            }
        }

        empty_condition.notify_one();

        return future;
    }
};


class AsyncActivity {
  public:
    template <class F, typename... Args>
    auto run_async(F&& f, Args&&... args) {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            parallel_jobs_++;
        }
        auto result = f(std::forward<Args>(args)...);
        {
            std::unique_lock<std::mutex> lock(mutex_);
            parallel_jobs_--;
        }
        cond_var_.notify_one();
        return result;
    }

    template <class F, typename... Args>
    auto run_unique(F&& f, Args&&... args) {
        std::unique_lock<std::mutex> lock(mutex_);
        while (parallel_jobs_ > 0) {
            cond_var_.wait(lock);
        }
        return f(std::forward<Args>(args)...);
    }

  private:
    size_t parallel_jobs_ = 0;

    std::mutex mutex_;
    std::condition_variable cond_var_;
};

#endif // __THREADING_HPP__
