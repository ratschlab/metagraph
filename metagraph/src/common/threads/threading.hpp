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
#include <cassert>
#include <chrono>


void set_num_threads(unsigned int num_threads);
unsigned int get_num_threads();

// Compute chunk size for splitting |total_size| items across |num_threads|.
// num_threads = 0 is treated as 1 (single-threaded).
inline size_t get_chunk_size(size_t total_size,
                             size_t max_chunk_size,
                             size_t num_threads,
                             bool one_chunk_if_single_thread = true,
                             size_t chunks_per_thread = 2) {
    num_threads = std::max<size_t>(1, num_threads);
    size_t num_chunks = num_threads * chunks_per_thread;
    return (one_chunk_if_single_thread && num_threads == 1)
            ? std::max<size_t>(1, total_size)
            : std::max<size_t>(1, std::min<size_t>(max_chunk_size,
                                                    (total_size + num_chunks - 1) / num_chunks));
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
        return emplace(false /*front*/, false /*force*/, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <class F, typename... Args>
    auto force_enqueue(F&& f, Args&&... args) {
        return emplace(false /*front*/, true /*force*/, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <class F, typename... Args>
    auto force_enqueue_front(F&& f, Args&&... args) {
        return emplace(true /*front*/, true /*force*/, std::forward<F>(f), std::forward<Args>(args)...);
    }

    // Exception behavior differs between the two execution paths:
    //   * Tasks run on a worker thread: if they throw, `future.get()` inside
    //     `wrapped_task` rethrows on the worker, which escapes the worker
    //     lambda and calls std::terminate (covered by MultiThreadException).
    //   * Tasks stolen by `help_while_waiting` (below) run on the caller, so
    //     an exception propagates out to the caller instead
    //     (covered by HelpWhileWaitingStealedException).
    template <class T>
    std::shared_future<T>& help_while_waiting(std::shared_future<T> &future) {
        auto ready = [&]() {
            return future.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        };
        while (workers.size() && !ready()) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                helper_wakeup.wait_for(lock, std::chrono::microseconds(100), [this, &ready]() {
                    return !tasks.empty() || ready();
                });
                if (ready())
                    break;
                if (tasks.empty())
                    continue;

                task = std::move(tasks.front());
                tasks.pop_front();
            }

            not_full.notify_one();
            task();
        }

        return future;
    }

    void join();

    void remove_waiting_tasks();

    ~ThreadPool();

    size_t num_threads() const { return num_threads_; }

  private:
    void initialize(size_t num_threads);

    size_t num_threads_;
    std::vector<std::thread> workers;
    size_t num_waiting_;
    std::deque<std::function<void()>> tasks;
    size_t max_num_tasks_;

    std::mutex queue_mutex;
    std::condition_variable has_work; // task available or joining
    std::condition_variable not_full; // queue has space
    std::condition_variable all_waiting; // all workers are idle
    std::condition_variable helper_wakeup; // task available for the helper or future is ready

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
            assert(!joining_);
            not_full.wait(lock, [this,force]() {
                return this->tasks.size() < this->max_num_tasks_ || force;
            });
            if (front) {
                tasks.emplace_front(std::move(wrapped_task));
            } else {
                tasks.emplace_back(std::move(wrapped_task));
            }
        }

        has_work.notify_one();
        helper_wakeup.notify_one();

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
