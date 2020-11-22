#ifndef __THREADING_HPP__
#define __THREADING_HPP__

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <atomic>


void set_num_threads(unsigned int num_threads);
unsigned int get_num_threads();


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
        return emplace(false, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <class F, typename... Args>
    auto force_enqueue(F&& f, Args&&... args) {
        return emplace(true, std::forward<F>(f), std::forward<Args>(args)...);
    }

    void join();

    void remove_waiting_tasks();

    ~ThreadPool();

  private:
    void initialize(size_t num_threads);

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    size_t max_num_tasks_;

    std::mutex queue_mutex;
    std::condition_variable empty_condition;
    std::condition_variable full_condition;

    bool joining_;
    bool stop_;

    template <class F, typename... Args>
    auto emplace(bool force, F&& f, Args&&... args) {
        using return_type = decltype(f(args...));
        auto wrap = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
        auto task_exception = std::make_shared<std::exception_ptr>(nullptr);
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            [f=std::move(wrap),task_exception]() mutable -> return_type {
                try {
                    return f();
                } catch (...) {
                    // catch and store the exception
                    *task_exception = std::current_exception();
                    std::rethrow_exception(*task_exception);
                }
            }
        );

        auto wrapped_task = [task,task_exception]() {
            (*task)();

            // if an exception was caught, throw it here
            if (*task_exception)
                std::rethrow_exception(*task_exception);
        };

        if (!workers.size()) {
            wrapped_task();
            return task->get_future();
        } else {
            std::unique_lock<std::mutex> lock(queue_mutex);
            full_condition.wait(lock, [this,force]() {
                return this->tasks.size() < this->max_num_tasks_ || force;
            });
            tasks.emplace(std::move(wrapped_task));
        }

        empty_condition.notify_one();

        return task->get_future();
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
