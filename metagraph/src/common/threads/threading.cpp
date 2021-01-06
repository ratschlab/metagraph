#include "threading.hpp"

#include <cassert>

static unsigned int NUM_THREADS_METAGRAPH_GLOBAL = 1;


void set_num_threads(unsigned int num_threads) {
    NUM_THREADS_METAGRAPH_GLOBAL = std::max(1u, num_threads);
}

unsigned int get_num_threads() {
    return NUM_THREADS_METAGRAPH_GLOBAL;
}


ThreadPool::ThreadPool(size_t num_workers, size_t max_num_tasks)
      : max_num_tasks_(std::max(max_num_tasks, size_t(1))), stop_(false) {
    initialize(num_workers);
}

ThreadPool::~ThreadPool() {
    stop_ = true;
    join();
}

void ThreadPool::join() {
    size_t num_workers = workers.size();

    if (!num_workers) {
        return;
    } else {
        std::lock_guard<std::mutex> lock(queue_mutex);
        assert(!joining_);
        joining_ = true;
    }
    empty_condition.notify_all();

    for (std::thread &worker : workers) {
        worker.join();
    }
    workers.clear();

    if (!stop_)
        initialize(num_workers);
}

void ThreadPool::remove_waiting_tasks() {
    std::unique_lock<std::mutex> lock(this->queue_mutex);
    this->tasks = std::queue<std::function<void()>>();
    this->empty_condition.notify_all();
}

void ThreadPool::initialize(size_t num_workers) {
    assert(!stop_);
    assert(workers.size() == 0);
    joining_ = false;

    if (!num_workers)
        return;

    for (size_t i = 0; i < num_workers; ++i) {
        workers.emplace_back([this]() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->empty_condition.wait(lock, [this]() {
                        return this->joining_ || !this->tasks.empty();
                    });
                    if (this->tasks.empty())
                        return;

                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }

                full_condition.notify_one();
                task();
            }
        });
    }
}
