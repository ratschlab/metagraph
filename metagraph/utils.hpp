#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <queue>
#include <atomic>

class DBG_succ;
class KMer;


namespace utils {

    uint64_t kFromFile(const std::string &infbase);

    /**
    * This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
    * as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
    * returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
    * second value set to true if G2(k2_node) < G1(k1_node).
    */
    std::pair<bool, bool> compare_nodes(const DBG_succ *G1, uint64_t k1_node,
                                        const DBG_succ *G2, uint64_t k2_node);

    /**
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    std::string get_filetype(const std::string &fname);

    std::deque<std::string> generate_strings(const std::string &alphabet,
                                             size_t length);

    void radix_sort(std::vector<KMer> &data, size_t k);
    void bucket_sort(std::vector<KMer> &data, size_t k);


    /** This returns the currently used memory by the process.
     *
     * The code was copied and has been modified from:
     * https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
     */
    class ThreadPool {
      public:
        ThreadPool(size_t num_workers);

        template <class F, typename... Args>
        auto enqueue(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
            using return_type = typename std::result_of<F(Args...)>::type;
            auto task = std::make_shared<std::packaged_task<return_type()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

            if (!workers.size()) {
                (*task)();
                return task->get_future();
            } else {
                std::unique_lock<std::mutex> lock(queue_mutex);
                tasks.emplace([task](){ (*task)(); });
            }
            condition.notify_one();

            return task->get_future();
        }

        void join();

        ~ThreadPool();

      private:
        void initialize(size_t num_threads);

        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;

        std::mutex queue_mutex;
        std::condition_variable condition;

        bool joining_;
        bool stop_;
    };


    class AsyncActivity {
      public:
        template <class F, typename... Args>
        auto run_async(F&& f, Args&&... args) -> decltype(f(args...)) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                parallel_jobs_++;
            }
            auto result = f(args...);
            {
                std::unique_lock<std::mutex> lock(mutex_);
                parallel_jobs_--;
            }
            cond_var_.notify_one();
            return std::move(result);
        }

        template <class F, typename... Args>
        auto run_unique(F&& f, Args&&... args) -> decltype(f(args...)) {
            std::unique_lock<std::mutex> lock(mutex_);
            while (parallel_jobs_ > 0) {
                cond_var_.wait(lock);
            }
            return f(args...);
        }

      private:
        size_t parallel_jobs_ = 0;

        std::mutex mutex_;
        std::condition_variable cond_var_;
    };

} // namespace utils

#endif // __UTILS_HPP__
