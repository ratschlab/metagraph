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
        ThreadPool(size_t num_workers, size_t max_num_tasks = 0);

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
                full_condition.wait(lock, [this]() {
                    return this->tasks.size() < this->max_num_tasks_;
                });
                tasks.emplace([task](){ (*task)(); });
            }
            empty_condition.notify_one();

            return task->get_future();
        }

        void join();

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


    /** A faster alternative to std::allocator<T>
     *
     * The code was copied and has been modified from:
     * https://probablydance.com/2014/11/09/plalloc-a-simple-stateful-allocator-for-node-based-containers/
     */
    template <typename T>
    class plalloc {
      public:
        typedef T value_type;

        plalloc() = default;
        template <typename U>
        plalloc(const plalloc<U>&) {}
        plalloc(const plalloc&) {}
        plalloc& operator=(const plalloc&) { return *this; }
        plalloc(plalloc&&) = default;
        plalloc& operator=(plalloc&&) = default;

        typedef std::true_type propagate_on_container_copy_assignment;
        typedef std::true_type propagate_on_container_move_assignment;
        typedef std::true_type propagate_on_container_swap;

        bool operator==(const plalloc &other) const { return this == &other; }
        bool operator!=(const plalloc &other) const { return !(*this == other); }

        T* allocate(size_t num_to_allocate) {
            if (num_to_allocate != 1)
                return static_cast<T*>(::operator new(sizeof(T) * num_to_allocate));

            if (available.size()) {
                T *result = available.back();
                available.pop_back();
                return result;
            }

            // first allocate 8, then double whenever
            // we run out of memory
            size_t to_allocate = 8 << memory.size();
            available.reserve(to_allocate);
            std::unique_ptr<value_holder[]> allocated(new value_holder[to_allocate]);
            value_holder *first_new = allocated.get();
            memory.emplace_back(std::move(allocated));
            size_t to_return = to_allocate - 1;
            for (size_t i = 0; i < to_return; ++i) {
                available.push_back(std::addressof(first_new[i].value));
            }
            return std::addressof(first_new[to_return].value);
        }
        void deallocate(T *ptr, size_t num_to_free) {
            if (num_to_free == 1) {
                available.push_back(ptr);
            } else {
                ::operator delete(ptr);
            }
        }

        // boilerplate that shouldn't be needed, except
        // libstdc++ doesn't use allocator_traits yet
        template<typename U>
        struct rebind {
            typedef plalloc<U> other;
        };

        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

        template<typename U, typename... Args>
        void construct(U *object, Args&&... args) {
            new (object) U(std::forward<Args>(args)...);
        }
        template<typename U, typename... Args>
        void construct(const U *object, Args &&... args) = delete;

        template<typename U>
        void destroy(U *object) { object->~U(); }

      private:
        union value_holder {
            value_holder() {}
            ~value_holder() {}
            T value;
        };

        std::vector<std::unique_ptr<value_holder[]>> memory;
        std::vector<T*> available;
    };

} // namespace utils

#endif // __UTILS_HPP__
