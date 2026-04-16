#include "test_helpers.hpp"

#include <algorithm>
#include <atomic>
#include <future>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <sdsl/int_vector_buffer.hpp>

#include "common/vector.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/algorithms.hpp"
#include "common/vectors/transpose.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "seq_io/formats.hpp"


namespace {

using namespace mtg;
using namespace mtg::seq_io;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

const std::vector<sdsl::bit_vector> vectors {
    { 0, 1, 1, 0, 0, 1, 0 },
    { 0, 1, 1, 1, 1, 1, 0 }
};

const std::vector<sdsl::bit_vector> matrix {
    { 0, 0 },
    { 1, 1 },
    { 1, 1 },
    { 0, 1 },
    { 0, 1 },
    { 1, 1 },
    { 0, 0 }
};

const std::vector<std::vector<uint64_t>> indices {
    { 1, 0 }, { 1, 1 },
    { 2, 0 }, { 2, 1 },
    { 3, 1 },
    { 4, 1 },
    { 5, 0 }, { 5, 1 }
};


std::vector<bit_vector_small> generate_rows() {
    std::vector<bit_vector_small> rows;
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(rows),
        [](const auto &a) {
            return bit_vector_small(a);
        });

    return rows;
}

TEST(Utils, call_rows) {
    auto rows = generate_rows();

    std::vector<const bit_vector_small*> rows_ptr;
    for (const auto &row : rows) {
        rows_ptr.push_back(&row);
    }

    uint64_t i = 0;
    std::function<void(const Vector<uint64_t> &)> call_row = [&](const Vector<uint64_t> &row_indices) {
        sdsl::bit_vector bv(rows_ptr.size());
        for (auto j : row_indices) {
            bv[j] = 1;
        }
        EXPECT_EQ(matrix[i++], bv) << i;
    };
    utils::call_rows(rows_ptr, call_row);

    ASSERT_EQ(7u, i);
}

TEST(Utils, TempFileInitialize) {
    utils::TempFile tmp;
}

TEST(Utils, TempFileOpenWrite) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
    }
}

TEST(Utils, TempFileOpenWriteRead) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
    }
}

TEST(Utils, TempFileCheckStateFlow) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEBUG_DEATH((void)tmp.ofstream().good(), "Can't write after reading");
    }
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEBUG_DEATH((void)tmp.ofstream().good(), "Can't write after reading");
    }
}

TEST(Utils, TempFileReadWritePairs) {
    {
        utils::TempFile tmp[2];

        ASSERT_TRUE(tmp[0].ofstream().good());
        ASSERT_TRUE(tmp[0].ifstream().good());

        ASSERT_TRUE(tmp[1].ofstream().good());
        ASSERT_TRUE(tmp[1].ifstream().good());

        EXPECT_TRUE(tmp[0].ifstream().peek() == std::ifstream::traits_type::eof());
        EXPECT_TRUE(tmp[1].ifstream().peek() == std::ifstream::traits_type::eof());
    }

    utils::TempFile tmp[2];

    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;
    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;

    std::ifstream& in0 = tmp[0].ifstream();
    std::ifstream& in1 = tmp[1].ifstream();

    std::string temp;
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    EXPECT_FALSE(in0 >> temp);
    EXPECT_FALSE(in1 >> temp);
    EXPECT_TRUE(in0.peek() == std::ifstream::traits_type::eof());
    EXPECT_TRUE(in1.peek() == std::ifstream::traits_type::eof());
}

template <typename T>
std::vector<T> uniqueize(const std::vector<T> &input) {
    auto vector = input;
    vector.erase(std::unique(vector.begin(), vector.end()), vector.end());
    return vector;
}

using utils::sample_indexes;

TEST(Utils, sample_indexes) {
    std::mt19937 gen;
    gen.seed(14);

    ASSERT_EQ(0u, sample_indexes(10, 0, gen).size());
    ASSERT_EQ(5u, sample_indexes(10, 5, gen).size());
    ASSERT_EQ(6u, sample_indexes(10, 6, gen).size());
    ASSERT_EQ(9u, sample_indexes(10, 9, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 10, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 11, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 1000, gen).size());

    ASSERT_EQ(0u, sample_indexes(0, 0, gen).size());
    ASSERT_EQ(0u, sample_indexes(0, 1000, gen).size());

    auto generated = sample_indexes(10, 0, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 5, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 6, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 9, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 10, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 11, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 1000, gen);
    ASSERT_EQ(generated, uniqueize(generated));

    generated = sample_indexes(100'000, 1'000, gen);
    EXPECT_TRUE(*std::min_element(generated.begin(), generated.end()) < 100'000u / 5);
    EXPECT_TRUE(*std::max_element(generated.begin(), generated.end()) > 4 * 100'000u / 5)
        << *std::max_element(generated.begin(), generated.end());

    ASSERT_EQ(2'000'000u, sample_indexes(100'000'000, 2'000'000, gen).size());
    ASSERT_EQ(12'000'000u, sample_indexes(100'000'000, 12'000'000, gen).size());
}


TEST(get_filetype, VCF) {
    EXPECT_EQ("VCF", file_format("file.VCF"));
    EXPECT_EQ("VCF", file_format("file.vcf"));
    EXPECT_EQ("VCF", file_format("file.VCF.gz"));
    EXPECT_EQ("VCF", file_format("file.vcf.gz"));
}

TEST(get_filetype, FASTA) {
    EXPECT_EQ("FASTA", file_format("file.FASTA"));
    EXPECT_EQ("FASTA", file_format("file.fasta"));
    EXPECT_EQ("FASTA", file_format("file.FASTA.gz"));
    EXPECT_EQ("FASTA", file_format("file.fasta.gz"));
}

TEST(get_filetype, FASTQ) {
    EXPECT_EQ("FASTQ", file_format("file.fq"));
    EXPECT_EQ("FASTQ", file_format("file.FQ"));
    EXPECT_EQ("FASTQ", file_format("file.fastq"));
    EXPECT_EQ("FASTQ", file_format("file.fq.gz"));
    EXPECT_EQ("FASTQ", file_format("file.FQ.gz"));
}

TEST(get_filetype, IncorrectFiletype) {
    EXPECT_EQ("", file_format("fq"));
    EXPECT_EQ("", file_format("fasta"));
    EXPECT_EQ("", file_format("bz2"));
    EXPECT_EQ("", file_format("file.gz"));
    EXPECT_EQ("", file_format("fasta.gz"));
}

bool colexicographically_greater(const std::string &s1, const std::string &s2) {
    return utils::colexicographically_greater(s1, s2);
}


TEST(colexicographically_greater, basics) {
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAA"));

    EXPECT_TRUE(colexicographically_greater("BAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("AAAB", "AAAA"));

    EXPECT_FALSE(colexicographically_greater("AAAA", "BAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAB"));
}

TEST(colexicographically_greater, different_lengths) {
    EXPECT_TRUE(colexicographically_greater("AAAAA", "AAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAAA", "BBBB"));

    EXPECT_TRUE(colexicographically_greater("AAAAB", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("ABAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("BAAAA", "AAAA"));
}


TEST(seq_equal, basics) {
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("ABAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 1));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 2));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 3));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 4));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 5));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 6));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 7));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 100));
}

TEST(ThreadPool, EmptyConstructor) {
    ThreadPool pool(1);
    ThreadPool pool2(2);
    ThreadPool pool4(4);
    ThreadPool pool20(20);
}

TEST(ThreadPool, EmptyTasks) {
    for (size_t i = 1; i < 20; ++i) {
        ThreadPool pool(i);
        for (size_t t = 0; t < 1000; ++t) {
            pool.enqueue([]() {});
        }
    }
}

TEST(ThreadPool, EmptyJoin) {
    ThreadPool pool(1);
    ThreadPool pool2(2);
    ThreadPool pool4(4);
    ThreadPool pool20(20);
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
}

TEST(ThreadPool, SingleThread) {
    ThreadPool pool(1);
    std::vector<size_t> result;
    std::mutex mu;

    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);

    pool.join();

    ASSERT_EQ(5u, result.size());
    for (int value : result) {
        ASSERT_EQ(1, value);
    }
}

TEST(ThreadPool, MultiThreadTwo) {
    ThreadPool pool(2);
    size_t num_tasks = 10000;

    std::vector<size_t> result;
    std::mutex mu;
    for (size_t t = 0; t < num_tasks; ++t) {
        pool.enqueue([&](size_t i) {
            size_t count = 3;
            for (size_t k = 0; k < num_tasks + i; ++k) {
                count += count * 3 - count / 2 + 1;
            }
            mu.lock();
            result.push_back(count);
            mu.unlock();
        }, t);
    }

    pool.join();

    ASSERT_EQ(num_tasks, result.size());
}

TEST(ThreadPool, MultiThreadFour) {
    ThreadPool pool(4);
    size_t num_tasks = 10000;

    std::vector<size_t> result;
    std::mutex mu;
    for (size_t t = 0; t < num_tasks; ++t) {
        pool.enqueue([&](size_t i) {
            size_t count = 3;
            for (size_t k = 0; k < num_tasks + i; ++k) {
                count += count * 3 - count / 2 + 1;
            }
            mu.lock();
            result.push_back(count);
            mu.unlock();
        }, t);
    }

    pool.join();

    ASSERT_EQ(num_tasks, result.size());
}

TEST(ThreadPool, MultiThread) {
    for (size_t i = 2; i < 20; ++i) {
        ThreadPool pool(i);
        std::vector<size_t> result;
        std::mutex mu;
        for (size_t t = 0; t < 1000; ++t) {
            pool.enqueue([&](size_t i) {
                mu.lock();
                result.push_back(i);
                mu.unlock();
            }, 1);
        }

        pool.join();

        ASSERT_EQ(1000u, result.size());
        for (int value : result) {
            ASSERT_EQ(1, value);
        }
    }
}

TEST(ThreadPool, MultiThreadFuture) {
    for (size_t i = 2; i < 20; ++i) {
        ThreadPool pool(i);

        std::vector<std::shared_future<size_t>> result;
        for (size_t t = 0; t < 1000; ++t) {
            result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
        }

        ASSERT_EQ(1000u, result.size());
        for (auto &value : result) {
            ASSERT_EQ(1u, value.get());
        }
    }
}

#ifndef _NO_DEATH_TEST
// Exceptions thrown in worker threads must kill the process (via std::terminate)
// even if the caller never calls future.get().
TEST(ThreadPool, MultiThreadException) {
    for (size_t num_workers : { 1, 10 }) {
        ASSERT_DEATH({
            ThreadPool pool(num_workers);
            pool.enqueue([]() { throw std::runtime_error("test exception"); });
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }, "");
    }
}
#endif // ifndef _NO_DEATH_TEST

// With 0 workers (synchronous), exceptions propagate back to the caller.
TEST(ThreadPool, SynchronousException) {
    ThreadPool pool(0);
    EXPECT_THROW(
        pool.enqueue([]() -> int { throw std::runtime_error("sync throw"); }),
        std::runtime_error
    );
}

// When help_while_waiting steals a throwing task, the exception propagates
// to the caller rather than crashing a worker thread.
TEST(ThreadPool, HelpWhileWaitingStealedException) {
    ThreadPool pool(1);

    // Block the only worker so it can't pick up subsequent tasks.
    std::promise<void> unblock;
    pool.enqueue([&]() { unblock.get_future().wait(); });
    std::this_thread::sleep_for(std::chrono::milliseconds(5));

    // Enqueue a task that throws — it stays in the queue since the worker is busy.
    pool.force_enqueue([]() -> int { throw std::runtime_error("stolen throw"); });

    // Enqueue the task we actually wait for — also stays in the queue.
    auto future = pool.force_enqueue([]() { return 42; });

    // help_while_waiting will steal the throwing task first, propagating the exception.
    EXPECT_THROW(pool.help_while_waiting(future).get(), std::runtime_error);

    unblock.set_value();
}

TEST(ThreadPool, DummyEmptyJoin) {
    ThreadPool pool(0);
    pool.join();
    pool.join();
    pool.join();
}

TEST(ThreadPool, DummyPoolRun) {
    ThreadPool pool(0);
    std::vector<size_t> result;
    std::mutex mu;
    for (size_t t = 0; t < 1000; ++t) {
        pool.enqueue([&](size_t i) {
            mu.lock();
            result.push_back(i);
            mu.unlock();
        }, 1);
    }

    pool.join();

    ASSERT_EQ(1000u, result.size());
    for (int value : result) {
        ASSERT_EQ(1, value);
    }
}

TEST(ThreadPool, DummyFuture) {
    ThreadPool pool(0);

    std::vector<std::shared_future<size_t>> result;
    for (size_t t = 0; t < 1000; ++t) {
        result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
    }

    ASSERT_EQ(1000u, result.size());
    for (auto &value : result) {
        ASSERT_EQ(1u, value.get());
    }
}


// ---------------------------------------------------------------------------
//  ThreadPool: force_enqueue_front ordering
// ---------------------------------------------------------------------------

// With 1 worker blocked, queued tasks accumulate. force_enqueue_front should
// place tasks at the head, so the worker processes them first after unblocking.
TEST(ThreadPool, ForceEnqueueFrontOrdering) {
    ThreadPool pool(1);

    std::promise<void> worker_started;
    std::promise<void> unblock;

    pool.enqueue([&]() {
        worker_started.set_value();
        unblock.get_future().wait();
    });
    worker_started.get_future().wait();

    std::vector<int> order;
    std::mutex mu;
    auto make_task = [&](int id) {
        return [&mu, &order, id]() {
            std::lock_guard<std::mutex> lock(mu);
            order.push_back(id);
        };
    };

    // Queue: push 1 back, push 2 back, push 3 front → [3, 1, 2]
    pool.force_enqueue(make_task(1));
    pool.force_enqueue(make_task(2));
    pool.force_enqueue_front(make_task(3));

    unblock.set_value();
    pool.join();

    ASSERT_EQ(3u, order.size());
    EXPECT_EQ(3, order[0]);
    EXPECT_EQ(1, order[1]);
    EXPECT_EQ(2, order[2]);
}

// ---------------------------------------------------------------------------
//  ThreadPool: force_enqueue bypasses queue limit
// ---------------------------------------------------------------------------

TEST(ThreadPool, ForceEnqueueBypassesQueueLimit) {
    // max_num_tasks = 2
    ThreadPool pool(1, 2);

    std::promise<void> worker_started;
    std::promise<void> unblock;

    pool.enqueue([&]() {
        worker_started.set_value();
        unblock.get_future().wait();
    });
    worker_started.get_future().wait();

    // Fill the queue to max_num_tasks
    pool.force_enqueue([]() {});
    pool.force_enqueue([]() {});

    // Queue is full. force_enqueue should still succeed (no deadlock).
    std::atomic<bool> task_ran{false};
    pool.force_enqueue([&]() { task_ran = true; });

    unblock.set_value();
    pool.join();

    EXPECT_TRUE(task_ran.load());
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — basic correctness
// ---------------------------------------------------------------------------

TEST(ThreadPool, HelpWhileWaitingBasic) {
    for (size_t num_threads : {1, 2, 4, 8}) {
        ThreadPool pool(num_threads);
        auto future = pool.enqueue([]() { return 42; });
        EXPECT_EQ(42, pool.help_while_waiting(future).get());
    }
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — caller must steal work to avoid deadlock
// ---------------------------------------------------------------------------

// The only worker is blocked. The future's task sits in the queue.
// Without work-stealing in help_while_waiting, this would deadlock.
TEST(ThreadPool, HelpWhileWaitingCallerStealsWork) {
    ThreadPool pool(1);

    std::promise<void> worker_started;
    std::promise<void> unblock;

    pool.enqueue([&]() {
        worker_started.set_value();
        unblock.get_future().wait();
    });
    worker_started.get_future().wait();

    // Worker is blocked. These tasks can only run if the caller steals them.
    std::atomic<std::thread::id> executor_thread_id{};
    pool.force_enqueue([&]() {
        executor_thread_id = std::this_thread::get_id();
    });

    auto future = pool.force_enqueue([]() { return 99; });

    // Must not deadlock — caller processes queued tasks.
    EXPECT_EQ(99, pool.help_while_waiting(future).get());

    // The intermediate task was executed on the calling thread, not a worker.
    EXPECT_EQ(std::this_thread::get_id(), executor_thread_id.load());

    unblock.set_value();
    pool.join();
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — no lost wakeup with many idle workers
// ---------------------------------------------------------------------------

// Tests that help_while_waiting wakes up promptly when the future completes.
// Queue is empty, many workers are idle, one worker is executing the future's
// task. When it completes, help_while_waiting must be woken via help_condition
// (the dedicated CV), not delayed by idle workers competing for notifications.
TEST(ThreadPool, HelpWhileWaitingNoLostWakeup) {
    for (int iter = 0; iter < 100; ++iter) {
        ThreadPool pool(8);

        std::promise<void> task_started;
        std::promise<void> release;

        // One task occupies one worker; the other 7 are idle.
        auto future = pool.enqueue([&]() {
            task_started.set_value();
            release.get_future().wait();
            return iter;
        });

        // Ensure the task is running (not still in the queue).
        task_started.get_future().wait();

        // Queue is now empty, 7 workers are idle. Release the task so
        // the future completes while help_while_waiting is waiting.
        release.set_value();

        // Must not deadlock.
        EXPECT_EQ(iter, pool.help_while_waiting(future).get());
    }
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — chunked work pattern
// ---------------------------------------------------------------------------

// Mimics the BRWT get_nonzero_rows pattern: enqueue many chunks, then
// call help_while_waiting for each future sequentially.
TEST(ThreadPool, HelpWhileWaitingChunkedWork) {
    ThreadPool pool(4);

    const int num_chunks = 40;
    std::vector<std::shared_future<int>> futures;
    for (int i = 0; i < num_chunks; ++i) {
        futures.push_back(pool.force_enqueue_front([i]() {
            int sum = 0;
            for (int j = 0; j < 1000; ++j) sum += i + j;
            return sum;
        }));
    }

    int total = 0;
    for (auto &f : futures) {
        total += pool.help_while_waiting(f).get();
    }

    int expected = 0;
    for (int i = 0; i < num_chunks; ++i) {
        for (int j = 0; j < 1000; ++j) expected += i + j;
    }
    EXPECT_EQ(expected, total);
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — dummy pool (0 workers)
// ---------------------------------------------------------------------------

TEST(ThreadPool, HelpWhileWaitingDummyPool) {
    ThreadPool pool(0);

    // 0 workers → task is executed inline in emplace(), future is already ready.
    auto future = pool.enqueue([]() { return 42; });
    EXPECT_EQ(std::future_status::ready,
              future.wait_for(std::chrono::seconds(0)));

    // help_while_waiting returns immediately (workers.size() == 0).
    EXPECT_EQ(42, pool.help_while_waiting(future).get());
}

// ---------------------------------------------------------------------------
//  ThreadPool: tasks that spawn subtasks — join must wait for all
// ---------------------------------------------------------------------------

TEST(ThreadPool, JoinWaitsForDynamicallySpawnedTasks) {
    ThreadPool pool(4);
    std::atomic<int> completed{0};

    for (int i = 0; i < 50; ++i) {
        pool.enqueue([&pool, &completed]() {
            // Each parent spawns 10 children.
            for (int j = 0; j < 10; ++j) {
                pool.force_enqueue([&completed]() { completed++; });
            }
            completed++;
        });
    }

    pool.join();
    EXPECT_EQ(550, completed.load());
}

// ---------------------------------------------------------------------------
//  ThreadPool: join-reinitialize cycle
// ---------------------------------------------------------------------------

TEST(ThreadPool, JoinReinitializeCycle) {
    ThreadPool pool(4);

    for (int cycle = 0; cycle < 5; ++cycle) {
        std::atomic<int> counter{0};
        for (int i = 0; i < 200; ++i) {
            pool.enqueue([&counter]() { counter++; });
        }
        pool.join();
        EXPECT_EQ(200, counter.load()) << "cycle " << cycle;
    }
}

// ---------------------------------------------------------------------------
//  ThreadPool: num_threads accessor
// ---------------------------------------------------------------------------

TEST(ThreadPool, NumThreads) {
    ThreadPool pool0(0);
    EXPECT_EQ(0u, pool0.num_threads());

    ThreadPool pool1(1);
    EXPECT_EQ(1u, pool1.num_threads());

    ThreadPool pool8(8);
    EXPECT_EQ(8u, pool8.num_threads());
}

// ---------------------------------------------------------------------------
//  ThreadPool: remove_waiting_tasks cancels pending work
// ---------------------------------------------------------------------------

TEST(ThreadPool, RemoveWaitingTasksCancelsPendingWork) {
    ThreadPool pool(1);

    std::promise<void> worker_started;
    std::promise<void> unblock;

    pool.enqueue([&]() {
        worker_started.set_value();
        unblock.get_future().wait();
    });
    worker_started.get_future().wait();

    // Worker is blocked. Queue up 100 tasks.
    std::atomic<int> counter{0};
    for (int i = 0; i < 100; ++i) {
        pool.force_enqueue([&counter]() { counter++; });
    }

    // Discard all pending tasks.
    pool.remove_waiting_tasks();

    unblock.set_value();
    pool.join();

    // The worker was busy, so none (or very few) of the 100 tasks ran.
    EXPECT_LT(counter.load(), 100);
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting under contention — many futures, few workers
// ---------------------------------------------------------------------------

// All workers are occupied with slow tasks. The caller must steal and execute
// fast tasks from the queue to complete quickly. Tests that the caller makes
// forward progress rather than starving.
TEST(ThreadPool, HelpWhileWaitingCallerMakesProgress) {
    ThreadPool pool(2);

    std::promise<void> workers_started;
    std::atomic<int> started_count{0};
    std::promise<void> unblock;
    auto unblock_future = unblock.get_future().share();

    // Occupy both workers with long tasks.
    for (int w = 0; w < 2; ++w) {
        pool.enqueue([&started_count, &workers_started, unblock_future]() {
            if (++started_count == 2)
                workers_started.set_value();
            unblock_future.wait();
        });
    }
    workers_started.get_future().wait();

    // Both workers are blocked. Enqueue work that only the caller can process.
    std::atomic<int> tasks_executed{0};
    for (int i = 0; i < 20; ++i) {
        pool.force_enqueue([&tasks_executed]() { tasks_executed++; });
    }
    auto future = pool.force_enqueue([]() { return true; });

    EXPECT_TRUE(pool.help_while_waiting(future).get());

    // The caller executed the future's task. It also processed tasks ahead of
    // the future in the queue.
    EXPECT_GT(tasks_executed.load(), 0);

    unblock.set_value();
    pool.join();

    // All 20 tasks eventually complete (some by caller, rest by workers after unblock).
    EXPECT_EQ(20, tasks_executed.load());
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — timeout with empty queue must not crash
// ---------------------------------------------------------------------------

// When help_while_waiting's wait_for times out with an empty queue and the
// future not yet ready, it must loop back — not access tasks.front() on an
// empty deque (undefined behavior / crash).
TEST(ThreadPool, HelpWhileWaitingTimeoutEmptyQueue) {
    ThreadPool pool(1);

    // The worker picks up this slow task immediately — queue becomes empty.
    // The task takes much longer than the 100μs wait_for timeout.
    auto future = pool.enqueue([]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        return 42;
    });

    // Small delay so the worker dequeues the task before we call HWW.
    std::this_thread::sleep_for(std::chrono::milliseconds(1));

    // Queue is empty, future not ready (worker still sleeping for 50ms).
    // help_while_waiting enters wait_for(100μs), times out, predicate false.
    // Without the tasks.empty() guard, this would access an empty deque.
    EXPECT_EQ(42, pool.help_while_waiting(future).get());
}

// Same scenario but with multiple workers — the future's task is picked up
// by one worker while others are idle. The queue is empty for extended time.
TEST(ThreadPool, HelpWhileWaitingTimeoutEmptyQueueMultiWorker) {
    ThreadPool pool(4);

    auto future = pool.enqueue([]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        return 99;
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(1));

    // 4 workers, 3 idle, 1 executing the slow task. Queue empty.
    // help_while_waiting will time out many times (~500x) before the future
    // is ready. Each timeout must not access the empty deque.
    EXPECT_EQ(99, pool.help_while_waiting(future).get());
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — cross-steal notification
// ---------------------------------------------------------------------------

// Multiple workers call help_while_waiting simultaneously (as happens in BRWT
// when workers execute slice_rows tasks that call get_nonzero_rows).
// Each enqueues a chunk task and waits for it. Cross-stealing occurs when
// worker A executes worker B's chunk. Worker B discovers its future is ready
// on the next wait_for timeout (100μs).
//
// The spin on all_done prevents workers from returning to the worker loop,
// so only the timeout (not the worker loop's notify) resolves the cross-steal.
TEST(ThreadPool, HelpWhileWaitingCrossStealNotification) {
    for (int iter = 0; iter < 50; ++iter) {
        const int N = 4;
        ThreadPool pool(N);

        std::atomic<int> chunks_enqueued{0};
        std::atomic<int> finished_hww{0};
        std::atomic<bool> all_done{false};

        for (int i = 0; i < N; ++i) {
            pool.force_enqueue([&]() {
                // Each worker enqueues one chunk task.
                auto f = pool.force_enqueue_front([]() { return 1; });

                // Barrier: wait until all workers have enqueued their chunks
                // before anyone starts stealing. This maximizes cross-steal
                // probability — all N chunks sit in the queue when HWW begins.
                chunks_enqueued.fetch_add(1, std::memory_order_release);
                while (chunks_enqueued.load(std::memory_order_acquire) < N) {
                    std::this_thread::yield();
                }

                pool.help_while_waiting(f).get();
                finished_hww.fetch_add(1, std::memory_order_release);

                // Spin here to prevent returning to the worker loop.
                // Without this, the worker loop's help_condition.notify_one()
                // would rescue stuck HWW callers, masking the bug.
                while (!all_done.load(std::memory_order_acquire)) {
                    std::this_thread::yield();
                }
            });
        }

        // Wait for all workers to finish help_while_waiting (with timeout).
        auto deadline = std::chrono::steady_clock::now()
                        + std::chrono::seconds(5);
        while (finished_hww.load(std::memory_order_acquire) < N) {
            ASSERT_LT(std::chrono::steady_clock::now(), deadline)
                << "iter " << iter << ": only "
                << finished_hww.load(std::memory_order_acquire) << "/" << N
                << " workers completed help_while_waiting (deadlock from "
                   "missing cross-steal notification)";
            std::this_thread::yield();
        }

        all_done.store(true, std::memory_order_release);
        pool.join();
    }
}

// ---------------------------------------------------------------------------
//  ThreadPool: remove_waiting_tasks unblocks enqueue
// ---------------------------------------------------------------------------

// If a thread is blocked in enqueue() because the queue is full,
// remove_waiting_tasks() should unblock it by notifying full_condition.
TEST(ThreadPool, RemoveWaitingTasksUnblocksEnqueue) {
    ThreadPool pool(1, 2);  // 1 worker, max 2 tasks in queue

    std::promise<void> worker_started;
    std::promise<void> unblock;

    pool.enqueue([&]() {
        worker_started.set_value();
        unblock.get_future().wait();
    });
    worker_started.get_future().wait();

    // Fill the queue to max_num_tasks.
    pool.force_enqueue([]() {});
    pool.force_enqueue([]() {});

    // A regular enqueue on another thread will block (queue full).
    std::atomic<bool> enqueue_completed{false};
    std::thread enqueuer([&]() {
        pool.enqueue([]() {});
        enqueue_completed.store(true, std::memory_order_release);
    });

    // Give the enqueuer time to block.
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    EXPECT_FALSE(enqueue_completed.load(std::memory_order_acquire));

    // Clear the queue — should unblock the enqueuer via full_condition.
    pool.remove_waiting_tasks();

    auto deadline = std::chrono::steady_clock::now()
                    + std::chrono::seconds(5);
    while (!enqueue_completed.load(std::memory_order_acquire)) {
        ASSERT_LT(std::chrono::steady_clock::now(), deadline)
            << "enqueue() still blocked after remove_waiting_tasks()";
        std::this_thread::yield();
    }

    enqueuer.join();
    unblock.set_value();
    pool.join();
}

// ---------------------------------------------------------------------------
//  ThreadPool: help_while_waiting — latency benchmark
// ---------------------------------------------------------------------------

// Measures the latency of help_while_waiting when the future is completed by
// a worker (not stolen). The worker loop's help_condition.notify_one() wakes
// HWW promptly; without it, HWW would rely on the 100μs poll timeout.
TEST(ThreadPool, HelpWhileWaitingLatencyBenchmark) {
    const int num_iters = 200;
    ThreadPool pool(1);

    // A releaser thread opens the gate after 200μs, ensuring HWW is already
    // blocking when the worker completes. The "above 200μs gate" number in
    // the output isolates wake-up overhead from the gate delay.
    long long total_us = 0;
    for (int i = 0; i < num_iters; ++i) {
        std::promise<void> gate;
        auto gate_future = gate.get_future().share();

        auto future = pool.enqueue([gate_future]() {
            gate_future.wait();  // block until released
            return 1;
        });

        // Give the worker time to pick up the task and block on the gate.
        std::this_thread::sleep_for(std::chrono::microseconds(100));

        // Release the worker from a separate thread after 200μs.
        // By then, HWW is definitely in wait_for.
        std::thread releaser([&gate]() {
            std::this_thread::sleep_for(std::chrono::microseconds(200));
            gate.set_value();
        });

        auto t0 = std::chrono::steady_clock::now();
        pool.help_while_waiting(future).get();
        auto t1 = std::chrono::steady_clock::now();

        total_us += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        releaser.join();
    }

    fprintf(stderr, "[BENCH] %d iterations: %lld μs total, %.1f μs/iter "
                    "(%.1f μs above 200μs gate)\n",
                    num_iters, total_us, (double)total_us / num_iters,
                    (double)total_us / num_iters - 200.0);
}

TEST(AsyncActivity, RunUniqueOnly) {
    AsyncActivity async;

    ASSERT_EQ(1, async.run_unique([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_unique([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_unique([](char a, char b, char c) {
                                         return std::string(1, a) + b + c;
                                      }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunAsyncOnly) {
    AsyncActivity async;

    ASSERT_EQ(1, async.run_async([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_async([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_async([](char a, char b, char c) {
                                        return std::string(1, a) + b + c;
                                     }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunBoth) {
    AsyncActivity async;

    ASSERT_EQ(1, async.run_async([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_async([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_async([](char a, char b, char c) {
                                        return std::string(1, a) + b + c;
                                     }, 'a', 'b', 'c'));

    ASSERT_EQ(1, async.run_unique([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_unique([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_unique([](char a, char b, char c) {
                                         return std::string(1, a) + b + c;
                                      }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunBothParallel) {
    AsyncActivity async;

    std::string result;

    std::thread t1([&](std::string *result) {
                    async.run_async(
                        [](std::string *result_) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(200));
                            result_->push_back('a');
                            return true;
                        }, result);
                    }, &result);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    async.run_unique([](std::string *result_) { result_->push_back('b'); }, &result);
    ASSERT_EQ(2u, result.size()) << result;
    EXPECT_EQ('b', result.back()) << result;

    std::thread t2([&](std::string *result) {
                    async.run_async(
                        [](std::string *result_) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(200));
                            result_->push_back('a');
                            return true;
                        }, result);
                    }, &result);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    async.run_unique([](std::string *result_) { result_->push_back('c'); }, &result);
    ASSERT_EQ(4u, result.size()) << result;
    EXPECT_EQ('c', result.back()) << result;

    t1.join();
    t2.join();
}

TEST(Misc, SplitString) {
    EXPECT_EQ(std::vector<std::string>(),
              utils::split_string("", ""));

    EXPECT_EQ(std::vector<std::string>({ "123" }),
              utils::split_string("123", ""));

    EXPECT_EQ(std::vector<std::string>({ "23" }),
              utils::split_string("123", "1"));

    EXPECT_EQ(std::vector<std::string>({ "1", "3" }),
              utils::split_string("123", "2"));

    EXPECT_EQ(std::vector<std::string>({ "1", "34", "67" }),
              utils::split_string("12342672", "2"));

    EXPECT_EQ(std::vector<std::string>({ "1", "34", "72" }),
              utils::split_string("126342672", "26"));

    EXPECT_EQ(std::vector<std::string>({ "126342672" }),
              utils::split_string("126342672", ""));
}

TEST(Misc, JoinStrings) {
    EXPECT_EQ("", utils::join_strings({ "" }, ""));
    EXPECT_EQ("23", utils::join_strings({ "23" }, ":"));
    EXPECT_EQ("1:3", utils::join_strings(
                std::vector<std::string>({ "1", "3" }), ":"));
    EXPECT_EQ("1:34:67", utils::join_strings(
                std::vector<std::string>({ "1", "34", "67" }), ":"));
    EXPECT_EQ("1::34::72", utils::join_strings(
                std::vector<std::string>({ "1", "34", "72" }), "::"));
    EXPECT_EQ("13472", utils::join_strings(
                std::vector<std::string>({ "1", "34", "72" }), ""));

    EXPECT_EQ(std::vector<std::string>({ "123" }),
              utils::split_string("123", ""));

    EXPECT_EQ(std::vector<std::string>({ "23" }),
              utils::split_string("123", "1"));

    EXPECT_EQ(".1", utils::join_strings({ "", "1" }, "."));
    EXPECT_EQ("1.", utils::join_strings({ "1", "" }, "."));
    EXPECT_EQ(".1", utils::join_strings({ "", "1" }, ".", false));
    EXPECT_EQ("1.", utils::join_strings({ "1", "" }, ".", false));
    EXPECT_EQ("1", utils::join_strings({ "", "1" }, ".", true));
    EXPECT_EQ("1", utils::join_strings({ "1", "" }, ".", true));

    EXPECT_EQ("..", utils::join_strings({ "", "", "" }, "."));
    EXPECT_EQ("..3", utils::join_strings({ "", "", "3" }, "."));
    EXPECT_EQ("1..", utils::join_strings({ "1", "", "" }, "."));
    EXPECT_EQ(".2.", utils::join_strings({ "", "2", "" }, "."));
    EXPECT_EQ("1.2.", utils::join_strings({ "1", "2", "" }, "."));
    EXPECT_EQ("1..3", utils::join_strings({ "1", "", "3" }, "."));
    EXPECT_EQ(".2.3", utils::join_strings({ "", "2", "3" }, "."));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, "."));

    EXPECT_EQ("..", utils::join_strings({ "", "", "" }, ".", false));
    EXPECT_EQ("..3", utils::join_strings({ "", "", "3" }, ".", false));
    EXPECT_EQ("1..", utils::join_strings({ "1", "", "" }, ".", false));
    EXPECT_EQ(".2.", utils::join_strings({ "", "2", "" }, ".", false));
    EXPECT_EQ("1.2.", utils::join_strings({ "1", "2", "" }, ".", false));
    EXPECT_EQ("1..3", utils::join_strings({ "1", "", "3" }, ".", false));
    EXPECT_EQ(".2.3", utils::join_strings({ "", "2", "3" }, ".", false));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, ".", false));

    EXPECT_EQ("", utils::join_strings({ "", "", "" }, ".", true));
    EXPECT_EQ("3", utils::join_strings({ "", "", "3" }, ".", true));
    EXPECT_EQ("1", utils::join_strings({ "1", "", "" }, ".", true));
    EXPECT_EQ("2", utils::join_strings({ "", "2", "" }, ".", true));
    EXPECT_EQ("1.2", utils::join_strings({ "1", "2", "" }, ".", true));
    EXPECT_EQ("1.3", utils::join_strings({ "1", "", "3" }, ".", true));
    EXPECT_EQ("2.3", utils::join_strings({ "", "2", "3" }, ".", true));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, ".", true));
}

TEST(Misc, RemoveSuffix) {
    EXPECT_EQ("", utils::remove_suffix("", ""));
    EXPECT_EQ("", utils::remove_suffix("", "A"));
    EXPECT_EQ("", utils::remove_suffix("", "1AC"));
    EXPECT_EQ("", utils::remove_suffix("1", "1"));
    EXPECT_EQ("1.", utils::remove_suffix("1.1AC", "1AC"));
    EXPECT_EQ("1", utils::remove_suffix("1.1AC", ".1AC"));
    EXPECT_EQ("acb", utils::remove_suffix("acb.cde", ".cde"));
    EXPECT_EQ("acb.cde", utils::remove_suffix("acb.cde.efg", ".efg"));
    EXPECT_EQ("acb", utils::remove_suffix("acb.cde.efg", ".cde.efg"));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", ""));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", ".cde"));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", "Aacb.cde.efg"));
}

TEST(Misc, CountIntersectionEmpty) {
    {
        std::vector<int> first;
        std::vector<int> second;
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second;
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first;
        std::vector<int> second { 1, 2, 3, 4, 5 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { -5, -4, -3, -2, -1 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { -1000 };
        std::vector<int> second { 1000 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
}

TEST(Misc, CountIntersection) {
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { 1, 2, 3, 4, 5 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { 0, 1, 2, 3, 4, 5, 6 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5, 6 };
        std::vector<int> second { 0, 1, 2, 3, 4, 5 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 0, 1, 2, 3, 4, 5 };
        std::vector<int> second { -5, -4, -3, -2, -1, 0 };
        EXPECT_EQ(1u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 0, 1000000 };
        std::vector<int> second { -1000000, 0, 1000000 };
        EXPECT_EQ(2u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
}

TEST(Vector, ReserveInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.reserve(1llu << 59), std::bad_alloc);
}

TEST(Vector, ResizeInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.resize(1llu << 59), std::bad_alloc);
}

TEST(Deque, ResizeInfinityCheckThrow) {
    std::deque<int> array;
    EXPECT_THROW(array.resize(1llu << 60), std::bad_alloc);
}

TEST(Misc, get_quantile) {
    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.75));

    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.75));

    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.75));
}


template <typename T>
std::vector<T> insert_reference_impl(std::vector<T> vector,
                                     const std::vector<uint64_t> &new_pos,
                                     T value) {
    assert(std::is_sorted(new_pos.begin(), new_pos.end()));
    for (auto i : new_pos) {
        vector.insert(vector.begin() + i, value);
    }
    return vector;
}

bitmap_vector to_bitmap(size_t size, const std::vector<uint64_t> &set_bits_pos) {
    assert(std::is_sorted(set_bits_pos.begin(), set_bits_pos.end()));
    sdsl::bit_vector mask(size, false);
    for (auto i : set_bits_pos) {
        mask[i] = true;
    }
    return bitmap_vector(std::move(mask));
}

TEST(Misc, insert_zeros_to_empty) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, 5, }) }) {
        std::vector<int> first;
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, insert_zeros_to_all_zeros) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {
        std::vector<int> first(50, 0);
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, insert_zeros_to_all_ones) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {
        std::vector<int> first(50, 1);
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, drag_and_mark_segments) {
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 0));
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 1));
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 100));

    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 0));
    EXPECT_EQ(std::vector<bool>({ 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 1));
    EXPECT_EQ(std::vector<bool>({ 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 100));

    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 0));
    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 1));
    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 100));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 4));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 2));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 1, 10, 4, 6 }), 10, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 1, 10, 4, 6 }), 6, 50));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 0, 0, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 1, 10, 6, 4 }), 6, 50));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 0, 0, 1, 0, 0, 0, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 0, 1, 1, 0, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 1, 1, 1, 1, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 3));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 100));
}


TEST(int_vector_buffer, small_buffer) {
    const size_t size = 10'000;
    const size_t width = 38;
    const size_t buffer_size = 100;
    const size_t num_threads = 4;

    sdsl::int_vector<> vector(size, 0, width);
    sdsl::util::set_random_bits(vector, 42);

    // initialize buffers
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < num_threads; ++i) {
        sdsl::int_vector_buffer<> buf(test_dump_basename + std::to_string(i),
                                      std::ios::out, buffer_size, width);
        for (auto it = vector.begin(); it != vector.end(); ++it) {
            buf.push_back(*it);
        }
    }

    // update buffers
    sdsl::int_vector<> final(size, 0, width);
    sdsl::util::set_random_bits(final, 12);

    std::vector<sdsl::int_vector_buffer<>> bufs(num_threads);

    std::vector<size_t> indexes(size);
    std::iota(indexes.begin(), indexes.end(), 0);
    std::mt19937 g(13);
    std::shuffle(indexes.begin(), indexes.end(), g);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < num_threads; ++i) {
        bufs[i] = sdsl::int_vector_buffer<>(test_dump_basename + std::to_string(i),
                                            std::ios::in|std::ios::out, buffer_size, width);
        for (size_t j : indexes) {
            bufs[i][j] -= vector[j];
            bufs[i][j] += final[j];
        }
    }

    for (auto &buf : bufs) {
        ASSERT_EQ(size, buf.size());
        for (size_t j = 0; j < final.size(); ++j) {
            ASSERT_EQ(final[j], buf[j]);
        }
        buf.close(true);
    }
}

TEST(int_vector_buffer, seek_to_front_65) {
    const size_t width = 1;
    const size_t buffer_size = 8; // buffer size 8 bytes
    const std::string &filename = test_dump_basename + "_ivb";

    // initialize buffer
    {
        sdsl::int_vector_buffer<> buf(filename, std::ios::out, buffer_size, width);
    }

    sdsl::int_vector_buffer<> buf(filename, std::ios::in | std::ios::out, buffer_size, width);
    for (size_t i = 0; i < 65; ++i) {
        buf.push_back(1);
    }

    ASSERT_EQ(1, buf[0]);
    ASSERT_EQ(1, buf[64]);

    buf.close(true);
}

TEST(int_vector_buffer, seek_to_front) {
    const size_t size = 1'000'000;
    const size_t width = 1;
    const size_t buffer_size = 10'000;
    const std::string &filename = test_dump_basename + "_ivb";

    sdsl::int_vector<> vector(size, 1, width);

    // initialize buffer
    {
        sdsl::int_vector_buffer<> buf(filename, std::ios::out, buffer_size, width);
    }

    sdsl::int_vector_buffer<> buf(filename, std::ios::in | std::ios::out, buffer_size, width);
    for (size_t i = 0; i < size; ++i) {
        buf.push_back(1);
    }

    for (size_t i = 0; i < size; ++i) {
        ASSERT_EQ(1, buf[i]);
    }

    buf.close(true);
}

TEST(int_vector_buffer, change_and_seek_to_front) {
    const size_t size = 1'000'000;
    const size_t width = 1;
    const size_t buffer_size = 10'000;
    const std::string &filename = test_dump_basename + "_ivb";

    sdsl::int_vector<> vector(size, 1, width);

    // initialize buffer
    {
        sdsl::int_vector_buffer<> buf(filename, std::ios::out, buffer_size, width);
        for (size_t i = 0; i < size; ++i) {
            buf.push_back(0);
        }
    }

    sdsl::int_vector_buffer<> buf(filename, std::ios::in | std::ios::out, buffer_size, width);
    for (size_t i = 0; i < size; ++i) {
        buf[i] = 1;
    }

    for (size_t i = 0; i < size; ++i) {
        ASSERT_EQ(1, buf[i]);
    }

    buf.close(true);
}

TEST(smooth_vector, empty) {
    std::vector<uint32_t> v;
    utils::smooth_vector(1, &v);
    EXPECT_THAT(v, testing::ElementsAre());
    utils::smooth_vector(2, &v);
    EXPECT_THAT(v, testing::ElementsAre());
    utils::smooth_vector(100, &v);
    EXPECT_THAT(v, testing::ElementsAre());
}

TEST(smooth_vector, window_size_1) {
    std::vector<uint32_t> v = { 1, 2, 3 };
    utils::smooth_vector(0, &v);
    EXPECT_THAT(v, testing::ElementsAre(1, 2, 3));
    utils::smooth_vector(1, &v);
    EXPECT_THAT(v, testing::ElementsAre(1, 2, 3));
}

TEST(smooth_vector, increasing) {
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(0, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, 2, 3, 4, 5, 6, 7, 8, 9));
        utils::smooth_vector(1, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, 2, 3, 4, 5, 6, 7, 8, 9));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(2, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, 2, 3, 4, 5, 6, 7, 8, 9));
    }
    {
        std::vector<uint32_t> v = { 1, 3, 5, 7, 9, 11, 13, 15, 17 };
        utils::smooth_vector(2, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, 2, 4, 6, 8, 10, 12, 14, 16));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(v.size(), &v);
        EXPECT_THAT(v, testing::ElementsAre(3, 4, 4, 5, 5, 6, 6, 7, 7));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(v.size() + 1, &v);
        EXPECT_THAT(v, testing::ElementsAre(3, 4, 4, 5, 5, 5, 6, 6, 7));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(v.size() * 2 - 2, &v);
        EXPECT_THAT(v, testing::ElementsAre(5, 5, 5, 5, 5, 5, 5, 5, 5));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(v.size() * 2 - 1, &v);
        EXPECT_THAT(v, testing::ElementsAre(5, 5, 5, 5, 5, 5, 5, 5, 5));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(v.size() * 2, &v);
        EXPECT_THAT(v, testing::ElementsAre(5, 5, 5, 5, 5, 5, 5, 5, 5));
    }
    {
        std::vector<uint32_t> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        utils::smooth_vector(100, &v);
        EXPECT_THAT(v, testing::ElementsAre(5, 5, 5, 5, 5, 5, 5, 5, 5));
    }
}

TEST(smooth_vector, zigzag) {
    {
        std::vector<int> v = { 1, -1, 1, -1, 1, -1, 1, -1, 1 };
        utils::smooth_vector(1, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, -1, 1, -1, 1, -1, 1, -1, 1));
        utils::smooth_vector(2, &v);
        EXPECT_THAT(v, testing::ElementsAre(1, 0, 0, 0, 0, 0, 0, 0, 0));
    }
    {
        std::vector<int> v = { 1, -1, 1, -1, 1, -1, 1, -1, 1 };
        utils::smooth_vector(3, &v);
        EXPECT_THAT(v, testing::ElementsAre(0, 0, 0, 0, 0, 0, 0, 0, 0));
    }
}

TEST(smooth_vector, check_asserts_odd) {
    std::vector<uint32_t> v(15, 0);
    for (size_t ws = 0; ws < 100; ++ws) {
        utils::smooth_vector(ws, &v);
    }
}

TEST(smooth_vector, check_asserts_even) {
    std::vector<uint32_t> v(16, 0);
    for (size_t ws = 0; ws < 100; ++ws) {
        utils::smooth_vector(ws, &v);
    }
}

// ============================================================
// ipow tests
// ============================================================

TEST(ipow, ZeroExpReturnsOne) {
    EXPECT_EQ(1u, utils::ipow(0, 0));
    EXPECT_EQ(1u, utils::ipow(1, 0));
    EXPECT_EQ(1u, utils::ipow(2, 0));
    EXPECT_EQ(1u, utils::ipow(100, 0));
    EXPECT_EQ(1u, utils::ipow(std::numeric_limits<uint64_t>::max(), 0));
}

TEST(ipow, ZeroBaseReturnsZero) {
    EXPECT_EQ(0u, utils::ipow(0, 1));
    EXPECT_EQ(0u, utils::ipow(0, 2));
    EXPECT_EQ(0u, utils::ipow(0, 100));
}

TEST(ipow, OneBaseAlwaysOne) {
    EXPECT_EQ(1u, utils::ipow(1, 0));
    EXPECT_EQ(1u, utils::ipow(1, 1));
    EXPECT_EQ(1u, utils::ipow(1, 100));
    EXPECT_EQ(1u, utils::ipow(1, std::numeric_limits<unsigned int>::max()));
}

TEST(ipow, ExpOneReturnsBase) {
    EXPECT_EQ(2u, utils::ipow(2, 1));
    EXPECT_EQ(5u, utils::ipow(5, 1));
    EXPECT_EQ(123456789u, utils::ipow(123456789, 1));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              utils::ipow(std::numeric_limits<uint64_t>::max(), 1));
}

TEST(ipow, SmallPowersOfTwo) {
    for (unsigned int e = 0; e < 64; ++e) {
        EXPECT_EQ(uint64_t(1) << e, utils::ipow(2, e)) << "2^" << e;
    }
}

TEST(ipow, PowersOfTen) {
    uint64_t expected = 1;
    for (unsigned int e = 0; e <= 19; ++e) {
        EXPECT_EQ(expected, utils::ipow(10, e)) << "10^" << e;
        if (e < 19)
            expected *= 10;
    }
}

TEST(ipow, SmallBaseSmallExp) {
    // brute-force check for small values
    for (uint64_t base = 2; base <= 20; ++base) {
        uint64_t result = 1;
        for (unsigned int exp = 0; exp <= 15; ++exp) {
            EXPECT_EQ(result, utils::ipow(base, exp)) << base << "^" << exp;
            if (result > std::numeric_limits<uint64_t>::max() / base)
                break;
            result *= base;
        }
    }
}

TEST(ipow, MatchesBOSSAlphabetSizes) {
    // DNA: sigma = alph_size - 1 = 5
    EXPECT_EQ(1u,   utils::ipow(5, 0));
    EXPECT_EQ(5u,   utils::ipow(5, 1));
    EXPECT_EQ(25u,  utils::ipow(5, 2));
    EXPECT_EQ(125u, utils::ipow(5, 3));
    EXPECT_EQ(uint64_t(5) * 5 * 5 * 5 * 5 * 5 * 5, utils::ipow(5, 7));

    // Protein: sigma = 21
    EXPECT_EQ(21u,    utils::ipow(21, 1));
    EXPECT_EQ(441u,   utils::ipow(21, 2));
    EXPECT_EQ(9261u,  utils::ipow(21, 3));
}

TEST(ipow, LargestNonOverflowing) {
    // 2^63 is the largest power of 2 that fits in uint64_t
    EXPECT_EQ(uint64_t(1) << 63, utils::ipow(2, 63));

    // 3^40 = 12157665459056928801, fits in uint64_t
    EXPECT_EQ(UINT64_C(12157665459056928801), utils::ipow(3, 40));

    // 5^27 = 7450580596923828125, fits
    EXPECT_EQ(UINT64_C(7450580596923828125), utils::ipow(5, 27));

    // 10^19 = 10000000000000000000, fits
    EXPECT_EQ(UINT64_C(10000000000000000000), utils::ipow(10, 19));

    // UINT64_MAX^1 fits
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              utils::ipow(std::numeric_limits<uint64_t>::max(), 1));
}

TEST(ipow, OverflowThrows) {
    // 2^64 overflows
    EXPECT_THROW(utils::ipow(2, 64), std::overflow_error);

    // 3^41 overflows (3^40 * 3 > UINT64_MAX)
    EXPECT_THROW(utils::ipow(3, 41), std::overflow_error);

    // 5^28 overflows
    EXPECT_THROW(utils::ipow(5, 28), std::overflow_error);

    // 10^20 overflows
    EXPECT_THROW(utils::ipow(10, 20), std::overflow_error);

    // large base with exp=2
    EXPECT_THROW(utils::ipow(uint64_t(1) << 32, 2), std::overflow_error);

    // UINT64_MAX^2 overflows
    EXPECT_THROW(utils::ipow(std::numeric_limits<uint64_t>::max(), 2),
                 std::overflow_error);

    // large base, large exp
    EXPECT_THROW(utils::ipow(1000000, 4), std::overflow_error);
}

TEST(ipow, OverflowBoundary) {
    // For base=2: 2^63 fits, 2^64 does not
    EXPECT_NO_THROW(utils::ipow(2, 63));
    EXPECT_THROW(utils::ipow(2, 64), std::overflow_error);

    // For base=10: 10^19 fits, 10^20 does not
    EXPECT_NO_THROW(utils::ipow(10, 19));
    EXPECT_THROW(utils::ipow(10, 20), std::overflow_error);

    // UINT64_MAX = 18446744073709551615
    // floor(UINT64_MAX^(1/2)) = 4294967295 = 2^32 - 1
    uint64_t sqrt_max = (uint64_t(1) << 32) - 1;
    EXPECT_NO_THROW(utils::ipow(sqrt_max, 2));
    EXPECT_THROW(utils::ipow(sqrt_max + 1, 2), std::overflow_error);
}

} // namespace
