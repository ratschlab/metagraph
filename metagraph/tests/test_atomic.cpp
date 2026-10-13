#include <atomic>
#include <thread>
#include <vector>

#include "gtest/gtest.h"
#include "common/vectors/vector_algorithm.hpp"

namespace {

TEST(BitVector, ConcurrentFetchAndModify) {
    for (auto order : { std::memory_order_relaxed, std::memory_order_acquire,
                        std::memory_order_release, std::memory_order_acq_rel,
                        std::memory_order_seq_cst }) {
        for (bool set : { false, true }) {
            sdsl::bit_vector bits(4096, !set);
            std::atomic<bool> start(false);
            std::atomic<uint64_t> changed(0);
            std::vector<std::thread> workers;
            for (size_t t = 0; t < 4; ++t) {
                workers.emplace_back([&]() {
                    while (!start.load(std::memory_order_acquire))
                        std::this_thread::yield();
                    uint64_t local_changed = 0;
                    for (uint64_t i = 0; i < bits.size(); ++i)
                        local_changed += set ? !fetch_and_set_bit(bits.data(), i, true, order)
                                             : fetch_and_unset_bit(bits.data(), i, true, order);
                    changed.fetch_add(local_changed, std::memory_order_relaxed);
                });
            }
            start.store(true, std::memory_order_release);
            for (auto &worker : workers)
                worker.join();

            EXPECT_EQ(bits.size(), changed.load());
            for (uint64_t i = 0; i < bits.size(); ++i)
                EXPECT_EQ(set, fetch_bit(bits.data(), i, true, std::memory_order_relaxed));
        }
    }
}

} // namespace
