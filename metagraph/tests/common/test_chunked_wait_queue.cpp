#include "common/chunked_wait_queue.hpp"

#include "gtest/gtest.h"

#include <chrono>
#include <functional>
#include <thread>

#include <cstdint>

namespace {

using namespace mg;
using common::ChunkedWaitQueue;

TEST(WaitQueue, Empty) {
    ChunkedWaitQueue<int32_t> under_test(20, 1);
    EXPECT_TRUE(under_test.empty());
    EXPECT_FALSE(under_test.full());
}

TEST(WaitQueue, PushPop) {
    ChunkedWaitQueue<int32_t> under_test(3, 1);
    under_test.set_out_file("/tmp/out");
    under_test.push_front(1);
    EXPECT_FALSE(under_test.full());
    under_test.push_front(2);
    under_test.push_front(3);
    EXPECT_TRUE(under_test.full());
    ChunkedWaitQueue<int32_t>::Iterator &iterator = under_test.begin();
    EXPECT_EQ(1, *iterator);
    EXPECT_EQ(2, *(++iterator));
    // the chunk wasn't yet cleaned up, so the queue should be full
    EXPECT_TRUE(under_test.full());

    EXPECT_EQ(3, *(++iterator));
    // at this point the first chunk was cleaned up, so the queue is not full any longer
    EXPECT_FALSE(under_test.full());

    under_test.push_front(4);
    EXPECT_EQ(4, *(++iterator));

    under_test.shutdown();
    ++iterator;
    EXPECT_TRUE(iterator == under_test.end());
    std::filesystem::remove("/tmp/out");
}

TEST(WaitQueue, Shutdown) {
    ChunkedWaitQueue<std::string> under_test(20, 2);
    under_test.shutdown();
    std::string v;
    EXPECT_TRUE(under_test.begin() == under_test.end());
}

void writeReadWaitQueue(uint32_t delay_read_ms, uint32_t delay_write_ms) {
    ChunkedWaitQueue<int32_t> under_test(20, 2);
    under_test.set_out_file("/tmp/out");
    struct Receiver {
        ChunkedWaitQueue<int32_t> *const under_test;
        uint32_t delay_read_ms;
        std::vector<int32_t> pop_result;

        void run() {
            for (auto &it = under_test->begin(); it != under_test->end(); ++it) {
                if (delay_read_ms > 0) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(delay_read_ms));
                }
                pop_result.push_back(*it);
            }
        }
    };

    Receiver receiver{&under_test, delay_read_ms, {}};

    // start 1 receiver thread eager to consume data
    std::thread receiverThread(std::bind(&Receiver::run, std::ref(receiver)));

    // start feeding the receiver with some data
    std::thread senderThread([&under_test, delay_write_ms] {
        for (uint32_t i = 0; i < 100; ++i) {
            if (delay_write_ms > 0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(delay_write_ms));
            }
            under_test.push_front(i);
        }
        under_test.shutdown();
    });

    receiverThread.join();
    senderThread.join();

    // make sure that the 1 receiver received all of the elements that were written
    for (int32_t i = 0; i < 100; ++i) {
        EXPECT_EQ(i, receiver.pop_result[i]);
    }
    std::filesystem::remove("/tmp/out");
}

TEST(WaitQueue, OneWriterOneReader) {
    writeReadWaitQueue(0, 0);
}

TEST(WaitQueue, OneWriterOneSlowReader) {
    writeReadWaitQueue(1, 0);
}

TEST(WaitQueue, OneSlowWriterOneReader) {
    writeReadWaitQueue(0, 1);
}

}  // namespace
