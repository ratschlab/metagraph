#include "gtest/gtest.h"

#include <vector>
#include <deque>

#include "common/ring_buffer.hpp"


namespace {

void test_buffer_back(RingBuffer<int> &buffer, size_t size) {
    EXPECT_EQ(size, buffer.capacity());
    std::deque<int> reference;

    for (size_t i = 0; i < 2 * buffer.capacity(); ++i) {
        buffer.push_back(i);
        reference.push_back(i);

        EXPECT_EQ(size, buffer.capacity());

        if (reference.size() == size)
            reference.pop_front();

        EXPECT_EQ(reference.back(), buffer.back());

        if (reference.size() == size) {
            EXPECT_EQ(reference.front(), buffer.front());
        }
    }
}

void test_buffer_front(RingBuffer<int> &buffer, size_t size) {
    EXPECT_EQ(size, buffer.capacity());
    std::deque<int> reference;

    for (size_t i = 0; i < 2 * buffer.capacity(); ++i) {
        buffer.push_front(i);
        reference.push_front(i);

        EXPECT_EQ(size, buffer.capacity());

        if (reference.size() == size)
            reference.pop_back();

        EXPECT_EQ(reference.front(), buffer.front());

        if (reference.size() == size) {
            EXPECT_EQ(reference.back(), buffer.back());
        }
    }
}

void test_buffer_front_and_back(RingBuffer<int> &buffer, size_t size) {
    EXPECT_EQ(size, buffer.capacity());
    std::deque<int> reference;

    for (size_t i = 0; i < 2 * buffer.capacity(); ++i) {
        if (i & 1) {
            buffer.push_front(i);
            reference.push_front(i);

            if (reference.size() == size)
                reference.pop_back();
        } else {
            buffer.push_back(i);
            reference.push_back(i);

            if (reference.size() == size)
                reference.pop_front();
        }

        EXPECT_EQ(size, buffer.capacity());

        if (reference.size() == size) {
            EXPECT_EQ(reference.front(), buffer.front());
            EXPECT_EQ(reference.back(), buffer.back());
        }
    }
}

TEST(RingBuffer, empty) {
    RingBuffer<int> buffer(0);
    EXPECT_EQ(0u, buffer.capacity());

    buffer.reset();

#ifndef NDEBUG
    ASSERT_DEATH(buffer.push_back(1), "");
    ASSERT_DEATH(buffer.push_front(1), "");
#endif
}

TEST(RingBuffer, PushBack) {
    for (size_t i = 1; i < 129; ++i) {
        RingBuffer<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_back(buffer, i);
        buffer.reset();
        test_buffer_back(buffer, i);
    }
}

TEST(RingBuffer, PushFront) {
    for (size_t i = 1; i < 129; ++i) {
        RingBuffer<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_front(buffer, i);
        buffer.reset();
        test_buffer_front(buffer, i);
    }
}

TEST(RingBuffer, PushFrontAndBack) {
    for (size_t i = 1; i < 129; ++i) {
        RingBuffer<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_front_and_back(buffer, i);
        buffer.reset();
        test_buffer_front_and_back(buffer, i);
    }
}

TEST(RingBuffer, PushAllCases) {
    for (size_t i = 1; i < 129; ++i) {
        RingBuffer<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_front_and_back(buffer, i);
        buffer.reset();
        test_buffer_front(buffer, i);
        buffer.reset();
        test_buffer_back(buffer, i);
        buffer.reset();
        test_buffer_front_and_back(buffer, i);
    }
}

} // namespace
