#include "gtest/gtest.h"

#include <vector>
#include <deque>

#include "rolling_window.hpp"


void test_buffer_back(RollingWindow<int> &buffer, size_t size) {
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

void test_buffer_front(RollingWindow<int> &buffer, size_t size) {
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

void test_buffer_front_and_back(RollingWindow<int> &buffer, size_t size) {
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

TEST(RollingWindow, empty) {
    RollingWindow<int> buffer(0);
    EXPECT_EQ(0u, buffer.capacity());

    buffer.reset();

    ASSERT_DEATH(buffer.push_back(1), "");
    ASSERT_DEATH(buffer.push_front(1), "");
}

TEST(RollingWindow, PushBack) {
    for (size_t i = 1; i < 129; ++i) {
        RollingWindow<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_back(buffer, i);
        buffer.reset();
        test_buffer_back(buffer, i);
    }
}

TEST(RollingWindow, PushFront) {
    for (size_t i = 1; i < 129; ++i) {
        RollingWindow<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_front(buffer, i);
        buffer.reset();
        test_buffer_front(buffer, i);
    }
}

TEST(RollingWindow, PushFrontAndBack) {
    for (size_t i = 1; i < 129; ++i) {
        RollingWindow<int> buffer(i);
        EXPECT_EQ(i, buffer.capacity());

        test_buffer_front_and_back(buffer, i);
        buffer.reset();
        test_buffer_front_and_back(buffer, i);
    }
}

TEST(RollingWindow, PushAllCases) {
    for (size_t i = 1; i < 129; ++i) {
        RollingWindow<int> buffer(i);
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
