#include "gtest/gtest.h"

#include "int_vector.hpp"
#include "test_helpers.hpp"


TEST(IntVector, call_nonzeros_all_zeros) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(std::lcm(64, w) + w * 2, 0, w);
        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                call_nonzeros(vector, begin, end, [](auto, auto) { EXPECT_TRUE(false); });
            }
        }
    }
}

TEST(IntVector, DISABLED_call_nonzeros_all_set_LONG_TEST) {
    for (size_t w = 1; w <= 64; ++w) {
        TEST_COUT << w;
        sdsl::int_vector<> vector((std::lcm(64, w) + w * 2) * w, 1, 1);
        vector.width(w);
        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                uint64_t it = begin;
                call_nonzeros(vector, begin, end, [&](auto i, auto count) {
                    EXPECT_LE(begin, i);
                    EXPECT_GT(end, i);
                    EXPECT_EQ(it++, i);
                    EXPECT_EQ(w < 64 ? (1llu << w) - 1 : static_cast<size_t>(-1),
                              count);
                    EXPECT_EQ(vector[i], count);
                });
                EXPECT_EQ(it, end);
            }
        }
    }
}


TEST(IntVector, DISABLED_call_nonzeros_all_ones_LONG_TEST) {
    for (size_t w = 1; w <= 64; ++w) {
        TEST_COUT << w;
        sdsl::int_vector<> vector(std::lcm(64, w) + w * 2, 1, w);
        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                uint64_t it = begin;
                call_nonzeros(vector, begin, end, [&](auto i, auto count) {
                    EXPECT_LE(begin, i);
                    EXPECT_GT(end, i);
                    EXPECT_EQ(it++, i);
                    EXPECT_EQ(1u, count);
                    EXPECT_EQ(vector[i], count);
                });
                EXPECT_EQ(it, end);
            }
        }
    }
}

TEST(IntVector, DISABLED_call_nonzeros_mostly_ones_LONG_TEST) {
    for (size_t w = 1; w <= 64; ++w) {
        TEST_COUT << w;
        sdsl::int_vector<> vector(std::lcm(64, w) + w * 2, 1, w);
        *vector.data() = 0;
        *(vector.data() + (vector.capacity() >> 6) - 1) = 0;
        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                uint64_t it = begin;
                call_nonzeros(vector, begin, end, [&](auto i, auto count) {
                    EXPECT_LE(begin, i);
                    EXPECT_GT(end, i);
                    while (vector[it] == 0) {
                        EXPECT_GT(i, it);
                        ++it;
                        ASSERT_GT(end, it);
                    }

                    EXPECT_EQ(1u, count);
                    EXPECT_EQ(vector[i], count);
                    EXPECT_EQ(it++, i);
                });
                for (; it < end; ++it) {
                    EXPECT_EQ(0u, vector[it]);
                }
            }
        }
    }
}

TEST(IntVector, DISABLED_call_nonzeros_sparse_LONG_TEST) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(std::lcm(64, w) + w * 2, 0, w);
        ASSERT_LT(65u, vector.size());
        size_t counter = 0;
        for (size_t i = 0; i < vector.size(); ++i) {
            if (i % 130 == 1)
                vector[i] = (++counter) % w;
        }
        TEST_COUT << w << " " << vector.size() << " " << counter;

        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                uint64_t it = begin;
                call_nonzeros(vector, begin, end, [&](auto i, auto count) {
                    EXPECT_LE(begin, i);
                    EXPECT_GT(end, i);
                    while (vector[it] == 0) {
                        EXPECT_GT(i, it);
                        ++it;
                        ASSERT_GT(end, it);
                    }

                    EXPECT_EQ(vector[i], count);
                    EXPECT_EQ(it++, i);
                });
                for (; it < end; ++it) {
                    EXPECT_EQ(0u, vector[it]);
                }
            }
        }
    }
}

TEST(IntVector, call_nonzeros_sparse) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(std::lcm(64, w) + w * 2, 0, w);
        ASSERT_LT(65u, vector.size());
        size_t counter = 0;
        for (size_t i = 0; i < vector.size(); ++i) {
            if (i % 130 == 1)
                vector[i] = (++counter) % w;
        }

        for (size_t begin = 0; begin < vector.size(); begin += vector.size() / 4) {
            for (size_t end = begin; end < vector.size(); ++end) {
                uint64_t it = begin;
                call_nonzeros(vector, begin, end, [&](auto i, auto count) {
                    EXPECT_LE(begin, i);
                    EXPECT_GT(end, i);
                    while (vector[it] == 0) {
                        EXPECT_GT(i, it);
                        ++it;
                        ASSERT_GT(end, it);
                    }

                    EXPECT_EQ(vector[i], count);
                    EXPECT_EQ(it++, i);
                });
                for (; it < end; ++it) {
                    EXPECT_EQ(0u, vector[it]) << vector.size()
                                              << "\n" << begin
                                              << " " << end
                                              << "\n" << w
                                              << "\n" << vector;
                }
            }
        }
    }
}
