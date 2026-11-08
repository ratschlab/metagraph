#include "gtest/gtest.h"

#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/threads/threading.hpp"
#include "test_helpers.hpp"


namespace {

using namespace mtg;

TEST(IntVector, call_nonzeros_all_zeros) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(300, 0, w);
        for (size_t begin = 0; begin < vector.size(); ++begin) {
            for (size_t end = begin; end < vector.size(); ++end) {
                call_nonzeros(vector, begin, end, [](auto, auto) { EXPECT_TRUE(false); });
            }
        }
    }
}

TEST(IntVector, call_nonzeros_all_set) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(300, 1, 1);
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

TEST(IntVector, call_nonzeros_all_ones) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(150, 1, w);
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

TEST(IntVector, call_nonzeros_mostly_ones) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(150, 1, w);
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

TEST(IntVector, call_nonzeros_sparse) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(200, 0, w);
        ASSERT_LT(65u, vector.size());
        size_t counter = 0;
        for (size_t i = 0; i < vector.size(); ++i) {
            if (i % 130 == 1)
                vector[i] = (++counter) % w;
        }

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

TEST(IntVector, call_nonzeros_sparse_every_4) {
    for (size_t w = 1; w <= 64; ++w) {
        sdsl::int_vector<> vector(1000, 0, w);
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

TEST(bit_vector, count_ones_empty) {
    sdsl::bit_vector v;
    EXPECT_EQ(0u, count_ones(v, 0, v.size()));
}

TEST(bit_vector, count_ones_all_zero) {
    {
        sdsl::bit_vector v(1, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(10, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(99999999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
}

TEST(bit_vector, count_ones_all_ones) {
    {
        sdsl::bit_vector v(1, 1);
        EXPECT_EQ(1u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(10, 1);
        EXPECT_EQ(10u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999, 1);
        EXPECT_EQ(999u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999999, 1);
        EXPECT_EQ(999999u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(99999999, 1);
        EXPECT_EQ(99999999u, count_ones(v, 0, v.size()));
    }
}

TEST(bit_vector, count_ones) {
    {
        sdsl::bit_vector v(1);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(10);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(999999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 100000 <= v.size(); i += 100000) {
            for (size_t j = i; j + 100000 <= v.size(); j += 100000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(99999999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 10000000 <= v.size(); i += 10000000) {
            for (size_t j = i; j + 10000000 <= v.size(); j += 10000000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector v(size);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    for (size_t size = 99999999; size < 9999999 + 200; ++size) {
        sdsl::bit_vector v(size);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 10000000 <= v.size(); i += 10000000) {
            for (size_t j = i; j + 10000000 <= v.size(); j += 10000000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
}

TEST(bit_vector, inner_prod_empty) {
    sdsl::bit_vector first;
    sdsl::bit_vector second;
    EXPECT_EQ(0u, inner_prod(first, second));
}

TEST(bit_vector, inner_prod_all_zero) {
    {
        sdsl::bit_vector first(1, 0);
        sdsl::bit_vector second(1, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10, 0);
        sdsl::bit_vector second(10, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999, 0);
        sdsl::bit_vector second(999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999, 0);
        sdsl::bit_vector second(999999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999, 0);
        sdsl::bit_vector second(99999999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_all_ones) {
    {
        sdsl::bit_vector first(1, 1);
        sdsl::bit_vector second(1, 1);
        EXPECT_EQ(1u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10, 1);
        sdsl::bit_vector second(10, 1);
        EXPECT_EQ(10u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999, 1);
        sdsl::bit_vector second(999, 1);
        EXPECT_EQ(999u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999, 1);
        sdsl::bit_vector second(999999, 1);
        EXPECT_EQ(999999u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999, 1);
        sdsl::bit_vector second(99999999, 1);
        EXPECT_EQ(99999999u, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_all_ones_offset) {
    for (size_t size = 0; size < 200; ++size) {
        sdsl::bit_vector first(size, 1);
        sdsl::bit_vector second(size, 1);
        EXPECT_EQ(size, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 99999999 + 200; ++size) {
        sdsl::bit_vector first(size, 1);
        sdsl::bit_vector second(size, 1);
        EXPECT_EQ(size, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_same) {
    {
        sdsl::bit_vector first(1);
        sdsl::bit_vector second(1);
        uint64_t prod = 0;
        uint64_t bits = 1;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10);
        sdsl::bit_vector second(10);
        uint64_t prod = 0;
        uint64_t bits = 10;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999);
        sdsl::bit_vector second(999);
        uint64_t prod = 0;
        uint64_t bits = 999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999);
        sdsl::bit_vector second(999999);
        uint64_t prod = 0;
        uint64_t bits = 999999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999);
        sdsl::bit_vector second(99999999);
        uint64_t prod = 0;
        uint64_t bits = 99999999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        uint64_t prod = 0;
        uint64_t bits = size;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 9999999 + 200; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        uint64_t prod = 0;
        uint64_t bits = size;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_disjoint) {
    {
        sdsl::bit_vector first(1);
        sdsl::bit_vector second(1);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10);
        sdsl::bit_vector second(10);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999);
        sdsl::bit_vector second(999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999);
        sdsl::bit_vector second(999999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999);
        sdsl::bit_vector second(99999999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 99999999 + 80; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
}

// the common set bits of |first| and |second| below their size, bit by bit
uint64_t inner_prod_naive(const sdsl::bit_vector &first, const sdsl::bit_vector &second) {
    uint64_t count = 0;
    for (size_t i = 0; i < first.size(); ++i) {
        count += first[i] && second[i];
    }
    return count;
}

// pseudo-random words over the logical words of |v|
void fill_words(sdsl::bit_vector *v, uint64_t seed) {
    for (size_t i = 0; i < (v->size() + 63) >> 6; ++i) {
        v->data()[i] = 18446744073709551557llu * (seed + i) + i * 32416189321llu;
    }
}

const std::vector<size_t> kSizesAroundWords = { 1, 2, 63, 64, 65, 127, 128, 129, 255, 256, 257,
                                                1000, 4095, 4096, 4097, 16384, 16385 };

// The pinned SDSL keeps a capacity above the rounded size after reserve(): inner_prod reads the
// words of the size, of both vectors, and masks the last of them, whatever either capacity
// (GPT review 5, finding 10: it read |first|'s capacity, over |second|'s end)
TEST(bit_vector, inner_prod_unequal_capacity_reserve) {
    for (size_t size : kSizesAroundWords) {
        for (bool reserve_first : { true, false }) {
            sdsl::bit_vector first(size);
            sdsl::bit_vector second(size);
            fill_words(&first, size);
            fill_words(&second, size + 7);
            sdsl::bit_vector &reserved = reserve_first ? first : second;
            reserved.reserve(4 * size + 4096);
            ASSERT_GT(reserved.capacity(), ((size + 63) >> 6) << 6) << size;
            // stale ones in every word above the size, and above the size in its last word
            for (size_t i = (size + 63) >> 6; i < reserved.capacity() >> 6; ++i) {
                reserved.data()[i] = ~uint64_t(0);
            }
            if (size & 0x3F)
                reserved.data()[size >> 6] |= ~sdsl::bits::lo_set[size & 0x3F];
            EXPECT_EQ(inner_prod_naive(first, second), inner_prod(first, second))
                    << size << (reserve_first ? " first" : " second") << " reserved";
        }
    }
    // the reviewer's case: two equal vectors of one word, the first with 64 words reserved
    sdsl::bit_vector a(64, 1);
    sdsl::bit_vector b(64, 1);
    a.reserve(4096);
    EXPECT_EQ(64u, inner_prod(a, b));
    EXPECT_EQ(64u, inner_prod(b, a));
}

// A shrinking resize keeps the capacity and the bits above the new size
TEST(bit_vector, inner_prod_shrinking_resize) {
    for (size_t size : kSizesAroundWords) {
        sdsl::bit_vector first(size + 1000, 1);
        first.resize(size);
        ASSERT_GT(first.capacity(), ((size + 63) >> 6) << 6) << size;
        sdsl::bit_vector second(size, 1);
        EXPECT_EQ(size, inner_prod(first, second)) << size;
        EXPECT_EQ(size, inner_prod(second, first)) << size;
        // and both shrunk from different sizes, with other bits
        sdsl::bit_vector other(size + 3000);
        fill_words(&other, size);
        other.resize(size);
        EXPECT_EQ(inner_prod_naive(first, other), inner_prod(first, other)) << size;
        EXPECT_EQ(inner_prod_naive(other, first), inner_prod(other, first)) << size;
    }
}

// push_back() grows the capacity by half its size
TEST(bit_vector, inner_prod_grown_by_push_back) {
    for (size_t size : kSizesAroundWords) {
        sdsl::bit_vector first;
        for (size_t i = 0; i < size; ++i) {
            first.push_back(i % 3 != 0);
        }
        sdsl::bit_vector second(size);
        fill_words(&second, size);
        EXPECT_EQ(inner_prod_naive(first, second), inner_prod(first, second)) << size;
        EXPECT_EQ(inner_prod_naive(second, first), inner_prod(second, first)) << size;
    }
}


TEST(bit_vector, to_sdsl_empty) {
    EXPECT_EQ(0u, to_sdsl(std::vector<bool>()).size());
    EXPECT_EQ(0u, to_sdsl(std::vector<uint8_t>()).size());
}

TEST(bit_vector, to_sdsl_zeros) {
    for (size_t i = 0; i < 1025; ++i) {
        auto vector_bool = to_sdsl(std::vector<bool>(i, false));
        EXPECT_EQ(i, vector_bool.size());
        EXPECT_EQ(0u, sdsl::util::cnt_one_bits(vector_bool));

        auto vector_uint8_t = to_sdsl(std::vector<uint8_t>(i, false));
        EXPECT_EQ(i, vector_uint8_t.size());
        EXPECT_EQ(0u, sdsl::util::cnt_one_bits(vector_uint8_t));
    }
}

TEST(bit_vector, to_sdsl_ones) {
    for (size_t i = 0; i < 1025; ++i) {
        auto vector_bool = to_sdsl(std::vector<bool>(i, true));
        ASSERT_EQ(i, vector_bool.size());
        ASSERT_EQ(i, sdsl::util::cnt_one_bits(vector_bool));

        auto vector_uint8_t = to_sdsl(std::vector<uint8_t>(i, true));
        ASSERT_EQ(i, vector_uint8_t.size());
        ASSERT_EQ(i, sdsl::util::cnt_one_bits(vector_uint8_t));

        auto vector_uint8_msb = to_sdsl(std::vector<uint8_t>(i, 4));
        ASSERT_EQ(i, vector_uint8_msb.size());
        ASSERT_EQ(i, sdsl::util::cnt_one_bits(vector_uint8_msb));

        auto vector_uint8_lsb = to_sdsl(std::vector<uint8_t>(i, 1));
        ASSERT_EQ(i, vector_uint8_lsb.size());
        ASSERT_EQ(i, sdsl::util::cnt_one_bits(vector_uint8_lsb));
    }
}

TEST(bit_vector, to_sdsl) {
    std::vector<bool> bits_bool;
    std::vector<uint8_t> bits_uint8_t_a;
    std::vector<uint8_t> bits_uint8_t_b;
    std::vector<uint64_t> ranks = { 0 };

    for (size_t i = 0; i < 10'000'000; ++i) {
        bits_bool.push_back((i + (i * i) % 31) % 2);
        bits_uint8_t_a.push_back(static_cast<bool>((i + (i * i) % 31) % 2));
        bits_uint8_t_b.push_back(((i + (i * i) % 31) % 2) ? 0xFF : 0);
        if (bits_uint8_t_a.back()) {
            ranks.back()++;
        }
        ranks.push_back(ranks.back());
    }

    // TODO: Don't use bit_vector_stat in this test.
    //       Just check the bits of the resulting sdsl::bit_vector.
    bit_vector_stat a(to_sdsl(bits_bool));
    bit_vector_stat b(to_sdsl(bits_uint8_t_a));
    bit_vector_stat c(to_sdsl(bits_uint8_t_b));

    for (size_t i = 0; i < bits_bool.size(); ++i) {
        EXPECT_EQ(ranks[i], a.rank1(i));
        EXPECT_EQ(ranks[i], b.rank1(i));
        EXPECT_EQ(ranks[i], c.rank1(i));
    }
}

TEST(bit_vector, autocorrelate) {
    std::vector<std::tuple<sdsl::bit_vector, sdsl::bit_vector, uint8_t>> results {
        { sdsl::bit_vector({}), sdsl::bit_vector({}), 5 },
        { sdsl::bit_vector({ 0, 0, 0, 0, 0 }), sdsl::bit_vector({ 0, 0, 0, 0, 0 }), 5 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 5 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 4 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 3 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 2 },
        { sdsl::bit_vector({ 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  }),
          sdsl::bit_vector({ 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  }), 3},
    };

    for (const auto &[input, output, offset] : results) {
        sdsl::bit_vector reference(input);

        if (input.size() >= offset) {
            for (size_t i = 0; i < reference.size(); ++i) {
                for (uint8_t j = 1; j < offset && i + j < reference.size(); ++j) {
                    reference[i] = reference[i] && input[i + j];
                }
            }
        }

        ASSERT_EQ(reference, output);
        EXPECT_EQ(reference, autocorrelate(input, offset));
    };
}


TEST(select_support_scan_offset, vector_zeros) {
    sdsl::bit_vector bv(156, false);
    sdsl::select_support_scan_offset<0> select(&bv);

    EXPECT_EQ(0UL, select.select_offset(1, 0));
    EXPECT_EQ(1UL, select.select_offset(2, 0));
    EXPECT_EQ(2UL, select.select_offset(3, 0));
    EXPECT_EQ(62UL, select.select_offset(63, 0));
    EXPECT_EQ(63UL, select.select_offset(64, 0));
    EXPECT_EQ(64UL, select.select_offset(65, 0));
    EXPECT_EQ(100UL, select.select_offset(101, 0));
    EXPECT_EQ(155UL, select.select_offset(156, 0));

    EXPECT_EQ(0UL + 1, select.select_offset(1, 1));
    EXPECT_EQ(1UL + 1, select.select_offset(2, 1));
    EXPECT_EQ(2UL + 1, select.select_offset(3, 1));
    EXPECT_EQ(62UL + 1, select.select_offset(63, 1));
    EXPECT_EQ(63UL + 1, select.select_offset(64, 1));
    EXPECT_EQ(64UL + 1, select.select_offset(65, 1));
    EXPECT_EQ(100UL + 1, select.select_offset(101, 1));
    EXPECT_EQ(155UL, select.select_offset(155, 1));

    EXPECT_EQ(0UL + 64, select.select_offset(1, 64));
    EXPECT_EQ(1UL + 64, select.select_offset(2, 64));
    EXPECT_EQ(2UL + 64, select.select_offset(3, 64));
    EXPECT_EQ(62UL + 64, select.select_offset(63, 64));
    EXPECT_EQ(63UL + 64, select.select_offset(64, 64));
    EXPECT_EQ(64UL + 64, select.select_offset(65, 64));
    EXPECT_EQ(155UL, select.select_offset(92, 64));

    EXPECT_EQ(100UL, select.select_offset(1, 100));
    EXPECT_EQ(101UL, select.select_offset(2, 100));
    EXPECT_EQ(102UL, select.select_offset(3, 100));
    EXPECT_EQ(155UL, select.select_offset(56, 100));
}

TEST(select_support_scan_offset, vector_ones) {
    sdsl::bit_vector bv(156, true);
    sdsl::select_support_scan_offset<> select(&bv);

    EXPECT_EQ(0UL, select.select_offset(1, 0));
    EXPECT_EQ(1UL, select.select_offset(2, 0));
    EXPECT_EQ(2UL, select.select_offset(3, 0));
    EXPECT_EQ(62UL, select.select_offset(63, 0));
    EXPECT_EQ(63UL, select.select_offset(64, 0));
    EXPECT_EQ(64UL, select.select_offset(65, 0));
    EXPECT_EQ(100UL, select.select_offset(101, 0));
    EXPECT_EQ(155UL, select.select_offset(156, 0));

    EXPECT_EQ(0UL + 1, select.select_offset(1, 1));
    EXPECT_EQ(1UL + 1, select.select_offset(2, 1));
    EXPECT_EQ(2UL + 1, select.select_offset(3, 1));
    EXPECT_EQ(62UL + 1, select.select_offset(63, 1));
    EXPECT_EQ(63UL + 1, select.select_offset(64, 1));
    EXPECT_EQ(64UL + 1, select.select_offset(65, 1));
    EXPECT_EQ(100UL + 1, select.select_offset(101, 1));
    EXPECT_EQ(155UL, select.select_offset(155, 1));

    EXPECT_EQ(0UL + 64, select.select_offset(1, 64));
    EXPECT_EQ(1UL + 64, select.select_offset(2, 64));
    EXPECT_EQ(2UL + 64, select.select_offset(3, 64));
    EXPECT_EQ(62UL + 64, select.select_offset(63, 64));
    EXPECT_EQ(63UL + 64, select.select_offset(64, 64));
    EXPECT_EQ(64UL + 64, select.select_offset(65, 64));
    EXPECT_EQ(155UL, select.select_offset(92, 64));

    EXPECT_EQ(100UL, select.select_offset(1, 100));
    EXPECT_EQ(101UL, select.select_offset(2, 100));
    EXPECT_EQ(102UL, select.select_offset(3, 100));
    EXPECT_EQ(155UL, select.select_offset(56, 100));
}

TEST(select_support_scan_offset, vector_every_seventh_set) {
    sdsl::bit_vector bv(10000, false);
    for (size_t i = 0; i < bv.size(); i += 7) {
        bv[i] = true;
    }

    sdsl::select_support_scan_offset<> select(&bv);

    bit_vector_stat ref(bv);

    const uint64_t num_set_bits = sdsl::util::cnt_one_bits(bv);

    for (size_t begin = 0; begin < bv.size(); ++begin) {
        uint64_t rank_offset = begin ? ref.rank1(begin - 1) : 0;

        for (size_t r = rank_offset + 1; r <= num_set_bits; ++r) {
            ASSERT_EQ(ref.select1(r), select.select_offset(r - rank_offset, begin))
                << "offset: " << begin
                << ", rank: " << r - rank_offset;
        }
    }
}

TEST(select_support_scan_offset, vector_random) {
    srand(41);

    for (size_t size = 1; size <= 1000; ++size) {
        sdsl::bit_vector bv(size);
        for (size_t i = 0; i < bv.size(); i += 7) {
            bv[i] = rand() & 1;
        }

        sdsl::select_support_scan_offset<> select(&bv);

        bit_vector_stat ref(bv);

        const uint64_t num_set_bits = sdsl::util::cnt_one_bits(bv);

        for (size_t begin = 0; begin < bv.size(); ++begin) {
            uint64_t rank_offset = begin ? ref.rank1(begin - 1) : 0;

            for (size_t r = rank_offset + 1; r <= num_set_bits; ++r) {
                ASSERT_EQ(ref.select1(r), select.select_offset(r - rank_offset, begin))
                    << "offset: " << begin
                    << ", rank: " << r - rank_offset;
            }
        }
    }
}

std::vector<uint64_t> get_test_size_range() {
    std::vector<uint64_t> result;
    for (size_t size = 1; size <= 1'000; ++size) {
        result.push_back(size);
    }
    for (size_t l = 10; l <= 21; ++l) {
        result.push_back(1llu << l);
        result.push_back((1llu << l) + 1);
        result.push_back((1llu << l) - 1);
    }
    return result;
}

TEST(generate_subindex, empty) {
    ThreadPool thread_pool(5);
    sdsl::bit_vector ref;
    bit_vector_stat col;
    ASSERT_EQ(sdsl::bit_vector(0),
        generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
}

TEST(generate_subindex, ref0_col0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector ref(size, 0);
        bit_vector_stat col(sdsl::bit_vector(size, 0));
        ASSERT_EQ(sdsl::bit_vector(0),
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1_col0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector ref(size, 1);
        bit_vector_stat col(sdsl::bit_vector(size, 0));
        ASSERT_EQ(sdsl::bit_vector(size, 0),
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1_col1) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector ref(size, 1);
        bit_vector_stat col(sdsl::bit_vector(size, 1));
        ASSERT_EQ(sdsl::bit_vector(size, 1),
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1first0_col1first0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bv[0] = 0;
        sdsl::bit_vector ref = bv;
        bit_vector_stat col(bv);
        ASSERT_EQ(sdsl::bit_vector(size - 1, 1),
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1last0_col1last0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bv[bv.size() - 1] = 0;
        sdsl::bit_vector ref = bv;
        bit_vector_stat col(bv);
        ASSERT_EQ(sdsl::bit_vector(size - 1, 1),
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1_col1first0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        sdsl::bit_vector ref = bv;
        bv[0] = 0;
        bit_vector_stat col(bv);
        ASSERT_EQ(bv,
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex, ref1_col1last0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        sdsl::bit_vector ref = bv;
        bv[bv.size() - 1] = 0;
        bit_vector_stat col(bv);
        ASSERT_EQ(bv,
            generate_subindex(col, ref, sdsl::util::cnt_one_bits(ref), thread_pool));
    }
}

TEST(generate_subindex_sparse, empty) {
    ThreadPool thread_pool(5);
    bit_vector_stat ref;
    bit_vector_stat col;
    ASSERT_EQ(sdsl::bit_vector(0), generate_subindex(col, ref, thread_pool));
}

TEST(generate_subindex_sparse, ref0_col0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        bit_vector_stat ref(sdsl::bit_vector(size, 0));
        bit_vector_stat col(sdsl::bit_vector(size, 0));
        ASSERT_EQ(sdsl::bit_vector(0), generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1_col0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        bit_vector_stat ref(sdsl::bit_vector(size, 1));
        bit_vector_stat col(sdsl::bit_vector(size, 0));
        ASSERT_EQ(sdsl::bit_vector(size, 0), generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1_col1) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        bit_vector_stat ref(sdsl::bit_vector(size, 1));
        bit_vector_stat col(sdsl::bit_vector(size, 1));
        ASSERT_EQ(sdsl::bit_vector(size, 1), generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1first0_col1first0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bv[0] = 0;
        bit_vector_stat ref(bv);
        bit_vector_stat col(bv);
        ASSERT_EQ(sdsl::bit_vector(size - 1, 1), generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1last0_col1last0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bv[bv.size() - 1] = 0;
        bit_vector_stat ref(bv);
        bit_vector_stat col(bv);
        ASSERT_EQ(sdsl::bit_vector(size - 1, 1), generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1_col1first0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bit_vector_stat ref(bv);
        bv[0] = 0;
        bit_vector_stat col(bv);
        ASSERT_EQ(bv, generate_subindex(col, ref, thread_pool));
    }
}

TEST(generate_subindex_sparse, ref1_col1last0) {
    ThreadPool thread_pool(5);
    for (uint64_t size : get_test_size_range()) {
        sdsl::bit_vector bv(size, 1);
        bit_vector_stat ref(bv);
        bv[bv.size() - 1] = 0;
        bit_vector_stat col(bv);
        ASSERT_EQ(bv, generate_subindex(col, ref, thread_pool));
    }
}

} // namespace
