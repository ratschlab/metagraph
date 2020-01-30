#include "int_vector_algorithm.hpp"

#include "common/utils/simd_utils.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include <sdsl/uint128_t.hpp>


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector) {
    sdsl::bit_vector result(vector.size(), 0);

    for (size_t i = 0; i < vector.size(); ++i) {
        if (vector[i])
            result[i] = 1;
    }

    return result;
}

sdsl::bit_vector to_sdsl(const std::vector<uint8_t> &vector) {
    sdsl::bit_vector result(vector.size(), 0);

    size_t i = 0;
#ifdef __AVX2__
    for (; i + 32 <= vector.size(); i += 32) {
        result.set_int(
            i,
            ~_mm256_movemask_epi8(_mm256_cmpeq_epi8(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&vector[i])),
                _mm256_setzero_si256()
            )),
            32
        );
    }
#endif

#ifdef __SSE2__
    for (; i + 16 <= vector.size(); i += 16) {
        result.set_int(
            i,
            ~_mm_movemask_epi8(_mm_cmpeq_epi8(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(&vector[i])),
                _mm_setzero_si128()
            )),
            16
        );
    }
#endif

    // at most 16 bits left
    uint16_t last_word = 0;
    size_t last_i = i;
    size_t j = 0;
    for (; i < vector.size(); ++i, ++j) {
        if (vector[i])
            last_word |= static_cast<uint16_t>(1) << j;
    }

    result.set_int(last_i, last_word, i - last_i);

    assert(static_cast<size_t>(std::count_if(vector.begin(), vector.end(),
                                             [](uint8_t a) { return a; }))
        == sdsl::util::cnt_one_bits(result));

    return result;
}

uint64_t count_ones(const sdsl::bit_vector &vector,
                    uint64_t begin, uint64_t end) {
    assert(begin <= end);
    assert(end <= vector.size());

    if (begin == end)
        return 0;

    if (end - begin <= 64)
        return sdsl::bits::cnt(vector.get_int(begin, end - begin));

    const uint64_t *data = vector.data() + (begin >> 6);
    const uint64_t *data_end = vector.data() + ((end + 63) >> 6);

    uint64_t count = 0;

    if (begin & 0x3F) {
        count += sdsl::bits::cnt((*data++) & (~sdsl::bits::lo_set[begin & 0x3F]));
    }

#ifdef __AVX2__
    size_t diff = ((data_end - data) >> 6) << 6;
    __m256i counts = popcnt_avx2_hs(data, diff);
    data += diff;

    for (; data + 4 <= data_end; data += 4) {
        counts = _mm256_add_epi64(
            counts,
            popcnt256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(data)))
        );
    }

    count += haddall_epi64(counts);
#endif

    while (data < data_end) {
        count += sdsl::bits::cnt(*data++);
    }

    if (end & 0x3F)
        count -= sdsl::bits::cnt((*(--data)) & (~sdsl::bits::lo_set[end & 0x3F]));

    return count;
}

uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second) {
    assert(first.size() == second.size());

    if (first.empty())
        return 0;

    const uint64_t *first_data = first.data();
    const uint64_t *first_end = first.data() + (first.capacity() >> 6);
    const uint64_t *second_data = second.data();

    uint64_t count = 0;

#ifdef __AVX2__
    size_t diff = ((first_end - first_data) >> 6) << 6;
    __m256i counts = inner_prod_avx2_hs(first_data, second_data, diff);
    first_data += diff;
    second_data += diff;

    for (; first_data + 4 <= first_end; first_data += 4, second_data += 4) {
        counts = _mm256_add_epi64(
            counts,
            popcnt256(_mm256_and_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(first_data)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(second_data))
            ))
        );
    }

    count += haddall_epi64(counts);
#endif

    while (first_data < first_end) {
        count += sdsl::bits::cnt(*(first_data++) & *(second_data++));
    }

    if (first.size() & 0x3F) {
        count -= sdsl::bits::cnt((*(--first_data))
            & *(--second_data)
            & (~sdsl::bits::lo_set[first.size() & 0x3F]));
    }

    return count;
}


inline sdsl::uint128_t pushback_epi64(const sdsl::uint128_t &v, uint64_t a) {
    return (v >> 64) | (sdsl::uint128_t(a) << 64);
}

sdsl::bit_vector autocorrelate(const sdsl::bit_vector &vector, uint8_t offset) {
    assert(offset < 64);

    if (vector.size() < offset)
        return vector;

    auto presence = vector;

    // process one word at a time
    // TODO: is it worth vectorizing this?
    size_t i = 0;
    auto dword = sdsl::uint128_t(vector.data()[0]) << 64;
    for (; i + 64 <= presence.size() - offset + 1; i += 64) {
        dword = pushback_epi64(dword, vector.data()[(i >> 6) + 1]);
        for (uint8_t j = 1; j < offset; ++j) {
            presence.data()[i >> 6] &= uint64_t(dword >> j);
        }
    }

    assert(presence.size() - i >= static_cast<size_t>(offset) - 1);
    assert(presence.size() - i <= 128 - static_cast<size_t>(offset) + 1);

    // handle last word
    if (vector.size() - i <= 64) {
        uint64_t word = vector.get_int(i, vector.size() - i);
        uint64_t word_masked = word;
        uint64_t mask = 0;
        for (uint8_t j = 1; j < offset; ++j) {
            mask |= uint64_t(1) << (vector.size() - i - j);
            word_masked &= (word >> j) | mask;
        }

        presence.set_int(i, word_masked, vector.size() - i);

    } else {
        dword = pushback_epi64(dword, vector.data()[(i >> 6) + 1])
            | (sdsl::uint128_t((1llu << offset) - 1) << (vector.size() - i));
        sdsl::uint128_t dword_masked = dword;
        for (uint8_t j = 1; j < offset; ++j) {
            dword_masked &= dword >> j;
        }
        presence.set_int(i, uint64_t(dword_masked));
        presence.set_int(i + 64, uint64_t(dword_masked >> 64), vector.size() - i - 64);
    }

#ifndef NDEBUG
    for (size_t i = 0; i < presence.size(); ++i) {
        bool b = vector[i];
        for (uint8_t j = 1; j < offset && i + j < presence.size(); ++j) {
            b &= vector[i + j];
        }
        assert(b == presence[i]);
    }
#endif

    return presence;
}
