#ifndef __SIMD_UTILS_HPP__
#define __SIMD_UTILS_HPP__

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include "common/algorithms.hpp"


/**
 * Helpers for fast map of h uniformly to the region [0, size)
 * (h * size) >> 64
 */

inline uint64_t restrict_to_fallback(uint64_t h, size_t size) {
    // adapted from:
    // https://stackoverflow.com/questions/28868367/getting-the-high-part-of-64-bit-integer-multiplication

    const uint64_t x0 = static_cast<uint32_t>(h);
    const uint64_t x1 = h >> 32;
    const uint64_t y0 = static_cast<uint32_t>(size);
    const uint64_t y1 = size >> 32;

    const uint64_t p01 = x0 * y1;

    const uint64_t middle_hi = (x1 * y0
        + ((x0 * y0) >> 32)
        + static_cast<uint32_t>(p01)) >> 32;

    assert(!size || x1 * y1 + middle_hi + (p01 >> 32) < size);
    return x1 * y1 + middle_hi + (p01 >> 32);
}

#ifdef __SIZEOF_INT128__
inline uint64_t restrict_to_int128_fallback(uint64_t h, size_t size) {
    assert(!size || ((static_cast<__uint128_t>(h) * size) >> 64) < size);
    assert(restrict_to_fallback(h, size) == ((static_cast<__uint128_t>(h) * size) >> 64));

    return (static_cast<__uint128_t>(h) * size) >> 64;
}
#endif

#ifdef __BMI2__
inline uint64_t restrict_to_bmi2(uint64_t h, size_t size) {

#ifndef NDEBUG
    uint64_t h_fallback1 = restrict_to_fallback(h, size);
#ifdef __SIZEOF_INT128__
    uint64_t h_fallback2 = restrict_to_int128_fallback(h, size);
#endif
#endif

    static_assert(sizeof(long long unsigned int) == sizeof(uint64_t));
    _mulx_u64(h, size, reinterpret_cast<long long unsigned int*>(&h));

    assert(!size || h < size);
    assert(h_fallback1 == h);

#ifdef __SIZEOF_INT128__
    assert(h_fallback2 == h);
#endif

    return h;
}
#endif


/**
 * A fast map of h uniformly to the region [0, size)
 * Returns floor((h * size) / 2^64)
 */
inline uint64_t restrict_to(uint64_t h, size_t size) {
#ifdef __BMI2__
    return restrict_to_bmi2(h, size);
#elif __SIZEOF_INT128__
    return restrict_to_int128_fallback(h, size);
#else
    return restrict_to_fallback(h, size);
#endif
}


#ifdef __AVX2__

inline __m256i restrict_to_mask_epi64(const uint64_t *hashes, size_t size, __m256i mask) {
    // TODO: is there a vectorized way of doing this?
    return _mm256_and_si256(
        _mm256_setr_epi64x(restrict_to(hashes[0], size),
                           restrict_to(hashes[1], size),
                           restrict_to(hashes[2], size),
                           restrict_to(hashes[3], size)),
        mask
    );
}


/**
 * Helpers for count_ones and inner_prod
 */

// from: https://github.com/WojciechMula/libalgebra/blob/master/libalgebra.h

// carry-save adder
inline void CSA256(__m256i *hi, __m256i *lo, __m256i a, __m256i b, __m256i c) {
    __m256i u = _mm256_xor_si256(a, b);
    *hi = _mm256_or_si256(_mm256_and_si256(a, b), _mm256_and_si256(u, c));
    *lo = _mm256_xor_si256(u, c);
}

// Lookup popcount
inline __m256i popcnt256(__m256i v) {
    __m256i lookup1 = _mm256_setr_epi8(
        4, 5, 5, 6, 5, 6, 6, 7,
        5, 6, 6, 7, 6, 7, 7, 8,
        4, 5, 5, 6, 5, 6, 6, 7,
        5, 6, 6, 7, 6, 7, 7, 8
    );

    __m256i lookup2 = _mm256_setr_epi8(
        4, 3, 3, 2, 3, 2, 2, 1,
        3, 2, 2, 1, 2, 1, 1, 0,
        4, 3, 3, 2, 3, 2, 2, 1,
        3, 2, 2, 1, 2, 1, 1, 0
    );

    __m256i low_mask = _mm256_set1_epi8(0x0f);
    __m256i lo = _mm256_and_si256(v, low_mask);
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
    __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
    __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);

    return _mm256_sad_epu8(popcnt1, popcnt2);
}

// Harley-Seal algorithm for popcount
inline __m256i popcnt_avx2_hs(const uint64_t *data, uint64_t size) {
    __m256i total = _mm256_setzero_si256();
    __m256i ones = _mm256_setzero_si256();
    __m256i twos = _mm256_setzero_si256();
    __m256i fours = _mm256_setzero_si256();
    __m256i eights = _mm256_setzero_si256();
    __m256i sixteens = _mm256_setzero_si256();
    __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

    #define LOAD(a) _mm256_loadu_si256((__m256i*)&data[i + (a * 4)])
    for (uint64_t i = 0; i + 64 <= size; i += 64) {
        CSA256(&twosA, &ones, ones, LOAD(0), LOAD(1));
        CSA256(&twosB, &ones, ones, LOAD(2), LOAD(3));
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, LOAD(4), LOAD(5));
        CSA256(&twosB, &ones, ones, LOAD(6), LOAD(7));
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsA, &fours, fours, foursA, foursB);
        CSA256(&twosA, &ones, ones, LOAD(8), LOAD(9));
        CSA256(&twosB, &ones, ones, LOAD(10), LOAD(11));
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, LOAD(12), LOAD(13));
        CSA256(&twosB, &ones, ones, LOAD(14), LOAD(15));
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsB, &fours, fours, foursA, foursB);
        CSA256(&sixteens, &eights, eights, eightsA, eightsB);
        total = _mm256_add_epi64(total, popcnt256(sixteens));
    }
    #undef LOAD

    total = _mm256_slli_epi64(total, 4);
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(eights), 3));
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(fours), 2));
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(twos), 1));
    total = _mm256_add_epi64(total, popcnt256(ones));

    return total;
}

// Harley-Seal algorithm for inner product of bit vectors
inline __m256i inner_prod_avx2_hs(const uint64_t *data1,
                                  const uint64_t *data2,
                                  uint64_t size) {
    __m256i total = _mm256_setzero_si256();
    __m256i ones = _mm256_setzero_si256();
    __m256i twos = _mm256_setzero_si256();
    __m256i fours = _mm256_setzero_si256();
    __m256i eights = _mm256_setzero_si256();
    __m256i sixteens = _mm256_setzero_si256();
    __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

    #define LOAD(a) _mm256_and_si256(_mm256_loadu_si256((__m256i*)&data1[i + (a * 4)]), \
                                     _mm256_loadu_si256((__m256i*)&data2[i + (a * 4)]))
    for (uint64_t i = 0; i + 64 <= size; i += 64) {
        CSA256(&twosA, &ones, ones, LOAD(0), LOAD(1));
        CSA256(&twosB, &ones, ones, LOAD(2), LOAD(3));
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, LOAD(4), LOAD(5));
        CSA256(&twosB, &ones, ones, LOAD(6), LOAD(7));
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsA, &fours, fours, foursA, foursB);
        CSA256(&twosA, &ones, ones, LOAD(8), LOAD(9));
        CSA256(&twosB, &ones, ones, LOAD(10), LOAD(11));
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, LOAD(12), LOAD(13));
        CSA256(&twosB, &ones, ones, LOAD(14), LOAD(15));
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsB, &fours, fours, foursA, foursB);
        CSA256(&sixteens, &eights, eights, eightsA, eightsB);
        total = _mm256_add_epi64(total, popcnt256(sixteens));
    }
    #undef LOAD

    total = _mm256_slli_epi64(total, 4);
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(eights), 3));
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(fours), 2));
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt256(twos), 1));
    total = _mm256_add_epi64(total, popcnt256(ones));

    return total;
}

inline uint64_t haddall_epi64(__m256i v) {
    // [ a, b, c, d ] -> [ a+c, b+d, c+a, d+b ]
    __m256i s1 = _mm256_add_epi64(v, _mm256_permute4x64_epi64(v, 0b01001110));

    // [ a+c, b+d, c+a, d+b ] -> a+c+b+d
    return _mm256_extract_epi64(
        _mm256_add_epi64(s1, _mm256_permute4x64_epi64(s1, 0b10001101)),
        0
    );
}


/**
 * Helpers for aligner
 */

// Bit reduce packed 64-bit numbers to 32 bits
inline __m128i cvtepi64_epi32(__m256i v) {
    return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(
        v,
        _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7)
    ));
}

inline __m256i rshiftpushback_epi32(__m256i v, uint32_t a) {
    return _mm256_insert_epi32(
        _mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 0)),
        a,
        7
    );
}

inline __m256i expandepi8_epi64(uint64_t a) {
    return _mm256_cvtepi8_epi32(_mm256_castsi256_si128(_mm256_insert_epi64(
        _mm256_undefined_si256(), a, 0
    )));
}


#endif

#endif // __SIMD_UTILS_HPP__
