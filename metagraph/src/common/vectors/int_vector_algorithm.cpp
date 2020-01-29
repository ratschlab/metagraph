#include "int_vector_algorithm.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif


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

#ifdef __AVX2__

// from: https://github.com/WojciechMula/libalgebra/blob/master/libalgebra.h

// carry-save added
void CSA256(__m256i *hi, __m256i *lo, __m256i a, __m256i b, __m256i c) {
    __m256i u = _mm256_xor_si256(a, b);
    *hi = _mm256_or_si256(_mm256_and_si256(a, b), _mm256_and_si256(u, c));
    *lo = _mm256_xor_si256(u, c);
}

__m256i popcnt256(__m256i v) {
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


__m256i popcnt_avx2_hs(const uint64_t *data, uint64_t size) {
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

__m256i inner_prod_avx2_hs(const uint64_t *data1,
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


#endif

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

    // [ a, b, c, d ] -> [ a+c, b+d, c+a, d+b ]
    __m256i s1 = _mm256_add_epi64(counts, _mm256_permute4x64_epi64(counts, 0b01001110));

    // [ a+c, b+d, c+a, d+b ] -> a+c+b+d
    count += _mm256_extract_epi64(
        _mm256_add_epi64(s1, _mm256_permute4x64_epi64(s1, 0b10001101)), 0
    );
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

    // [ a, b, c, d ] -> [ a+c, b+d, c+a, d+b ]
    __m256i s1 = _mm256_add_epi64(counts, _mm256_permute4x64_epi64(counts, 0b01001110));

    // [ a+c, b+d, c+a, d+b ] -> a+c+b+d
    count += _mm256_extract_epi64(
        _mm256_add_epi64(s1, _mm256_permute4x64_epi64(s1, 0b10001101)), 0
    );
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
