#include "int_vector_algorithm.hpp"

#include "common/utils/simd_utils.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector) {
    sdsl::bit_vector result(vector.size(), 0);
    for (size_t i = 0; i < vector.size(); ++i) {
        if (vector[i])
            result[i] = 1;
    }
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
