#include "vector_algorithm.hpp"

#include <sdsl/uint128_t.hpp>

#include "common/utils/simd_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/bit_vector_stat.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

using sdsl::uint128_t;

const uint64_t kBlockSize = 1'000'000 & ~0x3Full;
// Each block is a multiple of 64 bits for thread safety
static_assert((kBlockSize % 64) == 0);


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

    uint64_t i = 0;
#ifdef __AVX2__
    for ( ; i + 32 <= vector.size(); i += 32) {
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
    for ( ; i + 16 <= vector.size(); i += 16) {
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

    for ( ; i < vector.size(); i += 64) {
        uint64_t word = 0;
        uint8_t width = std::min(static_cast<uint64_t>(64), result.size() - i);
        for (int64_t j = i + width - 1; j >= static_cast<int64_t>(i); --j) {
            word = (word << 1) | static_cast<bool>(vector[j]);
        }

        result.set_int(i, word, width);
    }

    assert(static_cast<size_t>(std::count_if(vector.begin(), vector.end(),
                                             [](uint8_t a) { return a; }))
        == sdsl::util::cnt_one_bits(result));

    return result;
}

sdsl::bit_vector copy_to_bit_vector(const bit_vector &vector,
                                    double WORD_ACCESS_VS_SELECT_FACTOR) {
    sdsl::bit_vector result;

    if (vector.num_set_bits() <= vector.size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // sparse
        result = sdsl::bit_vector(vector.size(), false);

        uint64_t num_ones = vector.num_set_bits();
        for (uint64_t r = 1; r <= num_ones; ++r) {
            result[vector.select1(r)] = true;
        }

    } else if ((vector.size() - vector.num_set_bits())
                <= vector.size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // dense
        result = sdsl::bit_vector(vector.size(), true);

        uint64_t num_zeros = vector.size() - vector.num_set_bits();
        for (uint64_t r = 1; r <= num_zeros; ++r) {
            result[vector.select0(r)] = false;
        }

    } else {
        // moderate density
        result.resize(vector.size());

        uint64_t i;
        const uint64_t end = vector.size();
        uint64_t *data = result.data();
        for (i = 0; i + 64 <= end; i += 64, ++data) {
            *data = vector.get_int(i, 64);
        }
        if (i < vector.size()) {
            *data = vector.get_int(i, vector.size() - i);
        }
    }

    return result;
}

void add_to(const bit_vector &vector,
            sdsl::bit_vector *other,
            double WORD_ACCESS_VS_SELECT_FACTOR) {
    assert(other);
    assert(other->size() == vector.size());

    if (std::min(vector.num_set_bits(),
                 vector.size() - vector.num_set_bits())
            <= vector.size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // for very sparse or very dense vectors
        call_ones(vector, 0, vector.size(),
                  [other](auto i) { (*other)[i] = true; },
                  WORD_ACCESS_VS_SELECT_FACTOR);
    } else {
        // moderate density
        uint64_t i;
        const uint64_t end = vector.size();
        uint64_t *data = other->data();
        for (i = 0; i + 64 <= end; i += 64, ++data) {
            *data |= vector.get_int(i, 64);
        }
        if (i < vector.size()) {
            *data |= vector.get_int(i, vector.size() - i);
        }
    }
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

    for ( ; data + 4 <= data_end; data += 4) {
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

    for ( ; first_data + 4 <= first_end; first_data += 4, second_data += 4) {
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

// check if the bit vector is sparser than the threshold for its type
bool is_sparser(const bit_vector &bv,
                double threshold_stat,
                double threshold_sd,
                double threshold_rrr,
                double threshold_other) {
    if (bv.size() < 64u)
        return false;

    double density = static_cast<double>(bv.num_set_bits()) / bv.size();

    if (auto *adaptive = dynamic_cast<const bit_vector_adaptive *>(&bv)) {
        switch (adaptive->representation_tag()) {
            case bit_vector_adaptive::STAT_VECTOR:
                return density < threshold_stat;
            case bit_vector_adaptive::SD_VECTOR:
                return density < threshold_sd;
            case bit_vector_adaptive::RRR_VECTOR:
                return density < threshold_rrr;
            case bit_vector_adaptive::RANK_VECTOR:
                return density < threshold_stat;
        }
    }

    if (dynamic_cast<const bit_vector_stat *>(&bv))
        return density < threshold_stat;

    if (dynamic_cast<const bit_vector_sd *>(&bv))
        return density < threshold_sd;

    if (dynamic_cast<const bit_vector_rrr<> *>(&bv))
        return density < threshold_rrr;

    if (dynamic_cast<const bit_vector_rank *>(&bv))
        return density < threshold_stat;

    return density < threshold_other;
}

void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool) {
    uint64_t size = columns.at(0)->size();

    assert(result);
    assert(result->size() == size);

    const uint64_t block_size = std::max(kBlockSize, size / 100 / 64 * 64);

    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::future<void>> results;

    for (uint64_t i = 0; i < size; i += block_size) {
        results.push_back(thread_pool.enqueue([&](uint64_t begin, uint64_t end) {
            std::fill(result->data() + (begin >> 6),
                      result->data() + ((end + 63) >> 6),
                      0);

            for (const auto &col_ptr : columns) {
                assert(col_ptr);
                assert(col_ptr->size() == size);

                if (is_sparser(*col_ptr, 0, 0.2, 0.01, 0.01)) {
                    col_ptr->call_ones_in_range(begin, end,
                        [result](uint64_t k) { (*result)[k] = true; }
                    );
                } else {
                    assert((begin & 0x3F) == 0);

                    uint64_t i;
                    for (i = begin; i + 64 <= end; i += 64) {
                        result->data()[i / 64] |= col_ptr->get_int(i, 64);
                    }
                    if (i < end)
                        result->data()[i / 64] |= col_ptr->get_int(i, end - i);
                }
            }

        }, i, std::min(i + block_size, size)));
    }

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });
}

void compute_subindex(const bit_vector &column,
                      const sdsl::bit_vector &reference,
                      uint64_t begin, uint64_t end,
                      uint64_t reference_rank_begin,
                      sdsl::bit_vector *subindex);

sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference,
                                   uint64_t reference_num_set_bits,
                                   ThreadPool &thread_pool) {
    assert(column.size() == reference.size());

    sdsl::select_support_scan_offset<> reference_select(&reference);

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == reference_num_set_bits)
        return sdsl::bit_vector(reference_num_set_bits, true);

    sdsl::bit_vector subindex(reference_num_set_bits, false);

    const uint64_t block_size
      = std::max(reference_num_set_bits / 100,
                 static_cast<uint64_t>(1'000)) & ~0x3Full;
    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::future<void>> futures;

    uint64_t block_end = 0;

    for (uint64_t offset = 0; offset < reference_num_set_bits; offset += block_size) {
        // ........... [1     1       1    1  1] ............. //
        //              ^                     ^                //
        //          (rank = 1 mod 64)     (rank = 0 mod 64)    //
        //              |                     |                //
        //              |_______ BLOCK _______|                //

        uint64_t block_begin = block_end;

        if (offset + block_size >= reference_num_set_bits) {
            block_end = reference.size();
        } else {
            block_end = reference_select.select_offset(block_size + 1, block_begin);
            assert(count_ones(reference, block_begin, block_end) == block_size);
        }

        futures.push_back(thread_pool.enqueue(
            [&](uint64_t begin, uint64_t end, uint64_t offset_rank) {
                compute_subindex(column, reference, begin, end, offset_rank, &subindex);
            },
            block_begin, block_end, offset
        ));
    }

    std::for_each(futures.begin(), futures.end(), [](auto &f) { f.wait(); });

    return subindex;
}

sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const bit_vector_stat &reference) {
    assert(column.size() == reference.size());

    uint64_t shrinked_size = reference.num_set_bits();

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == shrinked_size)
        return sdsl::bit_vector(shrinked_size, true);

    sdsl::bit_vector subindex(shrinked_size, false);

    compute_subindex(column, reference.data(), 0, reference.size(), 0, &subindex);

    return subindex;
}

void compute_subindex(const bit_vector &column,
                      const sdsl::bit_vector &reference,
                      uint64_t begin, uint64_t end,
                      uint64_t offset,
                      sdsl::bit_vector *subindex) {
    assert(column.size() == reference.size());
    assert(begin <= end);
    assert(end <= reference.size());

    if (begin == end)
        return;

    uint64_t popcount = column.rank1(end - 1)
                        - (begin ? column.rank1(begin - 1) : 0);

    // check if all zeros
    if (!popcount)
        return;

    uint64_t i = begin;
    uint64_t rank = offset;

    // check if all ones
    if (popcount == end - begin) {
        for ( ; i + 64 <= end; i += 64, rank += 64) {
            subindex->set_int(rank, 0xFFFF, 64);
        }
        if (begin < end) {
            subindex->set_int(rank, sdsl::bits::lo_set[end - begin], end - begin);
        }
        return;
    }


#if __BMI2__
    if (is_sparser(column, 0, 0.02, 0.01, 0.01)) {
#endif
        // sparse or without __BMI2__
        column.call_ones_in_range(begin, end, [&](uint64_t j) {
            assert(j >= i);
            assert(reference[j]);

            rank += count_ones(reference, i, j);

            assert(rank < subindex->size());

            (*subindex)[rank++] = true;
            i = j + 1;
        });
#if __BMI2__
    } else {
        // dense
        for ( ; i + 64 <= end; i += 64) {
            uint64_t mask = reference.get_int(i, 64);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, 64))
                subindex->set_int(rank, _pext_u64(a, mask), popcount);
            rank += popcount;
        }
        if (i < end) {
            uint64_t mask = reference.get_int(i, end - i);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, end - i))
                subindex->set_int(rank, _pext_u64(a, mask), popcount);
        }
    }
#endif
}

// TODO: merge with compute_subindex
// indexes are distinct and sorted
sdsl::bit_vector subvector(const bit_vector &col,
                           const std::vector<uint64_t> &indexes) {
    assert(indexes.size() <= col.size());

    sdsl::bit_vector shrinked(indexes.size(), 0);

    uint64_t max_rank = col.num_set_bits();
    if (!max_rank)
        return shrinked;

    // the case of uncompressed vector
    if (dynamic_cast<const bit_vector_stat *>(&col)) {
        for (size_t j = 0; j < indexes.size(); ++j) {
            if (col[indexes[j]])
                shrinked[j] = true;
        }
        return shrinked;
    }

    uint64_t cur_rank = 1;
    uint64_t next_pos = col.select1(1);

    for (size_t j = 0; j < indexes.size(); ++j) {
        if (indexes[j] < next_pos)
            continue;

        if (indexes[j] == next_pos) {
            shrinked[j] = true;
            continue;
        }

        // indexes[j] > next_pos
        if (col[indexes[j]]) {
            shrinked[j] = true;
            continue;
        }

        // we found a zero, update next_pos
        cur_rank = col.rank1(indexes[j]) + 1;
        if (cur_rank > max_rank)
            return shrinked;

        next_pos = col.select1(cur_rank);
    }

    return shrinked;
}


inline uint128_t pushback_epi64(const uint128_t &v, const uint128_t &a) {
    return (v >> 64) | (a << 64);
}

sdsl::bit_vector autocorrelate(const sdsl::bit_vector &vector, uint8_t offset) {
    assert(offset < 64);

    if (vector.size() < offset)
        return vector;

    auto presence = vector;

    // process one word at a time
    // TODO: is it worth vectorizing this?
    size_t i = 0;
    auto dword = uint128_t(vector.data()[0]) << 64;
    for ( ; i + 64 <= presence.size() - offset + 1; i += 64) {
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
            | (uint128_t((1llu << offset) - 1) << (vector.size() - i));
        uint128_t dword_masked = dword;
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
