#include "vector_algorithm.hpp"

#include <cmath>
#include <mutex>

#include <sdsl/uint128_t.hpp>

#include "common/utils/simd_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

using sdsl::uint128_t;


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

sdsl::int_vector<> pack_vector(sdsl::int_vector<>&& vector, uint8_t width) {
    if (width == vector.width()) {
        return std::move(vector);
    } else {
        return pack_vector(vector, width);
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
            case bit_vector_adaptive::IL4096_VECTOR:
                return density < threshold_stat;
        }
    }

    if (dynamic_cast<const bit_vector_stat *>(&bv))
        return density < threshold_stat;

    if (dynamic_cast<const bit_vector_sd *>(&bv))
        return density < threshold_sd;

    if (dynamic_cast<const bit_vector_rrr<> *>(&bv))
        return density < threshold_rrr;

    if (dynamic_cast<const bit_vector_il<> *>(&bv))
        return density < threshold_stat;

    if (dynamic_cast<const bit_vector_il<4096> *>(&bv))
        return density < threshold_stat;

    return density < threshold_other;
}

void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool) {
    uint64_t size = columns.at(0)->size();

    assert(result);
    assert(result->size() == size);

    const uint64_t block_size
        = std::max(size / 100, static_cast<uint64_t>(1'000'000)) & ~0x3Full;
    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::shared_future<void>> results;

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

std::unique_ptr<bit_vector_adaptive>
compute_or(const std::vector<const bit_vector *> &columns,
           uint64_t *buffer,
           ThreadPool &thread_pool) {
    const uint64_t vector_size = columns.at(0)->size();
    const uint64_t block_size = std::max(vector_size / 100,
                                         static_cast<uint64_t>(100'000));

    std::vector<std::pair<uint64_t *, uint64_t *>>
            results((vector_size + block_size - 1) / block_size);

    std::vector<std::shared_future<void>> futures;
    std::atomic<uint64_t> num_set_bits = 0;

    for (uint64_t i = 0; i < vector_size; i += block_size) {
        futures.push_back(thread_pool.enqueue([&](uint64_t begin, uint64_t end) {
            // estimate the maximum number of set bits in the result
            uint64_t buffer_size = 0;
            for (size_t j = 0; j < columns.size(); ++j) {
                buffer_size += columns[j]->rank1(end - 1)
                                - (begin ? columns[j]->rank1(begin - 1) : 0);
            }

            std::pair<uint64_t *, uint64_t *> &result = results[begin / block_size];
            std::pair<uint64_t *, uint64_t *> buffer_pos;
            std::pair<uint64_t *, uint64_t *> buffer_merge;

            #pragma omp critical
            {
                // reserve space for 3 buffers of size |buffer_size| each
                result.first = result.second = buffer;
                buffer_pos.first = buffer_pos.second = buffer + buffer_size;
                buffer_merge.first = buffer_merge.second = buffer + 2 * buffer_size;
                buffer += 3 * buffer_size;
            }

            // write the positions of ones in the first vector to |result|
            columns[0]->call_ones_in_range(begin, end,
                [&](uint64_t k) { *(result.second++) = k; }
            );

            // iterate through all other columns and merge them sequentially
            //TODO: use multiway merge if there are more than two columns?
            for (size_t j = 1; j < columns.size(); ++j) {
                buffer_pos.second = buffer_pos.first;
                columns[j]->call_ones_in_range(begin, end,
                    [&](uint64_t k) { *(buffer_pos.second++) = k; }
                );

                std::merge(buffer_pos.first, buffer_pos.second,
                           result.first, result.second, buffer_merge.first);
                buffer_merge.second = buffer_merge.first
                                        + (buffer_pos.second - buffer_pos.first)
                                        + (result.second - result.first);
                result.swap(buffer_merge);
            }

            result.second = std::unique(result.first, result.second);

            num_set_bits += result.second - result.first;

        }, i, std::min(i + block_size, vector_size)));
    }

    std::for_each(futures.begin(), futures.end(), [](auto &f) { f.wait(); });

    return std::make_unique<bit_vector_smart>(
        [&](const auto &callback) {
            for (const auto &result : results) {
                std::for_each(result.first, result.second, callback);
            }
        },
        vector_size,
        num_set_bits
    );
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

    std::vector<std::shared_future<void>> futures;

    uint64_t block_end = 0;

    // Each block has |block_size| many bits in |reference| set to `1`
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
                                   const bit_vector &reference,
                                   ThreadPool &thread_pool) {
    assert(column.size() == reference.size());

    uint64_t reference_num_set_bits = reference.num_set_bits();

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == reference_num_set_bits)
        return sdsl::bit_vector(reference_num_set_bits, true);

    sdsl::bit_vector subindex(reference_num_set_bits, false);

    const uint64_t block_size
      = std::max(reference_num_set_bits / 100,
                 static_cast<uint64_t>(1'000)) & ~0x3Full;
    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::shared_future<void>> futures;

    uint64_t end = 0;

    // Each block has |block_size| many bits in |reference| set to `1`
    for (uint64_t offset = 0; offset < reference_num_set_bits; offset += block_size) {
        // ........... [1     1       1    1  1] ............. //
        //              ^                     ^                //
        //          (rank = 1 mod 64)     (rank = 0 mod 64)    //
        //              |                     |                //
        //              |_______ BLOCK _______|                //

        uint64_t begin = end;
        end = offset + block_size < reference_num_set_bits
                ? reference.select1(offset + block_size + 1)
                : reference.size();

        futures.push_back(thread_pool.enqueue([&,begin,end]() {
            column.call_ones_in_range(begin, end, [&](uint64_t j) {
                assert(reference[j]);
                subindex[reference.rank1(j) - 1] = true;
            });
        }));
    }

    std::for_each(futures.begin(), futures.end(), [](auto &f) { f.wait(); });

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

    // check if all ones
    if (popcount == end - begin) {
        for ( ; i + 64 <= end; i += 64, offset += 64) {
            subindex->set_int(offset, sdsl::bits::lo_set[64], 64);
        }
        if (i < end) {
            subindex->set_int(offset, sdsl::bits::lo_set[end - i], end - i);
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

            offset += count_ones(reference, i, j);

            assert(offset < subindex->size());

            (*subindex)[offset++] = true;
            i = j + 1;
        });
#if __BMI2__
    } else {
        // dense
        for ( ; i + 64 <= end; i += 64) {
            uint64_t mask = reference.get_int(i, 64);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, 64))
                subindex->set_int(offset, _pext_u64(a, mask), popcount);
            offset += popcount;
        }
        if (i < end) {
            uint64_t mask = reference.get_int(i, end - i);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, end - i))
                subindex->set_int(offset, _pext_u64(a, mask), popcount);
        }
    }
#endif
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


uint64_t footprint_sd_vector(uint64_t size, uint64_t num_set_bits) {
    uint8_t logn = sdsl::bits::hi(size) + 1;
    uint8_t logm = sdsl::bits::hi(num_set_bits) + 1;
    if (logm == logn)
        logm--;
    uint64_t low = num_set_bits * (logn - logm);
    uint64_t high = num_set_bits + (1ULL << logm);
    uint64_t high_select1 = footprint_select_support_mcl(high, num_set_bits);
    uint64_t high_select0 = footprint_select_support_mcl(high, high - num_set_bits);
    return low + high + high_select1 + high_select0;
}

uint64_t footprint_select_support_mcl(uint64_t size, uint64_t num_set_bits) {
    uint64_t sb = (num_set_bits + 4095) / 4096;
    double avg_diff = 1. * size / sb * 4095 / 4096;
    uint64_t blocks = ((sdsl::bits::hi((uint64_t)(avg_diff)) + 1) * 64 + 64 + 8)
                        * (num_set_bits / 4096);
    uint64_t miniblock_flags = 0;
    if (num_set_bits % 4096) {
        // the last block is large
        blocks += 4096 * (sdsl::bits::hi(size - 1) + 1) + 64 + 8;
        // bitmap with miniblock flags
        miniblock_flags += (sb + 63) / 64 * 64 + 64;
    } else {
        // empty bitmap with miniblock flags
        miniblock_flags += 64;
    }
    uint64_t offsets = (sb * (sdsl::bits::hi(size) + 1) + 63) / 64 * 64 + 64 + 8;
    return 64 + offsets + miniblock_flags + blocks;
}

uint64_t footprint_rank_support_v5(uint64_t size) {
    return ((((size + 63) >> 11) + 1) << 1) * 64 + 64; // ~0.062n
}
