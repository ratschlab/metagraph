#include "benchmark/benchmark.h"

#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/data_generation.hpp"


namespace {

template <uint64_t size, uint8_t density_percent>
static void BM_count_ones(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    sdsl::bit_vector bv
        = generator.generate_random_column(1'000'000, density_percent / 100.);

    for (auto _ : state) {
        count_ones(bv, 640, 640 + size);
    }
}

BENCHMARK_TEMPLATE(BM_count_ones, 0, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 1, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 2, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 3, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 5, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 10, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 64, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 100, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 200, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 1000, 50);
BENCHMARK_TEMPLATE(BM_count_ones, 10000, 50);



template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_or_call_ones(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    sdsl::bit_vector result(bv.size(), false);

    for (auto _ : state) {
        bv.call_ones_in_range(0, bv.size(),
            [&](uint64_t k) { result[k] = true; }
        );
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_stat, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_stat, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_stat, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_stat, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_stat, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_sd, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_rrr<>, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_rrr<>, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_rrr<>, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_rrr<>, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_call_ones, bit_vector_rrr<>, 100000, 50) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_or_get_int(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    sdsl::bit_vector result(bv.size(), false);

    for (auto _ : state) {
        uint64_t begin = 0;
        uint64_t end = bv.size();

        uint64_t i;
        for (i = begin; i + 64 <= end; i += 64) {
            result.data()[i / 64] |= bv.get_int(i, 64);
        }
        if (i < end)
            result.data()[i / 64] |= bv.get_int(i, end - i);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_stat, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_stat, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_stat, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_stat, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_stat, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_sd, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_rrr<>, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_rrr<>, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_rrr<>, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_rrr<>, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_or_get_int, bit_vector_rrr<>, 100000, 50) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_subvector_via_call_ones(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    sdsl::bit_vector bv
        = generator.generate_random_column(size, density_percent / 100.);

    sdsl::bit_vector reference
        = generator.generate_random_column(size, density_percent / 100.);

    reference |= bv;

    bit_vector_type vector(std::move(bv));

    sdsl::bit_vector subindex(sdsl::util::cnt_one_bits(reference), false);

    for (auto _ : state) {
        uint64_t rank = 0;
        uint64_t i = 0;

        vector.call_ones([&](uint64_t j) {
            assert(j >= i);
            assert(reference[j]);

            rank += count_ones(reference, i, j);

            assert(rank < subindex.size());

            subindex[rank++] = true;
            i = j + 1;
        });
    }
}

BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_stat, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_stat, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_stat, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_stat, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_stat, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_sd, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_sd, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_sd, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_sd, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_sd, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_rrr<>, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_rrr<>, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_rrr<>, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_rrr<>, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_call_ones, bit_vector_rrr<>, 100000, 50) -> Unit(benchmark::kMicrosecond);


#if __BMI2__

template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_subvector_via_pext(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    sdsl::bit_vector bv
        = generator.generate_random_column(size, density_percent / 100.);

    sdsl::bit_vector reference
        = generator.generate_random_column(size, density_percent / 100.);

    reference |= bv;

    bit_vector_type vector(std::move(bv));

    sdsl::bit_vector subindex(sdsl::util::cnt_one_bits(reference), false);

    for (auto _ : state) {
        uint64_t rank = 0;
        uint64_t i = 0;

        uint64_t end = reference.size();

        for (i = 0; i + 64 <= end; i += 64) {
            uint64_t mask = reference.data()[i / 64];
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = vector.get_int(i, 64))
                subindex.set_int(rank, _pext_u64(a, mask), popcount);
            rank += popcount;
        }
        if (i < end) {
            uint64_t mask = reference.get_int(i, end - i);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = vector.get_int(i, end - i))
                subindex.set_int(rank, _pext_u64(a, mask), popcount);
        }
    }
}

BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_stat, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_stat, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_stat, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_stat, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_stat, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_sd, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_sd, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_sd, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_sd, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_sd, 100000, 50) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_rrr<>, 100000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_rrr<>, 100000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_rrr<>, 100000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_rrr<>, 100000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_subvector_via_pext, bit_vector_rrr<>, 100000, 50) -> Unit(benchmark::kMicrosecond);

#endif // __BMI2__

} // namespace
