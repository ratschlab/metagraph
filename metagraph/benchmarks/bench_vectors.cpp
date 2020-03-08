#include "benchmark/benchmark.h"

#include "common/vectors/bit_vector.hpp"
#include "common/data_generation.hpp"


namespace mg {
namespace bm {

template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_query_rank(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv.rank1((i++ * 87'178'291'199) % size);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 15) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 40) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 60) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 70) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 90) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 95) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 98) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<63>, 10000000, 99) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_sequential_query_rank(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv.rank1((i++) % size);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 15) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 40) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 60) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 70) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 90) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 95) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 98) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_sequential_query_rank, bit_vector_rrr<63>, 10000000, 99) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_query_access(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv[(i++ * 87'178'291'199) % size];
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 15) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 40) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 60) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 70) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 90) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 95) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 98) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<63>, 10000000, 99) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_query_inverse_select(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv.inverse_select((i++ * 87'178'291'199) % size);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 15) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 40) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 60) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 70) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 90) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 95) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 98) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_inverse_select, bit_vector_rrr<63>, 10000000, 99) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_query_select(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)
    );

    uint64_t i = 0;
    uint64_t num_set_bits = bv.num_set_bits();
    for (auto _ : state) {
        bv.select1((i++ * 87'178'291'199) % num_set_bits);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 15) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 20) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 30) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 40) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 60) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 70) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 90) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 95) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 98) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_select, bit_vector_rrr<63>, 10000000, 99) -> Unit(benchmark::kMicrosecond);

} // namespace bm
} // namespace mg
