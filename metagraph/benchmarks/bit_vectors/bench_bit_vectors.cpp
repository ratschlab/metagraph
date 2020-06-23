#include "benchmark/benchmark.h"

#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/data_generation.hpp"


namespace {

#define DEFINE_BV_QUERY_BENCHMARK(NAME, OPERATION, MIN_INDEX, RANGE) \
template <class bit_vector_type, uint64_t size, uint8_t density_percent> \
static void BM_bv_query_##NAME(benchmark::State& state) { \
    DataGenerator gen; \
    gen.set_seed(42); \
 \
    bit_vector_type bv(gen.generate_random_column(size, density_percent / 100.)); \
    uint64_t max_index = bv.RANGE(); \
    if (max_index <= MIN_INDEX) \
        return; \
 \
    uint64_t i = 0; \
    for (auto _ : state) { \
        benchmark::DoNotOptimize(bv.OPERATION(MIN_INDEX + (i++ * 87'178'291'199) % (max_index - MIN_INDEX))); \
    } \
    state.counters["Density"] \
        = static_cast<double>(bv.num_set_bits()) / bv.size(); \
} \
 \
template <class bit_vector_type, uint64_t size, uint8_t density_percent> \
static void BM_bv_query_sequential_##NAME(benchmark::State& state) { \
    DataGenerator gen; \
    gen.set_seed(42); \
 \
    bit_vector_type bv(gen.generate_random_column(size, density_percent / 100.)); \
    uint64_t max_index = bv.RANGE(); \
    if (max_index <= MIN_INDEX) \
        return; \
 \
    uint64_t i = 0; \
    for (auto _ : state) { \
        benchmark::DoNotOptimize(bv.OPERATION(MIN_INDEX + (i++) % (max_index - MIN_INDEX))); \
    } \
    state.counters["Density"] \
        = static_cast<double>(bv.num_set_bits()) / bv.size(); \
}

DEFINE_BV_QUERY_BENCHMARK(rank,           rank1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(access,         operator[],     0, size);
DEFINE_BV_QUERY_BENCHMARK(inverse_select, inverse_select, 0, size);
DEFINE_BV_QUERY_BENCHMARK(select,         select1,        1, num_set_bits);
DEFINE_BV_QUERY_BENCHMARK(next,           next1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(prev,           prev1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(cond_rank,   conditional_rank1, 0, size);

#define INST_BV_QUERY_BENCHMARK(BV_TYPE, NAME) \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 0) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 1) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 2) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 5) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 9) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 10) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 15) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 20) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 30) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 40) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 50) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 60) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 70) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 80) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 90) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 95) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 98) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 99) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 100) -> Unit(benchmark::kMicrosecond); \
\
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 0) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 1) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 2) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 5) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 9) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 10) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 15) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 20) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 30) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 40) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 50) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 60) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 70) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 80) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 90) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 95) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 98) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 99) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 100) -> Unit(benchmark::kMicrosecond);

INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, rank);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, access);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, inverse_select);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, select);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, next);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, prev);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, cond_rank);

} // namespace
