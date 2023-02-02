#include "benchmark/benchmark.h"

#include <sdsl/sd_vector.hpp>

#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/sd_vector_builder_disk.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/data_generation.hpp"


namespace {

#define DEFINE_BV_QUERY_BENCHMARK(NAME, OPERATION, MIN_INDEX, RANGE) \
template <class bit_vector_type, uint64_t size, unsigned int density_promille> \
static void BM_bv_query_##NAME(benchmark::State& state) { \
    DataGenerator gen; \
    gen.set_seed(42); \
 \
    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.)); \
    uint64_t max_index = bv.RANGE(); \
    if (max_index <= MIN_INDEX) \
        return; \
 \
    uint64_t i = 0; \
    uint64_t result = 0; \
    for (auto _ : state) { \
        result += utils::get_first(bv.OPERATION(MIN_INDEX + (i++ * 87'178'291'199) % (max_index - MIN_INDEX))); \
    } \
    benchmark::DoNotOptimize(result); \
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size(); \
} \
 \
template <class bit_vector_type, uint64_t size, unsigned int density_promille> \
static void BM_bv_query_sequential_##NAME(benchmark::State& state) { \
    DataGenerator gen; \
    gen.set_seed(42); \
 \
    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.)); \
    uint64_t max_index = bv.RANGE(); \
    if (max_index <= MIN_INDEX) \
        return; \
 \
    uint64_t i = 0; \
    uint64_t result = 0; \
    for (auto _ : state) { \
        result += utils::get_first(bv.OPERATION(MIN_INDEX + (i++) % (max_index - MIN_INDEX))); \
    } \
    benchmark::DoNotOptimize(result); \
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size(); \
}

DEFINE_BV_QUERY_BENCHMARK(rank,           rank1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(access,         operator[],     0, size);
DEFINE_BV_QUERY_BENCHMARK(inverse_select, inverse_select, 0, size);
DEFINE_BV_QUERY_BENCHMARK(select,         select1,        1, num_set_bits);
DEFINE_BV_QUERY_BENCHMARK(next,           next1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(prev,           prev1,          0, size);
DEFINE_BV_QUERY_BENCHMARK(cond_rank,   conditional_rank1, 0, size);

template <class bit_vector_type, uint64_t size, unsigned int density_promille>
static void BM_bv_query_get_int(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.));
    const uint64_t max_index = bv.size() - 64;

    uint64_t i = 0;
    uint64_t result = 0;
    for (auto _ : state) {
        result += bv.get_int((i++ * 87'178'291'199) % max_index, 64);
    }
    benchmark::DoNotOptimize(result); \
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size();
}

template <class bit_vector_type, uint64_t size, unsigned int density_promille>
static void BM_bv_query_sequential_get_int(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.));
    const uint64_t max_index = bv.size() - 64;

    uint64_t i = 0;
    uint64_t result = 0;
    for (auto _ : state) {
        result += bv.get_int(i++ % max_index, 64);
    }
    benchmark::DoNotOptimize(result); \
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size();
}

template <class bit_vector_type, uint64_t size, unsigned int density_promille>
static void BM_bv_query_call_ones(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.));
    const size_t BS = 50'000;
    const uint64_t max_index = bv.size() - BS;

    uint64_t i = 0;
    uint64_t result = 0;
    for (auto _ : state) {
        uint64_t begin = (i++ * 87'178'291'199) % max_index;
        bv.call_ones_in_range(begin, begin + BS, [&](uint64_t x) { result += x; });
    }
    benchmark::DoNotOptimize(result);
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size();
}

template <class bit_vector_type, uint64_t size, unsigned int density_promille>
static void BM_bv_query_sequential_call_ones(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    bit_vector_type bv(gen.generate_random_column(size, density_promille / 1000.));
    const size_t BS = 50'000;
    const uint64_t max_index = bv.size() - BS;

    uint64_t i = 0;
    uint64_t result = 0;
    for (auto _ : state) {
        uint64_t begin = i % max_index;
        i += BS;
        bv.call_ones_in_range(begin, begin + BS, [&](uint64_t x) { result += x; });
    }
    benchmark::DoNotOptimize(result);
    state.counters["Density"] = static_cast<double>(bv.num_set_bits()) / bv.size();
}

#define INST_BV_QUERY_BENCHMARK(BV_TYPE, NAME) \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 0) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 5) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 10) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 20) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 50) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 90) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 100) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 150) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 200) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 300) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 400) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 500) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 600) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 700) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 800) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 900) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 950) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 980) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 990) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 995) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, BV_TYPE, 10000000, 1000) -> Unit(benchmark::kMicrosecond); \
\
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 0) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 5) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 10) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 20) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 50) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 90) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 100) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 150) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 200) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 300) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 400) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 500) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 600) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 700) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 800) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 900) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 950) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 980) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 990) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 995) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_sequential_##NAME, BV_TYPE, 10000000, 1000) -> Unit(benchmark::kMicrosecond);

INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, rank);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, access);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, inverse_select);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, select);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, next);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, prev);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, cond_rank);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, get_int);
INST_BV_QUERY_BENCHMARK(bit_vector_rrr<63>, call_ones);


template <class t_ifstream, uint64_t t>
static void BM_bv_query_random_sd_vector_access_every_nth_bit_set(benchmark::State& state) {
    const uint64_t size = 1e12;

    sdsl::sd_vector<> bv;
    {
        const std::string fname = "../tests/data/bit_vector_dump_test";
        sdsl::sd_vector_builder_disk<> builder(size, (size + t - 1) / t, fname);
        for (size_t i = 0; i < size; i += t) {
            builder.set(i);
        }
        builder.finish();
        t_ifstream in(fname);
        bv.load(in);
    }

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++ * 87'178'291'199) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 10000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 3000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 1000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 500) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 200) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, std::ifstream, 100) -> Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 10000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 3000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 1000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 500) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 200) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, sdsl::mmap_ifstream, 100) -> Unit(benchmark::kMicrosecond);

} // namespace
