#include "benchmark/benchmark.h"

#include <sdsl/sd_vector.hpp>

#include "common/unix_tools.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/sd_vector_disk.hpp"
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


template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_sequential_sd_vector_access(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector<> bv(gen.generate_random_column(size, density_percent / 100.));

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_sequential_sd_vector_disk_access(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector_disk<> bv;
    {
        auto vec = gen.generate_random_column(size, density_percent / 100.);
        const std::string fname = "../tests/data/bit_vector_dump_test";
        sdsl::sd_vector_disk_builder builder(vec.size(), sdsl::util::cnt_one_bits(vec), fname);
        call_ones(vec, [&](uint64_t i) { builder.set(i); });
        bv = sdsl::sd_vector_disk<>(builder);
    }

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_random_sd_vector_access(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector<> bv(gen.generate_random_column(size, density_percent / 100.));

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++ * 87'178'291'199) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_random_sd_vector_disk_access(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector_disk<> bv;
    {
        auto vec = gen.generate_random_column(size, density_percent / 100.);
        const std::string fname = "../tests/data/bit_vector_dump_test";
        sdsl::sd_vector_disk_builder builder(vec.size(), sdsl::util::cnt_one_bits(vec), fname);
        call_ones(vec, [&](uint64_t i) { builder.set(i); });
        bv = sdsl::sd_vector_disk<>(builder);
    }

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++ * 87'178'291'199) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_sequential_sd_vector_rank(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector<> bv(gen.generate_random_column(size, density_percent / 100.));
    sdsl::sd_vector<>::rank_1_type rank(&bv);

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(rank((1 + i++) % bv.size()));
    }

    state.counters["RAM"] = get_curr_RSS();
}

template <uint64_t size, uint8_t density_percent>
static void BM_bv_query_sequential_sd_vector_disk_rank(benchmark::State& state) {
    DataGenerator gen;
    gen.set_seed(42);

    sdsl::sd_vector_disk<> bv;
    {
        auto vec = gen.generate_random_column(size, density_percent / 100.);
        const std::string fname = "../tests/data/bit_vector_dump_test";
        sdsl::sd_vector_disk_builder builder(vec.size(), sdsl::util::cnt_one_bits(vec), fname);
        call_ones(vec, [&](uint64_t i) { builder.set(i); });
        bv = sdsl::sd_vector_disk<>(builder);
    }
    sdsl::sd_vector_disk<>::rank_1_type rank(&bv);

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(rank((1 + i++) % bv.size()));
    }

    state.counters["RAM"] = get_curr_RSS();
}


#define INST_BV_BENCHMARK(NAME) \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 0) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 1) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 2) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 5) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 10) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 20) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 30) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 50) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 70) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 80) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 90) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 95) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 98) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 99) -> Unit(benchmark::kMicrosecond); \
BENCHMARK_TEMPLATE(BM_bv_query_##NAME, 100'000'000, 100) -> Unit(benchmark::kMicrosecond); \

INST_BV_BENCHMARK(sequential_sd_vector_disk_access);
INST_BV_BENCHMARK(sequential_sd_vector_access);

INST_BV_BENCHMARK(random_sd_vector_disk_access);
INST_BV_BENCHMARK(random_sd_vector_access);

INST_BV_BENCHMARK(sequential_sd_vector_disk_rank);
INST_BV_BENCHMARK(sequential_sd_vector_rank);


template <uint64_t t>
static void BM_bv_query_random_sd_vector_access_every_nth_bit_set(benchmark::State& state) {
    const uint64_t size = 1e12;
    sdsl::sd_vector_builder builder(size, (size + t - 1) / t);
    for (size_t i = 0; i < size; i += t) {
        builder.set(i);
    }
    sdsl::sd_vector<> bv(builder);

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++ * 87'178'291'199) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, 10000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, 3000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, 1000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, 500) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_access_every_nth_bit_set, 200) -> Unit(benchmark::kMicrosecond);

template <uint64_t t>
static void BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set(benchmark::State& state) {
    const uint64_t size = 1e12;

    sdsl::sd_vector_disk<> bv;
    {
        const std::string fname = "../tests/data/bit_vector_dump_test";
        sdsl::sd_vector_disk_builder builder(size, (size + t - 1) / t, fname);
        for (size_t i = 0; i < size; i += t) {
            builder.set(i);
        }
        bv = sdsl::sd_vector_disk<>(builder);
    }

    uint64_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bv[(i++ * 87'178'291'199) % bv.size()]);
    }

    state.counters["RAM"] = get_curr_RSS();
}

BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 10000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 3000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 1000) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 500) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 200) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bv_query_random_sd_vector_disk_access_every_nth_bit_set, 100) -> Unit(benchmark::kMicrosecond);

} // namespace
