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
        generator.generate_random_column(size, density_percent / 100.)->to_vector()
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv.rank1((i++ * 87'178'291'199) % size);
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_rank, bit_vector_rrr<>, 10000000, 99) -> Unit(benchmark::kMicrosecond);


template <class bit_vector_type, uint64_t size, uint8_t density_percent>
static void BM_bit_vector_query_access(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    bit_vector_type bv(
        generator.generate_random_column(size, density_percent / 100.)->to_vector()
    );

    uint64_t i = 0;
    for (auto _ : state) {
        bv[(i++ * 87'178'291'199) % size];
    }
}

BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 5) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 9) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 10) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 50) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 80) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_bit_vector_query_access, bit_vector_rrr<>, 10000000, 99) -> Unit(benchmark::kMicrosecond);

} // namespace bm
} // namespace mg
