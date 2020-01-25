#include "benchmark/benchmark.h"

#include <vector>

#include "method_constructors.hpp"

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"


std::vector<double> get_densities(uint64_t num_cols,
                                  const std::vector<double> &vector) {
    auto densities = vector;
    if (densities.size() == 1) {
        densities.assign(num_cols, densities[0]);
    } else if (densities.size() != num_cols) {
        std::cout << "ERROR: wrong number of column counts" << std::endl;
        exit(1);
    }
    return densities;
}

template <size_t density_numerator,
          size_t density_denominator,
          size_t rows_arg = 1000000,
          size_t cols_arg = 1000,
          size_t unique_arg = 100>
static void BM_BRWTCompressSparse(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    for (auto _ : state) {
        auto density_arg = std::vector<double>(
            unique_arg, double(density_numerator) / double(density_denominator)
        );
        std::vector<std::unique_ptr<bit_vector>> generated_columns;
        generated_columns = generator.generate_random_columns(
            rows_arg,
            unique_arg,
            get_densities(unique_arg, density_arg),
            std::vector<uint32_t>(unique_arg, cols_arg / unique_arg)
        );

        std::unique_ptr<BinaryMatrix> matrix;

        matrix = generate_brwt_from_rows(
            std::move(generated_columns),
            2,
            true,
            2
        );
    }
}

BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 10)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 100)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 1000)->Unit(benchmark::kMillisecond);
