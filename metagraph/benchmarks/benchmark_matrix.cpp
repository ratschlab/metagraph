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

static void BM_BRWTCompress(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    for (auto _ : state) {
        std::vector<std::unique_ptr<bit_vector>> generated_columns;
        generated_columns = generator.generate_random_columns(
            1'000'000,
            100,
            get_densities(100, std::vector<double>(100, 0.01)),
            std::vector<uint32_t>(100, 1'000 / 100)
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

BENCHMARK(BM_BRWTCompress)->Unit(benchmark::kMillisecond);
