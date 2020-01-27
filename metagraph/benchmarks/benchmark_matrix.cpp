#include "benchmark/benchmark.h"

#include "benchmark_graph_helpers.hpp"
#include "method_constructors.hpp"

#include "annotation/annotation_converters.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/annotated_dbg.hpp"


std::vector<double> get_densities(uint64_t num_cols, std::vector<double> densities) {
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
          size_t unique_arg = 100,
          size_t arity_arg = 2,
          bool greedy_arg = true,
          size_t relax_arg = 2>
static void BM_BRWTCompressSparse(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

    auto density_arg = std::vector<double>(
        unique_arg,
        static_cast<double>(density_numerator) / density_denominator
    );
    std::vector<std::unique_ptr<bit_vector>> generated_columns;
    generated_columns = generator.generate_random_columns(
        rows_arg,
        unique_arg,
        get_densities(unique_arg, density_arg),
        std::vector<uint32_t>(unique_arg, cols_arg / unique_arg)
    );

    std::unique_ptr<BinaryMatrix> matrix;

    size_t i = 0;
    for (auto _ : state) {
        if (i)
            throw std::runtime_error("This benchmark will fail on the second iteration");

        matrix = generate_brwt_from_rows(
            std::move(generated_columns),
            arity_arg,
            greedy_arg,
            relax_arg
        );

        ++i;
    }
}

BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 10)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 100)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 1000)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);




static void BM_BRWTCompressTranscripts(benchmark::State& state) {
    auto anno_graph = build_anno_graph("../tests/data/transcripts_1000.fa");

    std::unique_ptr<annotate::MultiBRWTAnnotator> annotator;

    size_t i = 0;
    for (auto _ : state) {
        if (i)
            throw std::runtime_error("This benchmark will fail on the second iteration");

        const auto *column = dynamic_cast<const annotate::ColumnCompressed<>*>(
            &anno_graph->get_annotation()
        );

        if (!column)
            throw std::runtime_error("This shouldn't happen");

        annotator = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
            const_cast<annotate::ColumnCompressed<>&&>(*column),
            state.range(0),
            state.range(0)
        );

        ++i;
    }
}

BENCHMARK(BM_BRWTCompressTranscripts)->Unit(benchmark::kMillisecond)
                                     ->Iterations(1)
                                     ->Arg(1)
                                     ->Arg(4);
