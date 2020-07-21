#include <random>

#include <benchmark/benchmark.h>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/load/load_annotation.hpp"
#include "cli/config/config.hpp"


namespace {

using namespace mtg;


std::unique_ptr<annot::MultiLabelEncoded<std::string>> load_annotation() {
    if (!std::getenv("ANNO")) {
        std::cerr << "Set environment variable ANNO" << std::endl;
        exit(1);
    }

    auto annotation = cli::initialize_annotation(std::getenv("ANNO"));

    if (!annotation->load(std::getenv("ANNO"))) {
        std::cerr << "Can't load annotation from "
                  << std::getenv("ANNO") << std::endl;
        exit(1);
    }

    return annotation;
}

// generate a deterministic sequence of pseudo-random numbers
std::vector<uint64_t> random_numbers(size_t size, uint64_t min, uint64_t max) {
    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(min, max);

    std::vector<uint64_t> numbers(size);
    for (uint64_t &number : numbers) {
        number = dis(gen);
    }
    return numbers;
}


static void BM_anno_get_row(benchmark::State &state) {
    auto anno = load_annotation();
    const annot::binmat::BinaryMatrix &annotation_matrix = anno->get_matrix();

    auto rows = random_numbers(100'000, 0, annotation_matrix.num_rows() - 1);

    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(annotation_matrix.get_row(rows[i++ % rows.size()]));
    }
}
BENCHMARK(BM_anno_get_row) -> Unit(benchmark::kMicrosecond);

static void BM_anno_get_rows(benchmark::State &state) {
    auto anno = load_annotation();
    const annot::binmat::BinaryMatrix &annotation_matrix = anno->get_matrix();

    auto rows = random_numbers(100'000, 0, annotation_matrix.num_rows() - 1);
    std::sort(rows.begin(), rows.end());

    for (auto _ : state) {
        benchmark::DoNotOptimize(annotation_matrix.get_rows(rows));
    }
}
BENCHMARK(BM_anno_get_rows) -> Unit(benchmark::kMillisecond);

static void BM_anno_get_rows_unique(benchmark::State &state) {
    auto anno = load_annotation();
    if (!dynamic_cast<const annot::binmat::RainbowMatrix*>(&anno->get_matrix())) {
        state.SkipWithError("This is not a Rainbow type of matrix. Skipped.");
        return;
    }
    const annot::binmat::RainbowMatrix &rb_matrix
        = dynamic_cast<const annot::binmat::RainbowMatrix&>(anno->get_matrix());

    auto rows = random_numbers(100'000, 0, rb_matrix.num_rows() - 1);
    std::sort(rows.begin(), rows.end());

    for (auto _ : state) {
        benchmark::DoNotOptimize(rb_matrix.get_rows(&rows));
    }
}
BENCHMARK(BM_anno_get_rows_unique) -> Unit(benchmark::kMillisecond);

} // namespace
