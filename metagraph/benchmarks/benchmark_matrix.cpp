#include "benchmark/benchmark.h"

#include <vector>

#include "method_constructors.hpp"

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/annotation_converters.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "seq_io/sequence_io.hpp"
#include "common/algorithms.hpp"
#include "common/utils/string_utils.hpp"


std::vector<double> get_densities(uint64_t num_cols, const std::vector<double> &vector) {
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
          size_t unique_arg = 100,
          size_t arity_arg = 2,
          bool greedy_arg = true,
          size_t relax_arg = 2>
static void BM_BRWTCompressSparse(benchmark::State& state) {
    DataGenerator generator;
    generator.set_seed(42);

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

    for (auto _ : state) {
        matrix = generate_brwt_from_rows(
            std::move(generated_columns),
            arity_arg,
            greedy_arg,
            relax_arg
        );
    }
}

BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 10)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 100)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_BRWTCompressSparse, 1, 1000)->Unit(benchmark::kMillisecond);


std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    read_fasta_file_critical(filename,
                             [&](kseq_t *stream) {
                                 std::string name(stream->name.s);
                                 for (const auto &label : utils::split_string(name, "|")) {
                                     sequences.push_back(stream->seq.s);
                                     labels.push_back(label);
                                 }
                             },
                             true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(sequences);
    std::shared_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(&constructor)) };
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annotate::ColumnCompressed<>>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    return anno_graph;
}

static void BM_BRWTCompressTranscripts(benchmark::State& state) {
    auto anno_graph = build_anno_graph("../tests/data/transcripts_1000.fa");
    const auto &annotation = anno_graph->get_annotation();
    std::cout << annotation.num_objects() << " " << annotation.num_labels() << std::endl;

    std::unique_ptr<annotate::MultiBRWTAnnotator> annotator;
    for (auto _ : state) {
        const auto *column = dynamic_cast<const annotate::ColumnCompressed<>*>(
            &anno_graph->get_annotation()
        );

        if (!column)
            throw std::runtime_error("This shouldn't happen");

        utils::set_verbose(true);
        annotator = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
            const_cast<annotate::ColumnCompressed<>&&>(*column),
            state.range(0),
            state.range(0)
        );
    }
}

BENCHMARK(BM_BRWTCompressTranscripts)->Unit(benchmark::kMillisecond)
                                     ->Iterations(1)
                                     ->Arg(1)
                                     ->Arg(4);
