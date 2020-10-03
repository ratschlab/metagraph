#include "benchmark/benchmark.h"

#include <string>
#include <vector>

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/annotation_converters.hpp"
#include "cli/query.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "seq_io/sequence_io.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;

using mtg::seq_io::kseq_t;


template <class Annotation>
std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    seq_io::read_fasta_file_critical(filename, [&](kseq_t *stream) {
        for (const auto &label : utils::split_string(stream->name.s, "|")) {
            sequences.emplace_back(stream->seq.s);
            labels.push_back(label);
        }
    }, true);

    size_t k = 12;

    boss::BOSSConstructor constructor(k - 1);
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new boss::BOSS(&constructor));
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<Annotation>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    return anno_graph;
}


std::unique_ptr<AnnotatedDBG> build_query_graph(const AnnotatedDBG &anno_graph,
                                                const std::string &query_filename) {
    return cli::construct_query_graph(
        anno_graph,
        [&](std::function<void(const std::string&)> call_sequence) {
            seq_io::read_fasta_file_critical(
                query_filename,
                [&](kseq_t *stream) { call_sequence(stream->seq.s); }
            );
        },
        1
    );
}


std::vector<std::string> queries {
    "../tests/data/transcripts_100.fa",
    "../tests/data/transcripts_1000.fa",
};


enum QueryMode : bool { NORMAL, FAST };


template <QueryMode fast>
static void BM_GetTopLabels(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto anno_graph = build_anno_graph<annot::RowCompressed<>>(graph_file);

    for (auto _ : state) {
        std::unique_ptr<AnnotatedDBG> query_graph;
        if constexpr(fast)
            query_graph = build_query_graph(*anno_graph, queries[state.range(0)]);

        const auto &graph = fast ? *query_graph : *anno_graph;

        seq_io::read_fasta_file_critical(
            queries[state.range(0)],
            [&](kseq_t *stream) { graph.get_top_labels(stream->seq.s, -1); }
        );
    }
}

BENCHMARK_TEMPLATE(BM_GetTopLabels, QueryMode::NORMAL)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);
BENCHMARK_TEMPLATE(BM_GetTopLabels, QueryMode::FAST)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);


template <QueryMode fast>
static void BM_GetTopLabelSignatures(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto anno_graph = build_anno_graph<annot::RowCompressed<>>(graph_file);

    for (auto _ : state) {
        std::unique_ptr<AnnotatedDBG> query_graph;
        if constexpr(fast)
            query_graph = build_query_graph(*anno_graph, queries[state.range(0)]);

        const auto &graph = fast ? *query_graph : *anno_graph;

        seq_io::read_fasta_file_critical(
            queries[state.range(0)],
            [&](kseq_t *stream) { graph.get_top_label_signatures(stream->seq.s, -1); }
        );
    }
}

BENCHMARK_TEMPLATE(BM_GetTopLabelSignatures, QueryMode::NORMAL)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);
BENCHMARK_TEMPLATE(BM_GetTopLabelSignatures, QueryMode::FAST)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);


template <size_t file_index>
static void BM_BRWTCompressTranscripts(benchmark::State& state) {
    auto anno_graph = build_anno_graph<annot::ColumnCompressed<>>(queries[file_index]);

    std::unique_ptr<annot::MultiBRWTAnnotator> annotator;

    size_t i = 0;
    for (auto _ : state) {
        if (i++)
            throw std::runtime_error("This benchmark will fail on the second iteration");

        const auto *column = dynamic_cast<const annot::ColumnCompressed<>*>(
            &anno_graph->get_annotation()
        );

        if (!column)
            throw std::runtime_error("This shouldn't happen");

        annotator = annot::convert_to_greedy_BRWT<annot::MultiBRWTAnnotator>(
            const_cast<annot::ColumnCompressed<>&&>(*column),
            state.range(0),
            state.range(0)
        );
    }
}

BENCHMARK_TEMPLATE(BM_BRWTCompressTranscripts, 0)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1)
    ->Arg(1)
    ->Arg(4);

} // namespace
