#include "benchmark/benchmark.h"

#include <string>
#include <vector>

#include "benchmark_graph_helpers.hpp"

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"


namespace BENCHMARK_QUERY {

std::vector<std::string> queries {
    "../tests/data/transcripts_100.fa"
};

enum QueryMode : bool { NORMAL, FAST };

template <QueryMode fast>
static void BM_GetTopLabels(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto anno_graph = build_anno_graph<annotate::RowCompressed<>>(graph_file);

    for (auto _ : state) {
        read_fasta_file_critical(
            queries[state.range(0)],
            [&](kseq_t *stream) {
                std::string query(stream->seq.s);
                std::unique_ptr<AnnotatedDBG> query_graph;
                if constexpr(fast)
                    query_graph = build_query_graph(*anno_graph, query);

                const auto &graph = fast ? *query_graph : *anno_graph;
                graph.get_top_labels(query, -1);
            }
        );
    }
}

BENCHMARK_TEMPLATE(BM_GetTopLabels, QueryMode::NORMAL)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);
BENCHMARK_TEMPLATE(BM_GetTopLabels, QueryMode::FAST)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, queries.size() - 1, 1);


} // BENCHMARK_QUERY
