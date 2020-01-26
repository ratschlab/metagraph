#include "benchmark/benchmark.h"

#include "benchmark_graph_helpers.hpp"

#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"


static void BM_GetTopLabels(benchmark::State& state) {
    std::string query_file = "../tests/data/transcripts_100.fa";
    auto anno_graph = build_anno_graph("../tests/data/transcripts_1000.fa");

    for (auto _ : state) {
        read_fasta_file_critical(
            query_file,
            [&](kseq_t *stream) {
                anno_graph->get_top_labels(std::string(stream->seq.s), -1);
            }
        );
    }
}

BENCHMARK(BM_GetTopLabels)->Unit(benchmark::kMillisecond);
