#include <cmath>

#include <benchmark/benchmark.h>

#include "common/threads/threading.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "../benchmark_graph_helpers.hpp"

const std::string file_prefix = "/tmp/bm_mg_outfile.fasta.gz";


static void BM_ContigsWrite(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto graph = mg::bm::build_graph(graph_file);
    std::vector<std::string> contigs;
    graph->call_sequences([&](const auto &contig, auto&&) { contigs.push_back(contig); });

    set_num_threads(state.range(0));
    for (auto _ : state) {
        FastaWriter writer(file_prefix, "", false);
        for (const auto &contig : contigs) {
            writer.write(contig);
        }
    }
    set_num_threads(1);
}

BENCHMARK(BM_ContigsWrite)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);
