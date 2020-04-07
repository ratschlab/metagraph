#include <cmath>

#include <benchmark/benchmark.h>

#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "../benchmark_graph_helpers.hpp"

const std::string file_prefix = "/tmp/bm_mg_outfile.fasta.gz";


static void BM_UnitigsWrite(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto graph = mg::bm::build_graph(graph_file);
    std::vector<std::string> unitigs;
    graph->call_unitigs([&](const auto &unitig, auto&&) { unitigs.push_back(unitig); });

    for (auto _ : state) {
        FastaWriter writer(file_prefix, "", false);
        for (const auto &unitig : unitigs) {
            writer.write(unitig);
        }
    }
}

BENCHMARK(BM_UnitigsWrite)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(0, 10, 2);

static void BM_UnitigsExtractAndWrite(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto graph = mg::bm::build_graph(graph_file);

    for (auto _ : state) {
        FastaWriter writer(file_prefix, "", false);
        graph->call_unitigs([&](const auto &unitig, auto&&) { writer.write(unitig); });
    }
}

BENCHMARK(BM_UnitigsExtractAndWrite)
    ->Unit(benchmark::kMillisecond);
