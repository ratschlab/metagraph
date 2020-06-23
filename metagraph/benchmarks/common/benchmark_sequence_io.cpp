#include <cmath>
#include <random>

#include <benchmark/benchmark.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "common/threads/threading.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"


namespace {

using namespace mtg;

const std::string file_prefix = "/tmp/bm_mg_outfile.fasta.gz";


std::shared_ptr<DeBruijnGraph> build_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    seq_io::read_fasta_file_critical(filename,
                                     [&](seq_io::kseq_t *stream) {
                                         sequences.emplace_back(stream->seq.s);
                                     },
                                     true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(std::move(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor));
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    return graph;
}


template <size_t max_size>
static void BM_WriteRandomSequences(benchmark::State& state) {
    const std::string alphabet = "ATGCN";
    size_t count = 0;
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist4(0, 4);
    std::uniform_int_distribution<std::mt19937::result_type> dist1000(10, 1000);

    std::vector<std::string> sequences;
    while (count < max_size) {
        sequences.emplace_back(dist1000(rng), 'A');
        count += sequences.back().size() + sizeof(std::string);
        std::for_each(sequences.back().begin(),
                      sequences.back().end(),
                      [&](char &c) { c = alphabet[dist4(rng)]; });
    }

    set_num_threads(state.range(0));
    for (auto _ : state) {
        seq_io::FastaWriter writer(file_prefix, "", false, state.range(0) - 1);
        for (const auto &sequence : sequences) {
            writer.write(sequence);
        }
    }
    set_num_threads(1);
}

BENCHMARK_TEMPLATE(BM_WriteRandomSequences, 20000000)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);

BENCHMARK_TEMPLATE(BM_WriteRandomSequences, 200000000)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);


static void BM_WriteContigs(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto graph = build_graph(graph_file);
    std::vector<std::string> contigs;
    graph->call_sequences([&](const auto &contig, auto&&) { contigs.push_back(contig); });

    set_num_threads(state.range(0));
    for (auto _ : state) {
        seq_io::FastaWriter writer(file_prefix, "", false, state.range(0) - 1);
        for (const auto &contig : contigs) {
            writer.write(contig);
        }
    }
    set_num_threads(1);
}

BENCHMARK(BM_WriteContigs)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);

} // namespace
