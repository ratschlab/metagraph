#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#include "common/elias_fano_file_merger.hpp"
#include "common/utils/file_utils.hpp"

// Note: if testing with many chunks use 'ulimit -n <max_files>' to increase the maxium
// number of files the system allows you to open. On mac the default is only 256!

namespace {

using namespace mtg;


constexpr size_t ITEM_COUNT = 10'000;
const std::string chunk_prefix = "/tmp/bm_chunk_";
std::vector<std::string> sources;

std::vector<std::string> create_sources(size_t num_sources) {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

    std::vector<std::string> sources;
    sources.reserve(num_sources);
    for (uint32_t i = 0; i < num_sources; ++i) {
        std::vector<uint64_t> els(ITEM_COUNT);
        for (uint64_t j = 0; j < ITEM_COUNT; ++j) {
            els[i] = j * 20 + dist10(rng);
        }
        sources.push_back(chunk_prefix + std::to_string(i));
        std::ofstream f(sources.back(), std::ios::binary);
        f.write(reinterpret_cast<char *>(els.data()), els.size() * sizeof(uint64_t));
        f.close();
    }
    return sources;
}

static void BM_merge_files(benchmark::State &state) {
    if (state.thread_index == 0) {
        sources = create_sources(state.range(0));
    }
    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    const auto file_writer = [&out](const uint64_t &v) {
        out.write(reinterpret_cast<const char *>(&v), sizeof(uint64_t));
    };
    bool do_cleanup = false;
    for (auto _ : state) {
        common::merge_files<uint64_t>(sources, file_writer, do_cleanup);
    }
    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });
}

static void BM_merge_files_pairs(benchmark::State &state) {
    if (state.thread_index == 0) {
        sources = create_sources(state.range(0));
    }
    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    using Pair = std::pair<uint64_t, uint8_t>;
    const auto file_writer = [&out](const Pair &v) {
        out.write(reinterpret_cast<const char *>(&v), sizeof(Pair));
    };
    bool do_cleanup = false;
    for (auto _ : state) {
        common::merge_files<uint64_t, uint8_t>(sources, file_writer, do_cleanup);
    }
    std::for_each(sources.begin(), sources.end(),
                  [](const std::string &s) { std::filesystem::remove(s); });
}

BENCHMARK(BM_merge_files)->DenseRange(10, 100, 10);
BENCHMARK(BM_merge_files_pairs)->DenseRange(10, 100, 10);

} // namespace
