#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#include "common/file_merger.hpp"

constexpr size_t CHUNK_COUNT = 30;
constexpr size_t ITEM_COUNT = 100'000;
const std::string chunk_prefix = "/tmp/chunk_";

void create_sources() {
    std::mt19937 rng(123456);
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);


    for (uint32_t i = 0; i < CHUNK_COUNT; ++i) {
        std::vector<uint32_t> els(ITEM_COUNT);
        for (uint64_t j = 0; j < ITEM_COUNT; ++j) {
            els[i] = j + dist6(rng);
        }
        std::ofstream f(chunk_prefix + std::to_string(i), std::ios::binary);
        f.write(reinterpret_cast<char *>(els.data()), els.size() * sizeof(uint32_t));
        f.close();
    }
}

static void BM_merge_files(benchmark::State &state) {
    create_sources();
    std::vector<std::string> sources(CHUNK_COUNT);
    for (uint32_t i = 0; i < CHUNK_COUNT; ++i) {
        sources[i] = chunk_prefix + std::to_string(i);
    }
    std::ofstream out;
    char buffer[1024*1024];
    out.rdbuf()->pubsetbuf(buffer, 1024 * 1024);
    out.open("/tmp/out");
    const auto file_writer = [&out](const uint64_t &v) {
        out.write(reinterpret_cast<const char *>(&v), sizeof(uint64_t));
    };
    for (auto _ : state) {
        mg::common::merge_files<uint64_t>(sources, file_writer);
    }
}

BENCHMARK(BM_merge_files);

BENCHMARK_MAIN();
