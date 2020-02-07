#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>
#include <common/utils/file_utils.hpp>

#include "common/file_merger.hpp"

constexpr size_t ITEM_COUNT = 100'000;
const std::string chunk_prefix = "/tmp/chunk_";

void create_sources(int count) {
    std::mt19937 rng(123456);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

    for (uint32_t i = 0; i < count; ++i) {
        std::vector<uint32_t> els(ITEM_COUNT);
        for (uint64_t j = 0; j < ITEM_COUNT; ++j) {
            els[i] = j + dist10(rng);
        }
        std::ofstream f(chunk_prefix + std::to_string(i), std::ios::binary);
        f.write(reinterpret_cast<char *>(els.data()), els.size() * sizeof(uint32_t));
        f.close();
    }
}

static void BM_merge_files(benchmark::State &state) {
    create_sources(state.range(0));
    std::vector<std::string> sources(state.range(0));
    for (uint32_t i = 0; i < state.range(0); ++i) {
        sources[i] = chunk_prefix + std::to_string(i);
    }

    utils::TempFile tempfile("/tmp/");
    std::ofstream& out = tempfile.ofstream();
    const auto file_writer = [&out](const uint64_t &v) {
        out.write(reinterpret_cast<const char *>(&v), sizeof(uint64_t));
    };
    for (auto _ : state) {
        mg::common::merge_files<uint64_t>(sources, file_writer);
    }
}

static void BM_merge_files_pairs(benchmark::State &state) {
    create_sources(state.range(0));
    std::vector<std::string> sources(state.range(0));
    for (uint32_t i = 0; i < state.range(0); ++i) {
        sources[i] = chunk_prefix + std::to_string(i);
    }
    utils::TempFile tempfile("/tmp/");
    std::ofstream& out = tempfile.ofstream();
    using Pair = std::pair<uint64_t, uint8_t>;
    const auto file_writer = [&out](const Pair &v) {
        out.write(reinterpret_cast<const char *>(&v), sizeof(Pair));
    };
    for (auto _ : state) {
        mg::common::merge_files<uint64_t, uint8_t>(sources, file_writer);
    }
}

BENCHMARK(BM_merge_files_pairs)->DenseRange(10, 300, 10);;

BENCHMARK_MAIN();
