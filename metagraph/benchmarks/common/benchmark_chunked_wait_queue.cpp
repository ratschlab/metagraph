#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>
#include <common/threads/chunked_wait_queue.hpp>

#include "common/file_merger.hpp"

constexpr size_t CHUNK_COUNT = 10;
constexpr size_t ITEM_COUNT = 100'000;
const std::string chunk_prefix = "/tmp/chunk_";

void create_sources() {
    std::mt19937 rng(123456);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);


    for (uint32_t i = 0; i < CHUNK_COUNT; ++i) {
        std::vector<uint32_t> els(ITEM_COUNT);
        for (uint64_t j = 0; j < ITEM_COUNT; ++j) {
            els[i] = j + dist10(rng);
        }
        std::ofstream f(chunk_prefix + std::to_string(i), std::ios::binary);
        f.write(reinterpret_cast<char *>(els.data()), els.size() * sizeof(uint32_t));
        f.close();
    }
}

static void BM_queue_push_pop(benchmark::State &state) {
    mg::common::ChunkedWaitQueue<uint64_t> queue(100, 10);
    size_t sum = 0;
    for (auto _ : state) {
        for (uint32_t i = 0; i < 100; ++i) {
            queue.push(i);
        }
        queue.shutdown();
        auto &it = queue.begin();
        for (; it != queue.end(); ++it) {
            sum += *it;
        }
        for (uint32_t i = 0; i < 100; ++i) {
            --it;
        }
        queue.reset();
    }
    std::ofstream f("/tmp/dump");
    f << sum;
}

static void BM_queue_push_pop_back(benchmark::State &state) {
    mg::common::ChunkedWaitQueue<uint64_t> queue(100, 10);
    size_t sum = 0;
    for (auto _ : state) {
        for (uint32_t i = 0; i < 100; ++i) {
            queue.push(i);
        }
        queue.shutdown();
        auto &it = queue.begin();
        for (; it != queue.end(); ++it) {
            sum += *it;
        }
        for (uint32_t i = 0; i < 10; ++i) {
            --it;
            sum += *it;
        }
        queue.reset();
    }
    std::ofstream f("/tmp/dump");
    f << sum;
}


BENCHMARK(BM_queue_push_pop);
BENCHMARK(BM_queue_push_pop_back);
