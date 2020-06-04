#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/threads/chunked_wait_queue.hpp"


namespace {

using namespace mtg;

const std::string chunk_prefix = "/tmp/chunk_";


template <typename T>
static void BM_queue_push_pop(benchmark::State &state) {
    common::ChunkedWaitQueue<T> queue(10000);
    T sum = 0;
    for (auto _ : state) {
        for (uint32_t i = 0; i < 10000; ++i) {
            queue.push(std::move(i));
        }
        queue.shutdown();
        auto &it = queue.begin();
        for (; it != queue.end(); ++it) {
            sum += *it;
        }
        queue.reset();
    }
    std::ofstream f("/tmp/dump");
    f << uint64_t(sum);
}

BENCHMARK_TEMPLATE(BM_queue_push_pop, uint64_t);
BENCHMARK_TEMPLATE(BM_queue_push_pop, sdsl::uint128_t);
BENCHMARK_TEMPLATE(BM_queue_push_pop, sdsl::uint256_t);

} // namespace
