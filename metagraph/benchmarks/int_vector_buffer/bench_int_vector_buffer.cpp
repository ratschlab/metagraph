#include "benchmark/benchmark.h"

#include <filesystem>
#include <sdsl/int_vector_buffer.hpp>


namespace {

namespace fs = std::filesystem;

const size_t size = 1llu << 27;

std::string get_fname(benchmark::State &state) {
    if (!std::getenv("VEC")) {
        state.SkipWithError("Set environment variable VEC");
    }
    std::string fname = std::getenv("VEC");
    if (fs::exists(fname))
        throw std::runtime_error("File " + fname + " already exists!");

    return fname;
}


static void BM_int_vector_buffer_sequential_write_GiB(benchmark::State &state) {
    const size_t width = 6;
    const size_t buffer_size = state.range(0);

    const std::string fname = get_fname(state);
    for (auto _ : state) {
        sdsl::int_vector_buffer<> buf(fname, std::ios::out, buffer_size, width);
        for (size_t i = 0; i < size; ++i) {
            buf.push_back(i);
        }
        buf.close(true);
        state.counters["Speed"]
            = benchmark::Counter(size * width / 8,
                                 benchmark::Counter::kIsIterationInvariantRate,
                                 benchmark::Counter::OneK::kIs1024);
    }
}

BENCHMARK(BM_int_vector_buffer_sequential_write_GiB)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(4)->Range(512, 2 << 25);



static void BM_int_vector_buffer_sequential_read_GiB(benchmark::State &state) {
    const size_t width = 6;
    const size_t buffer_size = state.range(0);

    const std::string fname = get_fname(state);
    {
        sdsl::int_vector_buffer<> buf(fname, std::ios::out, 1024*1024, width);
        for (size_t i = 0; i < size; ++i) {
            buf.push_back(i);
        }
    }

    uint64_t result = 0;
    sdsl::int_vector_buffer<> buf(fname, std::ios::in, buffer_size, width);

    for (auto _ : state) {
        for (size_t i = 0; i < size; ++i) {
            result += buf[i];
        }
        state.counters["Speed"]
            = benchmark::Counter(size * width / 8,
                                 benchmark::Counter::kIsIterationInvariantRate,
                                 benchmark::Counter::OneK::kIs1024);
    }

    buf.close(true);
}

BENCHMARK(BM_int_vector_buffer_sequential_read_GiB)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(4)->Range(512, 2 << 25);


static void BM_int_vector_buffer_sequential_add_GiB(benchmark::State &state) {
    const size_t width = 6;
    const size_t buffer_size = state.range(0);

    const std::string fname = get_fname(state);
    {
        sdsl::int_vector_buffer<> buf(fname, std::ios::out, 1024*1024, width);
        for (size_t i = 0; i < size; ++i) {
            buf.push_back(i);
        }
    }

    sdsl::int_vector_buffer<> buf(fname, std::ios::in|std::ios::out, buffer_size, width);

    for (auto _ : state) {
        for (size_t i = 0; i < size; ++i) {
            buf[i] += 1;
        }
        state.counters["Speed"]
            = benchmark::Counter(size * width / 8,
                                 benchmark::Counter::kIsIterationInvariantRate,
                                 benchmark::Counter::OneK::kIs1024);
    }

    buf.close(true);
}

BENCHMARK(BM_int_vector_buffer_sequential_add_GiB)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(4)->Range(512, 2 << 25);

} // namespace
