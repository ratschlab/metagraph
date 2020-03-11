#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#include "common/elias_fano.hpp"
#include "common/vector.hpp"
#include "common/utils/file_utils.hpp"

namespace {
using namespace mg;
constexpr size_t ITEM_COUNT = 1e9 / 8; // test on 1GB worth of data
Vector<uint64_t> sorted(ITEM_COUNT);

static void init_sorted() {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

    uint64_t i = 0;
    std::for_each(sorted.begin(), sorted.end(), [&](uint64_t &v) {
        i += dist10(rng);
        v = i;
    });
}

static void BM_write_compressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        std::ofstream &out = tempfile.ofstream();
        common::EliasFanoEncoder<uint64_t> encoder(sorted, out);
        encoder.finish();
    }
}

static void BM_write_uncompressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        std::ofstream &out = tempfile.ofstream();
        out.write(reinterpret_cast<char *>(sorted.data()), sorted.size() * sizeof(uint64_t));
        out.close();
    }
}
uint64_t sum_compressed;
static void BM_read_compressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }
    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    common::EliasFanoEncoder<uint64_t> encoder(sorted, out);
    encoder.finish();
    for (auto _ : state) {
        common::EliasFanoDecoder<uint64_t> decoder(tempfile.ifstream());
        std::optional<uint64_t> value;
        sum_compressed = 0;
        while ((value = decoder.next()).has_value()) {
            sum_compressed += value.value();
        }
    }
}

uint64_t sum_uncompressed;
static void BM_read_uncompressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }
    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    out.write(reinterpret_cast<char *>(sorted.data()), sorted.size() * sizeof(uint64_t));
    out.close();
    for (auto _ : state) {
        std::ifstream &in = tempfile.ifstream();
        uint64_t value;
        sum_uncompressed = 0;
        while (in.read(reinterpret_cast<char *>(&value), sizeof(uint64_t))) {
            sum_uncompressed += value;
        }
        // making sure the compiler doesn't optimized away the reading and doing some
        // sanity check
        if (sum_compressed != sum_uncompressed) {
            std::cerr << "Error: Compressed and Non-compressed reads don't match. You "
                         "have a bug." << std::endl;
        }
    }
}

BENCHMARK(BM_write_compressed);
BENCHMARK(BM_write_uncompressed);
BENCHMARK(BM_read_compressed);
BENCHMARK(BM_read_uncompressed);
} // namespace
