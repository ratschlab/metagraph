#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#include "common/elias_fano.hpp"
#include "common/elias_fano_128.hpp"
#include "common/vector.hpp"
#include "common/utils/file_utils.hpp"

namespace {
using namespace mg;
constexpr size_t ITEM_COUNT = 1e8 / 8; // test on 1GB worth of data TODO:undo 8->9
Vector<uint64_t> sorted(ITEM_COUNT);
Vector<sdsl::uint128_t> sorted128(ITEM_COUNT / 2);

static void init_sorted() {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10000(0, 10000);

    uint64_t i = 0;
    std::for_each(sorted.begin(), sorted.end(), [&](uint64_t &v) {
        i += dist10000(rng);
        v = i;
    });
}

static void init_sorted128() {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist30(0, 30);
    std::uniform_int_distribution<std::mt19937::result_type> dist10000(0, 10000);
    sdsl::uint128_t i = 0;
    std::for_each(sorted128.begin(), sorted128.end(), [&](sdsl::uint128_t &v) {
        i += dist10000(rng);
        if (dist30(rng) < 1) { // increase the hi 64 bits every ~10th element
            i = ((sdsl::uint128_t)((uint64_t)(i >> 64) + 1) << 64) + (uint64_t)i;
        }
        v = i;
    });
}

static void BM_write_compressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        common::EliasFanoEncoder<uint64_t> encoder(sorted.size(), sorted.front(),
                                                   sorted.back(), tempfile.name());
        for (const auto &v : sorted) {
            encoder.add(v);
        }
        encoder.finish();
    }
}

static void BM_write_compressed128(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted128();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        common::EliasFanoEncoder<sdsl::uint128_t> encoder(
                sorted128.size(), sorted128.front(), sorted128.back(), tempfile.name());
        for (const auto &v : sorted128) {
            encoder.add(v);
        }
        encoder.finish();
    }
}

static void BM_write_compressed128_2(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted128();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        common::EliasFanoEncoder128<sdsl::uint128_t> encoder(
                sorted128.size(), sorted128.front(), sorted128.back(), tempfile.name());
        for (const auto &v : sorted128) {
            encoder.add(v);
        }
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

static void BM_write_uncompressed128(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted128();
    }

    for (auto _ : state) {
        utils::TempFile tempfile;
        std::ofstream &out = tempfile.ofstream();
        out.write(reinterpret_cast<char *>(sorted128.data()),
                  sorted128.size() * sizeof(sdsl::uint128_t));
        out.close();
    }
}

uint64_t sum_compressed;
static void BM_read_compressed(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted();
    }
    utils::TempFile tempfile;
    common::EliasFanoEncoder<uint64_t> encoder(sorted.size(), sorted.front(),
                                               sorted.back(), tempfile.name());
    for (const auto &v : sorted) {
        encoder.add(v);
    }
    encoder.finish();
    for (auto _ : state) {
        common::EliasFanoDecoder<uint64_t> decoder(tempfile.name());
        std::optional<uint64_t> value;
        sum_compressed = 0;
        while ((value = decoder.next()).has_value()) {
            sum_compressed += value.value();
        }
    }
}

sdsl::uint128_t sum_compressed128;
static void BM_read_compressed128_2(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted128();
    }
    utils::TempFile tempfile;
    common::EliasFanoEncoder128<sdsl::uint128_t> encoder(sorted128.size(), sorted128.front(),
                                               sorted128.back(), tempfile.name());
    for (const auto &v : sorted128) {
        encoder.add(v);
    }
    encoder.finish();
    for (auto _ : state) {
        common::EliasFanoDecoder128<sdsl::uint128_t> decoder(tempfile.name());
        std::optional<sdsl::uint128_t > value;
        sum_compressed128 = 0;
        while ((value = decoder.next()).has_value()) {
            sum_compressed128 += value.value();
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
        std::ifstream in = std::ifstream(tempfile.name(), std::ios::binary);
        uint64_t value;
        sum_uncompressed = 0;
        while (in.read(reinterpret_cast<char *>(&value), sizeof(uint64_t))) {
            sum_uncompressed += value;
        }
        // making sure the compiler doesn't optimized away the reading and doing some
        // sanity check
        if (sum_compressed != sum_uncompressed) {
            std::cerr << "Error: Compressed and Non-compressed reads don't match. You "
                         "have a bug. "
                      << sum_compressed << " vs. " << sum_uncompressed << std::endl;
        }
    }
}

sdsl::uint128_t sum_uncompressed128;
static void BM_read_uncompressed128(benchmark::State &state) {
    if (state.thread_index == 0) {
        init_sorted128();
    }
    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    out.write(reinterpret_cast<char *>(sorted128.data()), sorted128.size() * sizeof(sdsl::uint128_t));
    out.close();
    for (auto _ : state) {
        std::ifstream in = std::ifstream(tempfile.name(), std::ios::binary);
        sdsl::uint128_t value;
        sum_uncompressed128 = 0;
        while (in.read(reinterpret_cast<char *>(&value), sizeof(sdsl::uint128_t))) {
            sum_uncompressed128 += value;
        }
        // making sure the compiler doesn't optimized away the reading and doing some
        // sanity check
        if (sum_compressed128 != sum_uncompressed128) {
            std::cerr << "Error: Compressed and Non-compressed reads don't match. You "
                         "have a bug. " << std::endl;
        }
    }
}

BENCHMARK(BM_write_compressed);
BENCHMARK(BM_write_uncompressed);
BENCHMARK(BM_write_compressed128);
BENCHMARK(BM_write_compressed128_2);
BENCHMARK(BM_write_uncompressed128);
BENCHMARK(BM_read_compressed);
BENCHMARK(BM_read_compressed128_2);
BENCHMARK(BM_read_uncompressed128);
BENCHMARK(BM_read_uncompressed);
} // namespace
