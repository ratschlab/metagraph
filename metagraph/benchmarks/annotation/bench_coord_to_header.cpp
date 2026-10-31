// The header-to-record index of CoordToHeader (find_header) on a large synthetic .seqs: its bytes
// per header, its build time on one thread and on every thread, its lookups, and the same for
// the hash table it replaced (a tsl::hopscotch_map of views into the headers, about 40 bytes
// per bucket at load factor 0.9 rounded up to a power of two). Each runs once and prints its
// table to stderr; a lookup that answers wrong ends it with an error.
//
//     cd metagraph/build && make benchmarks
//     ./benchmarks --benchmark_filter=CoordToHeaderIndex
//
// $HEADER_INDEX_HEADERS sets the number of headers (default 20,000,000), $HEADER_INDEX_COLUMNS
// the number of columns they are spread over (default 1,000; AllTheBacteria has 1,027 batch
// columns, UniParc is one column), $HEADER_INDEX_THREADS the threads of the parallel build
// (default: every hardware thread; fewer on a shared machine).

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <benchmark/benchmark.h>
#include <tsl/hopscotch_map.h>

#include "annotation/coord_to_header.hpp"
#include "common/unix_tools.hpp"


namespace {

using namespace mtg;
using annot::CoordToHeader;

uint64_t env_number(const char *name, uint64_t fallback) {
    const char *value = std::getenv(name);
    return value ? std::strtoull(value, nullptr, 10) : fallback;
}

double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// |n| headers shaped like assembly contig names ("SAMN00012345.contig00017", 24 bytes), dealt
// round-robin over |num_columns| columns so that every column holds a slice of the samples
std::vector<std::vector<std::string>> synthetic_headers(uint64_t n, uint64_t num_columns) {
    std::vector<std::vector<std::string>> headers(num_columns);
    char buffer[64];
    for (uint64_t i = 0; i < n; ++i) {
        std::snprintf(buffer, sizeof(buffer), "SAMN%08llu.contig%05llu",
                      static_cast<unsigned long long>(i / 64),
                      static_cast<unsigned long long>(i % 64));
        headers[(i / 64) % num_columns].emplace_back(buffer);
    }
    return headers;
}

std::unique_ptr<CoordToHeader> mapping_of(std::vector<std::vector<std::string>> headers) {
    std::vector<std::vector<uint64_t>> num_kmers;
    for (const auto &column : headers) {
        num_kmers.emplace_back(column.size(), 1);
    }
    return std::make_unique<CoordToHeader>(std::move(headers), std::move(num_kmers));
}

void BM_CoordToHeaderIndex(benchmark::State &state) {
    const uint64_t n = env_number("HEADER_INDEX_HEADERS", 20'000'000);
    const uint64_t num_columns = std::max<uint64_t>(1, env_number("HEADER_INDEX_COLUMNS", 1000));
    const size_t threads = std::max<uint64_t>(
            1, env_number("HEADER_INDEX_THREADS", std::thread::hardware_concurrency()));
    for (auto _ : state) {
        const auto headers = synthetic_headers(n, num_columns);
        std::mt19937_64 rng(7);
        std::vector<std::pair<uint32_t, uint32_t>> probes;
        for (size_t i = 0; i < 1'000'000; ++i) {
            const uint32_t col = rng() % num_columns;
            if (headers[col].empty())
                continue;
            probes.emplace_back(col, rng() % headers[col].size());
        }

        std::cerr << "CoordToHeader header index, " << n << " headers in " << num_columns
                  << " columns (" << threads << " threads available)\n";
        // the compact index on one thread and on every thread (a fresh object each, so that
        // each builds)
        for (size_t t : { size_t(1), threads }) {
            auto cth = mapping_of(headers);
            const size_t rss_before = get_curr_RSS();
            auto start = std::chrono::steady_clock::now();
            const size_t distinct = cth->build_header_index(t);
            const double build_s = seconds_since(start);
            const size_t rss_after = get_curr_RSS();
            start = std::chrono::steady_clock::now();
            size_t wrong = 0;
            for (const auto &[col, s] : probes) {
                auto found = cth->find_header(headers[col][s]);
                wrong += !found || found->first != col || found->second != s;
            }
            const double hit_ns = seconds_since(start) * 1e9 / probes.size();
            start = std::chrono::steady_clock::now();
            for (const auto &[col, s] : probes) {
                wrong += cth->find_header(headers[col][s] + "x").has_value();
            }
            const double miss_ns = seconds_since(start) * 1e9 / probes.size();
            if (wrong || distinct != n) {
                state.SkipWithError("the compact index answered wrong");
                return;
            }
            std::fprintf(stderr, "  sorted array, %2zu thread(s): %6.2f bytes/header (%zu MiB; "
                         "RSS +%zu MiB), build %.2f s, lookup %.0f ns (present) %.0f ns "
                         "(absent)\n", t, double(cth->header_index_bytes()) / n,
                         cth->header_index_bytes() >> 20, (rss_after - rss_before) >> 20,
                         build_s, hit_ns, miss_ns);
        }
        // the hash table it replaced, built as coord_to_header.cpp built it
        {
            const size_t rss_before = get_curr_RSS();
            auto start = std::chrono::steady_clock::now();
            tsl::hopscotch_map<std::string_view, std::pair<uint64_t, size_t>> map;
            map.reserve(n);
            for (uint64_t col = 0; col < headers.size(); ++col) {
                for (size_t s = 0; s < headers[col].size(); ++s) {
                    map.try_emplace(std::string_view(headers[col][s]), col, s);
                }
            }
            const double build_s = seconds_since(start);
            const size_t rss_after = get_curr_RSS();
            start = std::chrono::steady_clock::now();
            size_t wrong = 0;
            for (const auto &[col, s] : probes) {
                auto it = map.find(headers[col][s]);
                wrong += it == map.end() || it->second.first != col || it->second.second != s;
            }
            const double hit_ns = seconds_since(start) * 1e9 / probes.size();
            if (wrong) {
                state.SkipWithError("the hash table answered wrong");
                return;
            }
            std::fprintf(stderr, "  hopscotch map (before E11), 1 thread: %zu buckets, "
                         "%6.2f bytes/header by RSS (+%zu MiB), build %.2f s, lookup %.0f ns "
                         "(present)\n", map.bucket_count(),
                         double(rss_after - rss_before) / n, (rss_after - rss_before) >> 20,
                         build_s, hit_ns);
        }
    }
}

BENCHMARK(BM_CoordToHeaderIndex)->Iterations(1)->Unit(benchmark::kSecond);

} // namespace
