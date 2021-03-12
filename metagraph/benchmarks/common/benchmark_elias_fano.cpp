#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#include "common/elias_fano/elias_fano.hpp"
#include "common/vector.hpp"
#include "common/utils/file_utils.hpp"


namespace {

using namespace mtg;
using mtg::elias_fano::EliasFanoEncoderBuffered;
using mtg::elias_fano::EliasFanoDecoder;

template <typename T>
class EliasFanoFixture : public benchmark::Fixture {
  public:
    std::vector<T> sorted;
    static T sum_compressed;
    static T sum_uncompressed;
    size_t size;

    EliasFanoFixture() { Unit(benchmark::TimeUnit::kMillisecond); init_sorted(); }

    void init_sorted() {

        if (!sorted.empty()) {
            return;
        }
        sorted.resize(1e8 / sizeof(T)); // test on 100MB worth of data
        std::mt19937 rng(123457);
        std::uniform_int_distribution<std::mt19937::result_type> dist30(0, 30);
        std::uniform_int_distribution<std::mt19937::result_type> dist10000(0, 10000);
        T i = 0;
        uint32_t half_size_bits = sizeof(T) * 4;
        std::for_each(sorted.begin(), sorted.end(), [&](T &v) {
            i += dist10000(rng);
            if (dist30(rng) < 1) { // increase the hi bits every ~30th element
                i = ((i >> half_size_bits) + 1) << half_size_bits;
            }
            v = i;
        });
        std::ifstream f;
    }

    void encode() {
        utils::TempFile tempfile;
        size = EliasFanoEncoderBuffered<T>::append_block(sorted, tempfile.name());
    }

    void write_compressed(benchmark::State &state) {
        for (auto _ : state) {
            encode();
        }
        std::cout << "Write compressed: compression factor for : " << sizeof(T)
                  << "-byte integers" << (double)sorted.size() * sizeof(T) / size
                  << std::endl;
    }

    void write_uncompressed(benchmark::State &state) {
        for (auto _ : state) {
            utils::TempFile tempfile;
            std::ofstream &out = tempfile.ofstream();
            out.write(reinterpret_cast<char *>(sorted.data()), sorted.size() * sizeof(T));
            out.close();
        }
    }

    void read_compressed(benchmark::State &state) {
        utils::TempFile tempfile;
        EliasFanoEncoderBuffered<T>::append_block(sorted, tempfile.name());
        for (auto _ : state) {
            EliasFanoDecoder<T> decoder(tempfile.name());
            std::optional<T> value;
            sum_compressed = 0;
            while ((value = decoder.next()).has_value()) {
                sum_compressed += value.value();
            }
        }
    }

    void read_uncompressed(benchmark::State &state) {
        utils::TempFile tempfile;
        std::ofstream &out = tempfile.ofstream();
        out.write(reinterpret_cast<char *>(sorted.data()), sorted.size() * sizeof(T));
        out.close();
        for (auto _ : state) {
            std::ifstream in = std::ifstream(tempfile.name(), std::ios::binary);
            T value;
            sum_uncompressed = 0;
            while (in.read(reinterpret_cast<char *>(&value), sizeof(T))) {
                sum_uncompressed += value;
            }
            // making sure the compiler doesn't optimized away the reading and doing some
            // sanity check
            if (sum_compressed != sum_uncompressed) {
                std::cerr << "Error: Compressed and Non-compressed reads don't match. "
                          << " for " << sizeof(T) << " bytes. You have a bug. "
                          << std::endl;
            }
        }
    }

    void read_uncompressed1024(benchmark::State &state) {
        utils::TempFile tempfile;
        std::ofstream &out = tempfile.ofstream();
        out.write(reinterpret_cast<char *>(sorted.data()), sorted.size() * sizeof(T));
        out.close();
        uint64_t res[128];
        for (auto _ : state) {
            std::ifstream in = std::ifstream(tempfile.name(), std::ios::binary);
            sum_uncompressed = 0;
            while (in.read(reinterpret_cast<char *>(res), 256)) {
                for (uint32_t i = 0; i < std::min(32L, in.gcount()/8); ++i) {
                    sum_uncompressed += res[i];
                }
            }
            // making sure the compiler doesn't optimized away the reading and doing some
            // sanity check
            if (sum_compressed != sum_uncompressed) {
                std::cerr << "Error: Compressed and Non-compressed reads don't match. "
                          << " for " << sizeof(T) << " bytes. You have a bug. "
                          << std::endl;
            }
        }
    }
};

template <typename T>
T EliasFanoFixture<T>::sum_compressed;
template <typename T>
T EliasFanoFixture<T>::sum_uncompressed;

BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_uncompressed64, uint64_t)
(benchmark::State &state) {
    write_uncompressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_compressed64, uint64_t)
(benchmark::State &state) {
    write_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_compressed64, uint64_t)
(benchmark::State &state) {
    read_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_uncompressed64, uint64_t)
(benchmark::State &state) {
    read_uncompressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_uncompressed1024, uint64_t)
(benchmark::State &state) {
    read_uncompressed1024(state);
}


BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_uncompressed128, sdsl::uint128_t)
(benchmark::State &state) {
    write_uncompressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_compressed128, sdsl::uint128_t)
(benchmark::State &state) {
    write_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_compressed128, sdsl::uint128_t)
(benchmark::State &state) {
    read_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_uncompressed128, sdsl::uint128_t)
(benchmark::State &state) {
    read_uncompressed(state);
}

BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_uncompressed256, sdsl::uint256_t)
(benchmark::State &state) {
    write_uncompressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_write_compressed256, sdsl::uint256_t)
(benchmark::State &state) {
    write_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_compressed256, sdsl::uint256_t)
(benchmark::State &state) {
    read_compressed(state);
}
BENCHMARK_TEMPLATE_F(EliasFanoFixture, BM_read_uncompressed256, sdsl::uint256_t)
(benchmark::State &state) {
    read_uncompressed(state);
}

} // namespace
