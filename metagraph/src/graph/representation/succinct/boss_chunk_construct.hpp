#ifndef __BOSS_CHUNK_CONSTRUCT_HPP__
#define __BOSS_CHUNK_CONSTRUCT_HPP__

#include <string>
#include <cstdint>
#include <filesystem>
#include <functional>

#include "kmer/kmer_collector_config.hpp"
#include "graph/representation/base/dbg_construct.hpp"
#include "boss_chunk.hpp"

namespace mtg {
namespace succinct {


class IBOSSChunkConstructor : public IGraphChunkConstructor<BOSS::Chunk> {
  public:
    virtual ~IBOSSChunkConstructor() {}

    static std::unique_ptr<IBOSSChunkConstructor>
    initialize(size_t k,
               bool canonical_mode = false,
               uint8_t bits_per_count = 0,
               const std::string &filter_suffix = "",
               size_t num_threads = 1,
               double memory_preallocated = 0,
               mtg::kmer::ContainerType container_type = mtg::kmer::ContainerType::VECTOR,
               const std::filesystem::path &swap_dir = "/tmp/",
               size_t max_disk_space_bytes = 1e9);

    virtual void add_sequence(std::string_view sequence, uint64_t count = 1) = 0;
    virtual void add_sequences(const std::function<void(CallString)> &generator) = 0;
    virtual void add_sequences(const std::function<void(CallStringCount)> &generator) = 0;

    virtual BOSS::Chunk* build_chunk() = 0;

    virtual uint64_t get_k() const = 0;
};

} // namespace succinct
} // namespace mtg

#endif // __BOSS_CHUNK_CONSTRUCT_HPP__
