#ifndef __BOSS_CHUNK_CONSTRUCT_HPP__
#define __BOSS_CHUNK_CONSTRUCT_HPP__

#include "graph/succinct/boss_chunk.hpp"

#include "graph/base/dbg_construct.hpp"

#include <string>
#include <cstdint>
#include <functional>

namespace mg {
namespace succinct {


class IBOSSChunkConstructor : public IGraphChunkConstructor<BOSS::Chunk> {
  public:
    virtual ~IBOSSChunkConstructor() {}

    static std::unique_ptr<IBOSSChunkConstructor>
    initialize(size_t k,
               bool use_sorted_set_disk = false,
               bool canonical_mode = false,
               bool count_kmers = false,
               const std::string &filter_suffix = "",
               size_t num_threads = 1,
               double memory_preallocated = 0);

    virtual void add_sequence(std::string &&sequence, uint64_t count = 1) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual BOSS::Chunk *build_chunk() = 0;

    virtual uint64_t get_k() const = 0;
};

} // namespace succinct
} // namespace mg

#endif // __BOSS_CHUNK_CONSTRUCT_HPP__
