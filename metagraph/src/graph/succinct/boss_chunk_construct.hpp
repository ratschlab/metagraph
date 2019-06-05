#ifndef __BOSS_CHUNK_CONSTRUCT_HPP__
#define __BOSS_CHUNK_CONSTRUCT_HPP__

#include "dbg_construct.hpp"
#include "boss_chunk.hpp"


class IBOSSChunkConstructor : public IGraphChunkConstructor<BOSS::Chunk> {
  public:
    virtual ~IBOSSChunkConstructor() {}

    static std::unique_ptr<IBOSSChunkConstructor>
    initialize(size_t k,
               bool canonical_mode = false,
               bool count_kmers = false,
               const std::string &filter_suffix = "",
               size_t num_threads = 1,
               double memory_preallocated = 0,
               bool verbose = false);

    virtual void add_sequence(std::string&& sequence) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual BOSS::Chunk* build_chunk() = 0;
};

#endif // __BOSS_CHUNK_CONSTRUCT_HPP__
