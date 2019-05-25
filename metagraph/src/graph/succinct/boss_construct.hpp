#ifndef __BOSS_CONSTRUCT_HPP__
#define __BOSS_CONSTRUCT_HPP__

#include "dbg_construct.hpp"
#include "boss_chunk.hpp"


class IBOSSChunkConstructor : public IGraphChunkConstructor<BOSS::Chunk> {
  public:
    virtual ~IBOSSChunkConstructor() {}

    static IBOSSChunkConstructor* initialize(size_t k,
                                             bool canonical_mode = false,
                                             const std::string &filter_suffix = "",
                                             size_t num_threads = 1,
                                             double memory_preallocated = 0,
                                             bool verbose = false);

    virtual void add_sequence(std::string&& sequence) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual BOSS::Chunk* build_chunk() = 0;
};


class BOSSConstructor : public IGraphConstructor<BOSS> {
  public:
    // see input arguments in IBOSSChunkConstructor::initialize
    template <typename... Args>
    BOSSConstructor(const Args&... args)
      : constructor_(IBOSSChunkConstructor::initialize(args...)) {}

    void add_sequence(std::string&& sequence) {
        constructor_->add_sequence(std::move(sequence));
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void add_sequences(const std::vector<std::string> &sequences) {
        constructor_->add_sequences([&sequences](const CallString &callback) {
            std::for_each(sequences.begin(), sequences.end(), callback);
        });
    }

    void build_graph(BOSS *graph) {
        auto chunk = constructor_->build_chunk();
        // initialize graph from the chunk built
        chunk->initialize_boss(graph);
        delete chunk;
    }

    static BOSS* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                         bool verbose = false) {
        return BOSS::Chunk::build_boss_from_chunks(chunk_filenames, verbose);
    }

  private:
    std::unique_ptr<IBOSSChunkConstructor> constructor_;
};

#endif // __BOSS_CONSTRUCT_HPP__
