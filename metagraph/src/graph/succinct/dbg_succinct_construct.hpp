#ifndef __DBG_SUCCINCT_CONSTRUCT_HPP__
#define __DBG_SUCCINCT_CONSTRUCT_HPP__

#include "dbg_construct.hpp"
#include "dbg_succinct.hpp"


class IBOSSChunkConstructor : public IGraphChunkConstructor<BOSS::Chunk> {
  public:
    virtual ~IBOSSChunkConstructor() {}

    static IBOSSChunkConstructor* initialize(size_t k,
                                             bool canonical_mode = false,
                                             const std::string &filter_suffix = "",
                                             size_t num_threads = 1,
                                             double memory_preallocated = 0,
                                             bool verbose = false);

    virtual void add_sequence(const std::string &sequence) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual BOSS::Chunk* build_chunk() = 0;
};


class BOSSConstructor : public IGraphConstructor<BOSS> {
  public:
    BOSSConstructor(size_t k,
                    bool canonical_mode = false,
                    const std::string &filter_suffix = "",
                    size_t num_threads = 1,
                    double memory_preallocated = 0,
                    bool verbose = false);

    void add_sequence(const std::string &sequence) {
        constructor_->add_sequence(sequence);
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void add_sequences(const std::vector<std::string> &sequences) {
        constructor_->add_sequences([&sequences](const CallString &callback) {
            std::for_each(sequences.begin(), sequences.end(), callback);
        });
    }

    void build_graph(BOSS *graph);

    static BOSS* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                         bool verbose = false);
  private:
    std::unique_ptr<IBOSSChunkConstructor> constructor_;
};

#endif // __DBG_SUCCINCT_CONSTRUCT_HPP__
