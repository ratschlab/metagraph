#ifndef __DBG_BITMAP_CONSTRUCT_HPP__
#define __DBG_BITMAP_CONSTRUCT_HPP__

#include "dbg_bitmap.hpp"
#include "dbg_construct.hpp"


class ISDChunkConstructor : public IChunkConstructor<DBGSD::Chunk> {
  public:
    virtual void add_sequence(const std::string &sequence) = 0;

    virtual DBGSD::Chunk* build_chunk() = 0;

    virtual size_t get_k() const = 0;

    virtual bool is_canonical_mode() const = 0;

    // TODO: replace with variadic template
    static ISDChunkConstructor* initialize(
        size_t k,
        bool canonical_mode = false,
        const std::string &filter_suffix = std::string(),
        size_t num_threads = 1,
        double memory_preallocated = 0,
        bool verbose = false
    );
};


template <typename KMER>
class SDChunkConstructor : public ISDChunkConstructor {
    friend ISDChunkConstructor;

  public:
    SDChunkConstructor(size_t k,
                       bool canonical_mode = false,
                       const std::string &filter_suffix = std::string(),
                       size_t num_threads = 1,
                       double memory_preallocated = 0,
                       bool verbose = false);

    inline void add_sequence(const std::string &sequence) {
        kmer_collector_.add_sequence(sequence);
    }

    inline void add_sequences(std::function<void(CallbackString)> generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    inline size_t get_k() const { return kmer_collector_.get_k(); }

    inline bool is_canonical_mode() const { return kmer_collector_.is_canonical_mode(); }

    DBGSD::Chunk* build_chunk();

  private:
    KmerCollector<KMER, KmerExtractor2Bit> kmer_collector_;
};


class DBGSDConstructor : public GraphConstructor {
  public:
    explicit DBGSDConstructor(size_t k,
                              bool canonical_mode = false,
                              const std::string &filter_suffix = std::string(),
                              size_t num_threads = 1,
                              double memory_preallocated = 0,
                              bool verbose = false)
          : constructor_(ISDChunkConstructor::initialize(
                k,
                canonical_mode,
                filter_suffix,
                num_threads,
                memory_preallocated,
                verbose)
          ) {}
    inline void add_sequence(const std::string &sequence) {
        constructor_->add_sequence(sequence);
    }

    inline void add_sequences(const std::vector<std::string> &sequences) {
        constructor_->add_sequences(
            [&sequences](const CallbackString &callback) {
                std::for_each(sequences.begin(), sequences.end(), callback);
            }
        );
    }

    inline void add_sequences(const std::function<void(CallbackString)> &callback) {
        constructor_->add_sequences(callback);
    }

    void build_graph(DBGSD *graph);

    DBGSD::Chunk* build_chunk();

    static DBGSD* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                          bool canonical_mode = false,
                                          bool verbose = false);
  private:
    std::unique_ptr<ISDChunkConstructor> constructor_;
};

#endif // __DBG_BITMAP_CONSTRUCT_HPP__
