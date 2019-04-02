#ifndef __DBG_SUCCINCT_CONSTRUCT_HPP__
#define __DBG_SUCCINCT_CONSTRUCT_HPP__

#include "dbg_construct.hpp"
#include "dbg_succinct.hpp"


class IDBGBOSSChunkConstructor : public IChunkConstructor<DBG_succ::Chunk> {
  public:
    virtual void add_sequence(const std::string &sequence) = 0;

    virtual DBG_succ::Chunk* build_chunk() = 0;

    // TODO: replace with variadic template
    static IDBGBOSSChunkConstructor* initialize(
        size_t k,
        bool canonical_mode = false,
        const std::string &filter_suffix = std::string(),
        size_t num_threads = 1,
        double memory_preallocated = 0,
        bool verbose = false
    );
};

template <typename KMER>
class DBGBOSSChunkConstructor : public IDBGBOSSChunkConstructor {
    friend IDBGBOSSChunkConstructor;

  public:
    DBGBOSSChunkConstructor(size_t k,
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

    DBG_succ::Chunk* build_chunk();

  private:
    KmerCollector<KMER, KmerExtractor> kmer_collector_;
};


class DBGSuccConstructor : public GraphConstructor {
  public:
    explicit DBGSuccConstructor(size_t k,
                                bool canonical_mode = false,
                                const std::string &filter_suffix = std::string(),
                                size_t num_threads = 1,
                                double memory_preallocated = 0,
                                bool verbose = false)
          : constructor_(IDBGBOSSChunkConstructor::initialize(
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

    void build_graph(DBG_succ *graph);

    static DBG_succ* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                             bool verbose = false);
  private:
    std::unique_ptr<IDBGBOSSChunkConstructor> constructor_;
};

#endif // __DBG_SUCCINCT_CONSTRUCT_HPP__
