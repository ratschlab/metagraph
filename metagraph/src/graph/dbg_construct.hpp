#ifndef __DBG_CONSTRUCT_HPP__
#define __DBG_CONSTRUCT_HPP__

#include <mutex>
#include <shared_mutex>

#include "dbg_succinct.hpp"
#include "dbg_sd.hpp"
#include "dbg_succinct_chunk.hpp"
#include "kmer_extractor.hpp"
#include "utils.hpp"


typedef std::function<void(const std::string&)> CallbackString;


template <typename KMER, class KmerExtractor>
class KmerCollector {
    using Extractor = KmerExtractor;
    using Sequence = std::vector<typename Extractor::TAlphabet>;
    Extractor kmer_extractor_;
  public:
    explicit KmerCollector(size_t k,
                           bool canonical_mode = false,
                           Sequence&& filter_suffix_encoded = {},
                           size_t num_threads = 1,
                           double memory_preallocated = 0,
                           bool verbose = false);

    inline size_t get_k() const { return k_; }

    inline size_t size() const { return kmers_.size(); }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    // TODO: another for std::string&& ?
    void add_sequence(const std::string &sequence);

    void add_sequences(const std::function<void(CallbackString)> &generate_sequences);

    template <typename... Args>
    inline void emplace_back(Args&&... args) {
        kmers_.emplace_back(std::forward<Args>(args)...);
    }

    inline Vector<KMER>& data() { return kmers_; }

    inline void
    call_kmers(const std::function<void(const KMER&)> &kmer_callback) const {
        std::for_each(kmers_.begin(), kmers_.end(), kmer_callback);
    }

    inline void clear() {
        kmers_.clear();
        sequences_storage_.clear();
        stored_sequences_size_ = 0;
        kmers_.shrink_to_fit();
    }

    void join();

    inline bool verbose() const { return verbose_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }

  private:
    void release_task_to_pool();

    size_t k_;
    Vector<KMER> kmers_;
    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;

    size_t num_threads_;
    utils::ThreadPool thread_pool_;

    std::vector<std::string> sequences_storage_;
    size_t stored_sequences_size_;

    bool verbose_;

    Sequence filter_suffix_encoded_;

    bool canonical_mode_;
};


template <class GraphChunk>
class IChunkConstructor {
  public:
    virtual ~IChunkConstructor() {}

    virtual void add_sequence(const std::string &sequence) = 0;

    virtual void add_sequences(std::function<void(CallbackString)> generate_sequences) = 0;

    virtual GraphChunk* build_chunk() = 0;
};

class IDBGBOSSChunkConstructor : public IChunkConstructor<DBG_succ::Chunk> {
  public:
    virtual void add_sequence(const std::string &sequence) = 0;

    virtual DBG_succ::Chunk* build_chunk() = 0;

    // TODO: replace with variadic template
    static IDBGBOSSChunkConstructor* initialize(
        size_t k,
        const std::string &filter_suffix = std::string(),
        size_t num_threads = 1,
        double memory_preallocated = 0,
        bool verbose = false
    );
};

template <typename KMER>
class DBGBOSSChunkConstructor : public IDBGBOSSChunkConstructor {
    friend IDBGBOSSChunkConstructor;
    using Extractor = KmerExtractor;
  public:
    DBGBOSSChunkConstructor(size_t k,
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
    KmerCollector<KMER, Extractor> kmer_collector_;
};



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
    using Extractor = KmerExtractor2Bit;
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
    KmerCollector<KMER, Extractor> kmer_collector_;
};



class GraphConstructor {
  public:
    virtual ~GraphConstructor() {}
};

class DBGSuccConstructor : public GraphConstructor {
  public:
    explicit DBGSuccConstructor(size_t k,
                                const std::string &filter_suffix = std::string(),
                                size_t num_threads = 1,
                                double memory_preallocated = 0,
                                bool verbose = false)
          : constructor_(IDBGBOSSChunkConstructor::initialize(
                k,
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

    inline void build_graph(DBG_succ *graph) {
        auto chunk = constructor_->build_chunk();
        // initialize graph from the chunk built
        chunk->initialize_graph(graph);
        delete chunk;
    }

    static DBG_succ* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                             bool verbose = false);
  private:
    std::unique_ptr<IDBGBOSSChunkConstructor> constructor_;
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

#endif // __DBG_CONSTRUCT_HPP__
