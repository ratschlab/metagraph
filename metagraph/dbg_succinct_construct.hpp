#ifndef __DBG_SUCCINCT_CONSTRUCT_HPP__
#define __DBG_SUCCINCT_CONSTRUCT_HPP__

#include "dbg_succinct.hpp"
#include "kmer.hpp"
#include "utils.hpp"


class DBGSuccConstructor {
  public:
    virtual void add_read(const std::string &read) = 0;

    virtual void add_reads(const std::vector<std::string> &reads) {
        for (const auto &sequence : reads) {
            add_read(sequence);
        }
    }

    virtual void build_graph(DBG_succ *graph) = 0;

    virtual size_t get_k() const = 0;
};


class SuffixArrayDBGSuccConstructor : public DBGSuccConstructor {
  public:
    explicit SuffixArrayDBGSuccConstructor(size_t k);

    void add_read(const std::string &read);

    void build_graph(DBG_succ *graph);

    size_t get_k() const { return k_; }

  private:
    size_t k_;
    std::string data_;
};


typedef std::function<void(const std::string&)> CallbackRead;


class KMerDBGSuccChunkConstructor {
  public:
    KMerDBGSuccChunkConstructor(size_t k,
                                const std::string &filter_suffix,
                                size_t num_threads = 1,
                                double memory_preallocated = 0,
                                bool verbose = false);

    void add_read(const std::string &sequence);

    void add_reads(std::function<void(CallbackRead)> generate_reads);

    DBG_succ::Chunk* build_chunk();

    size_t get_k() const { return k_; }

  private:
    void release_task_to_pool();

    size_t k_;
    std::vector<KMer> kmers_;
    size_t end_sorted_;
    std::mutex mutex_;

    size_t num_threads_;
    utils::ThreadPool thread_pool_;

    std::vector<std::string> reads_storage_;
    size_t stored_reads_size_;

    bool verbose_;

    std::vector<TAlphabet> filter_suffix_encoded_;
};


class KMerDBGSuccConstructor : public DBGSuccConstructor {
  public:
    explicit KMerDBGSuccConstructor(size_t k, size_t num_threads = 1)
      : constructor_(k, "", num_threads) {}

    void add_read(const std::string &read) { constructor_.add_read(read); }

    void build_graph(DBG_succ *graph);

    size_t get_k() const { return constructor_.get_k(); }

  private:
    KMerDBGSuccChunkConstructor constructor_;
};


#endif // __DBG_SUCCINCT_CONSTRUCT_HPP__
