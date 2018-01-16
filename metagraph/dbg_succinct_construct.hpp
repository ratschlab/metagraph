#ifndef __DBG_SUCCINCT_CONSTRUCT_HPP__
#define __DBG_SUCCINCT_CONSTRUCT_HPP__

#include "dbg_succinct.hpp"
#include "kmer.hpp"


class DBGSuccConstructor {
  public:
    virtual void add_read(const std::string &read) = 0;
    virtual void add_reads(const std::vector<std::string> &reads) = 0;

    virtual void build_graph(DBG_succ *graph) = 0;

    virtual size_t get_k() const = 0;
};


class KMerDBGSuccConstructor : public DBGSuccConstructor {
  public:
    explicit KMerDBGSuccConstructor(size_t k, size_t num_threads = 1);

    void add_read(const std::string &read);

    void add_reads(const std::vector<std::string> &reads);

    void build_graph(DBG_succ *graph);

    size_t get_k() const { return k_; }

  private:
    size_t k_;
    std::vector<KMer> kmers_;
    size_t num_threads_;
};


//TODO: implement
class SuffixArrayDBGSuccConstructor : public DBGSuccConstructor {
  public:
    explicit SuffixArrayDBGSuccConstructor(size_t k);

    void add_reads(const std::vector<std::string> &reads);

    void construct_graph(DBG_succ *graph);

    size_t get_k() const { return k_; }

  private:
    size_t k_;
    std::string data_;
};


class KMerDBGSuccChunkConstructor {
  public:
    KMerDBGSuccChunkConstructor(size_t k,
                                const std::string &filter_suffix,
                                size_t num_threads = 1);

    void add_read(const std::string &sequence);

    DBG_succ::Chunk* build_chunk();

  private:
    size_t k_;
    std::vector<KMer> kmers_;
    size_t num_threads_;
    std::vector<TAlphabet> filter_suffix_encoded_;
};


#endif // __DBG_SUCCINCT_CONSTRUCT_HPP__
