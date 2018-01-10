#ifndef __DBG_SUCCINCT_CHUNK_HPP__
#define __DBG_SUCCINCT_CHUNK_HPP__

#include <type_traits>

#include "dbg_succinct.hpp"

class KMer;


class DBG_succ::Chunk {
  public:
    virtual ~Chunk() {}

    virtual void push_back(TAlphabet W, TAlphabet F, bool last) = 0;
    virtual TAlphabet get_W_back() const = 0;
    virtual void alter_W_back(TAlphabet W) = 0;
    virtual void alter_last_back(bool last) = 0;

    virtual void extend(const Chunk &other) = 0;

    virtual uint64_t size() const = 0;

    virtual void initialize_graph(DBG_succ *graph) = 0;

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;
};


class DBG_succ::VectorChunk : public DBG_succ::Chunk {
  public:
    void push_back(TAlphabet W, TAlphabet F, bool last);
    TAlphabet get_W_back() const;
    void alter_W_back(TAlphabet W);
    void alter_last_back(bool last);

    void extend(const Chunk &other);
    void extend(const VectorChunk &other);

    uint64_t size() const;

    void initialize_graph(DBG_succ *graph);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    static VectorChunk* build_from_kmers(size_t k,
                                         std::vector<KMer> *kmers,
                                         unsigned int parallel = 1);

  private:
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_ = std::vector<uint64_t>(DBG_succ::alph_size, 0);
};


std::vector<TAlphabet> prepare_sequence(const std::string &sequence,
                                        size_t k = 0, bool add_bridge = false);


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
void sequence_to_kmers(const std::vector<TAlphabet> &seq,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix = {});


#endif // __DBG_SUCCINCT_CHUNK_HPP__
