#ifndef __DBG_SUCCINCT_CHUNK_HPP__
#define __DBG_SUCCINCT_CHUNK_HPP__

#include <type_traits>

#include "dbg_succinct.hpp"

class KMer;


class DBG_succ::Chunk {
  public:
    Chunk();

    void push_back(TAlphabet W, TAlphabet F, bool last);
    TAlphabet get_W_back() const;
    void alter_W_back(TAlphabet W);
    void alter_last_back(bool last);

    void extend(const Chunk &other);

    uint64_t size() const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    void initialize_graph(DBG_succ *graph) const;

    /**
     * Merge graph chunks from the vector
     * passed and release the chunks afterwards
     */
    static DBG_succ* build_graph_from_chunks(size_t k,
                                const std::vector<Chunk*> &graph_chunks,
                                bool verbose = false);

    /**
     * Initialize graph chunk from a list of sorted kmers.
     */
    static Chunk* build_from_kmers(size_t k, std::vector<KMer> *kmers);

  private:
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_;
};


#endif // __DBG_SUCCINCT_CHUNK_HPP__
