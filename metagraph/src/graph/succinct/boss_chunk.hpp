#ifndef __BOSS_CHUNK_HPP__
#define __BOSS_CHUNK_HPP__

#include <type_traits>

#include "boss.hpp"

/**
 * Represents a chunk of the BOSS representation of a deBrujin graph, i.e. the BOSS
 * table resulting from processing a subset of sorted kmers. Multiple Chunks can be
 * unified to form a BOSS table.
 * Unlike #BOSS, which supports rank and select operations on its #W_, #F_, last_
 * members, BOSS::Chunk simply stores #W_, #F_, last_ in an std::vector.
 */
class BOSS::Chunk {
  public:
    typedef uint8_t TAlphabet;

    /**
     * Creates an empty BOSS chunk with the given parameters
     * @param alph_size the alphabet size, without extra characters (i.e. sentinel).
     *        For DNA this will be 4.
     * @param k k-mer size
     * @param canonical if true, the BOSS table will be constructed with both a k-mer
     *        and its reverse complement
     */
    Chunk(uint64_t alph_size, size_t k, bool canonical);

    /**
     * Creates a BOSS Chunk with weights from the given kmers. Assumes that kmers are
     * distinct and sorted.
     * @param alph_size the alphabet size, without extra characters (i.e. sentinel).
     * For DNA this will be 4.
     * @param k k-mer size
     * @param canonical if true, build the chunk from canonicalized k-mers
     * @param kmers_with_counts the kmers to construct the chunk from
     * @param bits_per_count for weighted graphs, the number of bits used to store the
     * weight
     */
    template <typename Array>
    Chunk(uint64_t alph_size,
          size_t k,
          bool canonical,
          Array &kmers_with_counts,
          uint8_t bits_per_count);

    /**
     * Creates a BOSS Chunk (no weights) from the given kmers, which should be distinct
     * and sorted.
     */
    template <typename Array>
    Chunk(uint64_t alph_size, size_t k, bool canonical, Array &kmers);

    /**
     * Adds an entry into the BOSS table.
     * @param W the edge label
     * @param F the node label
     * @param last true if this is the last outgoing edge from F
     */
    void push_back(TAlphabet W, TAlphabet F, bool last);

    TAlphabet get_W_back() const { return W_.back(); }
    void alter_W_back(TAlphabet W) { W_.back() = W; }

    void alter_last_back(bool last) { last_.back() = last; }

    void extend(const Chunk &other);

    uint64_t size() const { return W_.size() - 1; }

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    void initialize_boss(BOSS *graph, sdsl::int_vector<> *weights = nullptr);

    /**
     * Merge BOSS chunks loaded from the files passed in #chunk_filenames into a DeBrujin
     * graph represented as a BOSS table.
     */
    static std::pair<BOSS* /* boss */, bool /* is_canonical */>
    build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                           bool verbose = false,
                           sdsl::int_vector<> *weights = nullptr);

    static constexpr auto kFileExtension = ".dbg.chunk";

  private:
    uint8_t get_W_width() const;

    size_t alph_size_;
    size_t k_;
    bool canonical_;
    // see the BOSS paper for the meaning of W_, last_ and F_
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_;
    sdsl::int_vector<> weights_;
};


#endif // __BOSS_CHUNK_HPP__
