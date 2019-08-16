#ifndef __BOSS_CHUNK_HPP__
#define __BOSS_CHUNK_HPP__

#include <type_traits>

#include "boss.hpp"


class BOSS::Chunk {
  public:
    typedef uint8_t TAlphabet;

    // Alphabet size without extra characters
    Chunk(uint64_t alph_size, size_t k);

    /**
     * Assumes that kmers are distinct and sorted
     */
    template <typename KMER, typename COUNT>
    Chunk(uint64_t alph_size,
          size_t k,
          const Vector<std::pair<KMER, COUNT>> &kmers,
          uint8_t bits_per_count = 8);

    template <typename KMER, typename COUNT>
    Chunk(uint64_t alph_size,
          size_t k,
          const utils::DequeStorage<std::pair<KMER, COUNT>> &kmers,
          uint8_t bits_per_count = 8);

    template <typename KMER>
    Chunk(uint64_t alph_size, size_t k, const Vector<KMER> &kmers);

    template <typename KMER>
    Chunk(uint64_t alph_size, size_t k, const utils::DequeStorage<KMER> &kmers);

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
     * Merge BOSS chunks loaded from the files passed
     */
    static BOSS*
    build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                           bool verbose = false,
                           sdsl::int_vector<> *weights = nullptr);

    static constexpr auto kFileExtension = ".dbg.chunk";

  private:
    uint8_t extended_alph_size() const;

    size_t alph_size_;
    size_t k_;
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_;
    sdsl::int_vector<> weights_;
};


#endif // __BOSS_CHUNK_HPP__
