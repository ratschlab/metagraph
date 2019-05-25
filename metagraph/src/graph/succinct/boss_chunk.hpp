#ifndef __BOSS_CHUNK_HPP__
#define __BOSS_CHUNK_HPP__

#include <type_traits>

#include "boss.hpp"


class BOSS::Chunk {
  public:
    typedef uint8_t TAlphabet;

    explicit Chunk(size_t k);

    /**
     * Assumes that |W| and |last| have dummy 0 at the first position
     */
    Chunk(size_t k,
          std::vector<TAlphabet>&& W,
          std::vector<bool>&& last,
          std::vector<uint64_t>&& F,
          sdsl::int_vector<>&& weights = sdsl::int_vector<>());

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
    const size_t alph_size_;
    const size_t bits_per_char_W_;
    size_t k_;
    std::vector<TAlphabet> W_;
    std::vector<bool> last_;
    std::vector<uint64_t> F_;
    sdsl::int_vector<> weights_;
};


#endif // __BOSS_CHUNK_HPP__
