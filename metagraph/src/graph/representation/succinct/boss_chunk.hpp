#ifndef __BOSS_CHUNK_HPP__
#define __BOSS_CHUNK_HPP__

#include <type_traits>

#include <sdsl/int_vector_buffer.hpp>

#include "boss.hpp"

/**
 * Structure storing data for the BOSS table (or part of it), without support for
 * 'rank' and 'select' operations.
 */
class BOSS::Chunk {
  public:
    typedef uint8_t TAlphabet;

    /**
     * Creates an empty BOSS chunk with the given parameters
     * @param alph_size the alphabet size. For DNA this will be 5 ($,A,C,G,T).
     * @param k k-mer size
     * @param canonical if true, the BOSS table will be constructed with both a k-mer
     *        and its reverse complement
     * @param swap_dir directory where to write vector buffers for construction with
     * streaming
     */
    Chunk(uint64_t alph_size, size_t k, bool canonical,
          const std::string &swap_dir = "");

    /**
     * Creates a BOSS Chunk with k-mer counts. Assumes that k-mers are distinct and
     * sorted.
     * @param alph_size the alphabet size. For DNA this will be 5 ($,A,C,G,T).
     * @param k k-mer size
     * @param canonical if true, the BOSS table will be constructed with both a k-mer
     *        and its reverse complement
     * @param kmers_with_counts the k-mers and their counts to construct the chunk from
     * @param bits_per_count for weighted graphs, the number of bits used to store the
     * weight counts
     * @param swap_dir directory where to write vector buffers for construction with
     * streaming
     */
    template <typename Array>
    Chunk(uint64_t alph_size,
          size_t k,
          bool canonical,
          const Array &kmers_with_counts,
          uint8_t bits_per_count = 0,
          const std::string &swap_dir = "");

    ~Chunk();

    /**
     * Adds an entry into the BOSS table.
     * @param W edge label
     * @param F last character
     * @param last true if this is the last outgoing edge from that node
     */
    void push_back(TAlphabet W, TAlphabet F, bool last);

    TAlphabet get_W_back() const { return W_[size_ - 1]; }
    void alter_W_back(TAlphabet W) { W_[size_ - 1] = W; }

    void alter_last_back(bool last) { last_[size_ - 1] = last; }

    void extend(const Chunk &other);

    uint64_t size() const { return size_; }

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    void initialize_boss(BOSS *graph, sdsl::int_vector<> *weights = nullptr);

    /**
     * Merge BOSS chunks loaded from the files passed in #chunk_filenames and construct
     * the full BOSS table.
     */
    static std::pair<BOSS* /* boss */, bool /* is_canonical */>
    build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                           bool verbose = false,
                           sdsl::int_vector<> *weights = nullptr,
                           const std::string &swap_dir = "");

    static constexpr auto kFileExtension = ".dbg.chunk";

  private:
    uint8_t get_W_width() const;

    size_t alph_size_;
    size_t k_;
    bool canonical_;
    // see the BOSS paper for the meaning of W_, last_ and F_
    std::string dir_;
    sdsl::int_vector_buffer<> W_;
    sdsl::int_vector_buffer<1> last_;
    std::vector<uint64_t> F_;
    sdsl::int_vector_buffer<> weights_;
    uint64_t size_;
};


#endif // __BOSS_CHUNK_HPP__
