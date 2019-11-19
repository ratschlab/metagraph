#ifndef __KMER_COLLECTOR_HPP__
#define __KMER_COLLECTOR_HPP__

#include "common/threading.hpp"

typedef std::function<void(const std::string&)> CallString;
typedef std::function<void(const std::string&, uint64_t)> CallStringCount;

/**
 * Collects k-mers extracted by a KmerExtractor into Container.
 * @tparam KMER
 * @tparam KmerExtractor
 * @tparam Container
 */
template <typename KMER, class KmerExtractor, class Container>
class KmerStorage {
    using Extractor = KmerExtractor;
    using Sequence = std::vector<typename Extractor::TAlphabet>;
    Extractor kmer_extractor_;

    static_assert(std::is_base_of<typename Container::key_type, KMER>::value);
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);

  public:
    using Key = typename Container::key_type;
    using Value = typename Container::value_type;
    using Data = typename Container::storage_type;

    KmerStorage(size_t k,
                bool both_strands_mode = false,
                Sequence&& filter_suffix_encoded = {},
                size_t num_threads = 1,
                double memory_preallocated = 0,
                bool verbose = false);

    inline size_t get_k() const { return k_; }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    void add_sequence(std::string&& sequence, uint64_t count = 1);

    void add_sequences(const std::function<void(CallString)> &generate_sequences);
    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences);

    void insert_dummy(const KMER &dummy_kmer);

    inline Data& data() { join(); return kmers_.data(); }

    void clear() { join(); kmers_.clear(); }

    inline bool verbose() const { return verbose_; }
    inline bool is_both_strands_mode() const { return both_strands_mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }

  private:
    void release_task_to_pool();
    void join();

    size_t k_;
    Container kmers_;

    size_t num_threads_;
    ThreadPool thread_pool_;

    std::vector<std::pair<std::string, uint64_t>> buffered_sequences_;
    size_t stored_sequences_size_;

    bool verbose_;

    Sequence filter_suffix_encoded_;

    bool both_strands_mode_;
};

#endif // __KMER_COLLECTOR_HPP__
