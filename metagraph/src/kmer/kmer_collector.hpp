#ifndef __KMER_COLLECTOR_HPP__
#define __KMER_COLLECTOR_HPP__

#include <cstddef>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "common/threads/threading.hpp"
#include "common/batch_accumulator.hpp"


namespace mtg {
namespace kmer {

typedef std::function<void(const std::string&)> CallString;
typedef std::function<void(const std::string&, uint64_t)> CallStringCount;

/**
 * Collects k-mers extracted by a KmerExtractor into Container. K-mers are
 * extracted from sequences (reads) that can be added one by one using
 * #add_sequence or they can be provided by a generator passed on to #add_sequences.
 * The resulting de-duped and sorted k-mers can then be accessed via #data().
 * @tparam KMER           The type of k-mers being accumulated by the class, one of
 * KmerExtractor::Kmer64/128/256.
 * @tparam KmerExtractor  Extracts k-mers from reads.
 * @tparam Container      Accumulates the resulting k-mers, can be #SortedSet,
 * #SortedSetDisk, or #SortedMultiset.
 */
template <typename KMER, class KmerExtractor, class Container>
class KmerCollector {
    KmerExtractor kmer_extractor_;

    static_assert(std::is_same_v<typename Container::key_type, typename KMER::WordType>);
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);

  public:
    using Extractor = KmerExtractor;
    using Key = typename Container::key_type;
    using Value = typename Container::value_type;
    using Data = typename Container::result_type;
    using Kmer = KMER;

    enum Mode { BASIC, CANONICAL_ONLY, BOTH };

    /**
     * @param  k The k-mer length
     * @param  mode Collection mode
     * @param  filter_suffix_encoded  Keep only k-mers with the given suffix. Useful for
     * sharding the collection process.
     * @param  num_threads The number of threads in the pool processing incoming sequences
     * @param  memory_preallocated The number of bytes to reserve in the container
     */
    KmerCollector(size_t k,
                  Mode mode = BASIC,
                  std::vector<typename Extractor::TAlphabet>&& filter_suffix_encoded = {},
                  size_t num_threads = 1,
                  double memory_preallocated = 0,
                  const std::filesystem::path &swap_dir = "/tmp/",
                  size_t disk_cap_bytes = 1e9);

    ~KmerCollector();

    inline size_t get_k() const { return k_; }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    void add_sequence(std::string_view sequence, uint64_t count = 1) {
        // push read to the processing queue
        if (sequence.size() >= k_)
            batcher_.push_and_pay(sequence.size(), sequence, count);
    }

    void add_sequence(std::string&& sequence, uint64_t count = 1) {
        // push read to the processing queue
        if (sequence.size() >= k_)
            batcher_.push_and_pay(sequence.size(), std::move(sequence), count);
    }

    size_t buffer_size() const { return buffer_size_; }

    void add_sequences(const std::function<void(CallString)> &generate_sequences);
    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences);
    void add_sequences(std::vector<std::string>&& sequences);
    void add_sequences(std::vector<std::pair<std::string, uint64_t>>&& sequences);

    // FYI: This function should be used only in special cases.
    //      In general, use `add_sequences` if possible, to make use of multiple threads.
    void add_kmer(const KMER &kmer) { kmers_->insert(&kmer.data(), &kmer.data() + 1); }

    // FYI: This returns a container with integer representation of k-mers.
    //      Use reinterpret_cast to cast them back to k-mers.
    inline Data& data() { join(); return kmers_->data(); }

    void clear() { join(); kmers_->clear(); }

    Container& container() { join(); return *kmers_; }

    inline Mode get_mode() const { return mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }
    inline std::filesystem::path tmp_dir() const { return tmp_dir_; }

  private:
    void join();

    size_t k_;
    std::unique_ptr<Container> kmers_;

    size_t num_threads_;
    ThreadPool thread_pool_;

    BatchAccumulator<std::pair<std::string, uint64_t>> batcher_;

    std::vector<typename Extractor::TAlphabet> filter_suffix_encoded_;

    Mode mode_;

    std::filesystem::path tmp_dir_;

    size_t buffer_size_;
};

/** Visible For Testing */
template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   typename KmerCollector<KMER, KmerExtractor, Container>::Mode mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix);

/** Visible For Testing */
template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 typename KmerCollector<KMER, KmerExtractor, Container>::Mode mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix);

} // namespace kmer
} // namespace mtg

#endif // __KMER_COLLECTOR_HPP__
