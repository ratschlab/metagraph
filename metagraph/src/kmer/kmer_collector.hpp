#ifndef __KMER_COLLECTOR_HPP__
#define __KMER_COLLECTOR_HPP__

#include "common/threads/threading.hpp"
#include "common/batch_accumulator.hpp"

namespace mg {
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
    using Extractor = KmerExtractor;
    using Sequence = std::vector<typename Extractor::TAlphabet>;
    Extractor kmer_extractor_;

    static_assert(std::is_base_of<typename Container::key_type, KMER>::value);
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);

  public:
    using Key = typename Container::key_type;
    using Value = typename Container::value_type;
    using Data = typename Container::result_type;

    /**
     * @param  k                      The k-mer length
     * @param  both_strands_mode      If true, both a sequence and its
     * reverse complement will be added (only makes sense for DNA sequences)
     * @param  filter_suffix_encoded  Keep only k-mers with the given
     * suffix. Useful for sharding the collection process.
     * @param  num_threads            The number of threads in the pool
     * processing incoming sequences
     * @param  memory_preallocated    The number of bytes to reserve in
     * the container
     */
    KmerCollector(size_t k,
                  bool both_strands_mode = false,
                  Sequence&& filter_suffix_encoded = {},
                  size_t num_threads = 1,
                  double memory_preallocated = 0);

    inline size_t get_k() const { return k_; }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    void add_sequence(std::string_view sequence, uint64_t count = 1) {
        // push read to the processing queue
        if (sequence.size() >= k_)
            batch_accumulator_.push_and_pay(sequence.size() - k_ + 1, sequence, count);
    }

    void add_sequence(std::string&& sequence, uint64_t count = 1) {
        // push read to the processing queue
        if (sequence.size() >= k_)
            batch_accumulator_.push_and_pay(sequence.size() - k_ + 1, std::move(sequence), count);
    }

    void add_sequences(const std::function<void(CallString)> &generate_sequences);
    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences);

    // FYI: This function should be used only in special cases.
    //      In general, use `add_sequences` if possible, to make use of multiple threads.
    void add_kmer(const KMER &kmer) { kmers_.insert(&kmer, &kmer + 1); }

    inline Data &data() { join(); return kmers_.data(); }

    void clear() { join(); kmers_.clear(); }

    inline bool is_both_strands_mode() const { return both_strands_mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }

    /**
     * Sends sequences accumulated in #batch_accumulator_ for processing
     * on the thread pool. */
    void add_batch(std::vector<std::pair<std::string, uint64_t>>&& sequences);

  private:
    void join();

    size_t k_;
    Container kmers_;

    size_t num_threads_;
    ThreadPool thread_pool_;

    BatchAccumulator<std::pair<std::string, uint64_t>, size_t> batch_accumulator_;

    Sequence filter_suffix_encoded_;

    bool both_strands_mode_;
};

/** Visible For Testing */
template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix,
                   bool remove_redundant = true);

/** Visible For Testing */
template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix);

} // namespace kmer
} // namespace mg

#endif // __KMER_COLLECTOR_HPP__
