#ifndef __KMER_COLLECTOR_HPP__
#define __KMER_COLLECTOR_HPP__

#include "kmer_extractor.hpp"
#include "kmer_collector.hpp"
#include "threading.hpp"
#include "sorted_set.hpp"

typedef std::function<void(const std::string&)> CallString;


template <typename KMER, class KmerExtractor>
class KmerCollector {
    using Extractor = KmerExtractor;
    using Sequence = std::vector<typename Extractor::TAlphabet>;
    Extractor kmer_extractor_;
  public:
    explicit KmerCollector(size_t k,
                           bool both_strands_mode = false,
                           Sequence&& filter_suffix_encoded = {},
                           size_t num_threads = 1,
                           double memory_preallocated = 0,
                           bool verbose = false);

    inline size_t get_k() const { return k_; }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    void add_sequence(std::string&& sequence);

    void add_sequences(const std::function<void(CallString)> &generate_sequences);

    inline Vector<KMER>& data() { join(); return kmers_.data(); }

    void clear() { join(); kmers_.clear(); }

    inline bool verbose() const { return verbose_; }
    inline bool is_both_strands_mode() const { return both_strands_mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }

  private:
    void release_task_to_pool();
    void join();

    size_t k_;
    SortedSet<KMER> kmers_;

    size_t num_threads_;
    ThreadPool thread_pool_;

    std::vector<std::string> sequences_storage_;
    size_t stored_sequences_size_;

    bool verbose_;

    Sequence filter_suffix_encoded_;

    bool both_strands_mode_;
};


template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads = 1,
                                size_t offset = 0);

#endif // __KMER_COLLECTOR_HPP__