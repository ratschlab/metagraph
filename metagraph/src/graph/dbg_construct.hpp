#ifndef __DBG_CONSTRUCT_HPP__
#define __DBG_CONSTRUCT_HPP__

#include <mutex>
#include <shared_mutex>

#include "kmer_extractor.hpp"


typedef std::function<void(const std::string&)> CallbackString;


template <class GraphChunk>
class IChunkConstructor {
  public:
    virtual ~IChunkConstructor() {}

    virtual void add_sequence(const std::string &sequence) = 0;

    virtual void add_sequences(std::function<void(CallbackString)> generate_sequences) = 0;

    virtual GraphChunk* build_chunk() = 0;
};


class GraphConstructor {
  public:
    virtual ~GraphConstructor() {}
};


template <typename KMER, class KmerExtractor>
class KmerCollector {
    using Extractor = KmerExtractor;
    using Sequence = std::vector<typename Extractor::TAlphabet>;
    Extractor kmer_extractor_;
  public:
    explicit KmerCollector(size_t k,
                           bool canonical_mode = false,
                           Sequence&& filter_suffix_encoded = {},
                           size_t num_threads = 1,
                           double memory_preallocated = 0,
                           bool verbose = false);

    inline size_t get_k() const { return k_; }

    inline size_t size() const { return kmers_.size(); }

    inline size_t suffix_length() const { return filter_suffix_encoded_.size(); }

    // TODO: another for std::string&& ?
    void add_sequence(const std::string &sequence);

    void add_sequences(const std::function<void(CallbackString)> &generate_sequences);

    template <typename... Args>
    inline void emplace_back(Args&&... args) {
        kmers_.emplace_back(std::forward<Args>(args)...);
    }

    inline Vector<KMER>& data() { return kmers_; }

    inline void
    call_kmers(const std::function<void(const KMER&)> &kmer_callback) const {
        std::for_each(kmers_.begin(), kmers_.end(), kmer_callback);
    }

    inline void clear() {
        kmers_.clear();
        sequences_storage_.clear();
        stored_sequences_size_ = 0;
        kmers_.shrink_to_fit();
    }

    void join();

    inline bool verbose() const { return verbose_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }
    inline size_t num_threads() const { return num_threads_; }
    inline size_t alphabet_size() const { return kmer_extractor_.alphabet.size(); }

  private:
    void release_task_to_pool();

    size_t k_;
    Vector<KMER> kmers_;
    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;

    size_t num_threads_;
    utils::ThreadPool thread_pool_;

    std::vector<std::string> sequences_storage_;
    size_t stored_sequences_size_;

    bool verbose_;

    Sequence filter_suffix_encoded_;

    bool canonical_mode_;
};


template <typename KMER>
void shrink_kmers(Vector<KMER> *kmers,
                  size_t num_threads,
                  bool verbose,
                  size_t offset = 0);

template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads = 1,
                                size_t offset = 0);

#endif // __DBG_CONSTRUCT_HPP__
