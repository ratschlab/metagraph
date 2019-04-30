#include "dbg_construct.hpp"

#include <type_traits>
#include <ips4o.hpp>

#include "kmer.hpp"
#include "unix_tools.hpp"
#include "reads_filtering.hpp"
#include "reverse_complement.hpp"

using TAlphabet = KmerExtractor::TAlphabet;

const size_t kMaxKmersChunkSize = 30'000'000;


template <class V>
void try_reserve(V *vector, size_t size, size_t min_size = 0) {
    size = std::max(size, min_size);

    while (size > min_size) {
        try {
            vector->reserve(size);
            return;
        } catch (const std::bad_alloc &exception) {
            size = min_size + (size - min_size) * 2 / 3;
        }
    }
    vector->reserve(min_size);
}

template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads,
                                size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          std::less<typename V::value_type>(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(array->begin() + offset, array->end());
    array->erase(unique_end, array->end());
}

template <typename KMER>
void shrink_kmers(Vector<KMER> *kmers,
                  size_t num_threads,
                  bool verbose,
                  size_t offset) {
    if (verbose) {
        std::cout << "Allocated capacity exceeded, filter out non-unique k-mers..."
                  << std::flush;
    }

    size_t prev_num_kmers = kmers->size();
    sort_and_remove_duplicates(kmers, num_threads, offset);

    if (verbose) {
        std::cout << " done. Number of kmers reduced from " << prev_num_kmers
                                                  << " to " << kmers->size() << ", "
                  << (kmers->size() * sizeof(KMER) >> 20) << "Mb" << std::endl;
    }
}

template <class Array, class Vector>
void extend_kmer_storage(const Array &temp_storage,
                         Vector *kmers,
                         size_t num_threads,
                         bool verbose,
                         std::mutex &mutex_resize,
                         std::shared_timed_mutex &mutex_copy) {
    // acquire the mutex to restrict the number of writing threads
    std::unique_lock<std::mutex> resize_lock(mutex_resize);

    if (kmers->size() + temp_storage.size() > kmers->capacity()) {
        std::unique_lock<std::shared_timed_mutex> reallocate_lock(mutex_copy);

        shrink_kmers(kmers, num_threads, verbose);

        try {
            try_reserve(kmers, kmers->size() + kmers->size() / 2,
                               kmers->size() + temp_storage.size());
        } catch (const std::bad_alloc &exception) {
            std::cerr << "ERROR: Can't reallocate. Not enough memory" << std::endl;
            exit(1);
        }
    }

    size_t offset = kmers->size();
    kmers->resize(kmers->size() + temp_storage.size());

    std::shared_lock<std::shared_timed_mutex> copy_lock(mutex_copy);

    resize_lock.unlock();

    std::copy(temp_storage.begin(),
              temp_storage.end(),
              kmers->begin() + offset);
}


template <typename KMER, class KmerExtractor>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Vector<KMER> *kmers,
                   const std::vector<TAlphabet> &suffix,
                   size_t num_threads,
                   bool verbose,
                   std::mutex &mutex_resize,
                   std::shared_timed_mutex &mutex_copy,
                   bool remove_redundant = true) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::kLogSigma);

    Vector<KMER> temp_storage;
    temp_storage.reserve(1.1 * kMaxKmersChunkSize);

    KmerExtractor kmer_extractor;

    generate_reads([&](const std::string &read) {
        kmer_extractor.sequence_to_kmers(read, k, suffix, &temp_storage);
        if (both_strands_mode) {
            auto rev_read = read;
            reverse_complement(rev_read.begin(), rev_read.end());
            kmer_extractor.sequence_to_kmers(rev_read, k, suffix, &temp_storage);
        }

        if (temp_storage.size() < kMaxKmersChunkSize)
            return;

        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage);
        }

        if (temp_storage.size() > 0.9 * kMaxKmersChunkSize) {
            extend_kmer_storage(temp_storage, kmers,
                                num_threads, verbose, mutex_resize, mutex_copy);
            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        if (remove_redundant) {
            sort_and_remove_duplicates(&temp_storage);
        }
        extend_kmer_storage(temp_storage, kmers,
                            num_threads, verbose, mutex_resize, mutex_copy);
    }
}

template <typename KMER, class KmerExtractor>
KmerCollector<KMER, KmerExtractor>
::KmerCollector(size_t k,
                bool both_strands_mode,
                Sequence&& filter_suffix_encoded,
                size_t num_threads,
                double memory_preallocated,
                bool verbose)
      : k_(k),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1,
                     std::max(static_cast<size_t>(1), num_threads_)),
        stored_sequences_size_(0),
        verbose_(verbose),
        filter_suffix_encoded_(std::move(filter_suffix_encoded)),
        both_strands_mode_(both_strands_mode) {
    assert(num_threads_ > 0);
    static_assert(KMER::kBitsPerChar == KmerExtractor::kLogSigma);

    try_reserve(&kmers_, memory_preallocated / sizeof(KMER));
    if (verbose_) {
        std::cout << "Preallocated "
                  << kmers_.capacity() * sizeof(KMER) / (1llu << 30)
                  << "Gb for the k-mer storage"
                  << ", capacity: " << kmers_.capacity() << " k-mers"
                  << std::endl;
    }
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>
::add_sequence(const std::string &sequence) {
    if (sequence.size() < k_)
        return;

    // put read into temporary storage
    stored_sequences_size_ += sequence.size();
    sequences_storage_.emplace_back(sequence);

    if (stored_sequences_size_ < kMaxKmersChunkSize)
        return;

    // extract all k-mers from sequences accumulated in the temporary storage
    release_task_to_pool();

    assert(!stored_sequences_size_);
    assert(!sequences_storage_.size());
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>
::add_sequences(const std::function<void(CallString)> &generate_sequences) {
    thread_pool_.enqueue(extract_kmers<KMER, Extractor>, generate_sequences,
                         k_, both_strands_mode_, &kmers_,
                         filter_suffix_encoded_,
                         num_threads_, verbose_,
                         std::ref(mutex_resize_), std::ref(mutex_copy_), true);
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>::release_task_to_pool() {
    auto *current_sequences_storage = new std::vector<std::string>();
    current_sequences_storage->swap(sequences_storage_);

    thread_pool_.enqueue(extract_kmers<KMER, Extractor>,
                         [current_sequences_storage](CallString callback) {
                             for (auto &&sequence : *current_sequences_storage) {
                                 callback(std::move(sequence));
                             }
                             delete current_sequences_storage;
                         },
                         k_, both_strands_mode_, &kmers_, filter_suffix_encoded_,
                         num_threads_, verbose_,
                         std::ref(mutex_resize_), std::ref(mutex_copy_), true);
    stored_sequences_size_ = 0;
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>::join() {
    release_task_to_pool();
    thread_pool_.join();

    if (verbose_) {
        std::cout << "Reading data has finished" << std::endl;
        std::cout << "Sorting k-mers and appending succinct"
                  << " representation from current bin...\t" << std::flush;
    }
    Timer timer;

    sort_and_remove_duplicates(&kmers_, num_threads_);

    if (verbose_)
        std::cout << timer.elapsed() << "sec" << std::endl;
}

template class KmerCollector<KmerExtractor::Kmer64, KmerExtractor>;
template class KmerCollector<KmerExtractor::Kmer128, KmerExtractor>;
template class KmerCollector<KmerExtractor::Kmer256, KmerExtractor>;
template class KmerCollector<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit>;
template class KmerCollector<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit>;
template class KmerCollector<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit>;
