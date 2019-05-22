#include "kmer_collector.hpp"

#include <type_traits>
#include <ips4o.hpp>

#include "kmer.hpp"
#include "unix_tools.hpp"
#include "reads_filtering.hpp"
#include "reverse_complement.hpp"

using TAlphabet = KmerExtractor::TAlphabet;

const size_t kMaxKmersChunkSize = 30'000'000;


template <typename KMER, class KmerExtractor>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   SortedSet<KMER> *kmers,
                   const std::vector<TAlphabet> &suffix,
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
            kmers->sort_and_remove_duplicates(&temp_storage, 1);
        }

        if (temp_storage.size() > 0.9 * kMaxKmersChunkSize) {
            kmers->insert(temp_storage.begin(), temp_storage.end());
            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        if (remove_redundant) {
            kmers->sort_and_remove_duplicates(&temp_storage, 1);
        }
        kmers->insert(temp_storage.begin(), temp_storage.end());
    }
}

template <typename KMER, class KmerExtractor>
KmerCollector<KMER, KmerExtractor>
::KmerCollector(size_t k,
                bool both_strands_mode,
                Sequence&& filter_suffix_encoded,
                size_t num_threads,
                double memory_preallocated,
                bool verbose,
                std::function<void(Vector<KMER>*)> cleanup)
      : k_(k),
        kmers_(num_threads, verbose, cleanup),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1,
                     std::max(static_cast<size_t>(1), num_threads_)),
        stored_sequences_size_(0),
        verbose_(verbose),
        filter_suffix_encoded_(std::move(filter_suffix_encoded)),
        both_strands_mode_(both_strands_mode) {
    assert(num_threads_ > 0);
    static_assert(KMER::kBitsPerChar == KmerExtractor::kLogSigma);

    kmers_.reserve(memory_preallocated / sizeof(KMER));
    if (verbose_) {
        std::cout << "Preallocated "
                  << kmers_.data().capacity() * sizeof(KMER) / (1llu << 30)
                  << "Gb for the k-mer storage"
                  << ", capacity: " << kmers_.data().capacity() << " k-mers"
                  << std::endl;
    }
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>
::add_sequence(std::string&& sequence) {
    if (sequence.size() < k_)
        return;

    // put read into temporary storage
    stored_sequences_size_ += sequence.size();
    sequences_storage_.push_back(std::move(sequence));

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
                         true);
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
                         true);
    stored_sequences_size_ = 0;
}

template <typename KMER, class KmerExtractor>
void KmerCollector<KMER, KmerExtractor>::join() {
    release_task_to_pool();
    thread_pool_.join();
}

template class KmerCollector<KmerExtractor::Kmer64, KmerExtractor>;
template class KmerCollector<KmerExtractor::Kmer128, KmerExtractor>;
template class KmerCollector<KmerExtractor::Kmer256, KmerExtractor>;
template class KmerCollector<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit>;
template class KmerCollector<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit>;
template class KmerCollector<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit>;
