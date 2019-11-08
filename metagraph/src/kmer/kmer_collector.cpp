#include "kmer_collector.hpp"

#include <type_traits>
#include <ips4o.hpp>

#include "kmer.hpp"
#include "unix_tools.hpp"
#include "reverse_complement.hpp"

const size_t kMaxKmersChunkSize = 30'000'000;


template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix,
                   bool remove_redundant = true) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(utils::is_instance<Container, SortedSet>{});
    static_assert(std::is_same_v<KMER, typename Container::key_type>);

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

            if (temp_storage.capacity() > 2 * kMaxKmersChunkSize)
                temp_storage = Vector<KMER>(1.1 * kMaxKmersChunkSize);

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

template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(utils::is_instance<Container, SortedMultiset>{});
    static_assert(std::is_same_v<KMER, typename Container::key_type>);

    using KmerCount = typename Container::count_type;

    Vector<KMER> temp_storage;
    Vector<std::pair<KMER, KmerCount>> temp_storage_with_counts;
    temp_storage_with_counts.reserve(1.1 * kMaxKmersChunkSize);

    KmerExtractor kmer_extractor;

    generate_reads([&](const std::string &read, uint64_t count) {
        count = std::min(count, kmers->max_count());

        kmer_extractor.sequence_to_kmers(read, k, suffix, &temp_storage);
        if (both_strands_mode) {
            auto rev_read = read;
            reverse_complement(rev_read.begin(), rev_read.end());
            kmer_extractor.sequence_to_kmers(rev_read, k, suffix, &temp_storage);
        }

        for (const KMER &kmer : temp_storage) {
            temp_storage_with_counts.emplace_back(kmer, count);
        }

        if (temp_storage.capacity() > 2 * kMaxKmersChunkSize)
            temp_storage = Vector<KMER>(1.1 * kMaxKmersChunkSize);

        temp_storage.resize(0);

        if (temp_storage_with_counts.size() > kMaxKmersChunkSize) {
            kmers->insert(temp_storage_with_counts.begin(),
                          temp_storage_with_counts.end());

            if (temp_storage_with_counts.capacity() > 2 * kMaxKmersChunkSize)
                temp_storage_with_counts
                    = Vector<std::pair<KMER, KmerCount>>(1.1 * kMaxKmersChunkSize);

            temp_storage_with_counts.resize(0);
        }
    });

    if (temp_storage_with_counts.size()) {
        kmers->insert(temp_storage_with_counts.begin(),
                      temp_storage_with_counts.end());
    }
}


template <typename KMER, class KmerExtractor, class Container>
KmerStorage<KMER, KmerExtractor, Container>
::KmerStorage(size_t k,
              bool both_strands_mode,
              Sequence&& filter_suffix_encoded,
              size_t num_threads,
              double memory_preallocated,
              bool verbose)
      : k_(k),
        kmers_(num_threads, verbose),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1,
                     std::max(static_cast<size_t>(1), num_threads_)),
        stored_sequences_size_(0),
        verbose_(verbose),
        filter_suffix_encoded_(std::move(filter_suffix_encoded)),
        both_strands_mode_(both_strands_mode) {
    assert(num_threads_ > 0);

    kmers_.reserve(memory_preallocated / sizeof(typename Container::value_type));
    if (verbose_) {
        std::cout << "Preallocated "
                  << (kmers_.data().capacity() * sizeof(typename Container::value_type) >> 30)
                  << "Gb for the k-mer storage"
                  << ", capacity: " << kmers_.data().capacity() << " k-mers"
                  << std::endl;
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
::add_sequence(std::string&& sequence, uint64_t count) {
    if (sequence.size() < k_)
        return;

    // put read into temporary storage
    stored_sequences_size_ += sequence.size();
    buffered_sequences_.emplace_back(std::move(sequence), count);

    if (stored_sequences_size_ < kMaxKmersChunkSize)
        return;

    // extract all k-mers from sequences accumulated in the temporary storage
    release_task_to_pool();

    assert(!stored_sequences_size_);
    assert(!buffered_sequences_.size());
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallString)> &generate_sequences) {
    if constexpr(utils::is_instance<Container, SortedSet>{}) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>, generate_sequences,
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_,
                             true);
    } else {
        thread_pool_.enqueue(count_kmers<KMER, Extractor, Container>,
                             [generate_sequences](CallStringCount callback) {
                                 generate_sequences([&](const std::string &seq) {
                                     callback(seq, 1);
                                 });
                             },
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
    if constexpr(utils::is_instance<Container, SortedSet>{}) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>,
                             [generate_sequences](CallString callback) {
                                 generate_sequences([&](const std::string &seq, uint64_t) {
                                     callback(seq);
                                 });
                             },
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_,
                             true);
    } else {
        thread_pool_.enqueue(count_kmers<KMER, Extractor, Container>,
                             generate_sequences,
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
::insert_dummy(const KMER &dummy_kmer) {
    kmers_.insert(&dummy_kmer, &dummy_kmer + 1);
};

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>::release_task_to_pool() {
    auto *buffered_sequences = new std::vector<std::pair<std::string, uint64_t>>();
    buffered_sequences->swap(buffered_sequences_);

    add_sequences([buffered_sequences](CallStringCount callback) {
        for (const auto& [sequence, count] : *buffered_sequences) {
            callback(sequence, count);
        }
        delete buffered_sequences;
    });

    stored_sequences_size_ = 0;
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>::join() {
    release_task_to_pool();
    thread_pool_.join();
}


#define INSTANTIATE_KMER_STORAGE(KMER_EXTRACTOR, KMER, CONTAINER, CLEANUP) \
    template class KmerStorage<KMER, \
                               KMER_EXTRACTOR, \
                               SortedSet<KMER, CONTAINER<KMER>, CLEANUP>>; \
    template class KmerStorage<KMER, \
                               KMER_EXTRACTOR, \
                               SortedMultiset<KMER, uint8_t, CONTAINER<std::pair<KMER, uint8_t>>, CLEANUP>>; \
    template class KmerStorage<KMER, \
                               KMER_EXTRACTOR, \
                               SortedMultiset<KMER, uint32_t, CONTAINER<std::pair<KMER, uint32_t>>, CLEANUP>>;


INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer64, Vector, utils::NoCleanup)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer128, Vector, utils::NoCleanup)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer256, Vector, utils::NoCleanup)

INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer64, Vector, utils::DummyKmersCleaner)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer128, Vector, utils::DummyKmersCleaner)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer256, Vector, utils::DummyKmersCleaner)

INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer64, Vector, utils::NoCleanup)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer128, Vector, utils::NoCleanup)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer256, Vector, utils::NoCleanup)


// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer64, DequeStorage, utils::NoCleanup)
// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer128, DequeStorage, utils::NoCleanup)
// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer256, DequeStorage, utils::NoCleanup)

// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer64, DequeStorage, utils::DummyKmersCleaner)
// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer128, DequeStorage, utils::DummyKmersCleaner)
// INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer256, DequeStorage, utils::DummyKmersCleaner)

// INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer64, DequeStorage, utils::NoCleanup)
// INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer128, DequeStorage, utils::NoCleanup)
// INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer256, DequeStorage, utils::NoCleanup)
