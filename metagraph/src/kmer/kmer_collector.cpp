#include "kmer_collector.hpp"

#include <type_traits>
#include <ips4o.hpp>

#include "kmer.hpp"
#include "unix_tools.hpp"
#include "reverse_complement.hpp"

const size_t kMaxKmersChunkSize = 30'000'000;


template <typename KMER, class KmerExtractor>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   SortedSet<KMER> *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix,
                   bool remove_redundant = true) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);

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

template <typename KMER, class KmerExtractor, typename KmerCount>
void count_kmers(std::function<void(CallString)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 SortedMultiset<KMER, KmerCount> *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);

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

        if (temp_storage.size() > kMaxKmersChunkSize) {
            kmers->insert(temp_storage.begin(), temp_storage.end());
            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        kmers->insert(temp_storage.begin(), temp_storage.end());
    }
}


template <typename KMER, class KmerExtractor, class Container>
KmerStorage<KMER, KmerExtractor, Container>
::KmerStorage(size_t k,
              bool both_strands_mode,
              Sequence&& filter_suffix_encoded,
              size_t num_threads,
              double memory_preallocated,
              bool verbose,
              std::function<void(Data*)> cleanup)
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
::add_kmer(std::string&& kmer, uint32_t count) {
    assert(kmer.size() == k_);
    KmerExtractor kmer_extractor;
    Vector<KMER> temp_storage;
    kmer_extractor.sequence_to_kmers(kmer, k_, filter_suffix_encoded_, &temp_storage);
    assert(temp_storage.size() == 1);
    if constexpr(std::is_base_of<SortedSet<KMER>, Container>::value) {
        (void)count;
        kmers_.insert(temp_storage.begin(), temp_storage.end());
    } else {
        typename Container::count_type c = count;
        kmers_.insert(temp_storage.begin(), temp_storage.end(), &c);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
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

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallString)> &generate_sequences) {
    if constexpr(std::is_base_of<SortedSet<KMER>, Container>::value) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor>, generate_sequences,
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_,
                             true);
    } else {
        thread_pool_.enqueue(count_kmers<KMER, Extractor, typename Container::count_type>,
                             generate_sequences,
                             k_, both_strands_mode_, &kmers_,
                             filter_suffix_encoded_);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>::insert_dummy(KMER dummy_kmer) {
    Vector<KMER> dummy_kmers = { dummy_kmer };
    kmers_.insert(dummy_kmers.begin(), dummy_kmers.end());
};

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>::release_task_to_pool() {
    auto *current_sequences_storage = new std::vector<std::string>();
    current_sequences_storage->swap(sequences_storage_);

    add_sequences([current_sequences_storage](CallString callback) {
        for (auto &&sequence : *current_sequences_storage) {
            callback(std::move(sequence));
        }
        delete current_sequences_storage;
    });

    stored_sequences_size_ = 0;
}

template <typename KMER, class KmerExtractor, class Container>
void KmerStorage<KMER, KmerExtractor, Container>::join() {
    release_task_to_pool();
    thread_pool_.join();
}

template class KmerStorage<KmerExtractor::Kmer64, KmerExtractor, SortedSet<KmerExtractor::Kmer64>>;
template class KmerStorage<KmerExtractor::Kmer128, KmerExtractor, SortedSet<KmerExtractor::Kmer128>>;
template class KmerStorage<KmerExtractor::Kmer256, KmerExtractor, SortedSet<KmerExtractor::Kmer256>>;
template class KmerStorage<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit, SortedSet<KmerExtractor2Bit::Kmer64>>;
template class KmerStorage<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit, SortedSet<KmerExtractor2Bit::Kmer128>>;
template class KmerStorage<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit, SortedSet<KmerExtractor2Bit::Kmer256>>;

template class KmerStorage<KmerExtractor::Kmer64, KmerExtractor, SortedMultiset<KmerExtractor::Kmer64, uint8_t>>;
template class KmerStorage<KmerExtractor::Kmer128, KmerExtractor, SortedMultiset<KmerExtractor::Kmer128, uint8_t>>;
template class KmerStorage<KmerExtractor::Kmer256, KmerExtractor, SortedMultiset<KmerExtractor::Kmer256, uint8_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer64, uint8_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer128, uint8_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer256, uint8_t>>;

template class KmerStorage<KmerExtractor::Kmer64, KmerExtractor, SortedMultiset<KmerExtractor::Kmer64, uint16_t>>;
template class KmerStorage<KmerExtractor::Kmer128, KmerExtractor, SortedMultiset<KmerExtractor::Kmer128, uint16_t>>;
template class KmerStorage<KmerExtractor::Kmer256, KmerExtractor, SortedMultiset<KmerExtractor::Kmer256, uint16_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer64, uint16_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer128, uint16_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer256, uint16_t>>;

template class KmerStorage<KmerExtractor::Kmer64, KmerExtractor, SortedMultiset<KmerExtractor::Kmer64, uint32_t>>;
template class KmerStorage<KmerExtractor::Kmer128, KmerExtractor, SortedMultiset<KmerExtractor::Kmer128, uint32_t>>;
template class KmerStorage<KmerExtractor::Kmer256, KmerExtractor, SortedMultiset<KmerExtractor::Kmer256, uint32_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer64, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer64, uint32_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer128, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer128, uint32_t>>;
template class KmerStorage<KmerExtractor2Bit::Kmer256, KmerExtractor2Bit, SortedMultiset<KmerExtractor2Bit::Kmer256, uint32_t>>;
