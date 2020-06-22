#include "kmer_collector.hpp"

#include <type_traits>

#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_set_disk.hpp"
#include "common/sorted_multiset_disk.hpp"
#include "common/unix_tools.hpp"
#include "kmer.hpp"
#include "kmer_extractor.hpp"
#include "kmer_to_int_converter.hpp"


namespace mtg {
namespace kmer {

const size_t kLargeBufferSize = 1'000'000;
const size_t kBufferSize = 100'000;


template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::value_type>);
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::key_type>);

    Vector<typename KMER::WordType> buffer;
    buffer.reserve(kBufferSize);

    KmerExtractor kmer_extractor;

    generate_reads([&](const std::string &read) {
        kmer_extractor.sequence_to_kmers(read, k, suffix,
                                         reinterpret_cast<Vector<KMER> *>(&buffer));
        if (both_strands_mode) {
            auto rev_read = read;
            reverse_complement(rev_read.begin(), rev_read.end());
            kmer_extractor.sequence_to_kmers(rev_read, k, suffix,
                                             reinterpret_cast<Vector<KMER> *>(&buffer));
        }

        if (buffer.size() > 0.9 * kBufferSize) {
            kmers->insert(buffer.begin(), buffer.end());

            if (buffer.capacity() > 2 * kBufferSize) {
                buffer = Vector<typename KMER::WordType>(kBufferSize);
            }

            buffer.resize(0);
        }
    });

    if (buffer.size())
        kmers->insert(buffer.begin(), buffer.end());
}

template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(utils::is_instance_v<Container, common::SortedMultiset>
                  || utils::is_instance_v<Container, common::SortedMultisetDisk>);
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::key_type>);

    using KmerCount = typename Container::count_type;

    Vector<KMER> buffer;
    Vector<std::pair<typename KMER::WordType, KmerCount>> buffer_with_counts;
    buffer.reserve(kBufferSize);
    buffer_with_counts.reserve(kBufferSize);

    KmerExtractor kmer_extractor;

    generate_reads([&](const std::string &read, uint64_t count) {
        count = std::min(count, kmers->max_count());

        kmer_extractor.sequence_to_kmers(read, k, suffix, &buffer);
        if (both_strands_mode) {
            auto rev_read = read;
            reverse_complement(rev_read.begin(), rev_read.end());
            kmer_extractor.sequence_to_kmers(rev_read, k, suffix, &buffer);
        }

        for (const KMER &kmer : buffer) {
            buffer_with_counts.emplace_back(kmer.data(), count);
        }

        if (buffer.capacity() > 2 * kBufferSize)
            buffer = Vector<KMER>(kBufferSize);

        buffer.resize(0);

        if (buffer_with_counts.size() > 0.9 * kBufferSize) {
            kmers->insert(buffer_with_counts.begin(),
                          buffer_with_counts.end());

            if (buffer_with_counts.capacity() > 2 * kBufferSize) {
                buffer_with_counts
                    = decltype(buffer_with_counts)(kBufferSize);
            }

            buffer_with_counts.resize(0);
        }
    });

    if (buffer_with_counts.size()) {
        kmers->insert(buffer_with_counts.begin(),
                      buffer_with_counts.end());
    }
}

template <typename KMER, class KmerExtractor, class Container>
KmerCollector<KMER, KmerExtractor, Container>
::KmerCollector(size_t k,
                bool both_strands_mode,
                std::vector<typename Extractor::TAlphabet>&& filter_suffix_encoded,
                size_t num_threads,
                double memory_preallocated,
                const std::filesystem::path &tmp_dir,
                size_t __attribute__((unused)) max_disk_space)
      : k_(k),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_), 1),
        batcher_([this](auto&& sequences) { add_batch(std::move(sequences)); },
                 kLargeBufferSize / sizeof(typename decltype(batcher_)::value_type),
                 kLargeBufferSize),
        filter_suffix_encoded_(std::move(filter_suffix_encoded)),
        both_strands_mode_(both_strands_mode) {
    assert(num_threads_ > 0);

    buffer_size_ = memory_preallocated / sizeof(typename Container::value_type);

    if constexpr(utils::is_instance_v<Data, common::ChunkedWaitQueue>) {
        tmp_dir_ = utils::create_temp_dir(tmp_dir, "kmers");
        kmers_ = std::make_unique<Container>(num_threads, buffer_size_,
                                             tmp_dir_, max_disk_space);
    } else {
        kmers_ = std::make_unique<Container>(num_threads, buffer_size_);
    }
    common::logger->trace(
            "Preallocated {} MiB for the k-mer storage, capacity: {} k-mers",
            kmers_->buffer_size() * sizeof(typename Container::value_type) >> 20,
            kmers_->buffer_size());
}

template <typename KMER, class KmerExtractor, class Container>
KmerCollector<KMER, KmerExtractor, Container>
::~KmerCollector() {
    if (!tmp_dir_.empty())
        std::filesystem::remove_all(tmp_dir_);
}

template <typename KMER, class KmerExtractor, class Container>
void KmerCollector<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallString)> &generate_sequences) {
    if constexpr(std::is_same_v<typename KMER::WordType, typename Container::value_type>) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>,
                             generate_sequences,
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_);
    } else {
        thread_pool_.enqueue(count_kmers<KMER, Extractor, Container>,
                             [generate_sequences](CallStringCount callback) {
                                 generate_sequences([&](const std::string &seq) { callback(seq, 1); });
                             },
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerCollector<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
    if constexpr(std::is_same_v<typename KMER::WordType, typename Container::value_type>) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>,
                             [generate_sequences](CallString callback) {
                                 generate_sequences([&](const std::string &seq, uint64_t) {
                                     callback(seq);
                                 });
                             },
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_);
    } else {
        thread_pool_.enqueue(count_kmers<KMER, Extractor, Container>,
                             generate_sequences,
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_);
    }
}

template <typename KMER, class KmerExtractor, class Container>
void KmerCollector<KMER, KmerExtractor, Container>
::add_batch(std::vector<std::pair<std::string, uint64_t>>&& sequences) {
    // we capture only a pointer to #sequences in the lambda expression to avoid
    // copying the sequences when the lambda is transformed to a callback which is passed
    // to std::bind by const reference in #add_sequences and hence is copied with all
    // its parameters captured.
    auto seqs = std::make_shared<std::vector<std::pair<std::string, uint64_t>>>(std::move(sequences));
    add_sequences([seqs](CallStringCount callback) {
        for (const auto &[seq, count] : *seqs) {
            callback(seq, count);
        }
    });
}

template <typename KMER, class KmerExtractor, class Container>
void KmerCollector<KMER, KmerExtractor, Container>::join() {
    batcher_.process_all_buffered();
    thread_pool_.join();
}

#define INSTANTIATE_KMER_STORAGE(KMER_EXTRACTOR, KMER) \
    template class KmerCollector<KMER, KMER_EXTRACTOR, \
            common::SortedSet<KMER::WordType>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, \
            common::SortedMultiset<KMER::WordType, uint8_t>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, \
            common::SortedMultiset<KMER::WordType, uint16_t>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, \
            common::SortedMultiset<KMER::WordType, uint32_t>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, common::SortedSetDisk<KMER::WordType>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, common::SortedMultisetDisk<KMER::WordType, uint8_t>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, common::SortedMultisetDisk<KMER::WordType, uint16_t>>; \
    template class KmerCollector<KMER, KMER_EXTRACTOR, common::SortedMultisetDisk<KMER::WordType, uint32_t>>;


INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer64)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer128)
INSTANTIATE_KMER_STORAGE(KmerExtractorBOSS, KmerExtractorBOSS::Kmer256)

INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::KmerBOSS64)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::KmerBOSS128)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::KmerBOSS256)

INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer64)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer128)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer256)

} // namespace kmer
} // namespace mtg
