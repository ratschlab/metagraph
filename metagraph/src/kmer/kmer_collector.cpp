#include "kmer_collector.hpp"

#include <type_traits>

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

namespace mg {
namespace kmer {

using namespace mg;

const size_t kMaxKmersChunkSize = 30'000'000;


template <typename KMER, class KmerExtractor, class Container>
void extract_kmers(std::function<void(CallString)> generate_reads,
                   size_t k,
                   bool both_strands_mode,
                   Container *kmers,
                   const std::vector<typename KmerExtractor::TAlphabet> &suffix,
                   bool remove_redundant) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::value_type>);
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::key_type>);

    Vector<KMER> temp_storage;
    Vector<typename KMER::WordType> *temp_storage_int
            = reinterpret_cast<Vector<typename KMER::WordType> *>(&temp_storage);
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
            kmers->sort_and_remove_duplicates(temp_storage_int, 1);
        }

        if (temp_storage.size() > 0.9 * kMaxKmersChunkSize) {
            kmers->insert(temp_storage_int->begin(), temp_storage_int->end());

            if (temp_storage.capacity() > 2 * kMaxKmersChunkSize)
                temp_storage = Vector<KMER>(1.1 * kMaxKmersChunkSize);

            temp_storage.resize(0);
        }
    });

    if (temp_storage.size()) {
        if (remove_redundant) {
            kmers->sort_and_remove_duplicates(temp_storage_int, 1);
        }
        kmers->insert(temp_storage_int->begin(), temp_storage_int->end());
    }
}

template <typename KMER, class KmerExtractor, class Container>
void count_kmers(std::function<void(CallStringCount)> generate_reads,
                 size_t k,
                 bool both_strands_mode,
                 Container *kmers,
                 const std::vector<typename KmerExtractor::TAlphabet> &suffix) {
    static_assert(KMER::kBitsPerChar == KmerExtractor::bits_per_char);
    static_assert(utils::is_instance<Container, common::SortedMultiset> {}
                  || utils::is_instance<Container, common::SortedMultisetDisk> {});
    static_assert(std::is_same_v<typename KMER::WordType, typename Container::key_type>);

    using KmerCount = typename Container::count_type;

    Vector<KMER> temp_storage;
    Vector<std::pair<KMER, KmerCount>> temp_storage_with_counts;
    using int_pair = std::pair<typename KMER::WordType, KmerCount>;
    Vector<int_pair> *temp_storage_with_counts_int
            = reinterpret_cast<Vector<int_pair> *>(&temp_storage);
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
            kmers->insert(temp_storage_with_counts_int->begin(),
                          temp_storage_with_counts_int->end());

            if (temp_storage_with_counts.capacity() > 2 * kMaxKmersChunkSize)
                temp_storage_with_counts
                    = Vector<std::pair<KMER, KmerCount>>(1.1 * kMaxKmersChunkSize);

            temp_storage_with_counts.resize(0);
        }
    });

    if (temp_storage_with_counts.size()) {
        kmers->insert(temp_storage_with_counts_int->begin(), temp_storage_with_counts_int->end());
    }
}

// removes redundant dummy BOSS k-mers from a sorted list
template <class KMER, class KMER_INT>
void cleanup_boss_kmers(Vector<KMER_INT> *kmers_int) {
    Vector<KMER> *kmers = reinterpret_cast<Vector<KMER> *> (kmers_int);

    assert(std::is_sorted(kmers->begin(), kmers->end(), utils::LessFirst()));
    assert(std::unique(kmers->begin(), kmers->end(), utils::EqualFirst()) == kmers->end());

    if (kmers->size() < 2)
        return;

    // last k-mer is never redundant. Start with the next one.
    uint64_t last = kmers->size() - 1;

    typename KMER::CharType edge_label, node_last_char;

    std::vector<uint64_t> last_kmer(1llu << KMER::kBitsPerChar, kmers->size());

    last_kmer[utils::get_first(kmers->at(last))[0]] = last;

    for (int64_t i = last - 1; i >= 0; --i) {
        const KMER &kmer = utils::get_first(kmers->at(i));
        node_last_char = kmer[1];
        edge_label = kmer[0];

        // assert((edge_label || node_last_char)
        //             && "dummy k-mer cannot be both source and sink dummy");

        if (!edge_label) {
            // sink dummy k-mer

            // skip if redundant
            if (node_last_char && KMER::compare_suffix(kmer, utils::get_first(kmers->at(last)), 0))
                continue;

        } else if (!node_last_char) {
            // source dummy k-mer

            // skip if redundant
            if (last_kmer[edge_label] < kmers->size()
                    && KMER::compare_suffix(kmer, utils::get_first(kmers->at(last_kmer[edge_label])), 1))
                continue;
        }

        // the k-mer is either not dummy, or not redundant -> keep the k-mer
        kmers->at(--last) = kmers->at(i);
        last_kmer[edge_label] = last;
    }

    kmers->erase(kmers->begin(), kmers->begin() + last);
}

template <class KmerExtractor, class KMER, class KMER_INT>
std::function<void(Vector<KMER_INT> *)> get_cleanup(bool clean_dummy_boss_kmers) {
    if constexpr(std::is_same_v<KmerExtractor, KmerExtractorBOSS>) {
        if (clean_dummy_boss_kmers) {
            return cleanup_boss_kmers<KMER, KMER_INT>;
        } else {
            return [](Vector<KMER_INT> *) {};
        }
    } else {
        std::ignore = clean_dummy_boss_kmers;
        return [](Vector<KMER_INT> *) {};
    }
}

template <typename KMER, class KmerExtractor, class Container>
KmerCollector<KMER, KmerExtractor, Container>
::KmerCollector(size_t k,
                bool both_strands_mode,
                Sequence&& filter_suffix_encoded,
                size_t num_threads,
                double memory_preallocated,
                const std::filesystem::path &tmp_dir,
                size_t __attribute__((unused)) max_disk_space)
      : k_(k),
        num_threads_(num_threads),
        thread_pool_(std::max(static_cast<size_t>(1), num_threads_) - 1,
                     std::max(static_cast<size_t>(1), num_threads_)),
        batch_accumulator_([this](auto&& sequences) { add_batch(std::move(sequences)); },
                           kMaxKmersChunkSize, kMaxKmersChunkSize, kMaxKmersChunkSize),
        filter_suffix_encoded_(std::move(filter_suffix_encoded)),
        both_strands_mode_(both_strands_mode),
        tmp_dir_(tmp_dir) {
    assert(num_threads_ > 0);
    using KMER_INT = typename Container::value_type;
    std::function<void(Vector<KMER_INT> *)> cleanup = get_cleanup<Extractor, KMER, KMER_INT>(
        filter_suffix_encoded_.empty()
    );
    buffer_size_ = memory_preallocated / sizeof(typename Container::value_type);
    if constexpr ((utils::is_instance<typename KmerCollector::Data, common::ChunkedWaitQueue> {})) {
        if (!filter_suffix_encoded_.empty()) {
            common::logger->error("Disk based sorting does not support chunking");
            exit(1);
        }
        kmers_ = std::make_unique<Container>(cleanup, num_threads, buffer_size_, tmp_dir,
                                             max_disk_space);
    } else {
        kmers_ = std::make_unique<Container>(cleanup, num_threads, buffer_size_);
    }
    common::logger->trace(
            "Preallocated {} MiB for the k-mer storage, capacity: {} k-mers",
            kmers_->buffer_size() * sizeof(typename Container::value_type) >> 20,
            kmers_->buffer_size());
}

template <typename KMER, class KmerExtractor, class Container>
size_t KmerCollector<KMER, KmerExtractor, Container>::buffer_size() const {
    return buffer_size_;
}

template <typename KMER, class KmerExtractor, class Container>
void KmerCollector<KMER, KmerExtractor, Container>
::add_sequences(const std::function<void(CallString)> &generate_sequences) {
    if constexpr(std::is_same_v<typename KMER::WordType, typename Container::value_type>) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>,
                             generate_sequences,
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_, true);
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
    if constexpr (std::is_same_v<typename KMER::WordType, typename Container::value_type>) {
        thread_pool_.enqueue(extract_kmers<KMER, Extractor, Container>,
                             [generate_sequences](CallString callback) {
                                 generate_sequences([&](const std::string &seq, uint64_t) {
                                     callback(seq);
                                 });
                             },
                             k_, both_strands_mode_, kmers_.get(),
                             filter_suffix_encoded_, true);
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
    batch_accumulator_.process_all_buffered();
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

INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer64)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer128)
INSTANTIATE_KMER_STORAGE(KmerExtractor2Bit, KmerExtractor2Bit::Kmer256)

} // namespace kmer
} // namespace mg
