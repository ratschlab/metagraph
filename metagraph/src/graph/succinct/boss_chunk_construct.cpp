#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "boss_chunk.hpp"
#include "kmer_collector.hpp"
#include "unix_tools.hpp"

/**
 * What type of data structure to use in the #KmerCollector.
 */
enum class ExtractorContainer {
  DEQUE,
  VECTOR,
  /**
   * Uses several vectors that are written to disk and then merged, as defined
   * in #SortedSetDisk
   */
  VECTOR_DISK
};
constexpr static ExtractorContainer kExtractorContainer =
    ExtractorContainer::VECTOR;
const static uint8_t kBitsPerCount = 8;


template <class Array>
void sort_and_remove_duplicates(Array *array,
                                size_t num_threads,
                                size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          utils::LessFirst<typename Array::value_type>(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(array->begin() + offset, array->end(),
                                  utils::EqualFirst<typename Array::value_type>());
    array->erase(unique_end, array->end());
}

template <typename Array>
void shrink_kmers(Array *kmers,
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
                  << (kmers->size() * sizeof(typename Array::value_type) >> 20) << "Mb" << std::endl;
    }
}

template <class Container, typename KMER>
inline KMER& push_back(Container &kmers, const KMER &kmer) {
    if constexpr(utils::is_pair<typename Container::value_type>::value) {
        kmers.emplace_back(kmer, 0);
        return kmers.back().first;
    } else {
        kmers.push_back(kmer);
        return kmers.back();
    }
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
// k is node length
template <typename Array>
void recover_source_dummy_nodes(size_t k,
                                Array *kmers,
                                size_t num_threads,
                                bool verbose) {
    using KMER = std::remove_reference_t<decltype(utils::get_first((*kmers)[0]))>;

    size_t dummy_begin = kmers->size();
    size_t num_dummy_parent_kmers = 0;

    for (size_t i = 0; i < dummy_begin; ++i) {
        const KMER &kmer = utils::get_first((*kmers)[i]);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        auto node_last_char = kmer[1];
        auto edge_label = kmer[0];
        // check if it's not a source dummy kmer
        if (node_last_char || !edge_label)
            continue;

        num_dummy_parent_kmers++;

        if (kmers->size() + 1 > kmers->capacity())
            shrink_kmers(kmers, num_threads, verbose, dummy_begin);

        push_back(*kmers, kmer).to_prev(k + 1, BOSS::kSentinelCode);
    }
    if (verbose) {
        std::cout << "Number of dummy k-mers with dummy prefix of length 1: "
                  << num_dummy_parent_kmers << std::endl;
    }
    sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

    if (verbose) {
        std::cout << "Number of dummy k-mers with dummy prefix of length 2: "
                  << kmers->size() - dummy_begin << std::endl;
    }

    for (size_t c = 3; c < k + 1; ++c) {
        size_t succ_dummy_begin = dummy_begin;
        dummy_begin = kmers->size();

        for (size_t i = succ_dummy_begin; i < dummy_begin; ++i) {
            if (kmers->size() + 1 > kmers->capacity())
                shrink_kmers(kmers, num_threads, verbose, dummy_begin);

            push_back(*kmers, utils::get_first((*kmers)[i])).to_prev(k + 1, BOSS::kSentinelCode);
        }
        sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

        if (verbose) {
            std::cout << "Number of dummy k-mers with dummy prefix of length " << c
                      << ": " << kmers->size() - dummy_begin << std::endl;
        }
    }
    ips4o::parallel::sort(kmers->begin(), kmers->end(),
                          utils::LessFirst<typename Array::value_type>(),
                          num_threads);
}

inline std::vector<KmerExtractorBOSS::TAlphabet>
encode_filter_suffix_boss(const std::string &filter_suffix) {
    KmerExtractorBOSS kmer_extractor;
    std::vector<typename KmerExtractorBOSS::TAlphabet> filter_suffix_encoded;
    std::transform(
        filter_suffix.begin(), filter_suffix.end(),
        std::back_inserter(filter_suffix_encoded),
        [&kmer_extractor](char c) {
            return c == BOSS::kSentinel
                            ? BOSS::kSentinelCode
                            : kmer_extractor.encode(c);
        }
    );
    return filter_suffix_encoded;
}

template <typename KmerStorage>
class BOSSChunkConstructor : public IBOSSChunkConstructor {
    friend IBOSSChunkConstructor;

    template <template <typename KMER> class KmerContainer, typename... Args>
    friend std::unique_ptr<IBOSSChunkConstructor>
    initialize_boss_chunk_constructor(size_t k, const Args& ...args);

  private:
    BOSSChunkConstructor(size_t k,
                         bool canonical_mode = false,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0,
                         bool verbose = false)
          : kmer_storage_(k + 1,
                          canonical_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          verbose) {
        if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
            kmer_storage_.insert_dummy(std::vector<KmerExtractorBOSS::TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_sequence(std::string&& sequence, uint64_t count) {
        kmer_storage_.add_sequence(std::move(sequence), count);
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_storage_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk() {
        auto &kmers = kmer_storage_.data();

        if (!kmer_storage_.suffix_length()) {
            if (kmer_storage_.verbose()) {
                std::cout << "Reconstructing all required dummy source k-mers..."
                          << std::endl;
            }
            Timer timer;

            // kmer_collector stores (BOSS::k_ + 1)-mers
            recover_source_dummy_nodes(kmer_storage_.get_k() - 1,
                                       &kmers,
                                       kmer_storage_.num_threads(),
                                       kmer_storage_.verbose());

            if (kmer_storage_.verbose())
                std::cout << "Dummy source k-mers were reconstructed in "
                          << timer.elapsed() << "sec" << std::endl;
        }

        if constexpr(std::is_same_v<typename KmerStorage::Data,
                                    utils::DequeStorage<typename KmerStorage::Value>>) {
            kmers.shrink_to_fit();
        }

        BOSS::Chunk *result;

        if constexpr(utils::is_pair<typename KmerStorage::Value>::value) {
            // kmer_collector stores (BOSS::k_ + 1)-mers
            result = new BOSS::Chunk(kmer_storage_.alphabet_size(),
                                     kmer_storage_.get_k() - 1,
                                     kmer_storage_.is_both_strands_mode(),
                                     kmers,
                                     kBitsPerCount);
        } else {
            // kmer_collector stores (BOSS::k_ + 1)-mers
            result = new BOSS::Chunk(kmer_storage_.alphabet_size(),
                                     kmer_storage_.get_k() - 1,
                                     kmer_storage_.is_both_strands_mode(),
                                     kmers);
        }

        kmer_storage_.clear();

        return result;
    }

    uint64_t get_k() const { return kmer_storage_.get_k() - 1; }

    KmerStorage kmer_storage_;
};

template <template <typename KMER> class KmerContainer, typename... Args>
static std::unique_ptr<IBOSSChunkConstructor>
initialize_boss_chunk_constructor(size_t k, const Args& ...args) {
    if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 64) {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer64>>(k, args...)
        );
    } else if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 128) {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer128>>(k, args...)
        );
    } else {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer256>>(k, args...)
        );
    }
}

template <typename KMER>
using KmerCounterVector = KmerCounter<KMER, KmerExtractorBOSS, uint8_t,
                                      Vector<std::pair<KMER, uint8_t>>,
                                      utils::NoCleanup>;
template <typename KMER>
using KmerCounterVectorClean = KmerCounter<KMER, KmerExtractorBOSS, uint8_t,
                                           Vector<std::pair<KMER, uint8_t>>,
                                           utils::DummyKmersCleaner>;
template <typename KMER>
using KmerCollectorVector = KmerCollector<KMER, KmerExtractorBOSS,
                                          Vector<KMER>,
                                          utils::NoCleanup>;
template <typename KMER>
using KmerCollectorVectorClean = KmerCollector<KMER, KmerExtractorBOSS,
                                               Vector<KMER>,
                                               utils::DummyKmersCleaner>;

template <typename KMER>
using KmerCollectorVectorDisk = KmerCollectorDisk<KMER, KmerExtractorBOSS>;

template <typename KMER>
using KmerCounterDeque = KmerCounter<KMER, KmerExtractorBOSS, uint8_t,
                                     utils::DequeStorage<std::pair<KMER, uint8_t>>,
                                     utils::NoCleanup>;
template <typename KMER>
using KmerCounterDequeClean = KmerCounter<KMER, KmerExtractorBOSS, uint8_t,
                                          utils::DequeStorage<std::pair<KMER, uint8_t>>,
                                          utils::DummyKmersCleaner>;
template <typename KMER>
using KmerCollectorDeque = KmerCollector<KMER, KmerExtractorBOSS,
                                         utils::DequeStorage<KMER>,
                                         utils::NoCleanup>;
template <typename KMER>
using KmerCollectorDequeClean = KmerCollector<KMER, KmerExtractorBOSS,
                                              utils::DequeStorage<KMER>,
                                              utils::DummyKmersCleaner>;

std::unique_ptr<IBOSSChunkConstructor>
IBOSSChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             bool count_kmers,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {

    #define OTHER_ARGS k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose

    if (count_kmers) {
        if (filter_suffix.size()) {
            switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerCounterVector>(OTHER_ARGS);
            case ExtractorContainer::DEQUE:
                return initialize_boss_chunk_constructor<KmerCounterDeque>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unsupported extractor container specified for "
                        "counter. Only VECTOR and DEQUE are supported");
            }
        } else {
            switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerCounterVectorClean>(OTHER_ARGS);
            case ExtractorContainer::DEQUE:
                return initialize_boss_chunk_constructor<KmerCounterDequeClean>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unsupported extractor container specified for "
                        "counter. Only VECTOR and DEQUE are supported");
            }
        }
    } else {
        if (filter_suffix.size()) {
            switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerCollectorVector>(OTHER_ARGS);
            case ExtractorContainer::DEQUE:
                return initialize_boss_chunk_constructor<KmerCollectorDeque>(OTHER_ARGS);
            case ExtractorContainer::VECTOR_DISK:
                return initialize_boss_chunk_constructor<KmerCollectorVectorDisk>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unknown extractor container: " +
                        to_string(static_cast<uint32_t>(kExtractorContainer)));
            }
        } else {
            switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerCollectorVectorClean>(OTHER_ARGS);
            case ExtractorContainer::DEQUE:
                return initialize_boss_chunk_constructor<KmerCollectorDequeClean>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unsupported extractor container specified for "
                        "collector. Only VECTOR and DEQUE are supported");
            }
        }
    }
}

