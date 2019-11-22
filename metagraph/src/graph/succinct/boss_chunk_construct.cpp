#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "utils/template_utils.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_set_disk.hpp"
#include "common/unix_tools.hpp"
#include "kmer/kmer_collector.hpp"
#include "boss_chunk.hpp"

namespace mg {
namespace succinct {

template <typename KMER>
using KmerMultsetVector = kmer::KmerCollector<
        KMER,
        KmerExtractorBOSS,
        common::SortedMultiset<KMER, uint8_t, Vector<std::pair<KMER, uint8_t>>>>;

template <typename KMER>
using KmerSetVector
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedSet<KMER, Vector<KMER>>>;

template <typename KMER>
using KmerSetDisk
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedSetDisk<KMER>>;

template <typename KMER,
          typename KmerCount = uint8_t,
          class Container = Vector<std::pair<KMER, KmerCount>>>
using KmerMultset
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedMultiset<KMER, KmerCount, Container>>;

/**
 * What type of data structure to use in the #KmerSet.
 */
enum class ExtractorContainer {
    VECTOR,
    /**
     * Uses several vectors that are written to disk and then merged, as defined
     * in #SortedSetDisk
     */
    VECTOR_DISK
};

constexpr static ExtractorContainer kExtractorContainer = ExtractorContainer::VECTOR;

const static uint8_t kBitsPerCount = 8;


template <class Array>
void sort_and_remove_duplicates(Array *array,
                                size_t num_threads,
                                size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          utils::LessFirst(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(array->begin() + offset, array->end(),
                                  utils::EqualFirst());
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
/**
 * Adds dummy source nodes for the given kmers.
 * The assumption is that only dummy sources with sentinels of of length 1 were added in
 * the  previous phases and this method will gradually add dummy sources with sentinels
 * of length 2,3,... up to k-1.
 * For example, if the input k-mers are {ACG, $TA, $CG}, the first k-mer will be
 * ignored (not a dummy source) and the last 2 k-mers will be expanded to sentinels of
 * length 2, by appending {$$T, $$C}
 * @tparam Container the data structure in which the k-mers were merged (e.g. a
 * ChunkedWaitQueue if using a SortedSetDisk or a Vector if using SortedSet).
 */
template < class Container>
static void recover_source_dummy_nodes(size_t k,
                                Container *kmers,
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
        // nothing to do if it's not a source dummy kmer
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
                          utils::LessFirst(),
                          num_threads);
}

template < typename U>
void recover_source_dummy_nodes(size_t k,
                                common::ChunkedWaitQueue<U> *kmers,
                                size_t num_threads,
                                bool verbose) {
    using Iterator = typename common::ChunkedWaitQueue<U>::iterator;
  Iterator iter =  kmers->begin();
    U first = *iter;

    using KMER = std::remove_reference_t<decltype(utils::get_first(first))>;

    size_t dummy_begin = kmers->size();
    size_t num_dummy_parent_kmers = 0;

    for (size_t i = 0; i < dummy_begin; ++i) {
        typename Array::value_type current = *iter;
        const KMER &kmer = utils::get_first(current);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        auto node_last_char = kmer[1];
        auto edge_label = kmer[0];
        // nothing to do if it's not a source dummy kmer
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
                          utils::LessFirst(),
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

template <typename KmerCollector>
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

        BOSS::Chunk *result;

        if constexpr(utils::is_pair<typename KmerCollector::Value>::value) {
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

    KmerCollector kmer_storage_;
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
        switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerMultsetVector>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unsupported extractor container specified for "
                        "counter. Only VECTOR is supported");
        }
    } else {
        switch (kExtractorContainer) {
            case ExtractorContainer::VECTOR:
                return initialize_boss_chunk_constructor<KmerSetVector>(OTHER_ARGS);
            case ExtractorContainer::VECTOR_DISK:
                return initialize_boss_chunk_constructor<KmerSetDisk>(OTHER_ARGS);
            default:
                throw std::logic_error(
                        "Unknown extractor container: " +
                        to_string(static_cast<uint32_t>(kExtractorContainer)));
        }
    }
}


} // namespace succinct
} // namespace mg
