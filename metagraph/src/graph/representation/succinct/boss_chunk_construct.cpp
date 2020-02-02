#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "common/circular_buffer.hpp"
#include "common/file_merger.hpp"
#include "common/logger.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_multiset_disk.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_set_disk_new.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "kmer/kmer_collector.hpp"
#include "boss_chunk.hpp"

namespace mg {
namespace succinct {

using namespace mg;
using common::logger;


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
                  size_t offset) {
    logger->trace("Allocated capacity exceeded, filter out non-unique k-mers...");

    size_t prev_num_kmers = kmers->size();
    sort_and_remove_duplicates(kmers, num_threads, offset);

    logger->trace("...done. Number of kmers reduced from {} to {}, {}Mb", prev_num_kmers,
                  kmers->size(), (kmers->size() * sizeof(typename Array::value_type) >> 20));
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
 * The assumption is that only dummy sources with sentinels of length 1 were added in
 * the  previous phases and this method will gradually add dummy sources with sentinels
 * of length 2, 3, ... up to k-1.
 * For example, if the input k-mers are {ACG, $TA, $CG}, the first k-mer will be
 * ignored (not a dummy source) and the last 2 k-mers will be expanded to sentinels of
 * length 2, by appending {$$T, $$C}
 * @tparam Container the data structure in which the k-mers were merged (e.g. a
 * ChunkedWaitQueue if using a SortedSetDisk or a Vector if using SortedSet).
 */
template <typename Container>
void recover_source_dummy_nodes(size_t k,
                                Container *kmers,
                                size_t num_threads,
                                ThreadPool & /* async_worker */) {
    using KMER = std::remove_reference_t<decltype(utils::get_first((*kmers)[0]))>;

    size_t dummy_begin = kmers->size();
    size_t num_dummy_parent_kmers = 0;

    for (size_t i = 0; i < dummy_begin; ++i) {
        const KMER &kmer = utils::get_first((*kmers)[i]);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        auto node_last_char = kmer[1];
        auto edge_label = kmer[0];
        // skip if it's not a source dummy kmer
        if (node_last_char || !edge_label)
            continue;

        num_dummy_parent_kmers++;

        if (kmers->size() + 1 > kmers->capacity())
            shrink_kmers(kmers, num_threads, dummy_begin);

        push_back(*kmers, kmer).to_prev(k + 1, BOSS::kSentinelCode);
    }
    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  num_dummy_parent_kmers);
    sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

    logger->trace("Number of dummy k-mers with dummy prefix of length 2: {}",
                  kmers->size() - dummy_begin);

    for (size_t c = 3; c < k + 1; ++c) {
        size_t succ_dummy_begin = dummy_begin;
        dummy_begin = kmers->size();

        for (size_t i = succ_dummy_begin; i < dummy_begin; ++i) {
            if (kmers->size() + 1 > kmers->capacity())
                shrink_kmers(kmers, num_threads, dummy_begin);

            push_back(*kmers, utils::get_first((*kmers)[i])).to_prev(k + 1, BOSS::kSentinelCode);
        }
        sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      kmers->size() - dummy_begin);
    }
    ips4o::parallel::sort(kmers->begin(), kmers->end(),
                          utils::LessFirst(),
                          num_threads);
}

/**
 * An element of #RecentKmers. Contains a k-mer and a flag that indicates if the k-mer
 * was identified as a redundant dummy source k-mer.
 * */
template <typename T>
struct Kmer {
    T kmer; // the actual k-mer or k-mer with count
    bool is_removed; // true if this is a redundant dummy source k-mer
};

template <typename T>
using RecentKmers = common::CircularBuffer<Kmer<T>>;
/**
 * Iterate backwards and check if there is a dummy incoming edge (a dummy edge with
 * identical suffix and label as the current one) into the same node.
 * If such an edge is found, then we remove it (because it's redundant).
 * @param kmer check for redundancy against this k-mer
 * @param buffer queue that stores the last  |Alphabet|^2 k-mers used for identifying
 * redundant dummy source k-mers
 */
template <typename T, typename KMER>
static void remove_redundant_dummy_source(const KMER &kmer, RecentKmers<T> *buffer) {
    using TAlphabet = typename KMER::CharType;
    TAlphabet curW = utils::get_first(kmer)[0];
    if (buffer->size() < 2 || curW == 0) {
        return;
    }
    using ReverseIterator = typename RecentKmers<T>::reverse_iterator;
    ReverseIterator it = buffer->rbegin();
    ++it;
    for (; KMER::compare_suffix(kmer, utils::get_first((*it).kmer), 1); ++it) {
        const KMER prev_kmer = utils::get_first((*it).kmer);
        assert((curW || prev_kmer[1]) && "Main dummy source k-mer must be unique");
        if (prev_kmer[0] == curW && !prev_kmer[1]) { // redundant dummy source k-mer
            (*it).is_removed = true;
            return;
        }
        if (it.at_begin())
            break;
    }
}


template <typename T>
void write_or_die(std::ofstream *f, const T &v) {
    if (!f->write(reinterpret_cast<const char *>(&v), sizeof(T))) {
        std::cerr << "Error: Writing of merged data failed." << std::endl;
        std::exit(EXIT_FAILURE);
    }
};


/**
 * Writes the front of #buffer to a file if it's not a redundant dummy kmer. If the k-mer
 * at the front is a (not-redundant) dummy k-mer it also writes its corresponding dummy
 * k-mer of prefix length 2 into #sorted_dummy_kmers.
 */
template <typename T>
uint8_t write_kmer(size_t k,
                   Vector<T> *dummy_kmers,
                   utils::BufferedAsyncWriter<T> *writer,
                   RecentKmers<T> *buffer,
                   common::SortedSetDisk<T> *sorted_dummy_kmers) {
    const Kmer<T> to_write = buffer->pop_front();
    if (to_write.is_removed) { // redundant dummy k-mer
        return 0;
    }
    writer->push(to_write.kmer);
    using KMER = std::remove_reference_t<decltype(utils::get_first(to_write.kmer))>;
    using TAlphabet = typename KMER::CharType;

    const KMER &kmer_to_write = utils::get_first(to_write.kmer);
    const TAlphabet node_last_char = kmer_to_write[1];
    const TAlphabet edge_label = kmer_to_write[0];
    if (node_last_char || !edge_label) { // not a dummy source kmer
        return 0;
    }

    if (dummy_kmers->size() == dummy_kmers->capacity()) {
        sorted_dummy_kmers->insert(dummy_kmers->begin(), dummy_kmers->end());
        dummy_kmers->resize(0);
    }

    push_back(*dummy_kmers, kmer_to_write).to_prev(k + 1, BOSS::kSentinelCode);
    return 1;
}

/**
 * Specialization of recover_dummy_nodes for a #common::ChunkedWaitQueue container
 * (used by #common::SortedSetDisk).
 * The method first removes redundant dummy source k-mers of prefix length 1, then
 * gradually constructs dummy source k-mers of prefix length 2..k and writes them into
 * separate files, de-duped and sorted. The final result is obtained by merging the
 * original #kmers (minus the redundant dummy source k-mers of prefix length 1) with  the
 * dummy source k-mers for prefix length 2..k which were generated and saved into files.
 */
template <typename T>
void recover_source_dummy_nodes(size_t k,
                                common::ChunkedWaitQueue<T> *kmers,
                                size_t num_threads,
                                ThreadPool &async_worker) {
    using KMER = std::remove_reference_t<decltype(utils::get_first(*(kmers->begin())))>;

    // name of the file containing dummy k-mers of given prefix length
    //TODO(ddanciu): make sure these path names are unique and can be changed
    const auto get_file_name
            = [](uint32_t pref_len) { return "/tmp/dummy" + std::to_string(pref_len); };

    const auto no_cleanup = [](typename common::SortedSetDisk<T>::storage_type *) {};

    const auto file_writer
            = [](std::ofstream &f) { return [&f](const T &v) { write_or_die<T>(&f, v); }; };

    std::vector<std::pair<std::string, std::ofstream>> files_to_merge;
    auto create_stream = [](const std::string &filename) {
        std::ofstream f(filename, std::ios::out | std::ios::binary);
        return std::make_pair(filename, std::move(f));
    };

    const std::string file_name = "/tmp/original_and_dummy_l1";
    files_to_merge.push_back(create_stream(file_name));
    files_to_merge.reserve(k + 1); // avoid re-allocations as we keep refs to elements
    std::ofstream *dummy_l1 = &files_to_merge.back().second;

    RecentKmers<T> recent_buffer((1llu << KMER::kBitsPerChar)
                                  * (1llu << KMER::kBitsPerChar));

    const std::string file_name_l2 = get_file_name(2);
    files_to_merge.push_back(create_stream(file_name_l2));
    std::ofstream *dummy_l2 = &files_to_merge.back().second;

    // this will contain dummy k-mers of prefix length 2
    common::SortedSetDisk<T> sorted_dummy_kmers(no_cleanup, num_threads,
                                                kmers->buffer_size(), "/tmp/chunk_",
                                                file_writer(*dummy_l2));
    Vector<T> dummy_kmers;
    dummy_kmers.reserve(sorted_dummy_kmers.buffer_size());

    // remove redundant dummy source k-mers of prefix length 1 and write them to a file
    // While traversing and removing redundant dummy source k-mers of prefix length 1,
    // we also  generate dummy k-mers of prefix length 2.
    size_t num_dummy_parent_kmers = 0;
    size_t num_parent_kmers = 0;
    // asynchronously writes a value of type T to a file stream
    utils::BufferedAsyncWriter<T> writer(file_name, dummy_l1);
    for (auto &it = kmers->begin(); it != kmers->end(); ++it) {
        num_parent_kmers++;
        const T el = *it;
        recent_buffer.push_back({ el, false });
        remove_redundant_dummy_source<T, KMER>(utils::get_first(el), &recent_buffer);
        if (recent_buffer.full()) {
            num_dummy_parent_kmers += write_kmer(k, &dummy_kmers, &writer, &recent_buffer,
                                                 &sorted_dummy_kmers);
        }
    }
    while (!recent_buffer.empty()) { // empty the buffer
        num_dummy_parent_kmers += write_kmer(k, &dummy_kmers, &writer, &recent_buffer,
                                             &sorted_dummy_kmers);
    }
    writer.flush();
    // push out the leftover dummy kmers
    sorted_dummy_kmers.insert(dummy_kmers.begin(), dummy_kmers.end());

    logger->trace("Total number of k-mers: {}", num_parent_kmers);
    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  num_dummy_parent_kmers);

    // generate dummy k-mers of prefix length 3..k
    common::SortedSetDisk<T> sorted_dummy_kmers2(no_cleanup, num_threads,
                                                 kmers->buffer_size(), "/tmp/chunk2_");
    common::SortedSetDisk<T> *source = &sorted_dummy_kmers;
    common::SortedSetDisk<T> *dest = &sorted_dummy_kmers2;
    for (size_t dummy_pref_len = 3; dummy_pref_len < k + 1; ++dummy_pref_len) {
        files_to_merge.push_back(create_stream(get_file_name(dummy_pref_len)));
        dest->clear(file_writer(files_to_merge.back().second));
        dummy_kmers.resize(0);
        size_t num_kmers = 0;
        for (auto &it = source->data().begin(); it != source->data().end(); ++it) {
            if (dummy_kmers.size() == dummy_kmers.capacity()) {
                dest->insert(dummy_kmers.begin(), dummy_kmers.end());
                dummy_kmers.resize(0);
            }

            push_back(dummy_kmers, utils::get_first(*it)).to_prev(k + 1, BOSS::kSentinelCode);
            num_kmers++;
        }
        // push out the leftover dummy kmers
        dest->insert(dummy_kmers.begin(), dummy_kmers.end());

        logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}",
                      dummy_pref_len - 1, num_kmers);

        std::swap(source, dest);
    }
    uint32_t num_kmers = 0;
    // iterate to merge the data and write it to disk
    for (auto &it = source->data().begin(); it != source->data().end(); ++it, ++num_kmers) {
    }
    logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}", k, num_kmers);

    // at this point, we have the original k-mers plus the  dummy k-mers with prefix
    // length x in /tmp/dummy_{x}, and we'll merge them all into a single stream
    std::vector<std::string> file_names;
    std::for_each(files_to_merge.begin(), files_to_merge.end(), [&file_names](auto &el) {
        el.second.flush();
        file_names.push_back(el.first);
    });

    kmers->reset();
    async_worker.enqueue([=]() {
        common::merge_files<T>(file_names, [&](const T &v) { kmers->push(v); });
        kmers->shutdown();
    });
}

inline std::vector<KmerExtractorBOSS::TAlphabet>
encode_filter_suffix_boss(const std::string &filter_suffix) {
    KmerExtractorBOSS kmer_extractor;
    std::vector<typename KmerExtractorBOSS::TAlphabet> filter_suffix_encoded;
    for (char c : filter_suffix) {
        filter_suffix_encoded.push_back(c == BOSS::kSentinel
                                        ? BOSS::kSentinelCode
                                        : kmer_extractor.encode(c));
    }
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
                         uint8_t bits_per_count = 0,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0)
          : kmer_collector_(k + 1,
                            canonical_mode,
                            encode_filter_suffix_boss(filter_suffix),
                            num_threads,
                            memory_preallocated),
            bits_per_count_(bits_per_count) {
        if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
            kmer_collector_.add_kmer(std::vector<KmerExtractorBOSS::TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_sequence(std::string_view sequence, uint64_t count) {
        kmer_collector_.add_sequence(sequence, count);
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk() {
        auto &kmers = kmer_collector_.data();

        if (!kmer_collector_.suffix_length()) {
            logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;

            // kmer_collector stores (BOSS::k_ + 1)-mers
            recover_source_dummy_nodes(kmer_collector_.get_k() - 1,
                                       &kmers,
                                       kmer_collector_.num_threads(),
                                       async_worker_);

            logger->trace("Dummy source k-mers were reconstructed in {} sec",
                          timer.elapsed());
        }

        BOSS::Chunk *result;

        if constexpr(utils::is_pair<typename KmerCollector::Value>::value) {
            // kmer_collector stores (BOSS::k_ + 1)-mers
            result = new BOSS::Chunk(kmer_collector_.alphabet_size(),
                                     kmer_collector_.get_k() - 1,
                                     kmer_collector_.is_both_strands_mode(),
                                     kmers,
                                     bits_per_count_);
        } else {
            // kmer_collector stores (BOSS::k_ + 1)-mers
            result = new BOSS::Chunk(kmer_collector_.alphabet_size(),
                                     kmer_collector_.get_k() - 1,
                                     kmer_collector_.is_both_strands_mode(),
                                     kmers);
        }

        kmer_collector_.clear();

        return result;
    }

    uint64_t get_k() const { return kmer_collector_.get_k() - 1; }

    KmerCollector kmer_collector_;
    uint8_t bits_per_count_;
    /** Used as an async executor for merging chunks from disk */
    ThreadPool async_worker_ = ThreadPool(1, 1);
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
using KmerSetVector
        = kmer::KmerCollector<KMER, KmerExtractorBOSS,
                              common::SortedSet<KMER, Vector<KMER>>>;

template <typename KMER>
using KmerMultsetVector8
        = kmer::KmerCollector<KMER, KmerExtractorBOSS,
                              common::SortedMultiset<KMER, uint8_t,
                                                     Vector<std::pair<KMER, uint8_t>>>>;

template <typename KMER>
using KmerMultsetVector16
        = kmer::KmerCollector<KMER, KmerExtractorBOSS,
                              common::SortedMultiset<KMER, uint16_t,
                                                     Vector<std::pair<KMER, uint16_t>>>>;

template <typename KMER>
using KmerMultsetVector32
        = kmer::KmerCollector<KMER, KmerExtractorBOSS,
                              common::SortedMultiset<KMER, uint32_t,
                                                     Vector<std::pair<KMER, uint32_t>>>>;

template <typename KMER>
using KmerSetDisk
        = kmer::KmerCollector<KMER, KmerExtractorBOSS,
                              common::SortedSetDisk<KMER>>;

template <typename KMER>
using KmerMultsetDiskVector8
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedMultisetDisk<KMER, uint8_t>>;

template <typename KMER>
using KmerMultsetDiskVector16
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedMultisetDisk<KMER, uint16_t>>;

template <typename KMER>
using KmerMultsetDiskVector32
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedMultisetDisk<KMER, uint32_t>>;

std::unique_ptr<IBOSSChunkConstructor>
IBOSSChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             uint8_t bits_per_count,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             kmer::ContainerType container_type) {
#define OTHER_ARGS k, canonical_mode, bits_per_count, filter_suffix, num_threads, memory_preallocated
    switch (container_type) {
        case kmer::ContainerType::VECTOR:
            if (!bits_per_count) {
                return initialize_boss_chunk_constructor<KmerSetVector>(OTHER_ARGS);
            } else if (bits_per_count <= 8) {
                return initialize_boss_chunk_constructor<KmerMultsetVector8>(OTHER_ARGS);
            } else if (bits_per_count <= 16) {
                return initialize_boss_chunk_constructor<KmerMultsetVector16>(OTHER_ARGS);
            } else if (bits_per_count <= 32) {
                return initialize_boss_chunk_constructor<KmerMultsetVector32>(OTHER_ARGS);
            } else {
                throw std::runtime_error(
                        "Error: trying to allocate too many bits per k-mer count");
            }
        case kmer::ContainerType::VECTOR_DISK:
            if (!bits_per_count) {
                return initialize_boss_chunk_constructor<KmerSetDisk>(OTHER_ARGS);
            } else if (bits_per_count <= 8) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector8>(OTHER_ARGS);
            } else if (bits_per_count <= 16) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector16>(OTHER_ARGS);
            } else if (bits_per_count <= 32) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector32>(OTHER_ARGS);
            } else {
                throw std::runtime_error(
                        "Error: trying to allocate too many bits per k-mer count");
            }
        default:
            logger->error("Invalid container type {}", container_type);
            std::exit(1);
    }
}

} // namespace succinct
} // namespace mg
