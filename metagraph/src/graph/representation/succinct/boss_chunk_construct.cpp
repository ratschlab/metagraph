#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "common/circular_buffer.hpp"
#include "common/elias_fano_file_merger.hpp"
#include "common/logger.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_multiset_disk.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_set_disk.hpp"
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
template <typename Data>
void recover_source_dummy_nodes(size_t k, size_t num_threads, Data *kmers) {
    using KMER = std::decay_t<decltype(utils::get_first((*kmers)[0]))>;

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
template <typename T, typename INT>
uint8_t write_kmer(size_t k,
                   std::function<INT(const T &v)> to_int,
                   Vector<T> *dummy_kmers,
                   common::EliasFanoEncoderBuffered<INT> *encoder,
                   RecentKmers<T> *buffer,
                   common::SortedSetDisk<T, INT> *sorted_dummy_kmers) {
    const Kmer<T> to_write = buffer->pop_front();
    if (to_write.is_removed) { // redundant dummy k-mer
        return 0;
    }
    encoder->add(to_int(to_write.kmer));
    using KMER = std::decay_t<decltype(utils::get_first(to_write.kmer))>;
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

// converters from T to its integer representation
template <typename T>
static typename T::WordType to_int(T v) {
    return v.data();
}

template <typename T, typename C>
static std::pair<typename T::WordType, C> to_int(std::pair<T, C> v) {
    return { v.first.data(), v.second };
}

/** Returns a lambda that writes compressed integers into encoder */
template <typename T, typename int_type>
std::function<void(const T &v)>
compressed_writer(common::EliasFanoEncoderBuffered<int_type> *encoder,
                  const std::function<int_type(const T &v)> &to_int) {
    return [encoder, &to_int](const T &v) { encoder->add(to_int(v)); };
};

/**
 * SFINAE structs to construct the integer type corresponding to T. Simply speaking, the
 * #type member of the structs will be T::Wortdype if T is not a pair, and
 * std::pair<T::WordType, C> if T is a pair.
 */
template <typename T, typename = void>
struct get_int_type {
    using type = typename T::WordType;
};

/**
 * Specializes get_int_type for an std::pair<T,C>. The #type member will be set to
 * std::pair<T::WordType, C>.
 */
template <typename T>
struct get_int_type<T, void_t<typename T::second_type>> {
    using type = std::pair<typename T::first_type::WordType, typename T::second_type>;
};

/**
 * Specialization of recover_dummy_nodes for a disk-based container, such as
 * #SortedSetDisk and #SortedMultisetDisk.
 * The method first removes redundant dummy source k-mers of prefix length 1, then
 * gradually constructs dummy source k-mers of prefix length 2..k and writes them into
 * separate files, de-duped and sorted. The final result is obtained by merging the
 * original #kmers (minus the redundant dummy source k-mers of prefix length 1) with  the
 * dummy source k-mers for prefix length 2..k which were generated and saved into files.
 */
template <typename KmerCollector>
void recover_source_dummy_nodes_disk(const KmerCollector &kmer_collector,
                                     typename KmerCollector::Data *kmers,
                                     ThreadPool &async_worker) {
    constexpr size_t CHUNK_QUEUE_BUFFER_SIZE = 10000;
    constexpr size_t CHUNK_QUEUE_FENCE_SIZE = 1; // no of elements to traverse backwards
    constexpr size_t ENCODER_BUFFER_SIZE = 100'000;

    std::filesystem::path tmp_dir = kmer_collector.tmp_dir();
    using T = typename KmerCollector::Value;
    using T_INT = typename get_int_type<T>::type;
    using KMER = typename utils::get_first_type<T>::type;

    // name of the file containing dummy k-mers of given prefix length
    const auto get_file_name = [&tmp_dir](uint32_t pref_len) {
        return tmp_dir / ("dummy_l" + std::to_string(pref_len));
    };

    const auto no_cleanup = [](typename common::SortedSetDisk<T>::storage_type *) {};

    std::vector<std::string> files_to_merge;

    const filesystem::path dummy_l1_path = tmp_dir / "dummy_source1";
    std::filesystem::create_directory(dummy_l1_path);
    const std::string orig_dummy_l1 = dummy_l1_path / "original_and_dummy_l1";
    files_to_merge.push_back(orig_dummy_l1);
    size_t k = kmer_collector.get_k() - 1;
    files_to_merge.reserve(k + 1);

    RecentKmers<T> recent_buffer((1llu << KMER::kBitsPerChar) * (1llu << KMER::kBitsPerChar));

    files_to_merge.push_back(get_file_name(2));

    const filesystem::path dummy_l2_path = tmp_dir / "dummy_source2";

    std::function<T_INT(const T &v)> to_intf = [](const T &v) { return to_int(v); };
    // this will contain dummy k-mers of prefix length 2
    common::EliasFanoEncoderBuffered<T_INT> dummy_l2(files_to_merge.back(), ENCODER_BUFFER_SIZE);
    common::SortedSetDisk<T, T_INT> sorted_dummy_kmers(
            no_cleanup, kmer_collector.num_threads(), kmer_collector.buffer_size(),
            dummy_l2_path, kmer_collector.max_disk_space(), [](const T &) {}, 100, to_intf);
    Vector<T> dummy_kmers;
    dummy_kmers.reserve(sorted_dummy_kmers.buffer_size());

    // traverse the input kmers and remove redundant dummy source k-mers of prefix length
    // 1. While traversing we also  generate dummy k-mers of prefix length 2.
    size_t num_dummy_parent_kmers = 0;
    size_t num_parent_kmers = 0;
    // contains original kmers and non-redundant source dummy k-mers with prefix length 1
    common::EliasFanoEncoderBuffered<T_INT> original_and_l1(orig_dummy_l1, ENCODER_BUFFER_SIZE);
    for (auto &it = kmers->begin(); it != kmers->end(); ++it) {
        num_parent_kmers++;
        const T el = *it;
        recent_buffer.push_back({ el, false });
        remove_redundant_dummy_source<T, KMER>(utils::get_first(el), &recent_buffer);
        if (recent_buffer.full()) {
            num_dummy_parent_kmers += write_kmer(k, to_intf, &dummy_kmers, &original_and_l1,
                                                 &recent_buffer, &sorted_dummy_kmers);
        }
    }
    while (!recent_buffer.empty()) { // empty the buffer
        num_dummy_parent_kmers += write_kmer(k, to_intf, &dummy_kmers, &original_and_l1,
                                             &recent_buffer, &sorted_dummy_kmers);
    }
    original_and_l1.finish();
    // push out the leftover dummy kmers
    sorted_dummy_kmers.insert(dummy_kmers.begin(), dummy_kmers.end());

    logger->trace("Total number of k-mers: {}", num_parent_kmers);
    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  num_dummy_parent_kmers);
    ThreadPool async_merge = ThreadPool(1, 1);
    // generate dummy k-mers of prefix length 3..k
    for (size_t dummy_pref_len = 3; dummy_pref_len < k + 1; ++dummy_pref_len) {
        const filesystem::path tmp_path
                = tmp_dir / ("dummy_source" + std::to_string(dummy_pref_len));
        const std::vector<string> chunk_files = sorted_dummy_kmers.files_to_merge();
        common::EliasFanoEncoderBuffered<T_INT> encoder(files_to_merge.back(),
                                                        ENCODER_BUFFER_SIZE);
        common::ChunkedWaitQueue<T> source(CHUNK_QUEUE_BUFFER_SIZE, CHUNK_QUEUE_FENCE_SIZE,
                                           compressed_writer(&encoder, to_intf));
        async_merge.enqueue([&chunk_files, &source]() {
            std::function<void(const T &)> on_new_item
                    = [&source](const T &v) { source.push(v); };
            common::merge_files<T, T_INT>(chunk_files, on_new_item);
            source.shutdown();
        });
        sorted_dummy_kmers.clear(tmp_path);
        dummy_kmers.resize(0);
        size_t num_kmers = 0;
        for (auto &it = source.begin(); it != source.end(); ++it) {
            if (dummy_kmers.size() == dummy_kmers.capacity()) {
                sorted_dummy_kmers.insert(dummy_kmers.begin(), dummy_kmers.end());
                dummy_kmers.resize(0);
            }

            push_back(dummy_kmers, utils::get_first(*it)).to_prev(k + 1, BOSS::kSentinelCode);
            num_kmers++;
        }
        // push out the leftover dummy kmers
        sorted_dummy_kmers.insert(dummy_kmers.begin(), dummy_kmers.end());
        encoder.finish();
        logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}",
                      dummy_pref_len - 1, num_kmers);
        files_to_merge.push_back(get_file_name(dummy_pref_len));
    }
    uint32_t num_kmers = 0;
    // iterate to merge the data and write it to disk
    const std::vector<string> chunk_files = sorted_dummy_kmers.files_to_merge();
    common::EliasFanoEncoderBuffered<T_INT> encoder(files_to_merge.back(),
                                                    ENCODER_BUFFER_SIZE);
    common::ChunkedWaitQueue<T> source(CHUNK_QUEUE_BUFFER_SIZE, CHUNK_QUEUE_FENCE_SIZE,
                                       compressed_writer(&encoder, to_intf));
    async_merge.enqueue([&chunk_files, &source]() {
        std::function<void(const T &)> on_new_item
                = [&source](const T &v) { source.push(v); };
        common::merge_files<T, T_INT>(chunk_files, on_new_item);
        source.shutdown();
    });
    for (auto &it = source.begin(); it != source.end(); ++it, ++num_kmers) {
    }
    encoder.finish();
    logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}", k, num_kmers);

    // at this point, we have the original k-mers plus the  dummy k-mers with prefix
    // length x in /tmp/dummy_{x}, and we'll merge them all into a single stream
    kmers->reset();
    async_worker.enqueue([kmers, files_to_merge, k, &tmp_dir]() {
        std::function<void(const T &)> on_new_item
                = [kmers](const T &v) { kmers->push(v); };
        common::merge_files<T, T_INT>(files_to_merge, on_new_item);
        kmers->shutdown();
        for (u_int32_t i = 1; i < k + 1; ++i) {
            std::filesystem::remove_all(tmp_dir / ("dummy_source" + std::to_string(i)));
        }
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
                         double memory_preallocated = 0,
                         const std::filesystem::path &tmp_dir = "/tmp",
                         size_t max_disk_space = 1e9)
        : kmer_collector_(k + 1,
                          canonical_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          tmp_dir,
                          max_disk_space),
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
        typename KmerCollector::Data &kmers = kmer_collector_.data();

        if (!kmer_collector_.suffix_length()) {
            logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;
            if constexpr ((utils::is_instance<typename KmerCollector::Data,
                                              common::ChunkedWaitQueue> {})) {
                recover_source_dummy_nodes_disk(kmer_collector_, &kmers, async_worker_);
            } else {
                // kmer_collector stores (BOSS::k_ + 1)-mers
                recover_source_dummy_nodes(kmer_collector_.get_k() - 1,
                                           kmer_collector_.num_threads(), &kmers);
            }
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
using KmerMultsetVector16 = kmer::KmerCollector<
        KMER,
        KmerExtractorBOSS,
        common::SortedMultiset<KMER, uint16_t, Vector<std::pair<KMER, uint16_t>>>>;

template <typename KMER>
using KmerMultsetVector32 = kmer::KmerCollector<
        KMER,
        KmerExtractorBOSS,
        common::SortedMultiset<KMER, uint32_t, Vector<std::pair<KMER, uint32_t>>>>;

template <typename KMER>
using KmerSetDisk
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, common::SortedSetDisk<KMER, typename KMER::WordType>>;

template <typename KMER>
using SortedMultisetDisk8
        = common::SortedMultisetDisk<KMER, typename KMER::WordType, uint8_t>;
template <typename KMER>
using KmerMultsetDiskVector8
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, SortedMultisetDisk8<KMER>>;

template <typename KMER>
using SortedMultisetDisk16
        = common::SortedMultisetDisk<KMER, typename KMER::WordType, uint16_t>;
template <typename KMER>
using KmerMultsetDiskVector16
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, SortedMultisetDisk16<KMER>>;

template <typename KMER>
using SortedMultisetDisk32
        = common::SortedMultisetDisk<KMER, typename KMER::WordType, uint32_t>;
template <typename KMER>
using KmerMultsetDiskVector32
        = kmer::KmerCollector<KMER, KmerExtractorBOSS, SortedMultisetDisk32<KMER>>;

std::unique_ptr<IBOSSChunkConstructor>
IBOSSChunkConstructor::initialize(size_t k,
                                  bool canonical_mode,
                                  uint8_t bits_per_count,
                                  const std::string &filter_suffix,
                                  size_t num_threads,
                                  double memory_preallocated,
                                  kmer::ContainerType container_type,
                                  const std::filesystem::path &tmp_dir,
                                  size_t max_disk_space_bytes) {
#define OTHER_ARGS k, canonical_mode, bits_per_count, filter_suffix, \
                   num_threads, memory_preallocated, tmp_dir, max_disk_space_bytes

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
