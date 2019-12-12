#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "boss_chunk.hpp"
#include "common/circular_buffer.hpp"
#include "common/file_merger.hpp"
#include "common/logger.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_set_disk.hpp"
#include "common/unix_tools.hpp"
#include "kmer/kmer_collector.hpp"
#include "utils/template_utils.hpp"

namespace mg {
namespace succinct {

using namespace mg;


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
 * What type of data structure to use in the #KmerSet for k-mer storage.
 * Note: for k-mer counting, only VECTOR is supported.
 */
enum class ExtractorContainer {
    VECTOR,
    /**
     * Uses several vectors that are written to disk and then merged, as defined
     * in #SortedSetDisk
     */
    VECTOR_DISK
};

constexpr static ExtractorContainer kExtractorContainer = ExtractorContainer::VECTOR_DISK;

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
                  size_t offset) {
    common::logger->trace("Allocated capacity exceeded, filter out non-unique k-mers...");

    size_t prev_num_kmers = kmers->size();
    sort_and_remove_duplicates(kmers, num_threads, offset);

    common::logger->trace(" done. Number of kmers reduced from {} to {}, {}Mb",
                          prev_num_kmers, kmers->size(),
                          (kmers->size() * sizeof(typename Array::value_type) >> 20));
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
                                size_t /* alphabet size */) {
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
            shrink_kmers(kmers, num_threads, dummy_begin);

        push_back(*kmers, kmer).to_prev(k + 1, BOSS::kSentinelCode);
    }
    common::logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                          num_dummy_parent_kmers);
    sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

    common::logger->trace("Number of dummy k-mers with dummy prefix of length 2: {}",
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

        common::logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}",
                              c, kmers->size() - dummy_begin);
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
template < typename T,typename TAlphabet>
static void remove_redundant_dummy_source(const T &kmer, RecentKmers<T> *buffer) {
    TAlphabet curW = kmer[0];
    if (buffer->size() < 2 || curW == 0) {
        return;
    }
    using Iterator = typename RecentKmers<T>::iterator;
    Iterator it = buffer->rbegin();
    --it;
    for (; T::compare_suffix(kmer, utils::get_first((*it).kmer), 1); --it) {
        const T prev_kmer = utils::get_first((*it).kmer);
        if (prev_kmer[0] == curW && !prev_kmer[1]) { // redundant dummy source k-mer
            (*it).is_removed = true;
            return;
        }
        if (it.at_begin())
            break;
    }
}


template <typename T>
void write_or_die(std::fstream *f, const T &v) {
    if (!f->write(reinterpret_cast<const char *>(&v), sizeof(T))) {
        std::cerr << "Error: Writing of merged data failed." << std::endl;
        std::exit(EXIT_FAILURE);
    }
};

/**
 * Writes the kmer at the front of #buffer to a file if it's not a redundant dummy
 * kmer. If the k-mer at the front is a (not-redundant) dummy k-mer it also writes its
 * corresponding dummy k-mer of prefix length 2 into #sorted_dummy_kmers.
 */
template <typename T>
uint8_t write_kmer(size_t k,
                   std::fstream *f,
                   Vector<T> *dummy_kmers,
                   ThreadPool *file_write_pool,
                   RecentKmers<T> *buffer,
                   common::SortedSetDisk<T> *sorted_dummy_kmers) {
    const Kmer<T> to_write = buffer->pop_front();
    if (to_write.is_removed) { // redundant dummy, &file_write_pool k-mer
        return 0;
    }
    file_write_pool->enqueue(write_or_die<T>, f, to_write.kmer);
    using KMER = std::remove_reference_t<decltype(utils::get_first(to_write.kmer))>;
    using TAlphabet = typename T::CharType;

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
                                size_t, /* num_threads unused */
                                size_t alphabet_size) {
    // name of the file containing dummy k-mers of given prefix length
    const auto get_file_name
            = [](uint32_t pref_len) { return "/tmp/dummy" + std::to_string(pref_len); };

    const auto no_cleanup = [](typename common::SortedSetDisk<T>::storage_type *) {};

    // asynchronously writes a value of type T to a file stream
    // TODO(ddanciu): profile this, as I suspect it's slow - use a WaitQueue instead
    ThreadPool file_write_pool = ThreadPool(1, 10000);
    const auto async_file_writer = [&file_write_pool](std::fstream &f) {
        return [&f, &file_write_pool](const T &v) {
            file_write_pool.enqueue(write_or_die<T>, &f, v);
        };
    };

    std::vector<std::pair<std::string, std::fstream *>> files_to_merge;
    auto create_stream = [](const std::string &name) {
        std::fstream *f = new std::fstream(name, std::ios::out | std::ios::binary);
        return std::pair<std::string, std::fstream *>({ name, f });
    };

    std::string file_name = "/tmp/original_and_dummy_l1";
    files_to_merge.push_back(create_stream(file_name));
    std::fstream *dummy_l1 = files_to_merge.back().second;

    using TAlphabet = typename T::CharType;
    RecentKmers<T> recent_buffer(alphabet_size * alphabet_size);

    std::string file_name_l2 = get_file_name(2);
    files_to_merge.push_back(create_stream(file_name_l2));
    std::fstream *dummy_l2 = files_to_merge.back().second;

    // this will contain dummy k-mers of prefix length 2
    common::SortedSetDisk<T> sorted_dummy_kmers(no_cleanup, async_file_writer(*dummy_l2));
    Vector<T> dummy_kmers;
    dummy_kmers.reserve(sorted_dummy_kmers.capacity());

    // remove redundant dummy source k-mers of prefix length 1 and write them to a file
    // While traversing and removing redundant dummy source k-mers of prefix length 1,
    // we also  generate dummy k-mers of prefix length 2.
    size_t num_dummy_parent_kmers = 0;
    for (auto &it = kmers->begin(); it != kmers->end(); ++it) {
        const T el = *it;
        recent_buffer.push_back({ el, false });
        remove_redundant_dummy_source<T, TAlphabet>(el, &recent_buffer);
        if (recent_buffer.full()) {
            num_dummy_parent_kmers += write_kmer(k, dummy_l1, &dummy_kmers, &file_write_pool,
                                                 &recent_buffer, &sorted_dummy_kmers);
        }
    }
    while (!recent_buffer.empty()) { // empty the buffer
        num_dummy_parent_kmers += write_kmer(k, dummy_l1, &dummy_kmers, &file_write_pool,
                                             &recent_buffer, &sorted_dummy_kmers);
    }

    // push out the leftover dummy kmers
    sorted_dummy_kmers.insert(dummy_kmers.begin(), dummy_kmers.end());

    common::logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                          num_dummy_parent_kmers);

    // generate dummy k-mers of prefix length 3..k
    common::SortedSetDisk<T> sorted_dummy_kmers2;
    common::SortedSetDisk<T> *source = &sorted_dummy_kmers;
    common::SortedSetDisk<T> *dest = &sorted_dummy_kmers2;
    for (size_t dummy_pref_len = 3; dummy_pref_len < k + 1; ++dummy_pref_len) {
        std::string file_name = get_file_name(dummy_pref_len);
        files_to_merge.push_back(create_stream(file_name));
        dest->clear(async_file_writer(*(files_to_merge.back().second)));
        dummy_kmers.resize(0);
        size_t num_kmers = 0;
        for (auto &it = source->data().begin(); it != source->data().end(); ++it) {
            const T &kmer = *it;
            if (dummy_kmers.size() == dummy_kmers.capacity()) {
                dest->insert(dummy_kmers.begin(), dummy_kmers.end());
                dummy_kmers.resize(0);
            }


            push_back(dummy_kmers, utils::get_first(kmer)).to_prev(k + 1, BOSS::kSentinelCode);
            num_kmers++;
        }
        // push out the leftover dummy kmers
        dest->insert(dummy_kmers.begin(), dummy_kmers.end());


        common::logger->trace("Number of dummy k-mers with dummy prefix of length {} :",
                              (dummy_pref_len - 1), num_kmers);
        std::swap(source, dest);
    }
    for (auto &it = source->data().begin(); it != source->data().end(); ++it) {
        // we iterate to merge the data and write it to disk
    }

    file_write_pool.join();

    // at this point, we have the original k-mers plus the  dummy k-mers with prefix
    // length x in /tmp/dummy{x}, and we'll merge them all into a single stream
    kmers->reset();
    std::vector<std::string> file_names;
    std::for_each(files_to_merge.begin(), files_to_merge.end(), [&file_names](auto &el) {
        delete el.second; // this will also close the stream
        file_names.push_back(el.first);
    });
    common::merge_files(file_names, kmers);
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
                         double memory_preallocated = 0)
          : kmer_storage_(k + 1,
                          canonical_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated) {
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
            common::logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;

            // kmer_collector stores (BOSS::k_ + 1)-mers
            recover_source_dummy_nodes(kmer_storage_.get_k() - 1,
                                       &kmers,
                                       kmer_storage_.num_threads(),
                                       kmer_storage_.alphabet_size());

            common::logger->trace("Dummy source k-mers were reconstructed in {} sec",
                                  timer.elapsed());
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
             double memory_preallocated) {

    #define OTHER_ARGS k, canonical_mode, filter_suffix, num_threads, memory_preallocated

    if (count_kmers) {
        return initialize_boss_chunk_constructor<KmerMultsetVector>(OTHER_ARGS);
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
