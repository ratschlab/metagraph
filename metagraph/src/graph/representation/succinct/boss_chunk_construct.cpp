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
#include "kmer/kmer_to_int_converter.hpp"
#include "boss_chunk.hpp"

namespace mg {
namespace succinct {

using namespace mg;
using common::logger;
using common::ChunkedWaitQueue;
using kmer::KmerCollector;
using TAlphabet = KmerExtractorBOSS::TAlphabet;
using utils::get_first;
using utils::get_first_type_t;


template <class Array>
void sort(Array *array, size_t num_threads, size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          utils::LessFirst(),
                          num_threads);
}

template <class Container, typename KMER>
void push_back(Container &kmers, const KMER &kmer) {
    if constexpr(utils::is_pair_v<typename Container::value_type>) {
        kmers.emplace_back(kmer, 0);
    } else {
        kmers.push_back(kmer);
    }
}

template <typename T>
void add_dummy_sink_kmers(size_t k, Vector<T> *kmers_p) {
    using KMER = get_first_type_t<T>;
    // TODO: this is a bit of a waste as we usually have less than 2^bits characters
    // Maybe define an ALPHABET_SIZE constant in KmerExtractorBOSS and KmerBOSS?
    constexpr uint8_t ALPHABET_LEN = 1 << KmerExtractorBOSS::bits_per_char;
    using INT = typename KMER::WordType;

    Vector<T> &kmers = *kmers_p;

    // points to the current k-mer with the given first character
    std::array<size_t, ALPHABET_LEN> first_char_it;
    std::array<size_t, ALPHABET_LEN> max_it;
    first_char_it[0] = 0;
    std::vector<INT> zeros(k + 1);
    for (uint32_t i = 1; i < ALPHABET_LEN; ++i) {
        zeros[k - 1] = i;
        KMER to_search = KMER(zeros, k + 1); // the $...$X->$ k-mer
        first_char_it[i] = std::lower_bound(kmers.begin(), kmers.end(), to_search,
                                            [](const T &a, const KMER &b) -> bool {
                                                return get_first(a) < b;
                                            })
                - kmers.begin();
        if (first_char_it[i] < kmers.size()) {
            std::cout << "\nTo search: " << first_char_it[i] << "\t"
                      << get_first(kmers[first_char_it[i]]).to_string(k + 1, "$ACGT")
                      << std::endl;
        }
        max_it[i - 1] = first_char_it[i];
    }
    max_it[ALPHABET_LEN - 1] = kmers.size();
    max_it[0] = kmers.size(); // only used when k==2 bc last char is $ for all dummy k-mers

    std::vector<KMER> last_dummy(ALPHABET_LEN, KMER(0));
    size_t size = kmers.size();
    for (size_t i = 1; i < size; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        KMER dummy_sink_kmer = kmer;
        dummy_sink_kmer.to_next(k + 1, BOSS::kSentinelCode);
        std::cout << "Dummy sink k-mer is: " << dummy_sink_kmer.to_string(k + 1, "$ACGT")
                  << std::endl;

        TAlphabet last_char = kmer[0];
        if (last_dummy[last_char] == dummy_sink_kmer) {
            continue;  // avoid generating duplicate dummy sink kmers
        }
        last_dummy[last_char] = dummy_sink_kmer;
        size_t &dummy_sink_it = first_char_it[last_char];

        while (dummy_sink_it < max_it[last_char]
               && KMER::less(get_first(kmers[dummy_sink_it]), dummy_sink_kmer)) {
            dummy_sink_it++;
        }
        if (dummy_sink_it < max_it[last_char]
            && !KMER::compare_suffix(get_first(kmers[dummy_sink_it]), dummy_sink_kmer)) {
            std::cout << "Not redundant\n";
            push_back(kmers, dummy_sink_kmer);
        }
        std::cout << "Redundant with "
                  << get_first(kmers[dummy_sink_it]).to_string(k + 1, "$ACGT") << std::endl;
    }
}

template <typename T>
void add_dummy_source_kmers(size_t k, Vector<T> *kmers_p, size_t end) {
    using KMER = get_first_type_t<T>;

    // points to the first k-mer that may be redundant with the current dummy source k-mer
    size_t dummy_it = 0;

    Vector<T> &kmers = *kmers_p;

    for (size_t i = 1; i < end; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        if (kmer[k] != get_first(kmers[i - 1])[k]) {
            // the last (most significant) character changed, need to start search from beginning
            dummy_it = 0;
        }

        if (KMER::compare_suffix(kmer, get_first(kmers[i -1 ]))) {
            continue; // i would generate the same dummy-1 k-mer as i-1, skip
        }

        KMER prev_kmer = kmer;
        prev_kmer.to_prev(k + 1, BOSS::kSentinelCode);

        while (dummy_it < end && KMER::less(get_first(kmers[dummy_it]), prev_kmer, 1)) {
            dummy_it++;
        }
        std::cout << "Checking  " << prev_kmer.to_string(k+1, "$ACGT") << " starting at " << get_first(kmers[dummy_it]).to_string(k+1, "$ACGT") << std::endl;
        size_t cur_pos = dummy_it;
        bool is_redundant = false;
        while (cur_pos < end
               && KMER::compare_suffix(get_first(kmers[cur_pos]), prev_kmer, 1)) {
            const KMER &current_kmer = get_first(kmers[cur_pos]);
            if (prev_kmer[0] == current_kmer[0]) {  // edge labels match
                is_redundant = true;
                break;
            }
            cur_pos++;
        }
        if (!is_redundant) {
            push_back(kmers, prev_kmer);
            std::cout << prev_kmer.to_string(k+1, "$ACGT") << " is not redundant!!!\n";
        }
    }
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
// k is node length
/**
 * Adds dummy source nodes for the given kmers.
 * The method first adds dummy sources with sentinels of length 1 (aka dummy-1 sources).
 * For each dummy-1 source node it checks if it's redundant, i.e. if there is another
 * non-dummy node that is identical except for the first character. For example, $ACGT is
 * redundant with TACGT. To do this efficiently, the method uses an iterator called
 * dummy_it that points to the first kmer with the suffix equal to or larger than the
 * current kmer. Since dummy-1 kmers for a given last character are actually ordered,
 * dummy_it only needs to move forward and must be reset when the last character changes.
 * For example, the k-mer ACG->T, generates the dummy-kmer $AC->G, while the k-mer AGG->T
 * generates the dummy k-mer $AG->T. Because ACG->T precedes AGG->T, so will their
 * generated dummy-1 kmers. To check if $AC->G is redundant, dummy_it is advanced until we
 * find or pass k-mers with the suffix AC. Then we check all kmers with the AC suffix, and
 * if we find one with the same edge label, such as TAC->G, then $AC->G is redundant and
 * will be skipped. The method will then gradually add dummy sources with sentinels of
 * length 2, 3, ... up to k-1.
 *
 * @tparam the type of kmers being processed, typically either a KMer64/128/256 or an
 * std::pair<KMer64/128/256, int8/16/32> if counting kmers.
 */
template <typename T>
void recover_dummy_nodes(size_t k, size_t num_threads, Vector<T> *kmers_p) {
    using KMER = get_first_type_t<T>;
    Vector<T> &kmers = *kmers_p;
    size_t original_end = kmers.size();

    add_dummy_sink_kmers(k, &kmers);
    size_t dummy_begin = kmers.size();

    add_dummy_source_kmers(k, &kmers, original_end);
    sort(&kmers, num_threads, dummy_begin);

    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  kmers.size() - dummy_begin);

    for (size_t c = 2; c < k + 1; ++c) {
        size_t dummy_end = kmers.size();

        for (size_t i = dummy_begin; i < dummy_end; ++i) {
            KMER kmer = get_first(kmers[i]);
            if (KMER::compare_suffix(kmer, get_first(kmers[i - 1]))) {
                continue; // i would generate the same dummy source k-mer as i-1, skip
            }

            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            push_back(kmers, kmer);
        }
        sort(&kmers, num_threads, dummy_end);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      kmers.size() - dummy_end);
        dummy_begin = dummy_end;
    }
    sort(&kmers, num_threads, 0);
    std::ofstream f("/tmp/kmers");
    for (size_t i = 0; i < kmers.size(); ++i) { // starting at 1 to skip the $$...$$ k-mer
        f << get_first(kmers[i]).to_string(k+1, "$ACGT") << std::endl;
    }
}

/**
 * An element of #RecentKmers. Contains a k-mer and a flag that indicates if the k-mer
 * was identified as a redundant dummy source k-mer.
 * */
template <typename T>
struct KmerBuffered {
    T data; // the actual k-mer or k-mer with count
    bool is_removed; // true if this is a redundant dummy source k-mer
};

template <typename T>
using RecentKmers = common::CircularBuffer<KmerBuffered<T>>;

template <typename T_INT, typename INT>
void encode(const INT &v, common::EliasFanoEncoderBuffered<T_INT> *encoder) {
    if constexpr (utils::is_pair_v<T_INT>) {
        encoder->add({ v, 0 });
    } else {
        encoder->add(v);
    }
}

template <typename KMER, typename T_INT>
static void push_dummy_sink(const KMER &kmer,
                            std::optional<KMER> *last_dummy_sink,
                            common::EliasFanoEncoderBuffered<T_INT> *encoder) {
    if (kmer[0] == 0) {
        *last_dummy_sink = kmer;
        return;
    }
    if (last_dummy_sink->has_value()) {
        if( !KMER::compare_suffix(last_dummy_sink->value(), kmer)) {
            encode(last_dummy_sink->value().data(), encoder);
        }
        last_dummy_sink->reset();
    }
}

/**
 * Push a k-mer into the buffer, then iterate backwards and check if there is a dummy
 * incoming edge (a dummy edge with identical suffix and label as the current one)
 * into the same node.
 * If such an edge is found, mark it for removal (because it's redundant).
 * @param kmer check for redundancy against this k-mer
 * @param buffer queue that stores the last  |Alphabet|^2 k-mers used for identifying
 * redundant dummy source k-mers
 */
template <typename T, typename T_INT>
static void push_and_remove_redundant_dummy_source(const T_INT &el,
                                                   RecentKmers<T> *buffer) {
    buffer->push_back({ el, false });

    using KMER = get_first_type_t<T>;
    const KMER &kmer = get_first(el);

    TAlphabet curW = kmer[0];
    if (buffer->size() < 2 || curW == 0) {
        return;
    }
    for (auto it = ++buffer->rbegin();
            KMER::compare_suffix(kmer, get_first((*it).data), 1); ++it) {
        KMER &prev_kmer = get_first((*it).data);
        assert((curW || prev_kmer[1]) && "Main dummy source k-mer must be unique");
        if (prev_kmer[0] == curW && !prev_kmer[1]) { // redundant dummy source k-mer
            (*it).is_removed = true;
            break;
        }
        if (it.at_begin())
            break;
    }
}

/**
 * Writes the front of #buffer to a file if it's not a redundant dummy kmer. If the k-mer
 * at the front is a (not-redundant) dummy k-mer it also writes its corresponding dummy
 * k-mer of prefix length 2 into #sorted_dummy_kmers.
 */
template <typename T, typename T_INT, typename INT>
uint8_t write_kmer(size_t k,
                   KmerBuffered<T> to_write,
                   common::EliasFanoEncoderBuffered<T_INT> *encoder,
                   std::vector<common::EliasFanoEncoderBuffered<INT>> *dummy_kmer_chunks) {
    static_assert(std::is_same_v<T_INT, get_int_t<T>>);

    if (to_write.is_removed) { // redundant dummy k-mer
        return 0;
    }
    if constexpr(utils::is_pair_v<T>) {
        encoder->add({ to_write.data.first.data(), to_write.data.second });
    } else {
        encoder->add(to_write.data.data());
    }
    auto &kmer_to_write = get_first(to_write.data);
    const TAlphabet node_last_char = kmer_to_write[1];
    const TAlphabet edge_label = kmer_to_write[0];
    if (node_last_char || !edge_label) { // not a dummy source kmer
        return 0;
    }

    kmer_to_write.to_prev(k + 1, BOSS::kSentinelCode);
    (*dummy_kmer_chunks)[kmer_to_write[0]].add(kmer_to_write.data());
    return 1;
}

/**
 * Specialization of recover_dummy_nodes for a disk-based container, such as
 * #SortedSetDisk and #SortedMultisetDisk.
 * The method first traverses the original kmers and generates dummy source k-mers of
 * prefix length 1 (dummy-1). At the next step, the non-redundant dummy-1 kmers are merged
 * with the original k-mers.
 * The method then gradually constructs dummy-i k-mers for i=2..k and writes them into
 * separate files, de-duped and sorted.
 * The final result is obtained by merging the original #kmers and dummy-1 kmers with
 * the dummy-k kmers, for k=2..k
 */
template <class KmerCollector, typename T>
void recover_source_dummy_nodes_disk(const KmerCollector &kmer_collector,
                                     ChunkedWaitQueue<T> *kmers,
                                     ThreadPool &async_worker) {
    constexpr size_t ENCODER_BUFFER_SIZE = 100'000;
    constexpr uint8_t ALPHABET_LEN = 1 << KmerExtractorBOSS::bits_per_char;

    const std::filesystem::path tmp_dir = kmer_collector.tmp_dir();
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using T_INT = get_int_t<T>; // either KMER_INT or <KMER_INT, count>
    using INT = typename KMER::WordType; // 64/128/256-bit integer

    size_t k = kmer_collector.get_k() - 1;

    std::string original_name = tmp_dir/"original_kmers";
    common::EliasFanoEncoderBuffered<T_INT> original_kmers(original_name,
                                                           ENCODER_BUFFER_SIZE);

    // traverse the input kmers and remove redundant dummy source k-mers of prefix length
    // 1. While traversing, we also  generate dummy k-mers of prefix length 2.
    size_t num_dummy_l1_kmers = 0;
    size_t num_parent_kmers = 0;

    RecentKmers<T> recent_buffer(std::pow(1llu << KMER::kBitsPerChar, 2));

    std::vector<common::EliasFanoEncoderBuffered<INT>> dummy_l1_chunks;
    std::vector<common::EliasFanoEncoderBuffered<INT>> dummy_l2_chunks;
    std::vector<common::EliasFanoEncoderBuffered<INT>> dummy_sink_chunks;
    std::vector<std::string> dummy_names(2*ALPHABET_LEN);
    std::vector<std::string> dummy_next_names(ALPHABET_LEN);
    std::vector<std::string> dummy_sink_names(ALPHABET_LEN);
    dummy_l1_chunks.reserve(ALPHABET_LEN);
    dummy_l2_chunks.reserve(ALPHABET_LEN);
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        dummy_names[i] = tmp_dir/("dummy_source_1_"+std::to_string(i));
        dummy_next_names[i]= tmp_dir/("dummy_source_2_"+std::to_string(i));
        dummy_names[i+ALPHABET_LEN]= tmp_dir/("dummy_sink_"+std::to_string(i));
        dummy_l1_chunks.emplace_back(dummy_names[i], ENCODER_BUFFER_SIZE);
        dummy_l2_chunks.emplace_back(dummy_next_names[i], ENCODER_BUFFER_SIZE);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i+ALPHABET_LEN], ENCODER_BUFFER_SIZE);
    }

    const KMER all_dummy = std::vector<TAlphabet>(k + 1, BOSS::kSentinelCode);
    for (auto &it = kmers->begin(); it != kmers->end(); ++it) {
        num_parent_kmers++;
        const T& v = *it;
        if constexpr(utils::is_pair_v<T>) {
            original_kmers.add({ v.first.data(), v.second });
        } else {
            original_kmers.add(v.data());
        }
        KMER kmer = get_first(v);
        kmer.to_prev(k + 1, BOSS::kSentinelCode);
        if (kmer != all_dummy) {
            TAlphabet curW = kmer[0];
            dummy_l1_chunks[curW].add(kmer.data());
        }
        KMER dummy_sink = get_first(v);
        dummy_sink.to_next(k+1, BOSS::kSentinelCode);
        if (dummy_sink != all_dummy) {
            TAlphabet first_char = dummy_sink[k];
            dummy_sink_chunks[first_char].add(kmer.data());
        }
    }
    original_kmers.finish();
    std::for_each(dummy_l1_chunks.begin(), dummy_l1_chunks.end(),
                  [](auto &v) { v.finish(); });
    std::for_each(dummy_sink_chunks.begin(), dummy_sink_chunks.end(),
                  [](auto &v) { v.finish(); });

    std::string original_and_dummy_l1_name = tmp_dir/"original_and_dummy_l1";
    common::EliasFanoEncoderBuffered<T_INT> original_and_l1(original_and_dummy_l1_name,
                                                            ENCODER_BUFFER_SIZE);

    std::optional<KMER> last_dummy_sink;

    // merge the original kmers with the dummy-1 kmers, while also removing the redundant
    // dummy-1 kmers
    const std::function<void(const T_INT &)> &on_new_item = [&](const T_INT &v) {
        const T &kmer_count = reinterpret_cast<const T &>(v);
        push_dummy_sink(utils::get_first(kmer_count), &last_dummy_sink, &original_and_l1);
        push_and_remove_redundant_dummy_source(kmer_count, &recent_buffer);
        if (recent_buffer.full()) {
            num_dummy_l1_kmers += write_kmer(k, recent_buffer.pop_front(),
                                             &original_and_l1, &dummy_l2_chunks);
        }
    };
    common::merge_dummy(original_name, dummy_names, on_new_item);
    if (last_dummy_sink.has_value()) {
        encode(last_dummy_sink.value().data(), &original_and_l1);
    }
    dummy_names.resize(ALPHABET_LEN);

    while (!recent_buffer.empty()) { // empty the buffer
        num_dummy_l1_kmers += write_kmer(k, recent_buffer.pop_front(), &original_and_l1,
                                         &dummy_l2_chunks);
    }
    original_and_l1.finish();
    std::for_each(dummy_l2_chunks.begin(), dummy_l2_chunks.end(),
                  [](auto &v) { v.finish(); });

    logger->trace("Total number of k-mers: {}", num_parent_kmers);
    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  num_dummy_l1_kmers);

    // stores the sorted original kmers and dummy-1 k-mers
    std::vector<std::string> files_to_merge;

    // generate dummy k-mers of prefix length 2..k
    for (size_t dummy_pref_len = 2; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will compress all sorted dummy k-mers of given prefix length
        files_to_merge.push_back(tmp_dir/("dummy_l" + std::to_string(dummy_pref_len)));
        common::EliasFanoEncoderBuffered<INT> encoder(files_to_merge.back(),
                                                      ENCODER_BUFFER_SIZE);

        std::vector<common::EliasFanoEncoderBuffered<INT>> dummy_next_chunks;
        std::swap(dummy_names, dummy_next_names);
        dummy_next_names.resize(0);
        dummy_next_chunks.reserve(ALPHABET_LEN);
        for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
            std::string name = tmp_dir/("dummy_source_"
                    + std::to_string(dummy_pref_len + 1) + "_" + std::to_string(i));
            dummy_next_names.push_back(name);
            dummy_next_chunks.emplace_back(dummy_next_names[i], ENCODER_BUFFER_SIZE);
        }
        size_t num_kmers = 0;
        const std::function<void(const INT &)> &write_dummy
                = [k, &num_kmers, &dummy_next_chunks, &encoder, &all_dummy](const INT &v) {
                      encoder.add(v);
                      KMER kmer(v);
                      kmer.to_prev(k + 1, BOSS::kSentinelCode);
                      if (kmer != all_dummy) {
                          const TAlphabet W = kmer[0];
                          dummy_next_chunks[W].add(kmer.data());
                          num_kmers++;
                      }
                  };
        common::merge_files(dummy_names, write_dummy);

        encoder.finish();
        std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                      [](auto &v) { v.finish(); });
        logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}",
                      dummy_pref_len, num_kmers);
    }
    std::for_each(dummy_next_names.begin(), dummy_next_names.end(),
                  [](const string &v) { std::filesystem::remove(v); });
    // at this point, we have the original k-mers plus the  dummy k-mers with prefix
    // length x in /tmp_dir/dummy_source_{x}, and we merge them all into a single stream
    kmers->reset();
    async_worker.enqueue([kmers, original_and_dummy_l1_name, files_to_merge]() {
        std::function<void(const T_INT &)> on_new_item
                = [kmers](const T_INT &v) { kmers->push(reinterpret_cast<const T &>(v)); };
        common::merge_dummy(original_and_dummy_l1_name, files_to_merge, on_new_item);
        kmers->shutdown();
    });
}

inline std::vector<TAlphabet>
encode_filter_suffix_boss(const std::string &filter_suffix) {
    KmerExtractorBOSS kmer_extractor;
    std::vector<TAlphabet> filter_suffix_encoded;
    for (char c : filter_suffix) {
        filter_suffix_encoded.push_back(c == BOSS::kSentinel
                                        ? BOSS::kSentinelCode
                                        : kmer_extractor.encode(c));
    }
    return filter_suffix_encoded;
}

template <typename T_TO, typename T_FROM>
ChunkedWaitQueue<T_TO>& reinterpret_container(ChunkedWaitQueue<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    return reinterpret_cast<ChunkedWaitQueue<T_TO>&>(container);
}

template <typename T_TO, typename T_FROM>
Vector<T_TO>& reinterpret_container(Vector<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    return reinterpret_cast<Vector<T_TO>&>(container);
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
            kmer_collector_.add_kmer(std::vector<TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_sequence(std::string_view sequence, uint64_t count) {
        kmer_collector_.add_sequence(sequence, count);
    }

    void add_sequences(const std::function<void(CallString)> &generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk() {
        typename KmerCollector::Data &kmer_ints = kmer_collector_.data();

        using T_INT = typename KmerCollector::Data::value_type;
        using T = get_kmer_t<typename KmerCollector::Kmer, T_INT>;
        auto &kmers = reinterpret_container<T>(kmer_ints);

        if (!kmer_collector_.suffix_length()) {
            logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;
            if constexpr ((std::is_same_v<typename KmerCollector::Data,
                                          ChunkedWaitQueue<T_INT>>)) {
                recover_source_dummy_nodes_disk(kmer_collector_, &kmers, async_worker_);
            } else {
                // kmer_collector stores (BOSS::k_ + 1)-mers
                static_assert(std::is_same_v<typename KmerCollector::Data, Vector<T_INT>>);
                recover_dummy_nodes(kmer_collector_.get_k() - 1,
                                           kmer_collector_.num_threads(), &kmers);
            }
            logger->trace("Dummy source k-mers were reconstructed in {} sec",
                          timer.elapsed());
        }

        BOSS::Chunk *result;

        if constexpr(utils::is_pair_v<typename KmerCollector::Value>) {
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
    if (k < 1 || k > 256 / KmerExtractorBOSS::bits_per_char - 1) {
        logger->error("For succinct graph, k must be between 2 and {}",
                      256 / KmerExtractorBOSS::bits_per_char - 1);
        exit(1);
    }

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
    = KmerCollector<KMER, KmerExtractorBOSS, common::SortedSet<typename KMER::WordType>>;

template <typename KMER>
using KmerMultsetVector8
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint8_t>>;

template <typename KMER>
using KmerMultsetVector16
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint16_t>>;

template <typename KMER>
using KmerMultsetVector32
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint32_t>>;

template <typename KMER>
using KmerSetDisk
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedSetDisk<typename KMER::WordType>>;

template <typename KMER>
using KmerMultsetDiskVector8
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint8_t>>;

template <typename KMER>
using KmerMultsetDiskVector16
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint16_t>>;

template <typename KMER>
using KmerMultsetDiskVector32
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint32_t>>;

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
