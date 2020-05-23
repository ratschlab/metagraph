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
inline void push_back(Container &kmers, const KMER &kmer) {
    if constexpr(utils::is_pair_v<typename Container::value_type>) {
        kmers.emplace_back(kmer, 0);
    } else {
        kmers.push_back(kmer);
    }
}

/**
 * Generates non-redundant dummy sink kmers (a1a2...ak->$) for the given #kmers_p.
 *
 * The method traverses #kmers_p and generates the corresponding dummy sink k-mer for each
 * kmer in #kmers_p. The dummy sink k-mer is added if it's not redundant, i.e. if there is
 * no kmer in #kmers_p with the same prefix. Since #kmers_p are ordered, the generated
 * dummy sink k-mers for a given edge label will also be ordered. The method keeps an
 * iterator for each edge label that points to the first k-mer with prefix greater or
 * equal to the current dummy k-mer. The current dummy k-mer is marked as redundant if its
 * corresponding iterator points to a kmer with identical prefix.
 * @tparam T the type of kmers being processed, typically either a KMer64/128/256 or an
 * std::pair<KMer64/128/256, int8/16/32> if counting kmers.
 * @param k node length (the actual k-mer has length k+1)
 * @param kmers_p list of sorted, non-redundant kmers of length k+1
 */
template <typename T>
void add_dummy_sink_kmers(size_t k, Vector<T> *kmers_p) {
    using KMER = get_first_type_t<T>;
    // TODO: this is a bit of a waste as we usually have less than 2^bits characters
    // Maybe define an ALPHABET_SIZE constant in KmerExtractorBOSS and KmerBOSS?
    constexpr uint8_t ALPHABET_LEN = 1 << KmerExtractorBOSS::bits_per_char;
    using INT = typename KMER::WordType;

    Vector<T> &kmers = *kmers_p;

    // points to the current k-mer with the given first character
    std::array<const T*, ALPHABET_LEN + 1> first_char_it;
    std::vector<INT> zeros(k + 1);
    for (uint32_t i = 1; i < ALPHABET_LEN; ++i) {
        zeros[k - 1] = i;
        first_char_it[i] = std::lower_bound(kmers.data(), kmers.data() + kmers.size(),
                                            KMER(zeros, k + 1), // the $$...i->$ k-mer
                                            [](const T &a, const KMER &b) -> bool {
                                                return get_first(a) < b;
                                            });
    }
    first_char_it[ALPHABET_LEN] = kmers.data() + kmers.size();

    std::vector<KMER> last_dummy(ALPHABET_LEN, KMER(0));
    size_t size = kmers.size();
    for (size_t i = 1; i < size; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // none of the original k-mers is a dummy k-mer
        assert(kmer[1] != 0 && kmer[0] != 0);

        KMER dummy_sink_kmer = kmer;
        dummy_sink_kmer.to_next(k + 1, BOSS::kSentinelCode);

        TAlphabet last_char = kmer[0];
        if (last_dummy[last_char] == dummy_sink_kmer) {
            continue; // avoid generating duplicate dummy sink kmers
        }
        last_dummy[last_char] = dummy_sink_kmer;

        while (first_char_it[last_char] < first_char_it[last_char + 1]
                && KMER::less(get_first(*first_char_it[last_char]), dummy_sink_kmer)) {
            first_char_it[last_char]++;
        }
        if (!KMER::compare_suffix(get_first(*first_char_it[last_char]), dummy_sink_kmer)
                || first_char_it[last_char] == first_char_it[last_char + 1]) {
            push_back(kmers, dummy_sink_kmer);
        }
    }
}

/**
 * Adds non-redundant dummy source nodes with prefix length 1 for the given kmers.
 *
 * For each kmer in #kmers_p, the method creates the corresponding dummy sources with
 * sentinels of length 1 (aka dummy-1 sources). The method then checks if the dummy-1
 * source node is redundant, i.e. if there is another non-dummy node that is identical
 * except for the first character. For example, $ACGT is redundant with TACGT. To do this
 * efficiently, the method uses an iterator called dummy_it that points to the first kmer
 * with the suffix equal to or larger than the current kmer. Since #kmers_p are ordered,
 * the generated dummy-1 kmers will also be ordered as long as the last character
 * in the original k-mers stays the same. This means that dummy_it only needs to move
 * forward and must be reset only when the last character in the original kmer changes.
 *
 * For example, the k-mer ACG->T, generates the dummy-kmer $AC->G, while the k-mer AGG->T
 * generates the dummy k-mer $AG->T. Because ACG->T precedes AGG->T, so will their
 * generated dummy-1 kmers. To check if $AC->G is redundant, dummy_it is advanced until we
 * find or pass k-mers with the suffix AC. Then we check all kmers with the AC suffix, and
 * if we find one with the same edge label, such as TAC->G, then $AC->G is redundant and
 * will be skipped.
 */
template <typename T>
void add_dummy_source_kmers(size_t k, Vector<T> *kmers_p, size_t end) {
    using KMER = get_first_type_t<T>;
    Vector<T> &kmers = *kmers_p;

    // points to the first k-mer that may be redundant with the current dummy source k-mer
    size_t real_it = 0;

    for (size_t i = 1; i < end; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // none of the original k-mers is a dummy k-mer
        assert(kmer[1] != 0 && kmer[0] != 0);

        if (kmer[k] != get_first(kmers[i - 1])[k]) {
            // the last (most significant) character changed -> reset the iterator
            real_it = 0;
        }

        if (KMER::compare_suffix(kmer, get_first(kmers[i - 1])))
            continue; // i would generate the same dummy-1 k-mer as i-1, skip

        KMER prev_kmer = kmer;
        prev_kmer.to_prev(k + 1, BOSS::kSentinelCode);

        while (real_it < end && KMER::less(get_first(kmers[real_it]), prev_kmer, 1)) {
            real_it++;
        }
        bool is_redundant = false;
        while (real_it < end
                && KMER::compare_suffix(get_first(kmers[real_it]), prev_kmer, 1)) {
            if (get_first(kmers[real_it])[0] == prev_kmer[0]) {  // edge labels match
                is_redundant = true;
                break;
            }
            real_it++;
        }
        if (!is_redundant) {
            push_back(kmers, prev_kmer);
        }
    }
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
/**
 * Adds dummy nodes for the given kmers.
 * The method first adds dummy sink kmers, then dummy sources with sentinels of length 1
 * (aka dummy-1 sources). The method will then gradually add dummy sources with sentinels
 * of length 2, 3, ... up to k-1.
 *
 * @param k the node length in the BOOS graph (so k-mer length is k+1)
 * @tparam T the type of kmers being processed, typically either a KMer64/128/256 or an
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
            if (KMER::compare_suffix(kmer, get_first(kmers[i - 1])))
                continue; // i would generate the same dummy source k-mer as i-1, skip

            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            push_back(kmers, kmer);
        }
        dummy_begin = dummy_end;
        sort(&kmers, num_threads, dummy_begin);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      kmers.size() - dummy_begin);
    }
    sort(&kmers, num_threads, 0);
}

/**
 * Contains recently added kmers and dummy-1 source kmers to facilitate checking for
 * redundancy. A dummy-1 source k-mer that is redundant is set to 0 and will be ignred.
 * */
template <typename KMER>
using RecentKmers = common::CircularBuffer<KMER>;
template <typename T>
using Encoder = common::EliasFanoEncoderBuffered<T>;
template <typename T>
using Decoder = common::EliasFanoDecoder<T>;

/**
 * Push a k-mer into the buffer, then iterate backwards and check if there is a dummy
 * incoming edge (a dummy edge with identical suffix and label as the current one)
 * into the same node.
 * If such an edge is found, mark it for removal (because it's redundant).
 * @param kmer check for redundancy against this k-mer
 * @param buffer queue that stores the last  |Alphabet|^2 k-mers used for identifying
 * redundant dummy source k-mers
 */
template <typename KMER>
static void remove_redundant_dummy_source(const KMER &kmer, RecentKmers<KMER> *buffer) {
    // kmer should not be a dummy k-mer
    assert(kmer[0] != BOSS::kSentinelCode && kmer[1] != BOSS::kSentinelCode);

    if (buffer->size() < 2) {
        return;
    }
    for (auto it = ++buffer->rbegin(); KMER::compare_suffix(kmer, *it, 1); ++it) {
        KMER &prev_kmer = *it;
        if (prev_kmer[1] == BOSS::kSentinelCode) { // redundant dummy source k-mer
            prev_kmer = KMER(0);
            break;
        }
        if (it.at_begin())
            break;
    }
}

constexpr uint8_t ALPHABET_LEN = 1 << KmerExtractorBOSS::bits_per_char;
constexpr uint32_t ENCODER_BUFFER_SIZE = 100'000;

/**
 * Splits #kmers by the first and last character into ALPHABET_LEN^2 blocks.
 */
template <typename T>
std::vector<std::string>
split(size_t k, const std::filesystem::path &tmp_dir, const ChunkedWaitQueue<T> &kmers) {
    logger->trace("Splitting k-mers into {} chunks...", ALPHABET_LEN);
    using T_INT = get_int_t<T>;
    std::vector<Encoder<T_INT>> sinks;
    std::vector<std::string> names(std::pow(ALPHABET_LEN, 2));
    for (uint32_t i = 0; i < names.size(); ++i) {
        names[i] = tmp_dir/("original_split_by01_" + std::to_string(i));
        sinks.emplace_back(names[i], ENCODER_BUFFER_SIZE);
    }

    size_t num_parent_kmers = 0;
    for (auto &it = kmers.begin(); it != kmers.end(); ++it) {
        num_parent_kmers++;
        const T &kmer = *it;
        TAlphabet first_char = get_first(kmer)[k];
        TAlphabet last_char = get_first(kmer)[0];
        uint32_t idx = first_char * ALPHABET_LEN + last_char;
        sinks[idx].add(reinterpret_cast<const T_INT&>(kmer));
    }
    std::for_each(sinks.begin(), sinks.end(), [](auto &f) { f.finish(); });
    logger->trace("Total number of k-mers: {}", num_parent_kmers);
    return names;
}


template <typename T_INT, typename KMER, typename INT>
void handle_dummy_sink(size_t k,
                       KMER kmer,
                       common::MergeDecoder<T_INT> &decoder,
                       Encoder<INT> *dummy_sink_enc,
                       INT *last_dummy_sink) {
    kmer.to_next(k + 1, BOSS::kSentinelCode);
    if (*last_dummy_sink == kmer.data()) { // generate only distinct dummy sink k-mers
        return;
    }
    *last_dummy_sink = kmer.data();

    while (decoder.top().has_value() && get_first(decoder.top().value()) < kmer.data()) {
        decoder.next();
    }
    if (!decoder.top().has_value()
        || !KMER::compare_suffix(
                reinterpret_cast<const KMER &>(get_first(decoder.top().value())), kmer)) {
        dummy_sink_enc->add(kmer.data());
    }
}

/**
 * Writes #to_write to #dummy_l1 if it's a dummy-1 kmer and it wasn't marked as redundant.
 */
template <typename T, typename INT>
void write_dummy_l1(const T &to_write, Encoder<INT> *dummy_l1) {
    auto kmer = get_first(to_write);
    if (kmer.data() == 0U) { // redundant dummy source k-mer
        return;
    }
    assert(kmer[0] != BOSS::kSentinelCode);

    const TAlphabet node_last_char = kmer[1];
    if (node_last_char == BOSS::kSentinelCode) { // only write dummy-1 source k-mers
        dummy_l1->add(kmer.data());
    }
}

template <typename KMER, typename INT>
void push_to_recent(const KMER &v, RecentKmers<KMER> *buffer, Encoder<INT> *dummy_l1) {
    assert(v[0] != BOSS::kSentinelCode); // no dummy sink k-mers allowed
    buffer->push_back(v);
    if (buffer->full()) {
        write_dummy_l1(buffer->pop_front(), dummy_l1);
    }
}

/**
 * Handles the next dummy source k-mer.
 * First, all non-dummy k-mers smaller than #dummy_source are written to #recent_buffer,
 * while also checking if any of the pushed k-mers is redundant with a previous dummy
 * source k-mer. For each of these k-mers, also generate the corresponding dummy sink k-mer
 * and check if it's redundant against #dummy_sink_it.
 * Finally, the dummy source k-mer is written to #recent_buffer.
 * @param k node length (k-mer length is k+1)
 * @param dummy_source_it traverses all k-mers with the same last char as #dummy_source
 * @param dummy_sink_it traverses all k-mers with the same *first* char as #dummy_source's last char
 * @param recent_buffer stores the recently seen non-dummy and dummy-1 k-mers (in order to check for redundancy)
 * @param dummy_l1 sink where the non-dummy1 source kmers are written
 * @param dummy_sink_chunk sink where the non-dummy sink k-mers are written
 * @param last_dummy_sink the last written dummy sink k-mer (to avoid writing duplicates)
 */
template <typename T_INT, typename KMER, typename INT>
void handle_dummy_source(size_t k,
                         const INT &dummy_source,
                         common::MergeDecoder<T_INT> &dummy_source_it,
                         common::MergeDecoder<T_INT> &dummy_sink_it,
                         RecentKmers<KMER> *recent_buffer,
                         Encoder<INT> *dummy_l1,
                         Encoder<INT> *dummy_sink_chunk,
                         INT *last_dummy_sink) {
    // Push all non-dummy kmers smaller than dummy source
    while (dummy_source_it.top().has_value()
           && get_first(dummy_source_it.top().value()) < dummy_source) {
        const KMER v(get_first(dummy_source_it.next().value()));
        push_to_recent(v, recent_buffer, dummy_l1);
        remove_redundant_dummy_source(v, recent_buffer);
        // push the dummy sink k-mer corresponding to v
        handle_dummy_sink(k, get_first(v), dummy_sink_it, dummy_sink_chunk, last_dummy_sink);
    }
}

/**
 * Generates non-redundant dummy-1 source k-mers and dummy sink kmers from #kmers.
 * @return a triplet containing the names of the original k-mer blocks, the dummy-1 source
 * k-mer blocks and the dummy sink k-mers
 */
template <typename T>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::string>
generate_dummy_1_kmers(size_t k,
                       const std::filesystem::path &tmp_dir,
                       ChunkedWaitQueue<T> *kmers) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using T_INT = get_int_t<T>; // either KMER::WordType or <KMER::WordType, count>
    using INT = typename KMER::WordType; // 64/128/256-bit integer

    // for a DNA alphabet, this will contain 16 chunks, split by kmer[0] and kmer[1]
    std::vector<std::string> original_names = split(k, tmp_dir, *kmers);

    std::vector<Encoder<INT>> dummy_l1_chunks;
    std::vector<Encoder<INT>> dummy_sink_chunks;
    std::vector<std::string> dummy_l1_names(ALPHABET_LEN);
    std::vector<std::string> dummy_sink_names(ALPHABET_LEN);
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        dummy_l1_names[i] = tmp_dir/("dummy_source_1_" + std::to_string(i));
        dummy_sink_names[i] = tmp_dir/("dummy_sink_" + std::to_string(i));
        dummy_l1_chunks.emplace_back(dummy_l1_names[i], ENCODER_BUFFER_SIZE);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i], ENCODER_BUFFER_SIZE);
    }

    logger->trace("Generating dummy-1 source kmers and dummy sink k-mers...");
    #pragma omp parallel for num_threads(4) schedule(static, 1)
    for (uint32_t first_ch = 1; first_ch < ALPHABET_LEN; ++first_ch) {  // skip $$..$
        RecentKmers<KMER> recent_buffer(1llu << 2 * KMER::kBitsPerChar);
        INT last_dummy_sink = 0;
        std::vector<std::string> first_ch_names(
                original_names.begin() + first_ch * ALPHABET_LEN,
                original_names.begin() + (first_ch + 1) * ALPHABET_LEN);
        std::vector<std::string> last_ch_names(ALPHABET_LEN);
        for(uint32_t i = 0; i<ALPHABET_LEN; ++i){
            last_ch_names[i] = original_names[i*ALPHABET_LEN + first_ch];
        }

        common::MergeDecoder<T_INT> it(first_ch_names, false);
        common::MergeDecoder<T_INT> dummy_sink_it(first_ch_names, false);
        common::MergeDecoder<T_INT> dummy_source_it(last_ch_names, false);
        INT prev_dummy_source(0);
        for (auto v = it.next(); v.has_value(); v = it.next()) {
            KMER dummy_source(get_first(v.value()));
            dummy_source.to_prev(k + 1, BOSS::kSentinelCode);
            if (dummy_source.data() != prev_dummy_source) {
                handle_dummy_source(k, dummy_source.data(), dummy_source_it, dummy_sink_it,
                                    &recent_buffer, &dummy_l1_chunks[first_ch],
                                    &dummy_sink_chunks[first_ch], &last_dummy_sink);
                push_to_recent(dummy_source, &recent_buffer, &dummy_l1_chunks[first_ch]);
                prev_dummy_source = dummy_source.data();
            }
        }
        // push leftover dummy_source_it
        handle_dummy_source(k, INT(-1), dummy_source_it, dummy_sink_it,
                            &recent_buffer, &dummy_l1_chunks[first_ch],
                            &dummy_sink_chunks[first_ch], &last_dummy_sink);
        while (!recent_buffer.empty()) { // add leftover elements from buffer
            write_dummy_l1(recent_buffer.pop_front(), &dummy_l1_chunks[first_ch]);
        }
    }
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        dummy_sink_chunks[i].finish();
        dummy_l1_chunks[i].finish();
    }

    // dummy sink k-mers are partitioned into blocks by kmer[1], so simply concatenating
    // the blocks will result in a single ordered block
    std::string dummy_sink_name = tmp_dir/"dummy_sink";
    common::concat(dummy_sink_names, dummy_sink_name);

    // similarly, the 16 blocks of the original k-mers can be concatenated in groups of
    // 4 without destroying the order
    std::vector<std::string> original_merged_names(ALPHABET_LEN);
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        original_merged_names[i] = tmp_dir/("original_split1_" + std::to_string(i));
        std::vector<std::string> source(ALPHABET_LEN);
        for (uint32_t j = 0; j < ALPHABET_LEN; ++j) {
            source[j] = original_names[j * ALPHABET_LEN + i];
        }
        common::concat(source, original_merged_names[i]);
    }
    return { original_merged_names, dummy_l1_names, dummy_sink_name };
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
void recover_dummy_nodes_disk(const KmerCollector &kmer_collector,
                              ChunkedWaitQueue<T> *kmers,
                              ThreadPool &async_worker) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using T_INT = get_int_t<T>; // either KMER_INT or <KMER_INT, count>
    using INT = typename KMER::WordType; // 64/128/256-bit integer

    size_t k = kmer_collector.get_k() - 1;
    const std::filesystem::path tmp_dir = kmer_collector.tmp_dir();

    std::string dummy_sink_name;
    std::vector<std::string> original_names;
    std::vector<std::string> dummy_names;
    std::tie(original_names, dummy_names, dummy_sink_name)
            = generate_dummy_1_kmers(k, tmp_dir, kmers);

    // stores the sorted original kmers and dummy-1 k-mers
    std::vector<std::string> files_to_merge;
    files_to_merge.push_back(dummy_sink_name);
    std::vector<std::string> dummy_next_names(ALPHABET_LEN);
    // generate dummy k-mers of prefix length 2..k
    logger->trace("Starting generating dummy-k source k-mers...");
    for (size_t dummy_pref_len = 1; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will compress all sorted dummy k-mers of given prefix length
        files_to_merge.push_back(tmp_dir/("dummy_l" + std::to_string(dummy_pref_len)));
        Encoder<INT> encoder(files_to_merge.back(), ENCODER_BUFFER_SIZE);

        std::vector<Encoder<INT>> dummy_next_chunks;
        for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
            dummy_next_names[i] = tmp_dir/("dummy_source_"
                    + std::to_string(dummy_pref_len + 1) + "_" + std::to_string(i));
            dummy_next_chunks.emplace_back(dummy_next_names[i], ENCODER_BUFFER_SIZE);
        }
        const std::function<void(const INT &)> &write_dummy = [&](const INT &v) {
            encoder.add(v);
            KMER kmer(v);
            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            dummy_next_chunks[kmer[0]].add(kmer.data());
        };
        size_t num_kmers = common::merge_files(dummy_names, write_dummy);

        encoder.finish();
        std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                      [](auto &v) { v.finish(); });
        std::swap(dummy_names, dummy_next_names);
        logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}",
                      dummy_pref_len, num_kmers);
    }
    const std::function<void(const INT &)> on_merge = [](const INT& ) {};
    common::merge_files(dummy_names, on_merge);
    // at this point, we have the original k-mers and dummy-1 k-mers in original_and_dummy_l1,
    // the dummy-x k-mers in dummy_source_{x}, and we merge them all into a single stream
    kmers->reset();
    async_worker.enqueue([kmers, original_names, files_to_merge]() {
        std::function<void(const T_INT &)> on_new_item
                = [kmers](const T_INT &v) { kmers->push(reinterpret_cast<const T &>(v)); };
        common::merge_dummy(original_names, files_to_merge, on_new_item);
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
                recover_dummy_nodes_disk(kmer_collector_, &kmers, async_worker_);
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
