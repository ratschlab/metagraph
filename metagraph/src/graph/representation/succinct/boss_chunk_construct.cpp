#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

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


namespace mtg {
namespace succinct {

using mtg::common::logger;
using mtg::common::ChunkedWaitQueue;
using mtg::kmer::get_int_t;
using mtg::kmer::get_kmer_t;
using mtg::kmer::KmerCollector;
using mtg::kmer::KmerExtractorBOSS;
using mtg::kmer::KmerExtractor2Bit;
using utils::get_first;
using utils::get_first_type_t;
using TAlphabet = KmerExtractorBOSS::TAlphabet;

constexpr size_t ENCODER_BUFFER_SIZE = 100'000;


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
    using KMER_INT = typename KMER::WordType;

    const size_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    Vector<T> &kmers = *kmers_p;

    // points to the current k-mer with the given first character
    std::vector<size_t> max_it(alphabet_size);
    std::vector<size_t> it(alphabet_size);
    for (TAlphabet c = 1; c < alphabet_size; ++c) {
        std::vector<KMER_INT> zeros(k + 1, 0);
        zeros[k - 1] = c;
        it[c] = std::lower_bound(kmers.data(), kmers.data() + kmers.size(),
                                 KMER(zeros, k + 1), // the $$...i->$ k-mer
                                 [](const T &a, const KMER &b) -> bool {
                                     return get_first(a) < b;
                                 }) - kmers.data();
        max_it[c - 1] = it[c];
    }
    max_it[alphabet_size - 1] = kmers.size();

    std::vector<KMER> last_dummy(alphabet_size, KMER(0));
    size_t size = kmers.size();
    for (size_t i = 1; i < size; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // none of the original k-mers is a dummy k-mer
        assert(kmer[1] != 0 && kmer[0] != 0);

        KMER dummy_sink = kmer;
        dummy_sink.to_next(k + 1, BOSS::kSentinelCode);

        TAlphabet last_char = kmer[0];

        while (it[last_char] < max_it[last_char]
                && KMER::less(get_first(kmers[it[last_char]]), dummy_sink)) {
            it[last_char]++;
        }
        if (last_dummy[last_char] != dummy_sink
            && (it[last_char] == max_it[last_char]
                || !KMER::compare_suffix(get_first(kmers[it[last_char]]), dummy_sink))) {
            push_back(kmers, dummy_sink);
            last_dummy[last_char] = dummy_sink;
        }
    }
}

/**
 * Adds non-redundant dummy source nodes with prefix length 1 for the given kmers.
 *
 * For each kmer in #kmers_p, the method creates the corresponding dummy sources with
 * sentinels of length 1 (aka dummy-1 sources). The method then checks if the dummy-1
 * source node is redundant, i.e. if there is another k-mer that is identical
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

    size_t dummy_source_begin = kmers.size();
    add_dummy_source_kmers(k, &kmers, original_end);

    ips4o::parallel::sort(kmers.begin() + dummy_source_begin, kmers.end(),
                          utils::LessFirst(), num_threads);

    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  kmers.size() - dummy_source_begin);

    for (size_t c = 2; c < k + 1; ++c) {
        size_t dummy_end = kmers.size();

        for (size_t i = dummy_source_begin; i < dummy_end; ++i) {
            KMER kmer = get_first(kmers[i]);
            if (KMER::compare_suffix(kmer, get_first(kmers[i - 1])))
                continue; // i would generate the same dummy source k-mer as i-1, skip

            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            push_back(kmers, kmer);
        }
        dummy_source_begin = dummy_end;
        ips4o::parallel::sort(kmers.begin() + dummy_source_begin, kmers.end(),
                              utils::LessFirst(), num_threads);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      kmers.size() - dummy_source_begin);
    }

    ips4o::parallel::sort(kmers.begin(), kmers.end(),
                          utils::LessFirst(), num_threads);
}

template <typename T>
using Encoder = common::EliasFanoEncoderBuffered<T>;
template <typename T>
using Decoder = common::EliasFanoDecoder<T>;

/**
 * Splits #kmers by W (kmer[0]) and F (kmer[k]) into |ALPHABET\{$}|^2 chunks.
 * T_REAL: type KmerExtractorT::KMerBOSS representing k-mers over the alphabet
 * without the sentinel character (e.g., ACGT).
 */
template <typename T_REAL>
std::vector<std::string> split(size_t k,
                               const std::filesystem::path &dir,
                               const ChunkedWaitQueue<T_REAL> &kmers) {
    using T_INT_REAL = get_int_t<T_REAL>;

    const uint8_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    size_t chunk_count = std::pow(alphabet_size, 2);

    logger->trace("Splitting k-mers into {} chunks...", chunk_count);

    std::vector<Encoder<T_INT_REAL>> sinks;
    std::vector<std::string> names(chunk_count);
    for (size_t i = 0; i < names.size(); ++i) {
        names[i] = dir/("real_F_W_" + std::to_string(i));
        sinks.emplace_back(names[i], ENCODER_BUFFER_SIZE);
    }

    size_t num_parent_kmers = 0;
    for (auto &it = kmers.begin(); it != kmers.end(); ++it) {
        const T_REAL &kmer = *it;
        TAlphabet F = get_first(kmer)[k];
        TAlphabet W = get_first(kmer)[0];
        size_t idx = F * alphabet_size + W;
        sinks[idx].add(reinterpret_cast<const T_INT_REAL &>(kmer));
        num_parent_kmers++;
    }
    std::for_each(sinks.begin(), sinks.end(), [](auto &f) { f.finish(); });
    logger->trace("Total number of non-dummy k-mers: {}", num_parent_kmers);
    return names;
}

template <typename Decoder, typename KMER>
void skip_same_suffix(const KMER &el, Decoder &decoder, size_t suf) {
    while (!decoder.empty()) {
        KMER kmer(decoder.top());
        if (!KMER::compare_suffix(kmer, el, suf)) {
            break;
        }
        decoder.pop();
    }
}

// construct word `...1001001001` with k ones for lifting k-mers
template <typename WordType>
inline WordType get_sentinel_delta(size_t char_width, size_t k) {
    assert(char_width * k <= sizeof(WordType) * 8);
    WordType word = 0;
    for (size_t i = 0; i < k; ++i) {
        word <<= char_width;
        word |= 1;
    }
    return word;
}

// transforms k-mer to the new character width
template <typename KMER_TO, typename KMER_FROM>
inline __attribute__((always_inline))
typename KMER_TO::WordType transform(const KMER_FROM &kmer, size_t k) {
    static constexpr size_t L1 = KMER_FROM::kBitsPerChar;
    static constexpr size_t L2 = KMER_TO::kBitsPerChar;
    static_assert(L2 >= L1);
    static_assert(L2 <= L1 + 1);
    assert(sizeof(typename KMER_TO::WordType)
            >= sizeof(typename KMER_FROM::WordType));

    if constexpr(L1 == L2) {
        return kmer.data();

    } else {
        typename KMER_TO::WordType word = 0;

        static constexpr uint64_t char_mask = (1ull << L1) - 1;

        for (int pos = L1 * (k - 1); pos >= 0; pos -= L1) {
            word <<= L2;
            assert(kmer[pos / L1] + 1 <= sdsl::bits::lo_set[L2]);
            word |= static_cast<uint64_t>(kmer.data() >> pos) & char_mask;
        }

        return word;
    }
}

static const uint16_t lookup_2_3[256] = {
    0,    1,    2,    3,    8,    9,    10,   11,   16,   17,   18,   19,   24,
    25,   26,   27,   64,   65,   66,   67,   72,   73,   74,   75,   80,   81,
    82,   83,   88,   89,   90,   91,   128,  129,  130,  131,  136,  137,  138,
    139,  144,  145,  146,  147,  152,  153,  154,  155,  192,  193,  194,  195,
    200,  201,  202,  203,  208,  209,  210,  211,  216,  217,  218,  219,  512,
    513,  514,  515,  520,  521,  522,  523,  528,  529,  530,  531,  536,  537,
    538,  539,  576,  577,  578,  579,  584,  585,  586,  587,  592,  593,  594,
    595,  600,  601,  602,  603,  640,  641,  642,  643,  648,  649,  650,  651,
    656,  657,  658,  659,  664,  665,  666,  667,  704,  705,  706,  707,  712,
    713,  714,  715,  720,  721,  722,  723,  728,  729,  730,  731,  1024, 1025,
    1026, 1027, 1032, 1033, 1034, 1035, 1040, 1041, 1042, 1043, 1048, 1049, 1050,
    1051, 1088, 1089, 1090, 1091, 1096, 1097, 1098, 1099, 1104, 1105, 1106, 1107,
    1112, 1113, 1114, 1115, 1152, 1153, 1154, 1155, 1160, 1161, 1162, 1163, 1168,
    1169, 1170, 1171, 1176, 1177, 1178, 1179, 1216, 1217, 1218, 1219, 1224, 1225,
    1226, 1227, 1232, 1233, 1234, 1235, 1240, 1241, 1242, 1243, 1536, 1537, 1538,
    1539, 1544, 1545, 1546, 1547, 1552, 1553, 1554, 1555, 1560, 1561, 1562, 1563,
    1600, 1601, 1602, 1603, 1608, 1609, 1610, 1611, 1616, 1617, 1618, 1619, 1624,
    1625, 1626, 1627, 1664, 1665, 1666, 1667, 1672, 1673, 1674, 1675, 1680, 1681,
    1682, 1683, 1688, 1689, 1690, 1691, 1728, 1729, 1730, 1731, 1736, 1737, 1738,
    1739, 1744, 1745, 1746, 1747, 1752, 1753, 1754, 1755
};

template <>
inline __attribute__((always_inline)) sdsl::uint128_t
transform<kmer::KMerBOSS<sdsl::uint128_t, 3>, kmer::KMerBOSS<uint64_t, 2>>(
        const kmer::KMerBOSS<uint64_t, 2> &kmer, size_t /*k*/) {
    const uint8_t *kmer_char = reinterpret_cast<const uint8_t *>(&kmer);
    // transform 64-bit kmer to 128 bits
    return static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[7]]) << 84
            | static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[6]]) << 72
            | static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[5]]) << 60
            | static_cast<uint64_t>(lookup_2_3[kmer_char[4]]) << 48
            | static_cast<uint64_t>(lookup_2_3[kmer_char[3]]) << 36
            | static_cast<uint64_t>(lookup_2_3[kmer_char[2]]) << 24
            | static_cast<uint64_t>(lookup_2_3[kmer_char[1]]) << 12
            | lookup_2_3[kmer_char[0]];
}

template <>
inline __attribute__((always_inline)) sdsl::uint256_t
transform<kmer::KMerBOSS<sdsl::uint256_t, 3>, kmer::KMerBOSS<sdsl::uint128_t, 2>>(
        const kmer::KMerBOSS<sdsl::uint128_t, 2> &kmer, size_t /*k*/) {
    const uint8_t *kmer_char = reinterpret_cast<const uint8_t *>(&kmer);
    // transform 128-bit kmer to 256 bits
    return sdsl::uint256_t(0, static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[15]]) << 52)
            | sdsl::uint256_t(static_cast<uint64_t>(lookup_2_3[kmer_char[14]]) << 48
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[13]]) << 36
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[12]]) << 24
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[11]]) << 12
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[10]])
                            ) << 120
            | sdsl::uint128_t(static_cast<uint64_t>(lookup_2_3[kmer_char[9]]) << 48
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[8]]) << 36
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[7]]) << 24
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[6]]) << 12
                                | static_cast<uint64_t>(lookup_2_3[kmer_char[5]])
                            ) << 60
            | (static_cast<uint64_t>(lookup_2_3[kmer_char[4]]) << 48
                        | static_cast<uint64_t>(lookup_2_3[kmer_char[3]]) << 36
                        | static_cast<uint64_t>(lookup_2_3[kmer_char[2]]) << 24
                        | static_cast<uint64_t>(lookup_2_3[kmer_char[1]]) << 12
                        | lookup_2_3[kmer_char[0]]);
}

// shift to the next dummy sink and add +1 to each character of the k-mer
template <typename KMER_TO, typename KMER_FROM>
inline typename KMER_TO::WordType get_sink_and_lift(const KMER_FROM &kmer, size_t k) {
    static constexpr int L1 = KMER_FROM::kBitsPerChar;
    static constexpr int L2 = KMER_TO::kBitsPerChar;
    static_assert(L2 >= L1);
    static_assert(L2 <= L1 + 1);
    assert(sizeof(typename KMER_TO::WordType)
            >= sizeof(typename KMER_FROM::WordType));

    static constexpr uint64_t first_char_mask_1 = (1ull << L1) - 1;

    typename KMER_TO::WordType word = (kmer.data() & first_char_mask_1) + 1;

    for (int pos = L1 * (k - 1); pos >= L1 * 2; pos -= L1) {
        word <<= L2;
        assert(kmer[pos / L1] + 1 <= sdsl::bits::lo_set[L2]);
        word |= (static_cast<uint64_t>(kmer.data() >> pos) & first_char_mask_1) + 1;
    }

    word <<= L2;

    return word;
}

/**
 * Generates non-redundant dummy-1 source k-mers and dummy sink kmers from #kmers.
 * @return a triplet containing the names of the original k-mer blocks, the dummy-1 source
 * k-mer blocks and the dummy sink k-mers
 */
template <typename T_REAL, typename T>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::string>
generate_dummy_1_kmers(size_t k,
                       size_t num_threads,
                       const std::filesystem::path &dir,
                       ChunkedWaitQueue<T_REAL> &kmers) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerExtractorBOSS::KmerBOSS
    using KMER_INT = typename KMER::WordType; // KmerExtractorBOSS::KmerBOSS::WordType

    using KMER_REAL = get_first_type_t<T_REAL>; // KmerExtractorT::KmerBOSS without sentinel
    using KMER_INT_REAL = typename KMER_REAL::WordType; // KmerExtractorT::KmerBOSS::WordType

    // for a DNA alphabet, this will contain 16 chunks, split by kmer[0] and kmer[1]
    std::vector<std::string> real_F_W = split(k, dir, kmers);

    const uint8_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    std::vector<Encoder<KMER_INT>> dummy_l1_chunks;
    std::vector<Encoder<KMER_INT>> dummy_sink_chunks;
    std::vector<std::string> dummy_l1_names(alphabet_size);
    std::vector<std::string> dummy_sink_names(alphabet_size);
    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        dummy_l1_names[i] = dir/("dummy_source_1_" + std::to_string(i));
        dummy_sink_names[i] = dir/("dummy_sink_" + std::to_string(i));
        dummy_l1_chunks.emplace_back(dummy_l1_names[i], ENCODER_BUFFER_SIZE);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i], ENCODER_BUFFER_SIZE);
    }

    logger->trace("Generating dummy-1 source k-mers and dummy sink k-mers...");
    uint64_t num_sink = 0;
    uint64_t num_source = 0;

    static constexpr size_t L = KMER::kBitsPerChar;
    KMER_INT kmer_delta = get_sentinel_delta<KMER_INT>(L, k + 1);
    // reset kmer[1] (the first character in k-mer, $ in dummy source) to zero
    kmer_delta &= ~KMER_INT(((1ull << L) - 1) << L);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (TAlphabet F = 0; F < alphabet_size; ++F) {

        // stream k-mers of pattern ***F*
        std::vector<std::string> F_chunks(real_F_W.begin() + F * alphabet_size,
                                          real_F_W.begin() + (F + 1) * alphabet_size);
        common::MergeDecoder<KMER_INT_REAL> it(F_chunks, false);

        std::vector<std::string> W_chunks;  // chunks with k-mers of the form ****F
        for (TAlphabet c = 0; c < alphabet_size; ++c) {
            W_chunks.push_back(real_F_W[c * alphabet_size + F]);
        }
        common::ConcatDecoder<KMER_INT_REAL> sink_gen_it(W_chunks);

        while (!it.empty()) {
            KMER_REAL dummy_source(it.pop());
            // skip k-mers that would generate identical source dummy k-mers
            skip_same_suffix(dummy_source, it, 0);
            dummy_source.to_prev(k + 1, 0);
            // generate dummy sink k-mers from all non-dummy kmers smaller than |dummy_source|
            while (!sink_gen_it.empty()
                    && sink_gen_it.top() < dummy_source.data()) {
                KMER_REAL v(sink_gen_it.pop());
                // skip k-mers with the same suffix as v, as they generate identical dummy
                // sink k-mers
                skip_same_suffix(v, sink_gen_it, 1);
                dummy_sink_chunks[F].add(get_sink_and_lift<KMER>(v, k + 1));
                num_sink++;
            }
            if (!sink_gen_it.empty()) {
                KMER_REAL top(sink_gen_it.top());
                if (KMER_REAL::compare_suffix(top, dummy_source, 1)) {
                    // The source dummy k-mer #dummy_source generated from #it is
                    // redundant iff it shares its suffix with another real k-mer (#top).
                    // In this case, #top generates a dummy sink k-mer redundant with #it.
                    // So if #dummy_source is redundant, the sink generated from #top is
                    // also redundant - so it's being skipped
                    skip_same_suffix(top, sink_gen_it, 1);
                    continue;
                }
            }
            // lift all and reset the first character to the sentinel 0 (apply mask)
            dummy_l1_chunks[F].add(transform<KMER>(dummy_source, k + 1) + kmer_delta);
            num_source++;
        }
        // handle leftover sink_gen_it
        while (!sink_gen_it.empty()) {
            KMER_REAL v(sink_gen_it.pop());
            skip_same_suffix(v, sink_gen_it, 1);
            dummy_sink_chunks[F].add(get_sink_and_lift<KMER>(v, k + 1));
            num_sink++;
        }
    }

    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        dummy_sink_chunks[i].finish();
        dummy_l1_chunks[i].finish();
    }

    logger->trace("Generated {} dummy sink and {} dummy source k-mers",
                  num_sink, num_source);

    // dummy sink k-mers are partitioned into blocks by F (kmer[1]), so simply
    // concatenating the blocks will result in a single ordered block
    logger->trace("Concatenating blocks of dummy sink k-mers ({} -> 1)...",
                  dummy_sink_names.size());
    std::string dummy_sink_name = dir/"dummy_sink";
    common::concat(dummy_sink_names, dummy_sink_name);

    // similarly, the 16 blocks of the original k-mers can be concatenated in groups of
    // 4 without destroying the order
    logger->trace("Concatenating blocks of original real k-mers ({} -> {})...",
                  real_F_W.size(), alphabet_size);
    std::vector<std::string> real_split_by_W;
    for (TAlphabet W = 0; W < alphabet_size; ++W) {
        std::vector<std::string> blocks;
        for (TAlphabet F = 0; F < alphabet_size; ++F) {
            blocks.push_back(real_F_W[F * alphabet_size + W]);
        }
        real_split_by_W.push_back(dir/("real_split_by_W_" + std::to_string(W)));
        common::concat(blocks, real_split_by_W.back());
    }
    return { real_split_by_W, dummy_l1_names, dummy_sink_name };
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
template <class KmerCollector, typename T_REAL, typename T>
void recover_dummy_nodes_disk(const KmerCollector &kmer_collector,
                              ChunkedWaitQueue<T_REAL> &kmers,
                              ChunkedWaitQueue<T> *kmers_out,
                              ThreadPool &async_worker) {
    using KMER_REAL = get_first_type_t<T_REAL>; // 64/128/256-bit KmerBOSS
    using T_INT_REAL = get_int_t<T_REAL>; // either KMER_INT or <KMER_INT, count>

    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS with sentinel $
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer

    const uint8_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    size_t k = kmer_collector.get_k() - 1;
    const std::filesystem::path dir = kmer_collector.tmp_dir();

    std::string dummy_sink_name;
    std::vector<std::string> real_split_by_W;
    std::vector<std::string> dummy_names;
    std::tie(real_split_by_W, dummy_names, dummy_sink_name)
            = generate_dummy_1_kmers<T_REAL, T>(k, kmer_collector.num_threads(), dir, kmers);

    // stores the sorted original kmers and dummy-1 k-mers
    std::vector<std::string> dummy_chunks = { dummy_sink_name };
    // generate dummy k-mers of prefix length 1..k
    logger->trace("Starting generating dummy-1..k source k-mers...");
    for (size_t dummy_pref_len = 1; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will compress all sorted dummy k-mers of given prefix length
        for (const std::string &f : dummy_names) {
            dummy_chunks.push_back(f);
        }

        std::vector<std::string> dummy_next_names(alphabet_size);
        std::vector<Encoder<KMER_INT>> dummy_next_chunks;
        for (TAlphabet i = 0; i < alphabet_size; ++i) {
            dummy_next_names[i] = dir/("dummy_source_"
                    + std::to_string(dummy_pref_len + 1) + "_" + std::to_string(i));
            dummy_next_chunks.emplace_back(dummy_next_names[i], ENCODER_BUFFER_SIZE);
        }

        KMER prev_kmer(0);
        uint64_t num_kmers = 0;
        const std::function<void(const KMER_INT &)> &write_dummy = [&](const KMER_INT &v) {
            KMER kmer(v);
            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            if (prev_kmer != kmer) {
                dummy_next_chunks[kmer[0]].add(kmer.data());
                prev_kmer = std::move(kmer);
            }
            num_kmers++;
        };
        common::merge_files(dummy_names, write_dummy, false);

        std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                      [](auto &v) { v.finish(); });
        dummy_names = std::move(dummy_next_names);
        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}",
                      dummy_pref_len, num_kmers);
    }
    // remove the last chunks with .up and .count
    const std::function<void(const KMER_INT &)> on_merge = [](const KMER_INT &) {};
    common::merge_files(dummy_names, on_merge);

    // at this point, we have the original k-mers and dummy-1 k-mers in original_and_dummy_l1,
    // the dummy-x k-mers in dummy_source_{x}, and we merge them all into a single stream
    kmers_out->reset();

    // add the main dummy source k-mer
    if constexpr (utils::is_pair_v<T>) {
        kmers_out->push({KMER(0), 0});
    } else {
        kmers_out->push(KMER(0));
    }

    const KMER_INT kmer_delta
            = get_sentinel_delta<KMER_INT>(KMER::kBitsPerChar, k + 1);

    // push all other dummy and non-dummy k-mers to |kmers_out|
    async_worker.enqueue([k, kmer_delta, kmers_out, real_split_by_W, dummy_chunks]() {
        common::Transformed<common::MergeDecoder<T_INT_REAL>, T> decoder(
            [&](const T_INT_REAL &v) {
                if constexpr (utils::is_pair_v<T>) {
                    return T(transform<KMER>(reinterpret_cast<const KMER_REAL &>(v.first), k + 1)
                                + kmer_delta,
                             v.second);
                } else {
                    return transform<KMER>(reinterpret_cast<const KMER_REAL &>(v), k + 1) + kmer_delta;
                }
            },
            real_split_by_W, true
        );

        common::Transformed<common::MergeDecoder<KMER_INT>, T> decoder_dummy(
            [](const KMER_INT &v) {
                if constexpr (utils::is_pair_v<T>) {
                    return T(reinterpret_cast<const KMER &>(v), 0);
                } else {
                    return reinterpret_cast<const KMER &>(v);
                }
            },
            dummy_chunks, true
        );

        while (!decoder.empty() && !decoder_dummy.empty()) {
            if (get_first(decoder.top()) < get_first(decoder_dummy.top())) {
                kmers_out->push(decoder.pop());
            } else {
                kmers_out->push(decoder_dummy.pop());
            }
        }
        while (!decoder.empty()) {
            kmers_out->push(decoder.pop());
        }
        while (!decoder_dummy.empty()) {
            kmers_out->push(decoder_dummy.pop());
        }

        kmers_out->shutdown();
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
  public:
    BOSSChunkConstructor(size_t k,
                         bool canonical_mode,
                         uint8_t bits_per_count,
                         const std::string &filter_suffix,
                         size_t num_threads,
                         double memory_preallocated,
                         const std::filesystem::path &tmp_dir,
                         size_t max_disk_space)
        : kmer_collector_(k + 1,
                          canonical_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          tmp_dir,
                          max_disk_space),
          bits_per_count_(bits_per_count) {
        if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)
            && (!utils::is_instance_v<typename KmerCollector::Data, ChunkedWaitQueue>
                || filter_suffix.size())) {
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
        BOSS::Chunk *result;

        if constexpr(std::is_same_v<typename KmerCollector::Extractor,
                                    KmerExtractorBOSS>) {
            typename KmerCollector::Data &kmer_ints = kmer_collector_.data();

            using T_INT = typename KmerCollector::Data::value_type;
            using T = get_kmer_t<typename KmerCollector::Kmer, T_INT>;
            auto &kmers = reinterpret_container<T>(kmer_ints);

            if (!kmer_collector_.suffix_length()) {
                logger->trace("Reconstructing all required dummy source k-mers...");
                Timer timer;
                if constexpr(std::is_same_v<typename KmerCollector::Data,
                                            ChunkedWaitQueue<T_INT>>) {
                    recover_dummy_nodes_disk(kmer_collector_, kmers, &kmers, async_worker_);
                } else {
                    // kmer_collector stores (BOSS::k_ + 1)-mers
                    static_assert(std::is_same_v<typename KmerCollector::Data, Vector<T_INT>>);
                    recover_dummy_nodes(kmer_collector_.get_k() - 1,
                                        kmer_collector_.num_threads(), &kmers);
                }
                logger->trace("Dummy source k-mers were reconstructed in {} sec",
                              timer.elapsed());
            }

            // kmer_collector stores (BOSS::k_ + 1)-mers
            result = new BOSS::Chunk(kmer_collector_.alphabet_size(),
                                     kmer_collector_.get_k() - 1,
                                     kmer_collector_.is_both_strands_mode(),
                                     kmers,
                                     bits_per_count_);
        } else {
            static_assert(std::is_same_v<typename KmerCollector::Extractor,
                                         KmerExtractor2Bit>);
            assert(!kmer_collector_.suffix_length());

            typename KmerCollector::Data &kmer_ints = kmer_collector_.data();

            using T_INT = typename KmerCollector::Data::value_type;
            using T = get_kmer_t<typename KmerCollector::Kmer, T_INT>;
            auto &kmers = reinterpret_container<T>(kmer_ints);

            logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;
            if constexpr(std::is_same_v<typename KmerCollector::Data,
                                        ChunkedWaitQueue<T_INT>>) {
#define INIT_CHUNK(KMER) \
    ChunkedWaitQueue<utils::replace_first_t<KMER, T>> queue(ENCODER_BUFFER_SIZE); \
    recover_dummy_nodes_disk(kmer_collector_, kmers, &queue, async_worker_); \
    logger->trace("Dummy source k-mers were reconstructed in {} sec", timer.elapsed()); \
    result = new BOSS::Chunk(KmerExtractorBOSS().alphabet.size(), \
                             kmer_collector_.get_k() - 1, \
                             kmer_collector_.is_both_strands_mode(), \
                             queue, \
                             bits_per_count_)

                if (kmer_collector_.get_k() * KmerExtractorBOSS::bits_per_char <= 64) {
                    INIT_CHUNK(KmerExtractorBOSS::Kmer64);
                } else if (kmer_collector_.get_k() * KmerExtractorBOSS::bits_per_char <= 128) {
                    INIT_CHUNK(KmerExtractorBOSS::Kmer128);
                } else {
                    INIT_CHUNK(KmerExtractorBOSS::Kmer256);
                }
            } else {
                throw std::runtime_error("Not implemented");
            }
        }

        kmer_collector_.clear();

        return result;
    }

    uint64_t get_k() const { return kmer_collector_.get_k() - 1; }

  private:
    KmerCollector kmer_collector_;
    uint8_t bits_per_count_;
    /** Used as an async executor for merging chunks from disk */
    ThreadPool async_worker_ = ThreadPool(1, 1);
};

template <template <typename, class> class KmerContainer, typename... Args>
static std::unique_ptr<IBOSSChunkConstructor>
initialize_boss_chunk_constructor(size_t k,
                                  bool canonical_mode,
                                  uint8_t bits_per_count,
                                  const std::string &filter_suffix,
                                  const Args& ...args) {
    if (k < 1 || k > 256 / KmerExtractorBOSS::bits_per_char - 1) {
        logger->error("For succinct graph, k must be between 2 and {}",
                      256 / KmerExtractorBOSS::bits_per_char - 1);
        exit(1);
    }

    // collect real k-mers in tight layout only with sorted set disk and if suffix is empty
    if (!filter_suffix.size()
            && utils::is_instance_v<typename KmerContainer<KmerExtractorBOSS::Kmer64,
                                                           KmerExtractorBOSS>::Data,
                                    ChunkedWaitQueue>) {
        if ((k + 1) * KmerExtractor2Bit::bits_per_char <= 64) {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractor2Bit::KmerBOSS64,
                                                                       KmerExtractor2Bit>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);

        } else if ((k + 1) * KmerExtractor2Bit::bits_per_char <= 128) {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractor2Bit::KmerBOSS128,
                                                                       KmerExtractor2Bit>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);

        } else {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractor2Bit::KmerBOSS256,
                                                                       KmerExtractor2Bit>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);
        }
    } else {
        if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 64) {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer64,
                                                                       KmerExtractorBOSS>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);

        } else if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 128) {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer128,
                                                                       KmerExtractorBOSS>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);

        } else {
            return std::make_unique<BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer256,
                                                                       KmerExtractorBOSS>>>(
                            k, canonical_mode, bits_per_count, filter_suffix, args...);
        }
    }
}

template <typename KMER, class KMER_EXTRACTOR>
using KmerSetVector
    = KmerCollector<KMER, KMER_EXTRACTOR, common::SortedSet<typename KMER::WordType>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetVector8
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedMultiset<typename KMER::WordType, uint8_t>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetVector16
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedMultiset<typename KMER::WordType, uint16_t>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetVector32
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedMultiset<typename KMER::WordType, uint32_t>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerSetDisk
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedSetDisk<typename KMER::WordType>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetDiskVector8
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedMultisetDisk<typename KMER::WordType, uint8_t>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetDiskVector16
    = KmerCollector<KMER, KMER_EXTRACTOR,
                    common::SortedMultisetDisk<typename KMER::WordType, uint16_t>>;

template <typename KMER, class KMER_EXTRACTOR>
using KmerMultsetDiskVector32
    = KmerCollector<KMER, KMER_EXTRACTOR,
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
} // namespace mtg
