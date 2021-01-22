#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "common/elias_fano_file_merger.hpp"
#include "common/logger.hpp"
#include "common/sorted_sets/sorted_multiset.hpp"
#include "common/sorted_sets/sorted_multiset_disk.hpp"
#include "common/sorted_sets/sorted_set.hpp"
#include "common/sorted_sets/sorted_set_disk.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "kmer/kmer_collector.hpp"
#include "kmer/kmer_to_int_converter.hpp"
#include "kmer/kmer_transform.hpp"
#include "boss_chunk.hpp"


namespace mtg {
namespace graph {
namespace boss {

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

const size_t ENCODER_BUFFER_SIZE = 100'000;

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
 * @tparam T_REAL the type of kmers being processed, typically either a KMer64/128/256 or
 * an std::pair<KMer64/128/256, int8/16/32> if counting kmers.
 * @tparam KMER "tight" KMer64/128/256 over the restricted alphabet (no sentinel)
 * @param k node length (the actual k-mer has length k+1)
 * @param kmers_p list of sorted, non-redundant kmers of length k+1
 */
template <typename KMER, typename T_REAL>
void add_dummy_sink_kmers(size_t k, const Vector<T_REAL> &kmers, Vector<KMER> *dummy_kmers) {
    using KMER_REAL = get_first_type_t<T_REAL>;
    using KMER_INT = typename KMER::WordType;

    KMER_INT kmer_delta = kmer::get_sentinel_delta<KMER_INT>(KMER::kBitsPerChar, k + 1);
    // reset kmer[0] (the last character in k-mer, $ in dummy sink) to zero
    kmer_delta &= ~KMER_INT(KMER::kFirstCharMask);

    const size_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    // points to the current k-mer with the given first character
    std::vector<size_t> max_it(alphabet_size);
    std::vector<size_t> it(alphabet_size);
    for (TAlphabet c = 1; c < alphabet_size; ++c) {
        std::vector<TAlphabet> zeros(k + 1, 0);
        zeros[k - 1] = c;
        it[c] = std::lower_bound(kmers.data(), kmers.data() + kmers.size(),
                                 KMER_REAL(zeros), // the AA...c->A k-mer
                                 [](const T_REAL &a, const KMER_REAL &b) -> bool {
                                     return get_first(a) < b;
                                 }) - kmers.data();
        max_it[c - 1] = it[c];
    }
    max_it[alphabet_size - 1] = kmers.size();

    std::vector<KMER_REAL> last_shifted(alphabet_size, KMER_REAL(1));
    size_t size = kmers.size();
    for (size_t i = 0; i < size; ++i) {
        TAlphabet last_char = get_first(kmers[i])[0];
        KMER_REAL shifted = get_first(kmers[i]);
        shifted.to_next(k + 1, 0);

        while (it[last_char] < max_it[last_char]
                && KMER_REAL::less(get_first(kmers[it[last_char]]), shifted)) {
            it[last_char]++;
        }
        if (last_shifted[last_char] != shifted
            && (it[last_char] == max_it[last_char]
                || !KMER_REAL::compare_suffix(get_first(kmers[it[last_char]]), shifted))) {
            dummy_kmers->emplace_back(kmer::transform<KMER>(shifted, k + 1) + kmer_delta);
            last_shifted[last_char] = shifted;
        }
    }
}

/**
 * Adds non-redundant dummy source nodes with prefix length 1 for the given kmers.
 *
 * For each kmer in #kmers, the method creates the corresponding dummy sources with
 * sentinels of length 1 (aka dummy-1 sources). The method then checks if the dummy-1
 * source node is redundant, i.e. if there is another k-mer that is identical
 * except for the first character. For example, $ACGT is redundant with TACGT. To do this
 * efficiently, the method uses an iterator called dummy_it that points to the first kmer
 * with the suffix equal to or larger than the current kmer. Since #kmers are ordered,
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
 *
 * Dummy-1 kmers are transformed to the alphabet that includes the sentinel ($) before
 * being added to #dummy_kmers.
 */
template <typename T_REAL, typename KMER>
void add_dummy_source_kmers(size_t k,
                            const Vector<T_REAL> &kmers,
                            Vector<KMER> *dummy_kmers) {
    using KMER_REAL = get_first_type_t<T_REAL>;
    using KMER_INT = typename KMER::WordType;

    KMER_INT kmer_delta = kmer::get_sentinel_delta<KMER_INT>(KMER::kBitsPerChar, k + 1);
    // reset kmer[1] (the first character in k-mer, $ in dummy source) to zero
    kmer_delta &= ~KMER_INT(KMER::kFirstCharMask << KMER::kBitsPerChar);

    // points to the first k-mer that may be redundant with the current dummy source k-mer
    size_t real_it = 0;

    for (size_t i = 0; i < kmers.size(); ++i) {
        const KMER_REAL &kmer = get_first(kmers[i]);

        if (i > 0 && kmer[k] != get_first(kmers[i - 1])[k]) {
            // the last (most significant) character changed -> reset the iterator
            real_it = 0;
        }

        if (i > 0 && KMER_REAL::compare_suffix(kmer, get_first(kmers[i - 1])))
            continue; // i would generate the same dummy-1 k-mer as i-1, skip

        KMER_REAL prev_kmer = kmer;
        prev_kmer.to_prev(k + 1, 0);

        while (real_it < kmers.size()
                && KMER_REAL::less(get_first(kmers[real_it]), prev_kmer, 1)) {
            real_it++;
        }
        bool is_redundant = false;
        while (real_it < kmers.size()
                && KMER_REAL::compare_suffix(get_first(kmers[real_it]), prev_kmer, 1)) {
            if (get_first(kmers[real_it])[0] == prev_kmer[0]) {  // edge labels match
                is_redundant = true;
                break;
            }
            real_it++;
        }
        if (!is_redundant) {
            dummy_kmers->emplace_back(transform<KMER>(prev_kmer, k + 1) + kmer_delta);
        }
    }
}

template <typename T>
inline T rev_comp(size_t k, T kmer, const std::vector<TAlphabet> &complement_code) {
    if constexpr(utils::is_pair_v<T>) {
        return T(kmer::reverse_complement(k, kmer.first, complement_code), kmer.second);
    } else {
        return kmer::reverse_complement(k, kmer, complement_code);
    }
}

template <typename T>
void add_reverse_complements(size_t k, size_t num_threads, Vector<T> *kmers) {
    size_t size = kmers->size();
    kmers->reserve(2 * size);

    logger->trace("Adding reverse-complement k-mers...");
    const std::vector<TAlphabet> complement_code = KmerExtractor2Bit().complement_code();
    const uint32_t BUF_SIZE = std::min(kmers->size() / num_threads, 10000UL);
    #pragma omp parallel num_threads(num_threads)
    {
        std::vector<T> buffer;
        buffer.reserve(BUF_SIZE);
        #pragma omp for schedule(static)
        for (T *kmer = kmers->data(); kmer < kmers->data() + size; ++kmer) {
            const T &rc = rev_comp(k + 1, *kmer, complement_code);
            if (get_first(rc) != get_first(*kmer)) {
                if (buffer.size() == buffer.capacity()) {
                    #pragma omp critical
                    {
                        // this is guaranteed to not reallocate
                        kmers->insert(kmers->end(), buffer.begin(), buffer.end());
                    }
                    buffer.resize(0);
                }
                buffer.push_back(std::move(rc));
            } else {
                if constexpr(utils::is_pair_v<T>) {
                    using C = typename T::second_type;
                    if (kmer->second >> (sizeof(C) * 8 - 1)) {
                        kmer->second = std::numeric_limits<C>::max();
                    } else {
                        kmer->second *= 2;
                    }
                }
            }
        }
        #pragma omp critical
        {
            kmers->insert(kmers->end(), buffer.begin(), buffer.end());
        }
    }
    logger->trace("Sorting all real kmers...");
    ips4o::parallel::sort(kmers->begin(), kmers->end(), utils::LessFirst(), num_threads);
}


template <typename T_TO, typename T_FROM>
ChunkedWaitQueue<T_TO>& reinterpret_container(ChunkedWaitQueue<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    return reinterpret_cast<ChunkedWaitQueue<T_TO>&>(container);
}

template <typename T_TO, typename T_FROM>
Vector<T_TO>& reinterpret_container(Vector<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    static_assert(!std::is_same_v<T_FROM, bool> && !std::is_same_v<T_TO, bool>);
    return reinterpret_cast<Vector<T_TO>&>(container);
}


/**
 * Adds dummy nodes for the given kmers and then builds a BOSS table (chunk).
 * The method first adds dummy sink kmers, then dummy sources with sentinels of length 1
 * (aka dummy-1 sources). The method will then gradually add dummy sources with sentinels
 * of length 2, 3, ... up to k-1.
 *
 * @tparam T a k-mer over the same alphabet as T_REAL plus the $ sentinel character (may
 * require one extra bit of storage per character)
 */
template <typename T, class KmerCollector>
BOSS::Chunk construct_boss_chunk(KmerCollector &kmer_collector,
                                 uint8_t bits_per_count,
                                 std::filesystem::path swap_dir) {
    using T_INT_REAL = typename KmerCollector::Data::value_type; // either KMER_INT or <KMER_INT, count>
    using T_REAL = get_kmer_t<typename KmerCollector::Kmer, T_INT_REAL>;

    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS with sentinel $
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer

    auto &kmers = reinterpret_container<T_REAL>(kmer_collector.data());

    size_t k = kmer_collector.get_k() - 1;
    size_t num_threads = kmer_collector.num_threads();
    bool both_strands_mode = kmer_collector.is_both_strands_mode();

    if (both_strands_mode) {
        add_reverse_complements(k, num_threads, &kmers);
    }

    logger->trace("Total number of real k-mers: {}", kmers.size());

    Vector<KMER> dummy_kmers;
    dummy_kmers.reserve(0.1 * kmers.size()); // assume at most 10% dummy sink k-mers
    add_dummy_sink_kmers(k, kmers, &dummy_kmers);
    // the total number of dummy source k-mers can be relatively accurately predicted
    // using the dummy sink kmers
    dummy_kmers.reserve(std::max(static_cast<int>(k) - 8, 2) * dummy_kmers.size());
    logger->trace("Added {} dummy sink k-mers", dummy_kmers.size());

    size_t dummy_source_begin = dummy_kmers.size();
    add_dummy_source_kmers(k, kmers, &dummy_kmers);

    ips4o::parallel::sort(dummy_kmers.begin() + dummy_source_begin, dummy_kmers.end(),
                          utils::LessFirst(), num_threads);

    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  dummy_kmers.size() - dummy_source_begin);

    for (size_t c = 2; c < k + 1; ++c) {
        size_t dummy_end = dummy_kmers.size();

        for (size_t i = dummy_source_begin; i < dummy_end; ++i) {
            KMER kmer = dummy_kmers[i];
            if (i > 0 && KMER::compare_suffix(kmer, dummy_kmers[i - 1]))
                continue; // i would generate the same dummy source k-mer as i-1, skip

            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            dummy_kmers.push_back(kmer);
        }
        dummy_source_begin = dummy_end;
        ips4o::parallel::sort(dummy_kmers.begin() + dummy_source_begin, dummy_kmers.end(),
                              utils::LessFirst(), num_threads);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      dummy_kmers.size() - dummy_source_begin);
    }

    ips4o::parallel::sort(dummy_kmers.begin(), dummy_kmers.end(),
                          utils::LessFirst(), num_threads);

    auto kmers_out = std::make_shared<ChunkedWaitQueue<T>>(ENCODER_BUFFER_SIZE);

    // add the main dummy source k-mer
    if constexpr(utils::is_pair_v<T>) {
        kmers_out->push({ KMER(0), 0 });
    } else {
        kmers_out->push(KMER(0));
    }

    ThreadPool async_worker(1, 1);
    async_worker.enqueue([k, &kmers, dummy_kmers{std::move(dummy_kmers)}, kmers_out]() {
        const KMER_INT kmer_delta
                = kmer::get_sentinel_delta<KMER_INT>(KMER::kBitsPerChar, k + 1);

        // merge #kmers and #dummy_kmers into #kmers_out
        size_t di = 0;
        for (const auto &kmer : kmers) {
            KMER lifted(transform<KMER>(get_first(kmer), k+1) + kmer_delta);
            if constexpr(utils::is_pair_v<T>) {
                while (di < dummy_kmers.size() && lifted > dummy_kmers[di]) {
                    kmers_out->push({ dummy_kmers[di], 0 });
                    di++;
                }
                kmers_out->push({ lifted, kmer.second });
            } else {
                while (di < dummy_kmers.size() && lifted > dummy_kmers[di]) {
                    kmers_out->push(dummy_kmers[di]);
                    di++;
                }
                kmers_out->push(lifted);
            }
        }
        for (size_t i = di; i < dummy_kmers.size(); ++i) {
            if constexpr(utils::is_pair_v<T>) {
                kmers_out->push({ dummy_kmers[i], 0 });
            } else {
                kmers_out->push(dummy_kmers[i]);
            }
        }
        kmers_out->shutdown();
    });

    BOSS::Chunk result(KmerExtractorBOSS().alphabet.size(),
                       k, both_strands_mode, *kmers_out,
                       bits_per_count, swap_dir);

    kmer_collector.clear();

    return result;
}

template <typename T>
using Encoder = common::EliasFanoEncoderBuffered<T>;
template <typename T>
using Decoder = common::EliasFanoDecoder<T>;

/**
 * Split k-mers from chunk |chunk_fname| into Sigma^2 chunks by F and W.
 */
template <typename T_REAL>
std::vector<std::string> split(size_t k, const std::string &chunk_fname) {
    using T_INT_REAL = get_int_t<T_REAL>;

    const uint8_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    size_t chunk_count = std::pow(alphabet_size, 2);

    logger->trace("Splitting k-mers from chunk {} into {} chunks...",
                  chunk_fname, chunk_count);

    std::vector<Encoder<T_INT_REAL>> sinks;
    std::vector<std::string> names(chunk_count);
    for (size_t i = 0; i < names.size(); ++i) {
        names[i] = chunk_fname + "_split_" + std::to_string(i);
        sinks.emplace_back(names[i], ENCODER_BUFFER_SIZE);
    }

    Decoder<T_INT_REAL> original_kmers(chunk_fname);
    size_t num_kmers = 0;
    while (std::optional<T_INT_REAL> orig = original_kmers.next()) {
        const T_REAL &kmer = reinterpret_cast<const T_REAL &>(orig.value());
        TAlphabet F = get_first(kmer)[k];
        TAlphabet W = get_first(kmer)[0];
        size_t idx = F * alphabet_size + W;
        sinks[idx].add(orig.value());
        num_kmers++;
    }
    std::for_each(sinks.begin(), sinks.end(), [](auto &f) { f.finish(); });
    logger->trace("Total number of real k-mers in {}: {}", chunk_fname, num_kmers);
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

/**
 * Generates ordered blocks of k-mers for the BOSS table.
 * Adds non-redundant dummy-1 source k-mers and dummy sink k-mers.
 * @return a triplet containing the names of the k-mer blocks.
 * 1) real k-mers (split by F),
 * 2) the dummy-1 source k-mers (split by F and W),
 * 3) the dummy sink k-mers (split by F).
 */
template <typename T_REAL, typename T>
std::tuple<std::vector<std::string>,
           std::vector<std::string>,
           std::vector<std::string>>
generate_dummy_1_kmers(size_t k,
                       size_t num_threads,
                       const std::filesystem::path &dir,
                       const std::vector<std::string> &real_F_W) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerExtractorBOSS::KmerBOSS
    using KMER_INT = typename KMER::WordType; // KmerExtractorBOSS::KmerBOSS::WordType

    using KMER_REAL = get_first_type_t<T_REAL>; // KmerExtractorT::KmerBOSS without sentinel
    using KMER_INT_REAL = typename KMER_REAL::WordType; // KmerExtractorT::KmerBOSS::WordType

    const uint8_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    using T_INT = get_int_t<T>;
    using T_INT_REAL = get_int_t<T_REAL>;
    std::vector<Encoder<T_INT>> real_chunks;
    std::vector<Encoder<KMER_INT>> dummy_sink_chunks;
    std::vector<std::string> real_names(alphabet_size);
    std::vector<std::string> dummy_sink_names(alphabet_size);
    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        real_names[i] = dir/("real_split_by_F_" + std::to_string(i));
        dummy_sink_names[i] = dir/("dummy_sink_" + std::to_string(i));
        real_chunks.emplace_back(real_names[i], ENCODER_BUFFER_SIZE);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i], ENCODER_BUFFER_SIZE);
    }
    std::vector<Encoder<KMER_INT>> dummy_l1_chunks;
    std::vector<std::string> dummy_l1_names;
    for (TAlphabet F = 0; F <= alphabet_size; ++F) {
        for (TAlphabet W = 0; W <= alphabet_size; ++W) {
            dummy_l1_names.push_back(dir/fmt::format("dummy_source_1_{}_{}", F, W));
            dummy_l1_chunks.emplace_back(dummy_l1_names.back(), ENCODER_BUFFER_SIZE);
        }
    }

    logger->trace("Generating dummy-1 source k-mers and dummy sink k-mers...");

    static constexpr size_t L = KMER::kBitsPerChar;
    KMER_INT kmer_delta = kmer::get_sentinel_delta<KMER_INT>(L, k + 1);
    // reset kmer[1] (the first character in k-mer, $ in dummy source) to zero
    KMER_INT kmer_delta_dummy = kmer_delta & ~KMER_INT(((1ull << L) - 1) << L);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (TAlphabet F = 0; F < alphabet_size; ++F) {
        // stream k-mers of pattern ***F*
        std::vector<std::string> F_chunks(real_F_W.begin() + F * alphabet_size,
                                          real_F_W.begin() + (F + 1) * alphabet_size);
        common::MergeDecoder<T_INT_REAL> it(F_chunks, false);
        std::vector<std::string> W_chunks;  // chunks with k-mers of the form ****F
        for (TAlphabet c = 0; c < alphabet_size; ++c) {
            W_chunks.push_back(real_F_W[c * alphabet_size + F]);
        }
        common::ConcatDecoder<KMER_INT_REAL> sink_gen_it(W_chunks);

        while (!it.empty()) {
            T_REAL v(it.pop());

            // write k-mer to its respective block |real_chunks[F]|
            if constexpr(utils::is_pair_v<T>) {
                real_chunks[F].add({ kmer::transform<KMER>(v.first, k + 1) + kmer_delta, v.second });
            } else {
                real_chunks[F].add(kmer::transform<KMER>(v, k + 1) + kmer_delta);
            }

            KMER_REAL dummy_source(utils::get_first(v));
            // skip k-mers that would generate identical source dummy k-mers
            if (!it.empty() && KMER_REAL::compare_suffix(KMER_REAL(utils::get_first(it.top())),
                                                         dummy_source, 0)) {
                continue;
            }

            dummy_source.to_prev(k + 1, 0);
            // generate dummy sink k-mers from all non-dummy kmers smaller than |dummy_source|
            while (!sink_gen_it.empty()
                    && sink_gen_it.top() < dummy_source.data()) {
                KMER_REAL v(sink_gen_it.pop());
                // skip k-mers with the same suffix as v, as they generate identical dummy
                // sink k-mers
                skip_same_suffix(v, sink_gen_it, 1);
                dummy_sink_chunks[F].add(kmer::get_sink_and_lift<KMER>(v, k + 1));
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
            assert(F == dummy_source[0]);
            TAlphabet next_F = dummy_source[k];
            // lift all and reset the first character to the sentinel 0 (apply mask)
            dummy_l1_chunks[(next_F+1) * (alphabet_size+1) + (F+1)].add(
                kmer::transform<KMER>(dummy_source, k + 1) + kmer_delta_dummy);
        }
        // handle leftover sink_gen_it
        while (!sink_gen_it.empty()) {
            KMER_REAL v(sink_gen_it.pop());
            skip_same_suffix(v, sink_gen_it, 1);
            dummy_sink_chunks[F].add(kmer::get_sink_and_lift<KMER>(v, k + 1));
        }
    }

    // remove all unmerged chunks
    common::remove_chunks(real_F_W);

    uint64_t num_sink = 0;
    uint64_t num_source = 0;
    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        real_chunks[i].finish();
        dummy_sink_chunks[i].finish();
        num_sink += dummy_sink_chunks[i].size();
    }
    for (auto &dummy_l1_chunk : dummy_l1_chunks) {
        dummy_l1_chunk.finish();
        num_source += dummy_l1_chunk.size();
    }

    logger->trace("Generated {} dummy sink and {} dummy source k-mers",
                  num_sink, num_source);

    return { real_names, dummy_l1_names, dummy_sink_names };
}

/**
 * Adds reverse complements
 */
template <typename T_REAL>
void add_reverse_complements(size_t k,
                             size_t num_threads,
                             size_t buffer_size,
                             const std::vector<std::string> &real_F_W) {
    using T_INT_REAL = get_int_t<T_REAL>; // either KMER_INT or <KMER_INT, count>

    const uint8_t alphabet_size = KmerExtractor2Bit().alphabet.size();

    assert(real_F_W.size() == std::pow(alphabet_size, 2));

    // initialize containers for each block of rc k-mers
    std::vector<std::unique_ptr<common::SortedSetDisk<T_INT_REAL>>> rc_set;
    for (size_t j = 0; j < real_F_W.size(); ++j) {
        std::string rc_dir = real_F_W[j] + "_rc";
        std::filesystem::create_directory(rc_dir);
        rc_set.emplace_back(
            new common::SortedSetDisk<T_INT_REAL>(
                    num_threads, buffer_size / real_F_W.size(),
                    rc_dir, std::numeric_limits<size_t>::max()));
    }

    logger->trace("Adding reverse complements...");

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < real_F_W.size(); ++j) {
        // read all k-mers from the block and keep it for the merge
        Decoder<T_INT_REAL> original_kmers(real_F_W[j], false);

        // initialize buffers for each output (F, W) chunk
        std::vector<Vector<T_INT_REAL>> buffers(real_F_W.size());
        for (auto &buf : buffers) {
            buf.reserve(10'000);
        }

        while (std::optional<T_INT_REAL> orig = original_kmers.next()) {
            const T_REAL &kmer = reinterpret_cast<const T_REAL &>(orig.value());
            const T_REAL &reverse
                = rev_comp(k + 1, kmer, KmerExtractor2Bit().complement_code());

            TAlphabet F = get_first(reverse)[k];
            TAlphabet W = get_first(reverse)[0];
            size_t idx = F * alphabet_size + W;
            auto &buffer = buffers[idx];

            buffer.push_back(reinterpret_cast<const T_INT_REAL &>(reverse));
            if (buffer.size() == buffer.capacity()) {
                rc_set[idx]->insert(buffer.begin(), buffer.end());
                buffer.resize(0);
            }
        }
        for (size_t idx = 0; idx < buffers.size(); ++idx) {
            rc_set[idx]->insert(buffers[idx].begin(), buffers[idx].end());
        }
    }

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < rc_set.size(); ++j) {
        // start merging original with reverse_complements kmers
        std::vector<std::string> chunks_to_merge = rc_set[j]->files_to_merge();
        chunks_to_merge.push_back(real_F_W[j]);

        Encoder<T_INT_REAL> out(real_F_W[j] + "_new", ENCODER_BUFFER_SIZE);
        const std::function<void(const T_INT_REAL &)> write
                = [&](const T_INT_REAL &v) { out.add(v); };
        common::merge_files(chunks_to_merge, write);
        out.finish();

        rc_set[j].reset();
        std::filesystem::remove_all(real_F_W[j] + "_rc");

        for (const auto &suffix : { "", ".up", ".count" }) {
            if (std::filesystem::exists(real_F_W[j] + "_new" + suffix))
                std::filesystem::rename(real_F_W[j] + "_new" + suffix,
                                        real_F_W[j] + suffix);
        }

        logger->trace("Augmented chunk {} with reverse complement k-mers."
                      " New size: {} bytes", real_F_W[j],
                      std::filesystem::file_size(real_F_W[j]));
    }
}

template <typename KMER>
std::vector<std::vector<std::string>>
reconstruct_dummy_source(const std::vector<std::string> &dummy_l1_names,
                         size_t k,
                         size_t num_threads,
                         std::filesystem::path dir);

template <typename T>
BOSS::Chunk build_boss(const std::vector<std::string> &real_names,
                       const std::vector<std::vector<std::string>> &dummy_source_names,
                       const std::vector<std::string> &dummy_sink_names,
                       std::filesystem::path dir,
                       size_t k,
                       bool both_strands_mode,
                       uint8_t bits_per_count,
                       std::filesystem::path swap_dir,
                       size_t num_threads);

/**
 * Reconstructs all the required dummy k-mers and builds a BOSS table (chunk).
 * The method first traverses the original kmers and generates dummy source k-mers of
 * prefix length 1 (dummy-1). At the next step, the non-redundant dummy-1 kmers are merged
 * with the original k-mers.
 * The method then gradually constructs dummy-i k-mers for i=2..k and writes them into
 * separate files, de-duped and sorted.
 * The final result is obtained by merging the original #kmers and dummy-1 kmers with
 * the dummy-k kmers, for k=2..k
 */
template <typename T, class KmerCollector>
BOSS::Chunk construct_boss_chunk_disk(KmerCollector &kmer_collector,
                                      uint8_t bits_per_count,
                                      std::filesystem::path swap_dir) {
    using T_INT_REAL = typename KmerCollector::Data::value_type; // either KMER_INT or <KMER_INT, count>
    using T_REAL = get_kmer_t<typename KmerCollector::Kmer, T_INT_REAL>;

    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS with sentinel $

    size_t k = kmer_collector.get_k() - 1;
    const std::filesystem::path dir = kmer_collector.tmp_dir();
    size_t num_threads = kmer_collector.num_threads();
    const uint64_t buffer_size = kmer_collector.buffer_size();
    bool both_strands_mode = kmer_collector.is_both_strands_mode();

    auto &container = kmer_collector.container();
    const std::vector<std::string> chunk_fnames = container.files_to_merge();

    // split each chunk by F and W
    std::vector<std::vector<std::string>> chunks_split(chunk_fnames.size());
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t i = 0; i < chunk_fnames.size(); ++i) {
        chunks_split[i] = split<T_REAL>(k, chunk_fnames[i]);
    }

    // release old chunks and the buffer
    container.clear();

    // for a DNA alphabet, this will contain 16 chunks, split by kmer[0] and kmer[1]
    std::vector<std::string> real_F_W(std::pow(KmerExtractor2Bit().alphabet.size(), 2));

    // merge each group of chunks into a single one for each F and W
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < real_F_W.size(); ++j) {
        std::vector<std::string> chunks_to_merge(chunks_split.size());
        for (size_t i = 0; i < chunks_split.size(); ++i) {
            chunks_to_merge[i] = chunks_split[i][j];
        }
        real_F_W[j] = dir/("real_F_W_" + std::to_string(j));
        Encoder<T_INT_REAL> out(real_F_W[j], ENCODER_BUFFER_SIZE);
        // merge chunks into a single one and remove them
        const std::function<void(const T_INT_REAL &)> write
                = [&](const T_INT_REAL &v) { out.add(v); };
        common::merge_files(chunks_to_merge, write, true);
        out.finish();
        logger->trace("Merged {} chunks into {} of size: {} bytes",
                      chunks_to_merge.size(), real_F_W[j],
                      std::filesystem::file_size(real_F_W[j]));
    }

    if (both_strands_mode) {
        // compute reverse complements k-mers and update the blocks #real_F_W
        add_reverse_complements<T_REAL>(k, num_threads, buffer_size, real_F_W);
    }

    // k-mer blocks split by F
    auto [real_names, dummy_l1_names, dummy_sink_names]
            = generate_dummy_1_kmers<T_REAL, T>(k, num_threads, dir, real_F_W);

    auto dummy_source_names
            = reconstruct_dummy_source<KMER>(dummy_l1_names, k, num_threads, dir);

    // merge all the original k-mers with the dummy k-mers and generate a BOSS table
    return build_boss<T>(real_names, dummy_source_names, dummy_sink_names,
                         dir, k, both_strands_mode, bits_per_count, swap_dir, num_threads);
}

template <typename KMER>
std::vector<std::vector<std::string>>
reconstruct_dummy_source(const std::vector<std::string> &dummy_l1_names,
                         size_t k,
                         size_t num_threads,
                         std::filesystem::path dir) {
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer
    const uint8_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    std::vector<std::string> names = dummy_l1_names;

    std::vector<std::vector<std::string>> result;

    // generate dummy k-mers of prefix length 1..k
    logger->trace("Starting generating dummy-1..k source k-mers...");
    for (size_t dummy_pref_len = 1; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will store blocks of source dummy k-mers of the current level
        uint64_t num_kmers = 0;
        result.emplace_back(alphabet_size);
        std::vector<std::string> next_names(alphabet_size * alphabet_size);
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (TAlphabet F = 0; F < alphabet_size; ++F) {
            // will store sorted dummy source k-mers of pattern ***F*
            result.back()[F] = dir/fmt::format("dummy_source_{}_{}", dummy_pref_len, F);
            Encoder<KMER_INT> dummy_chunk(result.back()[F], ENCODER_BUFFER_SIZE);

            std::vector<Encoder<KMER_INT>> dummy_next_chunks;
            for (TAlphabet i = 0; i < alphabet_size; ++i) {
                next_names[i * alphabet_size + F]
                    = dir/fmt::format("dummy_source_{}_{}_{}", dummy_pref_len + 1, i, F);
                dummy_next_chunks.emplace_back(next_names[i * alphabet_size + F],
                                               ENCODER_BUFFER_SIZE);
            }

            KMER prev_kmer(0);
            const std::function<void(const KMER_INT &)> &write_dummy = [&](const KMER_INT &v) {
                dummy_chunk.add(v);
                // ***F*
                KMER kmer(v);
                // $***F
                kmer.to_prev(k + 1, BOSS::kSentinelCode);
                if (prev_kmer != kmer) {
                    // $**XF
                    assert(kmer[0] == F);
                    dummy_next_chunks[kmer[k]].add(kmer.data());
                    prev_kmer = std::move(kmer);
                }
            };
            std::vector<std::string> F_chunk_names(names.begin() + F * alphabet_size,
                                                   names.begin() + (F + 1) * alphabet_size);
            common::merge_files(F_chunk_names, write_dummy, true);

            dummy_chunk.finish();
            std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                          [](auto &v) { v.finish(); });

            #pragma omp atomic
            num_kmers += dummy_chunk.size();
        }
        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}",
                      dummy_pref_len, num_kmers);
        names = std::move(next_names);
    }
    // remove the last chunks with .up and .count
    common::remove_chunks(names);

    return result;
}

template <typename T>
BOSS::Chunk build_boss_chunk(bool is_first_chunk,
                             const std::string &real_name,
                             const std::vector<std::string> &dummy_names,
                             size_t k,
                             bool both_strands_mode,
                             uint8_t bits_per_count,
                             std::filesystem::path swap_dir) {
    using T_INT = get_int_t<T>;
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS with sentinel $
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer

    auto kmers_out = std::make_shared<ChunkedWaitQueue<T_INT>>(ENCODER_BUFFER_SIZE);

    if (is_first_chunk) {
        // add the main dummy source k-mer
        if constexpr(utils::is_pair_v<T>) {
            kmers_out->push({ KMER_INT(0), 0 });
        } else {
            kmers_out->push(KMER_INT(0));
        }
    }

    // push all other dummy and non-dummy k-mers to |kmers_out|
    ThreadPool async_worker(1, 1);
    async_worker.enqueue([kmers_out, real_name, dummy_names]() {
        Decoder<T_INT> decoder(real_name, true);

        common::Transformed<common::MergeDecoder<KMER_INT>, T_INT> decoder_dummy(
            [](const KMER_INT &v) {
                if constexpr(utils::is_pair_v<T>) {
                    return std::make_pair(v, 0);
                } else {
                    return v;
                }
            },
            dummy_names, true
        );

        std::optional<T_INT> next = decoder.next();
        while (next && !decoder_dummy.empty()) {
            if (get_first(next.value()) < get_first(decoder_dummy.top())) {
                kmers_out->push(next.value());
                next = decoder.next();
            } else {
                kmers_out->push(decoder_dummy.pop());
            }
        }
        while (next) {
            kmers_out->push(next.value());
            next = decoder.next();
        }
        while (!decoder_dummy.empty()) {
            kmers_out->push(decoder_dummy.pop());
        }

        kmers_out->shutdown();
    });

    return BOSS::Chunk(KmerExtractorBOSS().alphabet.size(),
                       k, both_strands_mode,
                       reinterpret_container<T>(*kmers_out),
                       bits_per_count, swap_dir);
}

template <typename T>
BOSS::Chunk build_boss(const std::vector<std::string> &real_names,
                       const std::vector<std::vector<std::string>> &dummy_source_names,
                       const std::vector<std::string> &dummy_sink_names,
                       std::filesystem::path dir,
                       size_t k,
                       bool both_strands_mode,
                       uint8_t bits_per_count,
                       std::filesystem::path swap_dir,
                       size_t num_threads) {
    assert(real_names.size() == dummy_sink_names.size());
    assert(real_names.size() + 1 == dummy_source_names.front().size());
    assert(dummy_source_names.size() == k);

    if (k <= 1) {
        // For k <= 1 BOSS chunks can't be built independently.
        // Concatenate blocks and build a full BOSS table.
        logger->trace("Concatenating blocks of real k-mers ({} -> 1)...", real_names.size());
        std::string real_name = dir/"real_kmers";
        common::concat(real_names, real_name);

        logger->trace("Concatenating blocks of dummy sink k-mers ({} -> 1)...",
                      dummy_sink_names.size());
        std::vector<std::string> dummy_names { dir/"dummy_sink" };
        common::concat(dummy_sink_names, dummy_names.back());

        for (size_t i = 0; i < dummy_source_names.size(); ++i) {
            logger->trace("Concatenating blocks of dummy-{} source k-mers ({} -> 1)...",
                          i + 1, dummy_source_names[i].size());
            dummy_names.push_back(dir/fmt::format("dummy_source_{}", i));
            common::concat(dummy_source_names[i], dummy_names.back());
        }

        return build_boss_chunk<T>(true, real_name, dummy_names,
                                   k, both_strands_mode, bits_per_count, swap_dir);
    }

    // For k >= 2, BOSS chunks can be built independently, in parallel.

    // there are real k-mers that end with $, so we do the following trick
    const std::string empty_real_name = real_names.front() + "_empty";
    Encoder<get_int_t<T>> empty_real(empty_real_name, ENCODER_BUFFER_SIZE);
    empty_real.finish();

    // the first chunk has at most alphabet_size + 1 k-mers (where F = $)
    std::vector<std::string> dummy_names;
    for (size_t i = 0; i < dummy_source_names.size(); ++i) {
        dummy_names.push_back(dummy_source_names[i][0]);
    }
    BOSS::Chunk chunk = build_boss_chunk<T>(true, empty_real_name, dummy_names,
                                            k, both_strands_mode, bits_per_count, swap_dir);
    logger->trace("Chunk ..$. constructed");
    // construct all other chunks in parallel
    #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
    for (size_t F = 0; F < real_names.size(); ++F) {
        std::vector<std::string> dummy_names { dummy_sink_names[F] };
        for (size_t i = 0; i < dummy_source_names.size(); ++i) {
            dummy_names.push_back(dummy_source_names[i][F + 1]);
        }
        BOSS::Chunk next = build_boss_chunk<T>(false, // not the first chunk
                                               real_names[F], dummy_names,
                                               k, both_strands_mode, bits_per_count, swap_dir);
        logger->trace("Chunk ..{}. constructed", KmerExtractor2Bit().alphabet[F]);
        #pragma omp ordered
        chunk.extend(next);
        logger->trace("Chunk ..{}. appended", KmerExtractor2Bit().alphabet[F]);
    }
    return chunk;
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


template <typename KmerCollector>
class BOSSChunkConstructor : public IBOSSChunkConstructor {
  public:
    BOSSChunkConstructor(size_t k,
                         bool both_strands_mode,
                         uint8_t bits_per_count,
                         const std::string &filter_suffix,
                         size_t num_threads,
                         double memory_preallocated,
                         const std::filesystem::path &swap_dir,
                         size_t disk_cap_bytes)
        : swap_dir_(swap_dir),
          kmer_collector_(k + 1,
                          both_strands_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          swap_dir,
                          disk_cap_bytes,
                          both_strands_mode && filter_suffix.empty() /* keep only canonical k-mers */),
          bits_per_count_(bits_per_count) {
        if (filter_suffix.size()
                && filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
            kmer_collector_.add_kmer(std::vector<TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_sequence(std::string_view sequence, uint64_t count) {
        kmer_collector_.add_sequence(sequence, count);
    }

    void add_sequences(std::vector<std::string>&& sequences) {
        kmer_collector_.add_sequences(std::move(sequences));
    }

    void add_sequences(std::vector<std::pair<std::string, uint64_t>>&& sequences) {
        kmer_collector_.add_sequences(std::move(sequences));
    }

    BOSS::Chunk build_chunk() {
        BOSS::Chunk result;

        using T_INT = typename KmerCollector::Data::value_type;
        using T = get_kmer_t<typename KmerCollector::Kmer, T_INT>;

        if constexpr(std::is_same_v<typename KmerCollector::Extractor,
                                    KmerExtractorBOSS>) {
            assert(kmer_collector_.suffix_length());
            // kmer_collector stores (BOSS::k_ + 1)-mers

            typename KmerCollector::Data &kmer_ints = kmer_collector_.data();

            result = BOSS::Chunk(kmer_collector_.alphabet_size(),
                                 kmer_collector_.get_k() - 1,
                                 kmer_collector_.is_both_strands_mode(),
                                 reinterpret_container<T>(kmer_ints),
                                 bits_per_count_,
                                 swap_dir_);
        } else {
            static_assert(std::is_same_v<typename KmerCollector::Extractor,
                                         KmerExtractor2Bit>);
            assert(!kmer_collector_.suffix_length());

            logger->trace("Reconstructing all required dummy source k-mers...");
            Timer timer;

#define CONSTRUCT_CHUNK(KMER) \
    using T_FINAL = utils::replace_first_t<KMER, T>; \
    if constexpr(utils::is_instance_v<typename KmerCollector::Data, common::ChunkedWaitQueue>) { \
        result = construct_boss_chunk_disk<T_FINAL>(kmer_collector_, bits_per_count_, swap_dir_); \
    } else { \
        result = construct_boss_chunk<T_FINAL>(kmer_collector_, bits_per_count_, swap_dir_); \
    }

            if (kmer_collector_.get_k() * KmerExtractorBOSS::bits_per_char <= 64) {
                CONSTRUCT_CHUNK(KmerExtractorBOSS::Kmer64);
            } else if (kmer_collector_.get_k() * KmerExtractorBOSS::bits_per_char <= 128) {
                CONSTRUCT_CHUNK(KmerExtractorBOSS::Kmer128);
            } else {
                CONSTRUCT_CHUNK(KmerExtractorBOSS::Kmer256);
            }
        }

        return result;
    }

    uint64_t get_k() const { return kmer_collector_.get_k() - 1; }

  private:
    std::filesystem::path swap_dir_;
    KmerCollector kmer_collector_;
    uint8_t bits_per_count_;
};

template <template <typename, class> class KmerContainer, typename... Args>
static std::unique_ptr<IBOSSChunkConstructor>
initialize_boss_chunk_constructor(size_t k,
                                  bool canonical_mode,
                                  uint8_t bits_per_count,
                                  const std::string &filter_suffix,
                                  const Args& ...args) {
    if (k < 1 || k > 256 / KmerExtractorBOSS::bits_per_char - 1) {
        // DBGSuccinct::k = BOSS::k + 1
        logger->error("For succinct graph, k must be between 2 and {}",
                      256 / KmerExtractorBOSS::bits_per_char);
        exit(1);
    }

#define ARGS k, canonical_mode, bits_per_count, filter_suffix, args...

    // collect real k-mers in tight layout only if suffix is empty
    if (filter_suffix.empty()) {
        if ((k + 1) * KmerExtractor2Bit::bits_per_char <= 64) {
            using Container = KmerContainer<KmerExtractor2Bit::KmerBOSS64, KmerExtractor2Bit>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);

        } else if ((k + 1) * KmerExtractor2Bit::bits_per_char <= 128) {
            using Container = KmerContainer<KmerExtractor2Bit::KmerBOSS128, KmerExtractor2Bit>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);

        } else {
            using Container = KmerContainer<KmerExtractor2Bit::KmerBOSS256, KmerExtractor2Bit>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);
        }
    } else {
        if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 64) {
            using Container = KmerContainer<KmerExtractorBOSS::Kmer64, KmerExtractorBOSS>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);

        } else if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 128) {
            using Container = KmerContainer<KmerExtractorBOSS::Kmer128, KmerExtractorBOSS>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);
        } else {
            using Container = KmerContainer<KmerExtractorBOSS::Kmer256, KmerExtractorBOSS>;
            return std::make_unique<BOSSChunkConstructor<Container>>(ARGS);
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
                                  const std::filesystem::path &swap_dir,
                                  size_t disk_cap_bytes) {
#define OTHER_ARGS k, canonical_mode, bits_per_count, filter_suffix, \
                   num_threads, memory_preallocated, swap_dir, disk_cap_bytes

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
            logger->error("Invalid container type {}", (int)container_type);
            std::exit(1);
    }
}

} // namespace boss
} // namespace graph
} // namespace mtg
