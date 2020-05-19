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

template <typename T>
class Decoder {
  public:
    Decoder(const std::string& name) : name_(name), source_(name, ios::binary) {

    }
    std::optional<T> next() {
        T result;
        if (source_.read(reinterpret_cast<char *>(&result), sizeof(T))) {
            return result;
        }
        return std::nullopt;
    }
  private:
    std::string name_;
    std::ifstream source_;
};

template <typename T>
class ConcatDecoder {
  public:
    ConcatDecoder(const std::vector<std::string>& names) : names_(names), source_(names[0]) {

    }
    std::optional<T> next() {
        std::optional<T> next = source_.next();
        if (next.has_value()) {
            return next;
        } else {
            while (++idx_ < names_.size()) {
                source_ = Decoder<T>(names_[idx_]);
                next = source_.next();
                if (next.has_value()) {
                    return next;
                }
            }
            return std::nullopt;
        }
    }
  private:
    uint32_t idx_ = 0;
    std::vector<std::string> names_;
    Decoder<T> source_;
};

template <typename T>
class Encoder {
  public:
    Encoder(const std::string &file_name)
        : file_name_(file_name), sink_(file_name, ios::binary) {};

    void add(const T &value) {
        sink_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

    inline const std::string &name() { return file_name_; }

    void finish() { sink_.close(); }

  private:
    std::string file_name_;
    std::ofstream sink_;
};

/**
 * Push a k-mer into the buffer, then iterate backwards and check if there is a dummy
 * incoming edge (a dummy edge with identical suffix and label as the current one)
 * into the same node.
 * If such an edge is found, mark it for removal (because it's redundant).
 * @param kmer check for redundancy against this k-mer
 * @param buffer queue that stores the last  |Alphabet|^2 k-mers used for identifying
 * redundant dummy source k-mers
 */
template <typename T, typename INT>
static void push_and_remove_redundant_dummy_source(
        size_t k,
        const T &el,
        RecentKmers<T> *buffer,
        Encoder<INT> *encoder,
        std::vector<Encoder<INT>> *dummy_kmer_chunks,
        size_t *num_dummy_l1_kmers) {
    assert(get_first(el)[0] != BOSS::kSentinelCode);

    buffer->push_back({ el, false });
    if (buffer->full()) {
        *num_dummy_l1_kmers
                += write_kmer(k, buffer->pop_front(), encoder, dummy_kmer_chunks);
    }

    using KMER = get_first_type_t<T>;
    const KMER &kmer = get_first(el);
    if (buffer->size() < 2) {
        return;
    }

    TAlphabet curW = kmer[0];
    assert(curW != BOSS::kSentinelCode);

    for (auto it = ++buffer->rbegin();
         KMER::compare_suffix(kmer, get_first((*it).data), 1); ++it) {
        KMER &prev_kmer = get_first((*it).data);
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
template <typename T, typename INT>
uint8_t write_kmer(size_t k,
                   const KmerBuffered<T> &to_write,
                   Encoder<INT> *dummy_l1,
                   std::vector<Encoder<INT>> *dummy_kmer_chunks) {
    if (to_write.is_removed) { // redundant dummy k-mer
        return 0;
    }
    auto kmer_to_write = get_first(to_write.data);
    assert(kmer_to_write[0] != BOSS::kSentinelCode) ;
    const TAlphabet node_last_char = kmer_to_write[1];

    if (node_last_char != BOSS::kSentinelCode) { // not a dummy source kmer
        return 0;
    }
    dummy_l1->add(kmer_to_write.data());

    kmer_to_write.to_prev(k + 1, BOSS::kSentinelCode);
    (*dummy_kmer_chunks)[kmer_to_write[0]].add(kmer_to_write.data());
    return 1;
}

template <typename KMER, typename T_INT, typename INT>
std::optional<INT> next_dummy_source(size_t k, std::vector<Decoder<T_INT>> *dummy_source_it,
                                     std::vector<INT> *dummy_source_it_v) {
    if (dummy_source_it->empty()) {
        return {};
    }
    INT result = (*dummy_source_it_v)[0];
    uint32_t idx = 0; // start at 1, to skip the $$...$ k-mer
    for (uint32_t i = 1; i < dummy_source_it_v->size(); ++i) {
        if ((*dummy_source_it_v)[i] < result) {
            result = (*dummy_source_it_v)[i];
            idx = i;
        }
    }
    std::optional<T_INT> next = (*dummy_source_it)[idx].next();
    if (next.has_value()) {
        INT to_push = get_first(next.value());
        reinterpret_cast<KMER &>(to_push).to_prev(k + 1, BOSS::kSentinelCode);
        (*dummy_source_it_v)[idx] = to_push;
    } else {
        dummy_source_it->erase(dummy_source_it->begin() + idx);
        dummy_source_it_v->erase(dummy_source_it_v->begin() + idx);
    }
    return result;
}

constexpr size_t ENCODER_BUFFER_SIZE = 100'000;
constexpr uint8_t ALPHABET_LEN = 1 << KmerExtractorBOSS::bits_per_char;

/**
 * Splits #kmers by the first character into ALPHABET_LEN blocks.
 */
template <typename T>
std::vector<std::string> split_by_first(size_t k,
                                        const std::filesystem::path &tmp_dir,
                                        const ChunkedWaitQueue<T> &kmers) {
    std::vector<std::ofstream> original_chunks;
    std::vector<std::string> original_chunk_names(ALPHABET_LEN);
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        original_chunk_names[i] = tmp_dir / ("original_chunk_" + std::to_string(i));
        original_chunks.emplace_back(original_chunk_names[i], ios::binary);
    }

    logger->trace("Splitting k-mers into {} chunks...", ALPHABET_LEN);
    size_t num_parent_kmers = 0;
    for (auto &it = kmers.begin(); it != kmers.end(); ++it) {
        num_parent_kmers++;
        const T &kmer = *it;
        //push(kmer, &original_chunks[get_first(kmer)[k]]);
        original_chunks[get_first(kmer)[k]].write(reinterpret_cast<const char *>(&kmer),
                                                  sizeof(T));
    }
    std::for_each(original_chunks.begin(), original_chunks.end(),
                  [](std::ofstream &e) { e.close(); });
    logger->trace("Total number of k-mers: {}", num_parent_kmers);
    return original_chunk_names;
}

template <typename T, typename INT>
void push(const INT &v, RecentKmers<T> *buffer) {
    T to_push;
    if constexpr (utils::is_pair_v<T>) {
        to_push = { get_first_type_t<T>(v), 0 };
    } else {
        to_push = T(v);
    }
    assert(get_first(to_push)[0] != BOSS::kSentinelCode);
    buffer->push_back({ to_push, false });
}

/**
 * Generates dummy-1 source k-mers and dummy sink kmers from #kmers and merges them into
 * #merged_l1_name. The dummy-2 source kmers for each first character are returned.
 */
template <typename T>
std::pair<std::vector<std::string>, std::vector<std::string>>
generate_dummy_1_kmers(size_t k,
                       const std::filesystem::path &tmp_dir,
                       ChunkedWaitQueue<T> *kmers,
                       const std::string &dummy_l1_name,
                       const std::string &dummy_sink_name) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using T_INT = get_int_t<T>; // either KMER::WordType or <KMER::WordType, count>
    using INT = typename KMER::WordType; // 64/128/256-bit integer

    std::vector<std::string> original_chunk_names = split_by_first(k, tmp_dir, *kmers);

    std::vector<Encoder<INT>> dummy_l2_chunks;
    std::vector<Encoder<INT>> dummy_sink_chunks;
    std::vector<std::string> dummy_l2_names(ALPHABET_LEN);
    std::vector<std::string> dummy_sink_names(ALPHABET_LEN);
    std::vector<Decoder<T_INT>> dummy_source_it;
    std::vector<Decoder<T_INT>> dummy_sink_it;
    std::vector<INT> dummy_source_it_v;
    std::vector<std::optional<T_INT>> dummy_sink_it_v(ALPHABET_LEN);
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        original_chunk_names[i] = tmp_dir/("original_chunk_" + std::to_string(i));
        dummy_l2_names[i] = tmp_dir/("dummy_source_2_" + std::to_string(i));
        dummy_sink_names[i] = tmp_dir/("dummy_sink_" + std::to_string(i));
        dummy_l2_chunks.emplace_back(dummy_l2_names[i]);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i]);
        dummy_sink_it.emplace_back(original_chunk_names[i]);
        dummy_sink_it_v[i] = dummy_sink_it[i].next();
        dummy_source_it.emplace_back(original_chunk_names[i]);
        std::optional<T_INT> first = dummy_source_it.back().next();
        if (first.has_value() && get_first(first.value()) != INT(0)) {
            reinterpret_cast<KMER &>(first.value()).to_prev(k + 1, BOSS::kSentinelCode);
            dummy_source_it_v.push_back(get_first(first.value()));
        } else {
            dummy_source_it.pop_back();
        }
    }

    logger->trace("Generating dummy-1 source kmers and dummy sink k-mers...");

    std::vector<KMER> last_dummy(ALPHABET_LEN, KMER(0));
    RecentKmers<T> recent_buffer(1llu << 2 * KMER::kBitsPerChar);
    Encoder<INT> dummy_l1(dummy_l1_name);
    std::optional<INT> dummy_source = next_dummy_source<KMER, T_INT, INT>(k, &dummy_source_it, &dummy_source_it_v);
    size_t num_dummy_l1_kmers = 0;

    auto push_dummy_source = [&](const INT &max) {
        std::optional<INT> prev_dummy_source;
        while (dummy_source.has_value() && dummy_source.value() < max) {
            if (prev_dummy_source != dummy_source) { // skip identical dummy-1 k-mers
                push(dummy_source.value(), &recent_buffer);
                if (recent_buffer.full()) {
                    num_dummy_l1_kmers += write_kmer(k, recent_buffer.pop_front(),
                                                     &dummy_l1, &dummy_l2_chunks);
                }
                prev_dummy_source = dummy_source;
            }
            dummy_source
                    = next_dummy_source<KMER, T_INT, INT>(k, &dummy_source_it, &dummy_source_it_v);
        }
    };

    ConcatDecoder<T_INT> it(original_chunk_names);
    // skip the main dummy kmer
    it.next();
    for(auto v = it.next(); v.has_value(); v = it.next())
    {
        T kmer_with_count = reinterpret_cast<const T&>(v.value());
        const KMER &kmer = get_first(kmer_with_count);
        KMER dummy_sink = kmer;
        TAlphabet first_char = dummy_sink[0];
        dummy_sink.to_next(k+1, BOSS::kSentinelCode);

        if (last_dummy[first_char] != dummy_sink) { // distinct dummy sink k-mer
            last_dummy[first_char] = dummy_sink;

            while (dummy_sink_it_v[first_char].has_value()
                   && get_first(dummy_sink_it_v[first_char].value()) < dummy_sink.data()) {
                dummy_sink_it_v[first_char] = dummy_sink_it[first_char].next();
            }
            if (!dummy_sink_it_v[first_char].has_value()
                || !KMER::compare_suffix(reinterpret_cast<const KMER &>(
                                                 get_first(dummy_sink_it_v[first_char].value())),
                                         dummy_sink)) {
                dummy_sink_chunks[first_char].add(dummy_sink.data());
                std::cout << dummy_sink.to_string(k+1, "$ACGT") << " not redundant " << std::endl;
            }
        }
        push_dummy_source(kmer.data());
        push_and_remove_redundant_dummy_source(k, kmer_with_count, &recent_buffer, &dummy_l1,
                                               &dummy_l2_chunks, &num_dummy_l1_kmers);
    }
    push_dummy_source(INT(-1)); // leftover dummy sources
    while (!recent_buffer.empty()) { // add leftover elements from buffer
        num_dummy_l1_kmers += write_kmer(k, recent_buffer.pop_front(), &dummy_l1,
                                         &dummy_l2_chunks);
    }
    dummy_l1.finish();
    for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
        dummy_sink_chunks[i].finish();
        dummy_l2_chunks[i].finish();
    }
    common::concat(dummy_sink_names, dummy_sink_name);

    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  num_dummy_l1_kmers);

    return {original_chunk_names, dummy_l2_names};
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

    std::string dummy_l1_name = tmp_dir/"dummy_l1"; //TODO: all params or all return type
    std::string dummy_sink_name = tmp_dir/"dummy_sink";
    std::vector<std::string> original_names;
    std::vector<std::string> dummy_names;
    std::tie(original_names, dummy_names)
            = generate_dummy_1_kmers(k, tmp_dir, kmers, dummy_l1_name, dummy_sink_name);

    // stores the sorted original kmers and dummy-1 k-mers
    std::vector<std::string> files_to_merge;
    files_to_merge.push_back(dummy_sink_name);
    files_to_merge.push_back(dummy_l1_name);
    std::vector<std::string> dummy_next_names(ALPHABET_LEN);
    // generate dummy k-mers of prefix length 2..k
    logger->trace("Starting generating dummy-k source k-mers...");
    for (size_t dummy_pref_len = 2; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will compress all sorted dummy k-mers of given prefix length
        files_to_merge.push_back(tmp_dir/("dummy_l" + std::to_string(dummy_pref_len)));
        Encoder<INT> encoder(files_to_merge.back());

        std::vector<Encoder<INT>> dummy_next_chunks;
        for (uint32_t i = 0; i < ALPHABET_LEN; ++i) {
            dummy_next_names[i] = tmp_dir/("dummy_source_"
                    + std::to_string(dummy_pref_len + 1) + "_" + std::to_string(i));
            dummy_next_chunks.emplace_back(dummy_next_names[i]);
        }
        size_t num_kmers = 0;
        const std::function<void(const INT &)> &write_dummy = [&](const INT &v) {
            encoder.add(v);
            KMER kmer(v);
            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            dummy_next_chunks[kmer[0]].add(kmer.data());
            num_kmers++;
        };
        common::merge_files2(dummy_names, write_dummy);

        encoder.finish();
        std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                      [](auto &v) { v.finish(); });
        std::swap(dummy_names, dummy_next_names);
        logger->trace("Number of dummy k-mers with dummy prefix of length {} : {}",
                      dummy_pref_len, num_kmers);
    }
    std::for_each(dummy_names.begin(), dummy_names.end(),
                  [](const string &v) { std::filesystem::remove(v); });
    // at this point, we have the original k-mers and dummy-1 k-mers in original_and_dummy_l1,
    // the dummy-x k-mers in dummy_source_{x}, and we merge them all into a single stream
    kmers->reset();
    async_worker.enqueue([kmers, original_names, files_to_merge, k]() { //TODO: remove k
        std::function<void(const T_INT &)> on_new_item
                = [kmers, k](const T_INT &v) {
std::cout << get_first(reinterpret_cast<const T &>(v)).to_string(k+1, "$ACGT") << std::endl;
                      kmers->push(reinterpret_cast<const T &>(v)); };
        common::merge_dummy2(original_names, files_to_merge, on_new_item);
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
