#include "dbg_hash_fast.hpp"

#include <cassert>
#include <fstream>
#include <limits>

#include <tsl/ordered_set.h>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/serialization.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using mtg::kmer::KmerExtractor2Bit;


template <typename KMER = KmerExtractor2Bit::Kmer64>
class DBGHashFastImpl : public DBGHashFast::DBGHashFastInterface {
    using Kmer = KMER;
    using KmerWord = typename KMER::WordType;
    using TAlphabet = KmerExtractor2Bit::TAlphabet;

#if _PROTEIN_GRAPH
    using Flags = uint32_t;
#elif _DNA_CASE_SENSITIVE_GRAPH
    using Flags = uint16_t;
#elif _DNA5_GRAPH
    using Flags = uint16_t;
#elif _DNA_GRAPH
    using Flags = uint8_t;
#else
    static_assert(false, "invalid or undefined alphabet");
#endif

    using KmerIndex = tsl::ordered_set<KmerWord,
                                       utils::Hash<KmerWord>,
                                       std::equal_to<KmerWord>,
                                       std::allocator<KmerWord>,
                                       std::vector<KmerWord, std::allocator<KmerWord>>,
                                       std::uint64_t>;

    using KmerConstIterator = typename KmerIndex::const_iterator;

  public:
    DBGHashFastImpl(size_t k,
                    bool canonical_mode,
                    bool packed_serialization,
                    size_t reserve = 0)
          : k_(k),
            canonical_mode_(canonical_mode),
            packed_serialization_(packed_serialization),
            kIgnoreLastCharMask(~(KmerWord((1llu << kBitsPerChar) - 1)
                                        << static_cast<int>(kBitsPerChar * (k - 1)))) {
        kmers_.reserve(reserve);
        bits_.reserve(reserve);
    }

    void add_sequence(std::string_view sequence,
                      const std::function<void(node_index)> &on_insertion);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(std::string_view sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string_view sequence,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const;

    void call_nodes(const std::function<void(node_index)> &callback,
                    const std::function<bool()> &stop_early) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const {
        assert(in_graph(node));

        // TODO: use `next_kmer()`

        auto kmer = get_kmer(node);
        kmer.to_next(k_, seq_encoder_.encode(next_char));
        return get_node_index(kmer);
    }
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const {
        assert(in_graph(node));

        // TODO: check previous k-mer in vector similarly to `next_kmer()`

        auto kmer = get_kmer(node);
        kmer.to_prev(k_, seq_encoder_.encode(prev_char));
        return get_node_index(kmer);
    }

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        assert(in_graph(node));

        call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
    }
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        assert(in_graph(node));

        call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
    }

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index node) const {
        assert(in_graph(node));
        return outdegree(node) == 1;
    }
    bool has_multiple_outgoing(node_index node) const {
        assert(in_graph(node));
        return outdegree(node) > 1;
    }

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index node) const {
        assert(in_graph(node));
        return indegree(node) == 1;
    }

    node_index kmer_to_node(std::string_view kmer) const {
        assert(kmer.length() == k_);
        return get_node_index(seq_encoder_.encode(kmer));
    }

    std::string get_node_sequence(node_index node) const {
        assert(in_graph(node));
        return seq_encoder_.kmer_to_sequence(get_kmer(node), k_);
    }

    size_t get_k() const { return k_; }
    bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const {
        uint64_t nnodes = 0;
        call_nodes([&](auto) { nnodes++; }, [](){ return false; });
        return nnodes;
    }
    uint64_t max_index() const { return kmers_.size() * kAlphabetSize; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const {
        std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                          std::ios::binary);
        serialize(out);
    }

    bool load(std::istream &in);
    bool load(const std::string &filename) {
        std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                         std::ios::binary);
        return load(in);
    }

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const;

    const std::string& alphabet() const { return seq_encoder_.alphabet; }

  private:
    bool in_graph(node_index node) const {
        assert(node > 0 && node <= max_index());

        Flags flags = bits_[node_to_bucket(node)];

        return (flags >> ((node - 1) % kAlphabetSize)) & static_cast<Flags>(1);
    }

    Vector<std::pair<Kmer, bool>> sequence_to_kmers(std::string_view sequence,
                                                    bool canonical = false) const {
        return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, canonical);
    }

    size_t node_to_bucket(node_index node) const { return (node - 1) / kAlphabetSize; }
    node_index bucket_to_node(size_t bucket) const { return bucket * kAlphabetSize + 1; }

    node_index get_node_index(const Kmer &kmer) const;

    const Kmer get_kmer(node_index node) const {
        assert(in_graph(node));

        KmerWord c = (node - 1) % kAlphabetSize;

        return KMER(*(kmers_.begin() + node_to_bucket(node))
                        | (c << static_cast<int>(kBitsPerChar * (k_ - 1))));
    }

    KmerConstIterator next_kmer(node_index node) const;

    void print_internal_representation() const {
        for (auto it = kmers_.begin(); it != kmers_.end(); ++it) {
            size_t bucket = it - kmers_.begin();
            std::cout << bucket << ", " << bucket_to_node(bucket) << ": ";
            std::cout << KMER(*it).to_string(k_ - 1, seq_encoder_.alphabet) << ' ';
            std::cout << std::bitset<std::numeric_limits<Flags>::digits>(bits_[bucket]) << '\n';
        }
        std::cout << std::flush;
    }

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;
    std::vector<Flags> bits_;
    KmerExtractor2Bit seq_encoder_;

    bool packed_serialization_;

    const KmerWord kIgnoreLastCharMask;

    // Flags layout:
    //     <--- leading 0's --->
    //     < 1-bit flag for "may contain source kmer" >
    //     < alphabet_size length bit field for each possible k'th char in kmer>
    const uint8_t kAlphabetSize = seq_encoder_.alphabet.size();
    const Flags kLastCharMask = (Flags(1) << kAlphabetSize) - 1;
    const Flags kMayBeSourceKmer = Flags(1) << kAlphabetSize;

    static constexpr auto kExtension = DBGHashFast::kExtension;
    static constexpr uint8_t kBitsPerChar = KMER::kBitsPerChar;
};

// Insert sequence to graph and invoke callback |on_insertion| for each new
// node index augmenting the range [1,...,max_index], including those not
// pointing to any real node in graph. That is, the callback is invoked for
// all new real nodes and all new dummy node indexes allocated in graph.
// In short: max_index[after] = max_index[before] + {num_invocations}.
template <typename KMER>
void DBGHashFastImpl<KMER>::add_sequence(std::string_view sequence,
                                         const std::function<void(node_index)> &on_insertion) {
    auto add_seq = [&](const auto &sequence) {
        bool previous_valid = false;

        for (const auto &kmer_pair : sequence_to_kmers(sequence)) {
            const auto &[kmer, is_valid] = kmer_pair;
            if (!is_valid) {
                previous_valid = false;
                continue;
            }

            const KmerWord key = kmer.data() & kIgnoreLastCharMask;

            auto [iter, inserted] = kmers_.insert(key);

            // TODO: if previous k-mer wasn't inserted (and hence, had
            // been inserted earlier), compare the current k-mer with get_next.

            Flags char_flag = Flags(1) << kmer[k_ - 1];

            if (inserted) {
                bits_.push_back(previous_valid ? char_flag
                                               : char_flag | kMayBeSourceKmer);
                // call all indexes inserted (only one of them is a real node)
                uint64_t offset = 1 + (iter - kmers_.begin()) * kAlphabetSize;
                for (TAlphabet c = 0; c < kAlphabetSize; ++c) {
                    on_insertion(offset + c);
                }
            } else {
                auto &flags = bits_[iter - kmers_.begin()];
                flags |= char_flag;
                if (previous_valid)
                    flags &= ~kMayBeSourceKmer;
            }

            previous_valid = true;

            assert(iter != kmers_.end());
            assert(iter == kmers_.find(key));
        }
    };

    add_seq(sequence);

    if (canonical_mode_) {
        std::string rev_comp(sequence.begin(), sequence.end());
        reverse_complement(rev_comp.begin(), rev_comp.end());

        add_seq(rev_comp);
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
template <typename KMER>
void DBGHashFastImpl<KMER>::map_to_nodes_sequentially(
                                std::string_view sequence,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers(sequence)) {
        if (terminate())
            return;

        assert((kmers_.find(kmer.data() & kIgnoreLastCharMask) == kmers_.end())
               || get_node_index(kmer) == npos
               || kmer == get_kmer(get_node_index(kmer)));

        callback(is_valid ? get_node_index(kmer) : npos);
        // TODO: `next_kmer()` could speed this up
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
template <typename KMER>
void DBGHashFastImpl<KMER>::map_to_nodes(std::string_view sequence,
                                         const std::function<void(node_index)> &callback,
                                         const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers(sequence, canonical_mode_)) {
        if (terminate())
            return;

        assert(!get_node_index(kmer) || kmer == get_kmer(get_node_index(kmer)));

        callback(is_valid ? get_node_index(kmer) : npos);
        // TODO: `next_kmer()` could speed this up
    }
}


template <typename KMER>
typename DBGHashFastImpl<KMER>::KmerConstIterator
DBGHashFastImpl<KMER>::next_kmer(node_index node) const {
    assert(in_graph(node));

    KmerWord last_c = (node - 1) % kAlphabetSize;

    KmerWord next_kmer_prefix
        = (*(kmers_.begin() + node_to_bucket(node))
                | (last_c << static_cast<int>(kBitsPerChar * (k_ - 1))))
            >> kBitsPerChar;

    size_t bucket = node_to_bucket(node);
    KmerConstIterator next_it = kmers_.begin() + (bucket + 1);

    return next_it != kmers_.end() && next_it.key() == next_kmer_prefix
            ? next_it
            : kmers_.find(next_kmer_prefix);
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_outgoing_kmers(node_index node,
                                                const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    KmerConstIterator next_kmer_it = next_kmer(node);
    if (next_kmer_it == kmers_.end())
        return;

    size_t bucket = next_kmer_it - kmers_.begin();
    Flags flags = bits_[bucket];

    for (TAlphabet c = 0; c < kAlphabetSize; ++c, flags >>= 1) {
        if (flags & static_cast<Flags>(1)) {
            assert(in_graph(bucket_to_node(bucket) + c));

            callback(bucket_to_node(bucket) + c, seq_encoder_.decode(c));
        }
    }
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_incoming_kmers(node_index node,
                                                const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));

    const KMER &kmer = get_kmer(node);

    for (TAlphabet c = 0; c < kAlphabetSize; ++c) {
        KMER prev_kmer = kmer;
        prev_kmer.to_prev(k_, c);

        auto prev_kmer_index = get_node_index(prev_kmer);
        if (prev_kmer_index != npos)
            callback(prev_kmer_index, seq_encoder_.decode(c));
    }
}

template <typename KMER>
size_t DBGHashFastImpl<KMER>::outdegree(node_index node) const {
    assert(in_graph(node));

    KmerConstIterator next_kmer_it = next_kmer(node);

    return next_kmer_it != kmers_.end()
            ? sdsl::bits::cnt(bits_[next_kmer_it - kmers_.begin()] & kLastCharMask)
            : 0;
}

template <typename KMER>
size_t DBGHashFastImpl<KMER>::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;
    call_incoming_kmers(node, [&](auto, auto) { indegree++; });
    return indegree;
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    if (!(bits_[node_to_bucket(node)] & kMayBeSourceKmer))
        return false;

    const KMER &kmer = get_kmer(node);

    for (TAlphabet c = 0; c < kAlphabetSize; ++c) {
        KMER prev_kmer = kmer;
        prev_kmer.to_prev(k_, c);

        // return false if has at least one incoming
        if (get_node_index(prev_kmer) != npos)
            return false;
    }

    return true;
}

template <typename KMER>
void DBGHashFastImpl<KMER>::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    out.exceptions(out.badbit | out.failbit);

    serialize_number(out, k_);

    Serializer serializer(out);

    if (packed_serialization_) {
        serialize_number(out, kmers_.size());
        std::for_each(kmers_.begin(), kmers_.end(), serializer);
    } else {
        serialize_number(out, std::numeric_limits<uint64_t>::max());
        kmers_.serialize(serializer);
    }

    std::for_each(bits_.begin(), bits_.end(), serializer);

    serialize_number(out, canonical_mode_);
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::load(std::istream &in) {
    if (!in.good())
        return false;

    in.exceptions(in.badbit | in.failbit | in.eofbit);

    kmers_.clear();
    bits_.clear();

    try {
        k_ = load_number(in);

        Deserializer deserializer(in);

        uint64_t tag = load_number(in);

        if (tag < std::numeric_limits<uint64_t>::max()) {
            packed_serialization_ = true;

            const auto size = tag;
            kmers_.reserve(size + 1);
            for (uint64_t i = 0; i < size; ++i) {
                kmers_.insert(deserializer.operator()<KmerWord>());
            }

        } else {
            packed_serialization_ = false;
            kmers_ = KmerIndex::deserialize(deserializer, true);
        }

        bits_.resize(kmers_.size());
        for (auto &flags : bits_) {
            flags = deserializer.operator()<Flags>();
        }

        canonical_mode_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (!dynamic_cast<const DBGHashFastImpl*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashFastImpl*>(&other);

    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(canonical_mode_ == other_hash.canonical_mode_);
    assert(kmers_.size() == other_hash.kmers_.size());

    return kmers_ == other_hash.kmers_;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::get_node_index(const Kmer &kmer) const {
    const auto it = kmers_.find(kmer.data() & kIgnoreLastCharMask);
    if (it == kmers_.end())
        return npos;

    size_t bucket = it - kmers_.begin();

    if (!((bits_[bucket] >> kmer[k_ - 1]) & static_cast<Flags>(1)))
        return npos;

    assert(get_kmer(bucket_to_node(bucket) + kmer[k_ - 1]) == kmer);

    return bucket_to_node(bucket) + kmer[k_ - 1];
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_nodes(const std::function<void(node_index)> &callback,
                                       const std::function<bool()> &stop_early) const {
    for (size_t i = 0; i < kmers_.size(); ++i) {
        Flags flags = bits_[i];

        for (TAlphabet c = 0; c < kAlphabetSize; ++c, flags >>= 1) {
            if (stop_early())
                return;

            if (flags & static_cast<Flags>(1))
                callback(bucket_to_node(i) + c);
        }
    }
}

std::unique_ptr<DBGHashFast::DBGHashFastInterface>
DBGHashFast::initialize_graph(size_t k,
                              bool canonical_mode,
                              bool packed_serialization) {
    if (k < 1 || k > 256 / KmerExtractor2Bit::bits_per_char) {
        logger->error("For hash graph, k must be between 1 and {}",
                      256 / KmerExtractor2Bit::bits_per_char);
        exit(1);
    }

    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode, packed_serialization
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode, packed_serialization
        );
    } else {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode, packed_serialization
        );
    }
}

bool DBGHashFast::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        auto pos = in.tellg();
        auto k = load_number(in);
        in.seekg(pos, in.beg);

        // the actual value of |canonical| will be set in load
        hash_dbg_ = initialize_graph(k, false, false);
        return hash_dbg_->load(in) && in.good();
    } catch (...) {
        return false;
    }
}

bool DBGHashFast::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

} // namespace graph
} // namespace mtg
