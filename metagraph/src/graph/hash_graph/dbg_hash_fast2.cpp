#include "dbg_hash_fast2.hpp"

#include <cassert>

#include <tsl/robin_map.h>
#include <libmaus2/util/NumberSerialisation.hpp>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"


template <typename KMER = KmerExtractor2Bit::Kmer64>
class DBGHashFast2Impl : public DBGHashFast2::DBGHashFast2Interface {
    using Kmer = KMER;
    using TAlphabet = KmerExtractor2Bit::TAlphabet;
    using KmerPrefix = typename KMER::WordType;

#if _PROTEIN_GRAPH
    using Bits = uint32_t;
#elif _DNA_CASE_SENSITIVE_GRAPH
    using Bits = uint16_t;
#elif _DNA5_GRAPH
    using Bits = unsigned char;
#elif _DNA_GRAPH
    using Bits = unsigned char;
#else
    static_assert(false, "invalid or undefined alphabet");
#endif

    using KmerIndex = tsl::robin_map<Kmer,
                                     Bits,
                                     utils::Hash<Kmer>,
                                     std::equal_to<Kmer>,
                                     std::allocator<Kmer>>;

    using KmerIterator = typename KmerIndex::const_iterator;

  public:
    explicit DBGHashFast2Impl(size_t k,
                             bool canonical_mode,
                             bool packed_serialization);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const;

    void call_nodes(const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &stop_early) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index) const;
    bool has_multiple_outgoing(node_index) const;

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index) const;

    node_index kmer_to_node(const std::string &kmer) const;

    std::string get_node_sequence(node_index node) const;

    size_t get_k() const { return k_; }
    bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const { return kmers_.capacity() << KMER::kBitsPerChar; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const;

    const std::string& alphabet() const { return seq_encoder_.alphabet; }

    bool in_graph(node_index node) const;

  private:
    Vector<Kmer> sequence_to_kmers(const std::string &sequence, bool canonical = false) const {
        return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, canonical);
    }

    // TODO: sotre kmer_prefix as keys in kmers_
    auto find_kmer(const KmerPrefix &kmer_prefix) {
        return kmers_.find(KMER(kmer_prefix));
    }

    auto find_kmer(const KmerPrefix &kmer_prefix) const {
        return kmers_.find(KMER(kmer_prefix));
    }

    bool has_edge(typename KmerIndex::const_iterator out_iter, const TAlphabet c) const {
        if (out_iter == kmers_.end())
            return false;
        
        return (*out_iter).second & (Bits(1) << c);
    }

    bool has_edge(const Kmer &out_kmer, const TAlphabet c) const {
        const auto out_find = find_kmer(KmerPrefix(out_kmer.data()) & kIgnoreLastCharMask);
        return has_edge(out_find, c);
    }

    bool has_edge(const Kmer &out_kmer, const Kmer &in_kmer) const {
        return has_edge(out_kmer, in_kmer[1]);
    }

    node_index get_index(const KmerIterator &iter) const;
    node_index get_index(const Kmer &kmer) const;

    typename KmerIndex::iterator get_iter(node_index node);
    typename KmerIndex::const_iterator get_iter(node_index node) const;

    const Kmer get_kmer(node_index node,
                        std::optional<typename KmerIndex::const_iterator> kmer_iter = std::nullopt,
                        bool *has_no_incoming = NULL) const;

    void update_has_no_incoming_flag(node_index node);

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;
    KmerExtractor2Bit seq_encoder_;

    bool packed_serialization_;

    const KmerPrefix kIgnoreLastCharMask = ~(KmerPrefix(0)) >> KMER::kBitsPerChar;
    const Bits kHasNoIncomingFlag = Bits(1) << seq_encoder_.alphabet.size();
    static constexpr auto kExtension = DBGHashFast2::kExtension;
};

template <typename KMER>
DBGHashFast2Impl<KMER>::DBGHashFast2Impl(size_t k,
                                             bool canonical_mode,
                                             bool packed_serialization)
      : k_(k),
        canonical_mode_(canonical_mode),
        packed_serialization_(packed_serialization) {}

template <typename KMER>
void DBGHashFast2Impl<KMER>::add_sequence(const std::string &sequence,
                                            bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    bool first_kmer = true;

    for (const auto &kmer : sequence_to_kmers(sequence)) {
        Bits val = Bits(1) << kmer[0];
        // TODO: ">> kIgnoreLastCharMask" instead of "& kIgnoreLastCharMask"
        auto key = KmerPrefix(kmer.data()) & kIgnoreLastCharMask;

        auto [iter, inserted] = kmers_.insert(std::make_pair(key, val));
        if (!inserted)
            iter.value() |= val;

        if (first_kmer)
            update_has_no_incoming_flag(get_index(iter));

        first_kmer = false;

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(get_index(iter), true);
    }

    if (!canonical_mode_)
        return;

    first_kmer = true;

    for (const auto &kmer : sequence_to_kmers(seq_encoder_.reverse_complement(sequence))) {
        Bits val = Bits(1) << kmer[0];
        auto key = KmerPrefix(kmer.data()) & kIgnoreLastCharMask;

        auto [iter, inserted] = kmers_.insert(std::make_pair(key, val));
        if (!inserted)
            iter.value() |= val;

        if (first_kmer)
            update_has_no_incoming_flag(get_index(iter));

        first_kmer = false;

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(get_index(iter), true);
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
template <typename KMER>
void DBGHashFast2Impl<KMER>::map_to_nodes_sequentially(
                              std::string::const_iterator begin,
                              std::string::const_iterator end,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate) const {
    std::string sequence(begin, end);
    const auto &kmers = sequence_to_kmers(sequence);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {

        assert(it != kmers.end() || !is_valid);

        if (terminate())
            return;

        callback(is_valid ? get_index(*it++) : npos);
    }
    assert(it == kmers.end());
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
template <typename KMER>
void DBGHashFast2Impl<KMER>::map_to_nodes(const std::string &sequence,
                                            const std::function<void(node_index)> &callback,
                                            const std::function<bool()> &terminate) const {
    const auto &kmers = sequence_to_kmers(sequence, canonical_mode_);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {

        assert(it != kmers.end() || !is_valid);

        if (terminate())
            return;

        callback(is_valid ? get_index(*it++) : npos);
    }
    assert(it == kmers.end());
}

template <typename KMER>
void DBGHashFast2Impl<KMER>::call_outgoing_kmers(node_index node,
                                                 const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    KMER kmer = get_kmer(node);
    assert((*find_kmer(KmerPrefix(kmer.data()) & kIgnoreLastCharMask)).second
                && "val is not 0, there at least one k-mer");

    kmer.to_next(k_, 0);

    const auto next_kmer_prefix_it
        = find_kmer(KmerPrefix(kmer.data()) & kIgnoreLastCharMask);

    if (next_kmer_prefix_it == kmers_.end())
        return;

    const Bits val = (*next_kmer_prefix_it).second;

    for (size_t i = 0; i < seq_encoder_.alphabet.size(); ++i) {
        if ((val >> i) & 1) {
            auto next_kmer_index = get_index(next_kmer_prefix_it) + i;

            assert(in_graph(next_kmer_index));

            callback(next_kmer_index, i);
        }
    }
}

template <typename KMER>
void DBGHashFast2Impl<KMER>::call_incoming_kmers(node_index node,
                                                 const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto prev_kmer_iter = find_kmer(KmerPrefix(prev_kmer.data()) & kIgnoreLastCharMask);
        if (has_edge(prev_kmer_iter, c))
            callback(get_index(prev_kmer_iter), c);
    }
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::node_index
DBGHashFast2Impl<KMER>::traverse(node_index node, char next_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return get_index(kmer);
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::node_index
DBGHashFast2Impl<KMER>::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return get_index(kmer);
}

template <typename KMER>
void
DBGHashFast2Impl<KMER>
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

template <typename KMER>
void
DBGHashFast2Impl<KMER>
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

template <typename KMER>
size_t DBGHashFast2Impl<KMER>::outdegree(node_index node) const {
    assert(in_graph(node));

    const auto &kmer = get_kmer(node);
    const auto find = find_kmer(KmerPrefix(kmer.data()) & kIgnoreLastCharMask);
    const auto val = (*find).second;

#if _PROTEIN_GRAPH
    return __builtin_popcountl(val);
#elif _DNA_CASE_SENSITIVE_GRAPH
    return __builtin_popcount(val);
#elif _DNA5_GRAPH
    return __builtin_popcount(val);
#elif _DNA_GRAPH
    return __builtin_popcount(val);
#else
    static_assert(false, "invalid or undefined alphabet");
#endif
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) == 1;
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) > 1;
}

template <typename KMER>
size_t DBGHashFast2Impl<KMER>::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (has_edge(prev_kmer, kmer))
            indegree++;
    }

    return indegree;
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    bool has_no_incoming;
    const auto &kmer = get_kmer(node, std::nullopt, &has_no_incoming);
    std::ignore = kmer;
    return has_no_incoming;
}

template <typename KMER>
void DBGHashFast2Impl<KMER>::update_has_no_incoming_flag(node_index node) {
    assert(in_graph(node));

    auto kmer_iter = get_iter(node);
    const auto &kmer = get_kmer(node, kmer_iter);

    // TODO: make this another function, which is called in add_sequence
    bool has_no_incoming = true;

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (has_edge(prev_kmer, kmer))
            has_no_incoming = false;
    }

    auto val = kmer_iter.value();
    val &= ~kHasNoIncomingFlag;

    if (has_no_incoming)
        val |= kHasNoIncomingFlag;

    kmer_iter.value() = val;
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    return indegree(node) == 1;
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::node_index
DBGHashFast2Impl<KMER>::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    return get_index(seq_encoder_.encode(kmer));
}

template <typename KMER>
std::string DBGHashFast2Impl<KMER>::get_node_sequence(node_index node) const {
    assert(in_graph(node));

    return seq_encoder_.kmer_to_sequence(get_kmer(node), k_);
}

class Serializer {
  public:
    explicit Serializer(std::ostream &os) : os_(os) {}

    template <class T>
    void operator()(const T &value) {
        os_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

  private:
    std::ostream &os_;
};

class Deserializer {
  public:
    explicit Deserializer(std::istream &is) : is_(is) {}

    template <class T>
    T operator()() {
        T value;
        is_.read(reinterpret_cast<char *>(&value), sizeof(T));
        return value;
    }

  private:
    std::istream &is_;
};

template <typename KMER>
void DBGHashFast2Impl<KMER>::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    out.exceptions(out.badbit | out.failbit);

    serialize_number(out, k_);

    Serializer serializer(out);

    if (packed_serialization_) {
        serialize_number(out, kmers_.size());
        std::for_each(kmers_.begin(), kmers_.end(), serializer);
    } else {
        throw std::runtime_error("not implemented");
        //serialize_number(out, std::numeric_limits<uint64_t>::max());
        //kmers_.serialize(serializer);
    }

    serialize_number(out, canonical_mode_);
}

template <typename KMER>
void DBGHashFast2Impl<KMER>::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::load(std::istream &in) {
    if (!in.good())
        return false;

    in.exceptions(in.badbit | in.failbit | in.eofbit);

    kmers_.clear();

    try {
        k_ = load_number(in);

        Deserializer deserializer(in);

        uint64_t tag = load_number(in);

        if (tag < std::numeric_limits<uint64_t>::max()) {
            packed_serialization_ = true;

            const auto size = tag;
            kmers_.reserve(size + 1);
            for (uint64_t i = 0; i < size; ++i) {
                auto keyval = deserializer.operator()<std::pair<Kmer, Bits>>();
                kmers_.insert_or_assign(keyval.first, keyval.second);
            }

        } else {
            throw std::runtime_error("not implemented");
            //packed_serialization_ = false;
            //kmers_ = KmerIndex::deserialize(deserializer, true);
        }

        canonical_mode_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (!dynamic_cast<const DBGHashFast2Impl*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashFast2Impl*>(&other);

    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(canonical_mode_ == other_hash.canonical_mode_);
    assert(kmers_.size() == other_hash.kmers_.size());

    return kmers_ == other_hash.kmers_;
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::node_index
DBGHashFast2Impl<KMER>::get_index(const KmerIterator &kmer_iter) const {
    return ((node_index)(kmer_iter - kmers_.begin()) << KMER::kBitsPerChar) + 1;
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::node_index
DBGHashFast2Impl<KMER>::get_index(const Kmer &kmer) const {
    const auto find = find_kmer(KmerPrefix(kmer.data()) & kIgnoreLastCharMask);
    if (find == kmers_.end())
        return npos;

    return get_index(find);
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::KmerIndex::iterator
DBGHashFast2Impl<KMER>::get_iter(node_index node) {
    assert(in_graph(node));

    return (kmers_.begin() + ((node / KMER::kBitsPerChar) - 1));
}

template <typename KMER>
typename DBGHashFast2Impl<KMER>::KmerIndex::const_iterator
DBGHashFast2Impl<KMER>::get_iter(node_index node) const {
    assert(in_graph(node));

    return (kmers_.begin() + ((node / KMER::kBitsPerChar) - 1));
}

template <typename KMER>
const KMER DBGHashFast2Impl<KMER>
::get_kmer(node_index node,
           std::optional<typename KmerIndex::const_iterator> kmer_iter,
           bool *has_no_incoming) const {
    assert(in_graph(node));

    if (!kmer_iter)
        kmer_iter = get_iter(node);

    KMER prefix = kmer_iter->key();

    KmerPrefix c = node & ~kIgnoreLastCharMask;

    if (has_no_incoming)
        *has_no_incoming = ((kmer_iter->operator*()).second & kHasNoIncomingFlag) != 0u;

    return KMER(KmerPrefix(prefix.data()) | c);
}

template <typename KMER>
bool DBGHashFast2Impl<KMER>::in_graph(node_index node) const {
    //assert(node > 0 && node <= kmers_.size());
    std::ignore = node;
    return true;
}

template <typename KMER>
void DBGHashFast2Impl<KMER>::call_nodes(const std::function<void(node_index)> &callback,
                               const std::function<bool()> &stop_early) const {
    for (auto iter = kmers_.begin(); iter != kmers_.end() && !stop_early(); ++iter) {
        for (char c : seq_encoder_.alphabet) {
            if ((*iter).second & (Bits(1) << c))
                callback(get_index(iter) + c);
        }
    }
}

std::unique_ptr<DBGHashFast2::DBGHashFast2Interface>
DBGHashFast2::initialize_graph(size_t k,
                                 bool canonical_mode,
                                 bool packed_serialization) {
    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashFast2Impl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode, packed_serialization
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashFast2Impl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode, packed_serialization
        );
    } else {
        return std::make_unique<DBGHashFast2Impl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode, packed_serialization
        );
    }
}

DBGHashFast2::DBGHashFast2(size_t k,
                               bool canonical_mode,
                               bool packed_serialization) {
    hash_dbg_ = initialize_graph(k, canonical_mode, packed_serialization);
}

bool DBGHashFast2::load(std::istream &in) {
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

bool DBGHashFast2::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}
