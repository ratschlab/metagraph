#include "dbg_hash_ordered.hpp"

#include <cassert>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "hash_utils.hpp"


template <typename KMER = KmerExtractor2Bit::Kmer64>
class DBGHashOrderedImpl : public DBGHashOrdered::DBGHashOrderedInterface {
    using Kmer = KMER;
    using KmerIndex = tsl::ordered_set<Kmer,
                                       utils::Hash<Kmer>,
                                       std::equal_to<Kmer>,
                                       std::allocator<Kmer>,
                                       std::deque<Kmer, std::allocator<Kmer>>,
                                       std::uint64_t>;
  public:
    explicit DBGHashOrderedImpl(size_t k,
                                bool canonical_mode,
                                bool packed_serialization);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted) {
        add_sequence(sequence, [](){ return false; }, nodes_inserted);
    }

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    // `skip` is called before adding each k-mer into the graph and the k-mer
    // is skipped if `skip()` returns `false`.
    void add_sequence(const std::string &sequence,
                      const std::function<bool()> &skip,
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

    uint64_t num_nodes() const { return kmers_.size(); }

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

    node_index get_index(const Kmer &kmer) const;
    const Kmer& get_kmer(node_index node) const;

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;
    KmerExtractor2Bit seq_encoder_;

    bool packed_serialization_;

    static constexpr auto kExtension = DBGHashOrdered::kExtension;
};

template <typename KMER>
DBGHashOrderedImpl<KMER>::DBGHashOrderedImpl(size_t k,
                                             bool canonical_mode,
                                             bool packed_serialization)
      : k_(k),
        canonical_mode_(canonical_mode),
        packed_serialization_(packed_serialization) {}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::add_sequence(const std::string &sequence,
                                            const std::function<bool()> &skip,
                                            bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    if (sequence.size() < get_k())
        return;

    std::vector<bool> skipped;
    skipped.reserve(sequence.size() - get_k() + 1);

    for (const auto &kmer : sequence_to_kmers(sequence)) {
        skipped.push_back(skip());
        if (skipped.back())
            continue;

        auto index_insert = kmers_.insert(kmer);

        if (index_insert.second && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);
    }

    if (!canonical_mode_)
        return;

    auto it = skipped.end();
    for (const auto &kmer : sequence_to_kmers(seq_encoder_.reverse_complement(sequence))) {
        if (*(--it))
            continue;

        auto index_insert = kmers_.insert(kmer);

        if (index_insert.second && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
template <typename KMER>
void DBGHashOrderedImpl<KMER>::map_to_nodes_sequentially(
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
void DBGHashOrderedImpl<KMER>::map_to_nodes(const std::string &sequence,
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
void DBGHashOrderedImpl<KMER>::call_outgoing_kmers(node_index node,
                                                   const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        auto next = get_index(next_kmer);
        if (next != npos)
            callback(next, c);
    }
}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::call_incoming_kmers(node_index node,
                                                   const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto prev = get_index(prev_kmer);
        if (prev != npos)
            callback(prev, c);
    }
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::traverse(node_index node, char next_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return get_index(kmer);
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return get_index(kmer);
}

template <typename KMER>
void
DBGHashOrderedImpl<KMER>
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

template <typename KMER>
void
DBGHashOrderedImpl<KMER>
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

template <typename KMER>
size_t DBGHashOrderedImpl<KMER>::outdegree(node_index node) const {
    assert(in_graph(node));

    size_t outdegree = 0;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (get_index(next_kmer) != npos)
            outdegree++;
    }

    return outdegree;
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    bool outgoing_edge_detected = false;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (get_index(next_kmer) != npos) {
            if (outgoing_edge_detected)
                return false;

            outgoing_edge_detected = true;
        }
    }

    return outgoing_edge_detected;
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    bool outgoing_edge_detected = false;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (get_index(next_kmer) != npos) {
            if (outgoing_edge_detected)
                return true;

            outgoing_edge_detected = true;
        }
    }

    return false;
}

template <typename KMER>
size_t DBGHashOrderedImpl<KMER>::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (get_index(prev_kmer) != npos)
            indegree++;
    }

    return indegree;
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::has_no_incoming(node_index node) const {
    assert(in_graph(node));

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (get_index(prev_kmer) != npos)
            return false;
    }

    return true;
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    bool incoming_edge_detected = false;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (get_index(prev_kmer) != npos) {
            if (incoming_edge_detected)
                return false;

            incoming_edge_detected = true;
        }
    }

    return incoming_edge_detected;
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    return get_index(seq_encoder_.encode(kmer));
}

template <typename KMER>
std::string DBGHashOrderedImpl<KMER>::get_node_sequence(node_index node) const {
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
void DBGHashOrderedImpl<KMER>::serialize(std::ostream &out) const {
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

    serialize_number(out, canonical_mode_);
}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::load(std::istream &in) {
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
                kmers_.insert(deserializer.operator()<Kmer>());
            }

        } else {
            packed_serialization_ = false;

            kmers_ = KmerIndex::deserialize(deserializer, true);
        }

        canonical_mode_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (!dynamic_cast<const DBGHashOrderedImpl*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashOrderedImpl*>(&other);

    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(canonical_mode_ == other_hash.canonical_mode_);
    assert(kmers_.size() == other_hash.kmers_.size());

    return kmers_ == other_hash.kmers_;
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::get_index(const Kmer &kmer) const {
    auto find = kmers_.find(kmer);
    if (find == kmers_.end())
        return npos;

    return find - kmers_.begin() + 1;
}

template <typename KMER>
const KMER& DBGHashOrderedImpl<KMER>::get_kmer(node_index node) const {
    assert(in_graph(node));
    assert(node == get_index(*(kmers_.nth(node - 1))));

    return *(kmers_.nth(node - 1));
}

template <typename KMER>
bool DBGHashOrderedImpl<KMER>::in_graph(node_index node) const {
    assert(node > 0 && node <= kmers_.size());
    std::ignore = node;
    return true;
}


std::unique_ptr<DBGHashOrdered::DBGHashOrderedInterface>
DBGHashOrdered::initialize_graph(size_t k,
                                 bool canonical_mode,
                                 bool packed_serialization) {
    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode, packed_serialization
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode, packed_serialization
        );
    } else {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode, packed_serialization
        );
    }
}

DBGHashOrdered::DBGHashOrdered(size_t k,
                               bool canonical_mode,
                               bool packed_serialization) {
    hash_dbg_ = initialize_graph(k, canonical_mode, packed_serialization);
}

bool DBGHashOrdered::load(std::istream &in) {
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

bool DBGHashOrdered::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}
