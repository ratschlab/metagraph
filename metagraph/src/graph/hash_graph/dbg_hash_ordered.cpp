#include "dbg_hash_ordered.hpp"

#include <cassert>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"


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
    explicit DBGHashOrderedImpl(size_t k, bool canonical_mode = false);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const;

    size_t outdegree(node_index) const;
    size_t indegree(node_index) const;

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

    static constexpr auto kExtension = DBGHashOrdered::kExtension;
};

template <typename KMER>
DBGHashOrderedImpl<KMER>::DBGHashOrderedImpl(size_t k, bool canonical_mode)
      : k_(k), canonical_mode_(canonical_mode) {}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::add_sequence(const std::string &sequence,
                                            bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    for (const auto &kmer : sequence_to_kmers(sequence)) {
        auto index_insert = kmers_.insert(kmer);

        if (index_insert.second && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);
    }

    if (!canonical_mode_)
        return;

    for (const auto &kmer : sequence_to_kmers(seq_encoder_.reverse_complement(sequence))) {
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
    auto kmer = get_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return get_index(kmer);
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::traverse_back(node_index node, char prev_char) const {
    auto kmer = get_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return get_index(kmer);
}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::adjacent_outgoing_nodes(node_index node,
                                                       std::vector<node_index> *target_nodes) const {
    assert(target_nodes);

    call_outgoing_kmers(node, [&](auto i, char) { target_nodes->push_back(i); });
}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::adjacent_incoming_nodes(node_index node,
                                                       std::vector<node_index> *source_nodes) const {
    assert(source_nodes);

    call_incoming_kmers(node, [&](auto i, char) { source_nodes->push_back(i); });
}

template <typename KMER>
size_t DBGHashOrderedImpl<KMER>::outdegree(node_index node) const {
    assert(node);

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
size_t DBGHashOrderedImpl<KMER>::indegree(node_index node) const {
    assert(node);

    size_t indegree = 0;

    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto next_kmer = kmer;
        next_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (get_index(next_kmer) != npos)
            indegree++;
    }

    return indegree;
}

template <typename KMER>
typename DBGHashOrderedImpl<KMER>::node_index
DBGHashOrderedImpl<KMER>::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    return get_index(seq_encoder_.encode(kmer));
}

template <typename KMER>
std::string DBGHashOrderedImpl<KMER>::get_node_sequence(node_index node) const {
    assert(node > 0);
    assert(node <= kmers_.size());

    return seq_encoder_.kmer_to_sequence(get_kmer(node), k_);
}

template <typename KMER>
void DBGHashOrderedImpl<KMER>::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    serialize_number(out, kmers_.size());
    for (const auto &kmer : kmers_) {
        out.write(reinterpret_cast<const char *>(&kmer), sizeof(kmer));
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

    kmers_.clear();

    try {
        k_ = load_number(in);
        const auto size = load_number(in);
        kmers_.reserve(size + 1);
        Kmer kmer;
        for (uint64_t i = 0; i < size; ++i) {
            in.read(reinterpret_cast<char *>(&kmer), sizeof(kmer));
            kmers_.insert(kmer);
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
    assert(node > 0);
    assert(node <= kmers_.size());
    assert(node == get_index(*(kmers_.nth(node - 1))));

    return *(kmers_.nth(node - 1));
}


std::unique_ptr<DBGHashOrdered::DBGHashOrderedInterface>
DBGHashOrdered::initialize_graph(size_t k, bool canonical_mode) {
    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode
        );
    } else {
        return std::make_unique<DBGHashOrderedImpl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode
        );
    }
}

DBGHashOrdered::DBGHashOrdered(size_t k, bool canonical_mode) {
    hash_dbg_ = initialize_graph(k, canonical_mode);
}

bool DBGHashOrdered::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        auto pos = in.tellg();
        auto k = load_number(in);
        in.seekg(pos, in.beg);

        // the actual value of |canonical| will be set in load
        hash_dbg_ = initialize_graph(k, false);
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
