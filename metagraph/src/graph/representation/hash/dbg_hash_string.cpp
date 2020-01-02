#include "dbg_hash_string.hpp"

#include <cassert>

#include "common/serialization.hpp"
#include "common/vectors/bit_vector.hpp"
#include "alphabets.hpp"


#if _PROTEIN_GRAPH
    const std::string DBGHashString::alphabet_ = alphabets::kAlphabetProtein;
#elif _DNA_CASE_SENSITIVE_GRAPH
    const std::string DBGHashString::alphabet_ = alphabets::kAlphabetDNACaseSent;
#elif _DNA5_GRAPH
    const std::string DBGHashString::alphabet_ = alphabets::kAlphabetDNA5;
#elif _DNA_GRAPH
    const std::string DBGHashString::alphabet_ = alphabets::kAlphabetDNA;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

void DBGHashString::add_sequence(const std::string &sequence,
                                 bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    for (const auto &seq_encoded : encode_sequence(sequence)) {
        assert(sequence.size() >= k_);

        for (size_t i = 0; i + k_ - 1 < seq_encoded.size(); ++i) {
            auto index_insert = kmers_.insert(seq_encoded.substr(i, k_));

            if (index_insert.second && nodes_inserted)
                nodes_inserted->insert_bit(kmers_.size() - 1, true);
        }
    }
}

void DBGHashString::map_to_nodes(const std::string &sequence,
                                 const std::function<void(node_index)> &callback,
                                 const std::function<bool()> &terminate) const {
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
        callback(kmer_to_node(std::string(sequence.data() + i,
                                          sequence.data() + i + k_)));
    }
}

void DBGHashString
::map_to_nodes_sequentially(std::string::const_iterator begin,
                            std::string::const_iterator end,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    for (auto it = begin; it + k_ <= end && !terminate(); ++it) {
        callback(kmer_to_node(std::string(it, it + k_)));
    }
}

DBGHashString::node_index
DBGHashString::traverse(node_index node, char next_char) const {
    assert(node);
    auto kmer = node_to_kmer(node).substr(1) + next_char;
    return kmer_to_node(kmer);
}

DBGHashString::node_index
DBGHashString::traverse_back(node_index node, char prev_char) const {
    assert(node > 0 && node <= num_nodes());
    auto kmer = node_to_kmer(node);
    kmer.pop_back();
    return kmer_to_node(std::string(1, prev_char) + kmer);
}

void DBGHashString
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

void DBGHashString
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGHashString
::call_outgoing_kmers(node_index node,
                      const OutgoingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    auto prefix = node_to_kmer(node).substr(1);

    for (char c : alphabet_) {
        auto next = kmer_to_node(prefix + c);
        if (next != npos)
            callback(next, c);
    }
}

void DBGHashString
::call_incoming_kmers(node_index node,
                      const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    auto suffix = node_to_kmer(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto prev = kmer_to_node(std::string(1, c) + suffix);
        if (prev != npos)
            callback(prev, c);
    }
}

size_t DBGHashString::outdegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    size_t outdegree = 0;

    auto next = node_to_kmer(node).substr(1);
    next.push_back('\0');

    for (char c : alphabet_) {
        next.back() = c;

        if (kmer_to_node(next) != npos)
            outdegree++;
    }

    return outdegree;
}

bool DBGHashString::has_single_outgoing(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool outgoing_edge_detected = false;

    auto next = node_to_kmer(node).substr(1);
    next.push_back('\0');

    for (char c : alphabet_) {
        next.back() = c;

        if (kmer_to_node(next) != npos) {
            if (outgoing_edge_detected)
                return false;

            outgoing_edge_detected = true;
        }
    }

    return outgoing_edge_detected;
}

bool DBGHashString::has_multiple_outgoing(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool outgoing_edge_detected = false;

    auto next = node_to_kmer(node).substr(1);
    next.push_back('\0');

    for (char c : alphabet_) {
        next.back() = c;

        if (kmer_to_node(next) != npos) {
            if (outgoing_edge_detected)
                return true;

            outgoing_edge_detected = true;
        }
    }

    return false;
}

size_t DBGHashString::indegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    size_t indegree = 0;

    auto prev = node_to_kmer(node);
    prev.pop_back();
    prev.insert(0, 1, '\0');

    for (char c : alphabet_) {
        prev.front() = c;

        if (kmer_to_node(prev) != npos)
            indegree++;
    }

    return indegree;
}

bool DBGHashString::has_no_incoming(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    auto prev = node_to_kmer(node);
    prev.pop_back();
    prev.insert(0, 1, '\0');

    for (char c : alphabet_) {
        prev.front() = c;

        if (kmer_to_node(prev) != npos)
            return false;
    }

    return true;
}

bool DBGHashString::has_single_incoming(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool incoming_edge_detected = false;

    auto prev = node_to_kmer(node);
    prev.pop_back();
    prev.insert(0, 1, '\0');

    for (char c : alphabet_) {
        prev.front() = c;

        if (kmer_to_node(prev) != npos) {
            if (incoming_edge_detected)
                return false;

            incoming_edge_detected = true;
        }
    }

    return incoming_edge_detected;
}

void DBGHashString
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    for (auto it = kmers_.begin(); it != kmers_.end(); ++it) {
        callback(it - kmers_.begin() + 1, *it);
    }
}

DBGHashString::node_index
DBGHashString::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    auto find = kmers_.find(kmer);
    if (find == kmers_.end())
        return npos;

    return find - kmers_.begin() + 1;
}

std::string DBGHashString::node_to_kmer(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    assert((kmers_.begin() + (node - 1))->length() == k_);
    return *(kmers_.begin() + (node - 1));
}

class KmerSerializer {
  public:
    KmerSerializer(std::ostream &os, size_t k) : os_(os), k_(k) {}

    void operator()(const std::string &str) {
        assert(str.size() == k_);
        os_.write(str.data(), k_);
    }

    template <class T>
    void operator()(const T &value) {
        os_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

  private:
    std::ostream &os_;
    size_t k_;
};

class KmerDeserializer {
  public:
    KmerDeserializer(std::istream &is, size_t k) : is_(is), k_(k) {}

    template <class T>
    T operator()() {
        T value;
        if constexpr(std::is_same<T, std::string>::value) {
            value.resize(k_);
            is_.read(value.data(), k_);
        } else {
            is_.read(reinterpret_cast<char *>(&value), sizeof(T));
        }
        return value;
    }

  private:
    std::istream &is_;
    size_t k_;
};

void DBGHashString::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    out.exceptions(out.badbit | out.failbit);

    serialize_number(out, k_);
    KmerSerializer serializer(out, k_);
    kmers_.serialize(serializer);
}

void DBGHashString::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGHashString::load(std::istream &in) {
    if (!in.good())
        return false;

    in.exceptions(in.badbit | in.failbit | in.eofbit);

    try {
        k_ = load_number(in);
        KmerDeserializer deserializer(in, k_);
        kmers_ = KmerIndex::deserialize(deserializer, true);
        return in.good();
    } catch (...) {
        return false;
    }
}

bool DBGHashString::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

bool DBGHashString::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || num_nodes() != other.num_nodes()
            || is_canonical_mode() != other.is_canonical_mode())
        return false;

    if (!dynamic_cast<const DBGHashString*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashString*>(&other);
    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(kmers_.size() == other_hash.kmers_.size());
    assert(is_canonical_mode() == other_hash.is_canonical_mode());

    return kmers_ == other_hash.kmers_;
}

std::vector<std::string> DBGHashString::encode_sequence(const std::string &sequence) const {
    std::vector<std::string> results;

    #if _DNA_GRAPH
        auto it = sequence.begin();
        while (it != sequence.end()) {
            auto jt = std::find_if(it, sequence.end(),
                                   [&](char c) {
                                       return alphabet_.find(c) == std::string::npos;
                                   });
            if (jt - it >= static_cast<int>(k_))
                results.emplace_back(it, jt);

            if (jt == sequence.end())
                break;

            it = std::next(jt);
        }
    #elif _PROTEIN_GRAPH
        results.push_back(sequence);

        std::replace_if(results[0].begin(), results[0].end(),
                        [&](char c) {
                            return alphabet_.find(c) == std::string::npos;
                        },
                        'X');

        if (results[0].size() < k_)
            return {};
    #else
        results.push_back(sequence);

        std::replace_if(results[0].begin(), results[0].end(),
                        [&](char c) {
                            return alphabet_.find(c) == std::string::npos;
                        },
                        'N');

        if (results[0].size() < k_)
            return {};
    #endif

    return results;
}

const std::string& DBGHashString::alphabet() const { return alphabet_; }
