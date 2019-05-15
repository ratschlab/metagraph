#include "dbg_hash_string.hpp"

#include <cassert>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"


void DBGHashString::add_sequence(const std::string &sequence,
                                 bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes());

    auto seq_encoded = encode_sequence(sequence);

    if (sequence.size() < k_)
        return;

    for (size_t i = 0; i + k_ - 1 < seq_encoded.size(); ++i) {
        auto index_insert = indices_.emplace(seq_encoded.substr(i, k_),
                                             kmers_.size());
        if (index_insert.second) {
            kmers_.push_back(index_insert.first->first);
            if (nodes_inserted)
                nodes_inserted->insert_bit(kmers_.size() - 1, true);
        }
    }
}

void DBGHashString::map_to_nodes(const std::string &sequence,
                                 const std::function<void(node_index)> &callback,
                                 const std::function<bool()> &terminate) const {
    for (size_t i = 0; i + k_ <= sequence.size(); ++i) {
        auto node = kmer_to_node(std::string(&sequence[i], &sequence[i + k_]));

        callback(node);

        if (terminate())
            return;
    }
}

void DBGHashString
::map_to_nodes_sequentially(std::string::const_iterator begin,
                            std::string::const_iterator end,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    for (auto it = begin; it + k_ <= end; ++it) {
        auto node = kmer_to_node(std::string(it, it + k_));

        callback(node);

        if (terminate())
            return;
    }
}

void DBGHashString::extend_from_seed(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate,
                                           node_index seed) const {
    auto it = begin;
    if (it + k_ > end)
        return;
    if (seed == npos) {
        seed = kmer_to_node(std::string(it, it + k_));
        if (seed == npos || terminate())
            return;
        callback(seed);
        it ++;
    }
    // If the seed includes some mis-alignments, return.
    else if (node_to_kmer(seed).substr(1) != std::string(it, it + k_ - 1)) {
        return;
    }
    while (it + k_ <= end) {
        auto node = kmer_to_node(std::string(it, it + k_));
        if (node == npos || terminate())
            return;
        callback(node);
        it ++;
    }
}

DBGHashString::node_index DBGHashString::traverse(node_index node, char next_char) const {
    assert(node);
    auto kmer = node_to_kmer(node).substr(1) + next_char;
    return kmer_to_node(kmer);
}

DBGHashString::node_index
DBGHashString::traverse_back(node_index node, char prev_char) const {
    assert(node);
    auto kmer = node_to_kmer(node);
    kmer.pop_back();
    return kmer_to_node(std::string(1, prev_char) + kmer);
}

void DBGHashString
::adjacent_outgoing_nodes(node_index node,
                          std::vector<node_index> *target_nodes) const {
    assert(node);
    assert(target_nodes);

    auto prefix = node_to_kmer(node).substr(1);

    for (char c : alphabet_) {
        auto next = kmer_to_node(prefix + c);
        if (next != npos)
            target_nodes->push_back(next);
    }
}

void DBGHashString
::adjacent_incoming_nodes(node_index node,
                          std::vector<node_index> *source_nodes) const {
    assert(node);
    assert(source_nodes);

    auto suffix = node_to_kmer(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto next = kmer_to_node(std::string(1, c) + suffix);
        if (next != npos)
            source_nodes->push_back(next);
    }
}

void DBGHashString
::call_outgoing_kmers(node_index node,
                      const OutgoingEdgeCallback &callback) const {
    assert(node);

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
    assert(node);

    auto suffix = node_to_kmer(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto next = kmer_to_node(std::string(1, c) + suffix);
        if (next != npos)
            callback(next, c);
    }
}

size_t DBGHashString::outdegree(node_index node) const {
    assert(node);

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

size_t DBGHashString::indegree(node_index node) const {
    assert(node);

    size_t indegree = 0;

    auto next = node_to_kmer(node);
    next.pop_back();
    next.insert(0, 1, '\0');

    for (char c : alphabet_) {
        next.front() = c;

        if (kmer_to_node(next) != npos)
            indegree++;
    }
    return indegree;
}

void DBGHashString
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    for (const auto &kmer : indices_) {
        callback(kmer.second + 1, kmer.first);
    }
}

DBGHashString::node_index
DBGHashString::kmer_to_node(const std::string &kmer) const {
    if (kmer.length() != k_)
        throw std::runtime_error("Error: incompatible k-mer size");

    auto find = indices_.find(kmer);
    if (find == indices_.end())
        return npos;

    return find->second + 1;
}

std::string DBGHashString::node_to_kmer(node_index node) const {
    assert(node);
    assert(kmers_.at(node - 1).size() == k_);
    return std::string(kmers_.at(node - 1));
}

void DBGHashString::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, kmers_.size());
    serialize_number(out, k_);
    serialize_string_number_map(out, indices_);
}

void DBGHashString::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGHashString::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        size_t size = load_number(in);
        k_ = load_number(in);
        kmers_.resize(size);
        load_string_number_map(in, &indices_);
        //kmers_.resize(indices_.size());
        for (auto &kmer : indices_) {
            kmers_[kmer.second] = kmer.first;
        }
        return true;
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
    if (!dynamic_cast<const DBGHashString*>(&other)) {
        throw std::runtime_error("Not implemented");
        return false;
    }

    const auto& other_hash = *dynamic_cast<const DBGHashString*>(&other);
    if (k_ != other_hash.k_)
        return false;

    for (const auto &kmer : kmers_) {
        if (other_hash.indices_.find(kmer) == other_hash.indices_.end())
            return false;
    }

    return true;
}

std::string DBGHashString::encode_sequence(const std::string &sequence) const {
    std::string result = sequence;

    for (char &c : result) {
        if (alphabet_.find(c) == std::string::npos)
            c = 'N';
    }
    return result;
}
