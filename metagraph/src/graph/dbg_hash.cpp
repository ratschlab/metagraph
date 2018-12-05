#include "dbg_hash.hpp"

#include <cassert>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"


const std::string kAlphabet = "ACGTN";


void DBGHash::add_sequence(const std::string &sequence,
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
                nodes_inserted->insertBit(kmers_.size() - 1, true);
        }
    }
}

void DBGHash::map_to_nodes(const std::string &sequence,
                           const std::function<void(node_index)> &callback,
                           const std::function<bool()> &terminate) const {
    for (size_t i = 0; i + k_ - 1 < sequence.size(); ++i) {
        auto node = kmer_to_node(std::string(&sequence[i], &sequence[i + k_]));

        callback(node);

        if (terminate())
            return;
    }
}

DBGHash::node_index DBGHash::traverse(node_index node, char next_char) const {
    auto kmer = node_to_kmer(node).substr(1) + next_char;
    return kmer_to_node(kmer);
}

DBGHash::node_index DBGHash::traverse_back(node_index node, char prev_char) const {
    auto kmer = node_to_kmer(node);
    kmer.pop_back();
    return kmer_to_node(std::string(1, prev_char) + kmer);
}

DBGHash::node_index DBGHash::kmer_to_node(const std::string &kmer) const {
    if (kmer.length() != k_)
        throw std::runtime_error("Error: incompatible k-mer size");

    auto find = indices_.find(kmer);
    if (find == indices_.end())
        return npos;

    return find->second + 1;
}

std::string DBGHash::node_to_kmer(node_index i) const {
    assert(kmers_.at(i - 1).size() == k_);
    return std::string(kmers_.at(i - 1));
}

void DBGHash::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, kmers_.size());
    serialize_number(out, k_);
    serialize_string_number_map(out, indices_);
}

void DBGHash::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGHash::load(std::istream &in) {
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

bool DBGHash::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

std::string DBGHash::encode_sequence(const std::string &sequence) const {
    std::string result = sequence;

    for (char &c : result) {
        if (kAlphabet.find(c) == std::string::npos)
            c = 'N';
    }
    return result;
}
