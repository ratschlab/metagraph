#include "dbg_hash_ordered.hpp"

#include <cassert>

#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"
#include "helpers.hpp"


const std::string kAlphabet = "ACGTN";


DBGHashOrdered::DBGHashOrdered(size_t k, bool canonical_only)
      : k_(k), canonical_only_(canonical_only) {}

void DBGHashOrdered::add_sequence(const std::string &sequence,
                                  bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == num_nodes() + 1);

    for (const auto &kmer : sequence_to_kmers(sequence)) {
        auto index_insert = kmers_.insert(kmer);
        if (index_insert.second && nodes_inserted) {
            nodes_inserted->insertBit(kmers_.size() - 1, true);
        }
    }
}

void DBGHashOrdered::map_to_nodes(const std::string &sequence,
                                  const std::function<void(node_index)> &callback,
                                  const std::function<bool()> &terminate) const {
    for (const auto &kmer : sequence_to_kmers(sequence)) {
        auto node = get_index(kmer);

        callback(node);

        if (terminate())
            return;
    }
}

DBGHashOrdered::node_index
DBGHashOrdered::traverse(node_index node, char next_char) const {
    auto kmer = get_kmer(node);
    kmer.next_kmer(k_ - 1, seq_encoder_.encode(next_char));
    return get_index(kmer);
}

DBGHashOrdered::node_index
DBGHashOrdered::traverse_back(node_index node, char prev_char) const {
    auto kmer = get_kmer(node);
    kmer = kmer.prev_kmer(k_ - 1, seq_encoder_.encode(prev_char));
    return get_index(kmer);
}

DBGHashOrdered::node_index
DBGHashOrdered::kmer_to_node(const std::string &kmer) const {
    if (kmer.length() != k_)
        throw std::runtime_error("Error: incompatible k-mer size");

    assert(sequence_to_kmers(kmer).size() == 1);

    return get_index(sequence_to_kmers(kmer)[0]);
}

std::string DBGHashOrdered::node_to_kmer(node_index i) const {
    assert(i > 0);
    assert(i <= kmers_.size());
    return seq_encoder_.kmer_to_sequence(get_kmer(i));
}

void DBGHashOrdered::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    serialize_number(out, kmers_.size());
    for (const auto &kmer : kmers_) {
        out.write(reinterpret_cast<const char *>(&kmer), sizeof(kmer));
    }
    serialize_number(out, canonical_only_);
}

void DBGHashOrdered::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGHashOrdered::load(std::istream &in) {
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

        canonical_only_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

bool DBGHashOrdered::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

Vector<DBGHashOrdered::Kmer>
DBGHashOrdered::sequence_to_kmers(const std::string &sequence) const {
    if (sequence.size() < k_)
        return {};

    Vector<Kmer> kmers;
    kmers.reserve(sequence.size() - k_ + 1);

    seq_encoder_.sequence_to_kmers(seq_encoder_.encode(sequence), k_, {}, &kmers);
    if (!canonical_only_)
        return kmers;

    std::string rev_compl = sequence;
    reverse_complement(rev_compl.begin(), rev_compl.end());

    Vector<Kmer> rev_kmers;
    kmers.reserve(sequence.size() - k_ + 1);

    seq_encoder_.sequence_to_kmers(
        seq_encoder_.encode(rev_compl), k_, {}, &rev_kmers
    );

    assert(kmers.size() == rev_kmers.size());

    for (size_t i = 0; i < kmers.size(); ++i) {
        kmers[i] = std::min(kmers[i], rev_kmers[rev_kmers.size() - 1 - i]);
    }

    return kmers;
}

DBGHashOrdered::node_index DBGHashOrdered::get_index(const Kmer &kmer) const {
    auto find = kmers_.find(kmer);
    if (find == kmers_.end())
        return npos;

    return find - kmers_.begin() + 1;
}

DBGHashOrdered::Kmer DBGHashOrdered::get_kmer(node_index node) const {
    assert(node > 0);
    assert(node <= kmers_.size());
    return *(kmers_.begin() + (node - 1));
}
