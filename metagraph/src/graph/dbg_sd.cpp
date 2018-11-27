#include "dbg_sd.hpp"

#include <cassert>

#include "serialization.hpp"
#include "utils.hpp"
#include "helpers.hpp"


const std::string kAlphabet = "ACGT";


// TODO: for now, assume all kmers present
DBGSD::DBGSD(size_t k, bool canonical_only)
      : k_(k),
        canonical_only_(canonical_only),
        kmers_((1llu << (k_ << 1)) + 1, true),
        offset_((1llu << (k_ << 1)) - 1) {}

void DBGSD::add_sequence(const std::string &,
                         bit_vector_dyn *) {}

void DBGSD::map_to_nodes(const std::string &sequence,
                         const std::function<void(node_index)> &callback,
                         const std::function<bool()> &terminate) const {
    for (const auto &kmer : sequence_to_kmers(sequence)) {
        auto node = get_index(kmer);

        callback(node);

        if (terminate())
            return;
    }
}

bool DBGSD::find(const std::string &sequence,
                 double kmer_discovery_fraction) const {
    if (sequence.length() < k_)
        return false;

    const size_t num_kmers = sequence.length() - k_ + 1;
    const size_t max_kmers_missing = num_kmers * (1 - kmer_discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    map_to_nodes(sequence,
        [&](node_index node) {
            if (node) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; }
    );
    return num_kmers_missing <= max_kmers_missing;
}

DBGSD::node_index
DBGSD::traverse(node_index node, char next_char) const {
    auto kmer = get_kmer(node);
    kmer.next_kmer(k_ - 1, seq_encoder_.encode(next_char));
    return get_index(kmer);
}

DBGSD::node_index
DBGSD::traverse_back(node_index node, char prev_char) const {
    auto kmer = get_kmer(node);
    kmer = kmer.prev_kmer(k_ - 1, seq_encoder_.encode(prev_char));
    return get_index(kmer);
}

DBGSD::node_index
DBGSD::kmer_to_node(const std::string &kmer) const {
    if (kmer.length() != k_)
        throw std::runtime_error("Error: incompatible k-mer size");

    assert(sequence_to_kmers(kmer).size() == 1);

    return get_index(sequence_to_kmers(kmer)[0]);
}

std::string DBGSD::node_to_kmer(node_index i) const {
    assert(i > 0);
    assert(i <= kmers_.size());
    return seq_encoder_.kmer_to_sequence(get_kmer(i));
}

void DBGSD::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    kmers_.serialize(out);
    serialize_number(out, canonical_only_);
}

void DBGSD::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGSD::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        k_ = load_number(in);
        offset_ = (1llu << (k_ << 1)) - 1;
        kmers_.load(in);
        if (!in.good())
            return false;

        try {
            canonical_only_ = load_number(in);
            if (in.eof())
                canonical_only_ = false;
        } catch (...) {
            canonical_only_ = false;
        }

        return true;
    } catch (...) {
        return false;
    }
}

bool DBGSD::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

Vector<DBGSD::Kmer>
DBGSD::sequence_to_kmers(const std::string &sequence) const {
    auto kmers = seq_encoder_.sequence_to_kmers(sequence, k_);
    if (!canonical_only_)
        return kmers;

    std::string rev_compl = sequence;
    reverse_complement(rev_compl.begin(), rev_compl.end());
    auto rev_kmers = seq_encoder_.sequence_to_kmers(rev_compl, k_);

    assert(kmers.size() == rev_kmers.size());

    for (size_t i = 0; i < kmers.size(); ++i) {
        kmers[i] = std::min(kmers[i], rev_kmers[rev_kmers.size() - 1 - i]);
    }

    return kmers;
}

DBGSD::node_index DBGSD::get_index(const Kmer &kmer) const {
    assert(kmer.data() - offset_ <= kmers_.size());
    return kmer.data() - offset_;
}

DBGSD::Kmer DBGSD::get_kmer(node_index node) const {
    assert(node > 0);
    return Kmer(node + offset_);
}
