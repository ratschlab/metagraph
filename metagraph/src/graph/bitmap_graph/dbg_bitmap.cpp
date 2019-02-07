#include "dbg_bitmap.hpp"

#include <cassert>
#include <cmath>

#include "serialization.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "dbg_bitmap_construct.hpp"


// Assume all k-mers present
DBGSD::DBGSD(size_t k, bool canonical_mode)
      : alphabet(seq_encoder_.alphabet),
        k_(k),
        canonical_mode_(canonical_mode),
        kmers_(std::pow(static_cast<long double>(alphabet.size()), k_) + 1, true) {
    assert(k > 1);
    assert(kmers_.num_set_bits() == kmers_.size());
    if (k * std::log2(alphabet.size()) >= 64) {
        std::cerr << "ERROR: Too large k!"
                  << " Maximum allowed k with this alphabet is "
                  << static_cast<int>(64. / std::log2(alphabet.size())) - 1 << std::endl;
        exit(1);
    }
}

DBGSD::DBGSD(DBGSDConstructor *builder) : DBGSD(2) {
    assert(builder);

    builder->build_graph(this);
    assert(kmers_[0]);
}

void DBGSD::add_sequence(const std::string &, bit_vector_dyn *) {
    assert(false);
}

void DBGSD::map_to_nodes(const std::string &sequence,
                         const std::function<void(node_index)> &callback,
                         const std::function<bool()> &terminate) const {
    for (const auto &kmer : sequence_to_kmers(sequence, canonical_mode_)) {
        auto index = to_index(kmer);
        callback(kmers_[index] ? index : npos);

        if (terminate())
            return;
    }
}

DBGSD::node_index
DBGSD::traverse(node_index node, char next_char) const {
    assert(node);

    auto kmer = to_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return get_node(kmer);
}

DBGSD::node_index
DBGSD::traverse_back(node_index node, char prev_char) const {
    assert(node);

    auto kmer = to_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return get_node(kmer);
}

std::vector<DBGSD::node_index>
DBGSD::adjacent_outgoing_nodes(node_index node) const {
    std::vector<node_index> indices;
    for (char c : alphabet) {
        auto next_index = traverse(node, c);
        if (next_index != npos)
            indices.emplace_back(std::move(next_index));
    }
    return indices;
}

std::vector<DBGSD::node_index>
DBGSD::adjacent_incoming_nodes(node_index node) const {
    std::vector<node_index> indices;
    for (char c : alphabet) {
        auto next_index = traverse_back(node, c);
        if (next_index != npos)
            indices.emplace_back(std::move(next_index));
    }
    return indices;
}

DBGSD::node_index
DBGSD::kmer_to_node(const std::string &kmer) const {
    assert(kmer.size() == k_);
    assert(sequence_to_kmers(kmer).size() == 1);
    return get_node(sequence_to_kmers(kmer)[0]);
}

std::string DBGSD::node_to_kmer(node_index i) const {
    assert(i != npos);
    return seq_encoder_.kmer_to_sequence(to_kmer(i), k_);
}

void DBGSD::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    kmers_.serialize(out);
    serialize_number(out, canonical_mode_);
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
        kmers_.load(in);
        if (!in.good())
            return false;

        try {
            canonical_mode_ = load_number(in);
            if (in.eof())
                canonical_mode_ = false;
        } catch (...) {
            canonical_mode_ = false;
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

typedef uint8_t TAlphabet;
struct Edge {
    DBGSD::node_index id;
    std::vector<TAlphabet> source_kmer;
};

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void DBGSD::call_paths(Call<const std::vector<node_index>,
                            const std::vector<TAlphabet>&> callback,
                       bool split_to_contigs) const {
    uint64_t nnodes = num_nodes();
    // keep track of reached edges
    std::vector<bool> discovered(nnodes, false);
    // keep track of edges that are already included in covering paths
    std::vector<bool> visited(nnodes, false);
    // store all branch nodes on the way
    std::deque<Edge> nodes;
    std::vector<uint64_t> path;

    auto node_to_index = [&](const node_index &id) -> uint64_t {
        assert(id > 0);
        assert(id < kmers_.size());
        return kmers_.rank1(id) - 2;
    };

    auto index_to_node = [&](uint64_t id) -> node_index {
        assert(id < nnodes);
        return kmers_.select1(id + 2);
    };

    // start at the source node
    for (uint64_t j = 0; j < nnodes; ++j) {
        uint64_t i = index_to_node(j);
        assert(node_to_index(i) == j);
        if (visited[j])
            continue;

        //TODO: traverse backwards

        discovered[j] = true;
        nodes.push_back({ j, seq_encoder_.encode(node_to_kmer(i)) });
        nodes.back().source_kmer.pop_back();

        // keep traversing until we have worked off all branches from the queue
        while (!nodes.empty()) {
            uint64_t node = nodes.front().id;
            auto sequence = std::move(nodes.front().source_kmer);
            path.clear();
            nodes.pop_front();

            // traverse simple path until we reach its tail or
            // the first edge that has been already visited
            while (!visited[node]) {
                assert(node < nnodes);
                assert(discovered[node]);
                sequence.push_back(seq_encoder_.encode(node_to_kmer(index_to_node(node)).back()));
                path.push_back(index_to_node(node));
                visited[node] = true;

                auto out_kmers = adjacent_outgoing_nodes(index_to_node(node));
                std::vector<TAlphabet> kmer(sequence.end() - k_ + 1, sequence.end());

                if (out_kmers.size() == 1) {
                    node = node_to_index(out_kmers[0]);
                    discovered[node] = true;
                    continue;
                } else if (out_kmers.size() > 1) {
                    bool continue_traversal = false;
                    for (const auto &next : out_kmers) {
                        auto next_index = node_to_index(next);
                        if (!discovered[next_index]) {
                            continue_traversal = true;
                            discovered[next_index] = true;
                            nodes.push_back({ next_index, kmer });
                        }
                    }

                    if (split_to_contigs)
                        break;

                    if (continue_traversal) {
                        node = nodes.back().id;
                        nodes.pop_back();
                        continue;
                    } else {
                        break;
                    }
                }
            }

            if (path.size())
                callback(path, sequence);
        }
    }
}

void DBGSD::call_sequences(Call<const std::string&> callback,
                           bool split_to_contigs) const {
    std::string sequence;

    call_paths([&](const auto&, const auto &path) {
        if (path.size()) {
            sequence.clear();
            std::transform(path.begin(), path.end(),
                           std::back_inserter(sequence),
                           [&](char c) {
                               return seq_encoder_.decode(c);
                           });
            callback(sequence);
        }
    }, split_to_contigs);
}

/**
 * Traverse graph and iterate over all nodes
 */
void DBGSD::call_kmers(Call<node_index, const std::string&> callback) const {
    uint64_t nnodes = num_nodes();
    for (size_t i = 1; i <= nnodes; ++i) {
        callback(i, node_to_kmer(kmers_.select1(i + 1)));
    }
}

Vector<DBGSD::Kmer> DBGSD::sequence_to_kmers(const std::string &sequence,
                                             bool to_canonical) const {
    Vector<Kmer> kmers;
    seq_encoder_.sequence_to_kmers<uint64_t>(sequence, k_, {}, &kmers, to_canonical);
    return kmers;
}

DBGSD::node_index DBGSD::get_node(const Kmer &kmer) const {
    assert(kmer == sequence_to_kmers(seq_encoder_.kmer_to_sequence(kmer, k_))[0]);

    if (!canonical_mode_) {
        auto index = to_index(kmer);
        return kmers_[index] ? index : npos;
    } else {
        auto str = seq_encoder_.kmer_to_sequence(kmer, k_);
        reverse_complement(str.begin(), str.end());
        auto rev_kmer = sequence_to_kmers(str)[0];
        return kmers_[to_index(std::min(kmer, rev_kmer))] ? to_index(kmer) : npos;
    }
}

DBGSD::node_index DBGSD::to_index(const Kmer &kmer) const {
    assert(kmer == sequence_to_kmers(seq_encoder_.kmer_to_sequence(kmer, k_))[0]);
    return Kmer::KMerWordType(1u) + kmer.data();
}

DBGSD::Kmer DBGSD::to_kmer(node_index node) const {
    assert(node != npos);
    assert(node < kmers_.size());
    return Kmer(node - static_cast<Kmer::KMerWordType>(1u));
}

bool DBGSD::equals(const DBGSD &other, bool verbose) const {
    if (verbose) {
        if (k_ != other.k_) {
            std::cerr << "k: " << k_ << " != " << other.k_ << std::endl;
            return false;
        }

        if (canonical_mode_ != other.canonical_mode_) {
            std::cerr << "canonical: " << canonical_mode_
                      << " != " << other.canonical_mode_ << std::endl;
            return false;
        }

        if (kmers_.num_set_bits() != other.kmers_.num_set_bits()) {
            std::cerr << "setbits: " << kmers_.num_set_bits()
                      << " != " << other.kmers_.num_set_bits() << std::endl;
            return false;
        }

        uint64_t cur_one = 1;
        uint64_t mismatch = 0;
        kmers_.call_ones(
            [&](const auto &pos) {
                mismatch += (pos != other.kmers_.select1(cur_one++));
            }
        );
        return !mismatch;
    }

    if (k_ == other.k_
            && canonical_mode_ == other.canonical_mode_
            && kmers_.num_set_bits() == other.kmers_.num_set_bits()) {
        uint64_t cur_one = 1;
        uint64_t mismatch = 0;
        kmers_.call_ones(
            [&](const auto &pos) {
                mismatch += (pos != other.kmers_.select1(cur_one++));
            }
        );
        return !mismatch;
    }

    return false;
}

std::ostream& operator<<(std::ostream &out, const DBGSD &graph) {
    out << "k: " << graph.k_ << std::endl
        << "canonical: " << graph.canonical_mode_ << std::endl
        << "nodes:" << std::endl;
    uint64_t nnodes = graph.num_nodes();
    for (size_t i = 1; i <= nnodes; ++i) {
        out << graph.kmers_.select1(i + 1) << "\t"
            << graph.node_to_kmer(graph.kmers_.select1(i + 1)) << std::endl;
    }
    return out;
}
