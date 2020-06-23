#include "dbg_bitmap.hpp"

#include <cassert>
#include <fstream>
#include <cmath>

#include "common/serialization.hpp"
#include "dbg_bitmap_construct.hpp"

namespace mtg {
namespace bitmap_graph {

// Assume all k-mers present
DBGBitmap::DBGBitmap(size_t k, bool canonical_mode)
      : k_(k),
        canonical_mode_(canonical_mode),
        seq_encoder_(),
        kmers_(std::pow(static_cast<long double>(seq_encoder_.alphabet.size()), k_) + 1, true),
        complete_(true) {
    assert(k > 1);
    assert(kmers_.num_set_bits() == kmers_.size());
    if (k * std::log2(alphabet().size()) >= 64) {
        std::cerr << "ERROR: Too large k!"
                  << " Maximum allowed k with this alphabet is "
                  << static_cast<int>(64. / std::log2(alphabet().size())) - 1 << std::endl;
        exit(1);
    }
}

DBGBitmap::DBGBitmap(DBGBitmapConstructor *builder) : DBGBitmap(2) {
    assert(builder);

    builder->build_graph(this);
    assert(kmers_[0]);
}


void DBGBitmap::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers(sequence, canonical_mode_)) {
        if (terminate())
            return;

        callback(is_valid ? to_node(kmer) : npos);
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
void DBGBitmap::map_to_nodes_sequentially(std::string_view sequence,
                                          const std::function<void(node_index)> &callback,
                                          const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers(sequence)) {
        if (terminate())
            return;

        callback(is_valid ? to_node(kmer) : npos);
    }
}

DBGBitmap::node_index
DBGBitmap::traverse(node_index node, char next_char) const {
    assert(node > 0 && node <= num_nodes());

    auto kmer = node_to_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return to_node(kmer);
}

DBGBitmap::node_index
DBGBitmap::traverse_back(node_index node, char prev_char) const {
    assert(node > 0 && node <= num_nodes());

    auto kmer = node_to_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return to_node(kmer);
}

void DBGBitmap::call_outgoing_kmers(node_index node,
                                    const OutgoingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        auto next_index = to_node(next_kmer);
        if (next_index != npos)
            callback(next_index, c);
    }
}

size_t DBGBitmap::outdegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    size_t outdegree = 0;

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (to_node(next_kmer) != npos)
            outdegree++;
    }

    return outdegree;
}

bool DBGBitmap::has_single_outgoing(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool outgoing_edge_detected = false;

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (to_node(next_kmer) != npos) {
            if (outgoing_edge_detected)
                return false;

            outgoing_edge_detected = true;
        }
    }

    return outgoing_edge_detected;
}

bool DBGBitmap::has_multiple_outgoing(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool outgoing_edge_detected = false;

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        if (to_node(next_kmer) != npos) {
            if (outgoing_edge_detected)
                return true;

            outgoing_edge_detected = true;
        }
    }

    return false;
}

void DBGBitmap::call_incoming_kmers(node_index node,
                                    const OutgoingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto next_index = to_node(prev_kmer);
        if (next_index != npos)
            callback(next_index, c);
    }
}

size_t DBGBitmap::indegree(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    size_t indegree = 0;

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (to_node(prev_kmer) != npos)
            indegree++;
    }

    return indegree;
}

bool DBGBitmap::has_no_incoming(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (to_node(prev_kmer) != npos)
            return false;
    }

    return true;
}

bool DBGBitmap::has_single_incoming(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    bool incoming_edge_detected = false;

    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        if (to_node(prev_kmer) != npos) {
            if (incoming_edge_detected)
                return false;

            incoming_edge_detected = true;
        }
    }

    return incoming_edge_detected;
}

void DBGBitmap::adjacent_outgoing_nodes(node_index node,
                                        const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    call_outgoing_kmers(node, [&](node_index child, char) { callback(child); });
}

void DBGBitmap::adjacent_incoming_nodes(node_index node,
                                        const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());

    call_incoming_kmers(node, [&](node_index parent, char) { callback(parent); });
}

DBGBitmap::node_index DBGBitmap::to_node(const Kmer &kmer) const {
    auto index = kmer.data() + 1;
    assert(index < kmers_.size());
    assert(!complete_ || kmers_[index]);

    return complete_
        ? index
        : (kmers_[index] ? kmers_.rank1(index) - 1 : npos);
}

DBGBitmap::node_index DBGBitmap::kmer_to_node(std::string_view kmer) const {
    assert(kmer.size() == k_);
    return to_node(Kmer(seq_encoder_.encode(kmer)));
}

uint64_t DBGBitmap::node_to_index(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    return complete_ ? node : kmers_.select1(node + 1);
}

DBGBitmap::Kmer DBGBitmap::node_to_kmer(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    return Kmer { complete_ ? node - 1 : kmers_.select1(node + 1) - 1 };
}

std::string DBGBitmap::get_node_sequence(node_index node) const {
    assert(node > 0 && node <= num_nodes());
    assert(sequence_to_kmers(seq_encoder_.kmer_to_sequence(
        node_to_kmer(node), k_)).size() == 1);
    assert(node == to_node(sequence_to_kmers(seq_encoder_.kmer_to_sequence(
        node_to_kmer(node), k_))[0].first));

    return seq_encoder_.kmer_to_sequence(node_to_kmer(node), k_);
}

uint64_t DBGBitmap::num_nodes() const {
    assert(kmers_[0] && "The first bit must be always set to 1");
    return kmers_.num_set_bits() - 1;
}

void DBGBitmap::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    kmers_.serialize(out);
    serialize_number(out, canonical_mode_);
}

void DBGBitmap::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

bool DBGBitmap::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        k_ = load_number(in);

        auto pos = in.tellg();

        if (!kmers_.load(in)) {
            kmers_ = decltype(kmers_)();
            // backward compatibility for loading bit_vector_sd
            in.seekg(pos, in.beg);
            bit_vector_sd temp_vector;
            if (!temp_vector.load(in))
                return false;

            kmers_ = temp_vector.convert_to<bit_vector_smart>();
        }

        if (!in.good())
            return false;

        complete_ = (kmers_.size() == kmers_.num_set_bits());

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

bool DBGBitmap::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

Vector<std::pair<DBGBitmap::Kmer, bool>>
DBGBitmap::sequence_to_kmers(std::string_view sequence, bool to_canonical) const {
    return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, to_canonical);
}

bool DBGBitmap::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (dynamic_cast<const DBGBitmap*>(&other))
        return equals(*dynamic_cast<const DBGBitmap*>(&other), false);

    throw std::runtime_error("Not implemented");
}

bool DBGBitmap::equals(const DBGBitmap &other, bool verbose) const {
    if (!verbose) {
        return k_ == other.k_
                && canonical_mode_ == other.canonical_mode_
                && kmers_ == other.kmers_;
    }

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

    return kmers_ == other.kmers_;
}

void DBGBitmap::print(std::ostream &out) const {
    if (complete_) {
        out << "Complete graph" << std::endl;
        return;
    }

    std::string vertex_header = "Vertex";
    vertex_header.resize(get_k(), ' ');

    out << "Index"
        << "\t" << "Repr"
        << "\t" << vertex_header
        << std::endl;

    uint64_t nnodes = num_nodes();
    for (size_t node = 1; node <= nnodes; ++node) {
        out << node << "\t"
            << node_to_index(node) << "\t"
            << get_node_sequence(node) << std::endl;
    }
}

} // namespace bitmap_graph
} // namespace mtg
