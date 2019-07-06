#include "dbg_bitmap.hpp"

#include <cassert>
#include <cmath>

#include "serialization.hpp"
#include "utils.hpp"
#include "dbg_bitmap_construct.hpp"


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


void DBGBitmap::map_to_nodes(const std::string &sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    auto kmers = sequence_to_kmers(sequence, canonical_mode_);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {
        if (terminate())
            return;

        if (is_valid) {
            assert(it != kmers.end());
            callback(to_node(*it++));
        } else {
            callback(npos);
        }
    }
    assert(it == kmers.end());
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
void DBGBitmap::map_to_nodes_sequentially(std::string::const_iterator begin,
                                          std::string::const_iterator end,
                                          const std::function<void(node_index)> &callback,
                                          const std::function<bool()> &terminate) const {
    std::string sequence(begin, end);
    auto kmers = sequence_to_kmers(sequence);
    auto it = kmers.begin();
    for (bool is_valid : seq_encoder_.valid_kmers(sequence, k_)) {
        if (terminate())
            return;

        if (is_valid) {
            assert(it != kmers.end());
            callback(to_node(*it++));
        } else {
            callback(npos);
        }
    }
    assert(it == kmers.end());
}

DBGBitmap::node_index
DBGBitmap::traverse(node_index node, char next_char) const {
    assert(node);

    auto kmer = node_to_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return to_node(kmer);
}

DBGBitmap::node_index
DBGBitmap::traverse_back(node_index node, char prev_char) const {
    assert(node);

    auto kmer = node_to_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return to_node(kmer);
}

void DBGBitmap::call_outgoing_kmers(node_index node,
                                    const OutgoingEdgeCallback &callback) const {
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
    size_t outdegree = 0;
    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_next(k_, seq_encoder_.encode(c));

        auto next_index = to_node(next_kmer);
        if (next_index != npos)
            outdegree++;
    }

    return outdegree;
}

void DBGBitmap::call_incoming_kmers(node_index node,
                                    const OutgoingEdgeCallback &callback) const {
    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto next_index = to_node(next_kmer);
        if (next_index != npos)
            callback(next_index, c);
    }
}

size_t DBGBitmap::indegree(node_index node) const {
    size_t indegree = 0;
    const auto &kmer = node_to_kmer(node);

    for (char c : alphabet()) {
        auto next_kmer = kmer;
        next_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto next_index = to_node(next_kmer);
        if (next_index != npos)
            indegree++;
    }

    return indegree;
}

void DBGBitmap::adjacent_outgoing_nodes(node_index node,
                                        std::vector<node_index> *target_nodes) const {
    assert(target_nodes);

    call_outgoing_kmers(node, [target_nodes](node_index target, char) {
        target_nodes->push_back(target);
    });
}

void DBGBitmap::adjacent_incoming_nodes(node_index node,
                                        std::vector<node_index> *source_nodes) const {
    assert(source_nodes);

    call_incoming_kmers(node, [source_nodes](node_index source, char) {
        source_nodes->push_back(source);
    });
}

DBGBitmap::node_index
DBGBitmap::to_node(const Kmer &kmer) const {
    auto index = kmer.data() + 1;
    assert(index < kmers_.size());
    assert(!complete_ || kmers_[index]);

    return complete_
        ? index
        : (kmers_[index] ? kmers_.rank1(index) - 1 : npos);
}

DBGBitmap::node_index
DBGBitmap::kmer_to_node(const std::string &kmer) const {
    assert(kmer.size() == k_);
    return to_node(Kmer(seq_encoder_.encode(kmer)));
}

uint64_t DBGBitmap::node_to_index(node_index node) const {
    assert(node);
    assert(node < kmers_.num_set_bits());

    return complete_ ? node : kmers_.select1(node + 1);
}

DBGBitmap::Kmer DBGBitmap::node_to_kmer(node_index node) const {
    assert(node);
    assert(node < kmers_.num_set_bits());

    return Kmer { complete_ ? node - 1 : kmers_.select1(node + 1) - 1 };
}

std::string DBGBitmap::get_node_sequence(node_index node) const {
    assert(node);
    assert(sequence_to_kmers(seq_encoder_.kmer_to_sequence(
        node_to_kmer(node), k_)).size() == 1);
    assert(node == to_node(sequence_to_kmers(seq_encoder_.kmer_to_sequence(
        node_to_kmer(node), k_))[0]));

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

Vector<DBGBitmap::Kmer> DBGBitmap::sequence_to_kmers(const std::string &sequence,
                                                     bool to_canonical) const {
    return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, to_canonical);
}

bool DBGBitmap::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const DBGBitmap*>(&other))
        return equals(*dynamic_cast<const DBGBitmap*>(&other), false);

    throw std::runtime_error("Not implemented");
}

bool DBGBitmap::equals(const DBGBitmap &other, bool verbose) const {
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
        return cur_one == other.kmers_.num_set_bits() + 1 && !mismatch;
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
        return cur_one == other.kmers_.num_set_bits() + 1 && !mismatch;
    }

    return false;
}

void DBGBitmap::print(std::ostream &out) const {
    if (complete_) {
        out << "Complete graph" << std::endl;
        return;
    }

    auto vertex_header = std::string("Vertex");
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
