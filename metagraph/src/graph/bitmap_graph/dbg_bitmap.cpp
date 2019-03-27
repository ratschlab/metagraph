#include "dbg_bitmap.hpp"

#include <cassert>
#include <cmath>

#include "serialization.hpp"
#include "utils.hpp"
#include "dbg_bitmap_construct.hpp"


// Assume all k-mers present
DBGBitmap::DBGBitmap(size_t k, bool canonical_mode)
      : alphabet(seq_encoder_.alphabet),
        k_(k),
        canonical_mode_(canonical_mode),
        kmers_(std::pow(static_cast<long double>(alphabet.size()), k_) + 1, true),
        complete_(true) {
    assert(k > 1);
    assert(kmers_.num_set_bits() == kmers_.size());
    if (k * std::log2(alphabet.size()) >= 64) {
        std::cerr << "ERROR: Too large k!"
                  << " Maximum allowed k with this alphabet is "
                  << static_cast<int>(64. / std::log2(alphabet.size())) - 1 << std::endl;
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
    for (const auto &kmer : sequence_to_kmers(sequence, canonical_mode_)) {
        callback(to_node(kmer));

        if (terminate())
            return;
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
void DBGBitmap::map_to_nodes_sequentially(std::string::const_iterator begin,
                                          std::string::const_iterator end,
                                          const std::function<void(node_index)> &callback,
                                          const std::function<bool()> &terminate) const {
    for (const auto &kmer : sequence_to_kmers(std::string(begin, end))) {
        callback(to_node(kmer));

        if (terminate())
            return;
    }
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

    for (char c : alphabet) {
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

    for (char c : alphabet) {
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

    for (char c : alphabet) {
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

    for (char c : alphabet) {
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
        kmers_.load(in);
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

typedef uint8_t TAlphabet;
struct Edge {
    DBGBitmap::node_index id;
    std::vector<TAlphabet> source_kmer;
};

/**
 * Traverse graph and extract directed paths covering the graph
 * edge, edge -> edge, edge -> ... -> edge, ... (k+1 - mer, k+...+1 - mer, ...)
 */
void DBGBitmap::call_paths(Call<const std::vector<node_index>,
                                const std::vector<TAlphabet>&> callback,
                           bool split_to_contigs) const {
    uint64_t nnodes = num_nodes();
    // keep track of reached edges
    sdsl::bit_vector discovered(nnodes + 1, false);
    // keep track of edges that are already included in covering paths
    sdsl::bit_vector visited(nnodes + 1, false);
    // store all branch nodes on the way
    std::deque<Edge> nodes;
    std::vector<node_index> path;
    std::vector<node_index> target_nodes;

    // start at the source node
    for (node_index j = 1; j <= nnodes; ++j) {
        uint64_t i = node_to_index(j);
        if (!kmers_[i])
            continue;

        if (visited[j])
            continue;

        //TODO: traverse backwards

        discovered[j] = true;

        // TODO: replace this double conversion
        nodes.push_back({ j, seq_encoder_.encode(get_node_sequence(j)) });
        nodes.back().source_kmer.pop_back();

        // keep traversing until we have worked off all branches from the queue
        while (!nodes.empty()) {
            node_index node = nodes.front().id;
            auto sequence = std::move(nodes.front().source_kmer);
            path.clear();
            nodes.pop_front();
            assert(node);

            // traverse simple path until we reach its tail or
            // the first edge that has been already visited
            while (!visited[node]) {
                assert(node);
                assert(node <= nnodes);
                assert(discovered[node]);

                sequence.push_back(node_to_kmer(node)[k_ - 1]);
                path.push_back(node);
                visited[node] = true;

                target_nodes.clear();

                adjacent_outgoing_nodes(node, &target_nodes);

                if (!target_nodes.size())
                    continue;

                if (target_nodes.size() == 1) {
                    node = target_nodes[0];
                    discovered[node] = true;
                    continue;
                }

                std::vector<TAlphabet> kmer(sequence.end() - k_ + 1, sequence.end());

                bool continue_traversal = false;
                for (const auto &next : target_nodes) {
                    if (!discovered[next]) {
                        continue_traversal = true;
                        discovered[next] = true;
                        nodes.push_back({ next, kmer });
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

            if (path.size())
                callback(path, sequence);
        }
    }
}

void DBGBitmap::call_sequences(Call<const std::string&> callback,
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
void DBGBitmap::call_kmers(Call<node_index, const std::string&> callback) const {
    uint64_t nnodes = num_nodes();
    for (size_t node = 1; node <= nnodes; ++node) {
        if (kmers_[node_to_index(node)])
            callback(node, get_node_sequence(node));
    }
}

Vector<DBGBitmap::Kmer> DBGBitmap::sequence_to_kmers(const std::string &sequence,
                                                     bool to_canonical) const {
    return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, to_canonical);
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

std::ostream& operator<<(std::ostream &out, const DBGBitmap &graph) {
    out << "k: " << graph.k_ << std::endl
        << "complete: " << graph.complete_ << std::endl
        << "canonical: " << graph.canonical_mode_ << std::endl
        << "nodes:" << std::endl;

    uint64_t nnodes = graph.num_nodes();

    for (size_t node = 1; node <= nnodes; ++node) {
        auto i = graph.node_to_index(node);
        if (graph.kmers_[i])
            out << i << "\t" << graph.get_node_sequence(node) << std::endl;
    }
    return out;
}
