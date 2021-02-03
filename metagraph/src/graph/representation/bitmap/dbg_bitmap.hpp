#ifndef __DBG_BITMAP_HPP__
#define __DBG_BITMAP_HPP__

#include <iostream>

#include "graph/representation/base/sequence_graph.hpp"
#include "kmer/kmer_extractor.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"


namespace mtg {
namespace graph {

class DBGBitmapConstructor;


/**
 * Node-centric de Bruijn graph
 *
 * In canonical mode, for each k-mer in the graph, its
 * reverse complement is stored in the graph as well.
 */
class DBGBitmap : public DeBruijnGraph {
    friend DBGBitmapConstructor;

  public:
    // Initialize complete graph
    explicit DBGBitmap(size_t k, Mode mode = BASIC);

    // Initialize graph from builder
    explicit DBGBitmap(DBGBitmapConstructor *builder);

    void add_sequence(std::string_view,
                      const std::function<void(node_index)> & = [](node_index) {}) {
        throw std::runtime_error("Not implemented");
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(std::string_view sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string_view sequence,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const;
    void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index) const;
    bool has_multiple_outgoing(node_index) const;

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index) const;

    node_index kmer_to_node(std::string_view kmer) const;

    std::string get_node_sequence(node_index node) const;

    inline size_t get_k() const { return k_; }
    inline Mode get_mode() const { return mode_; }

    uint64_t num_nodes() const;

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DBGBitmap &other) const { return equals(other, false); }
    bool operator!=(const DBGBitmap &other) const { return !(*this == other); }

    bool operator==(const DeBruijnGraph &other) const;

    bool equals(const DBGBitmap &other, bool verbose = false) const;

    bool is_complete() const { return complete_; }

    const std::string& alphabet() const { return seq_encoder_.alphabet; }

    void print(std::ostream &out) const;

    using Chunk = bit_vector_smart;

    static constexpr auto kExtension = ".bitmapdbg";
    static constexpr auto kChunkFileExtension = ".dbgsdchunk";

  private:
    typedef mtg::kmer::KmerExtractor2Bit::Kmer64 Kmer;

    Vector<std::pair<Kmer, bool>> sequence_to_kmers(std::string_view sequence,
                                                    bool canonical = false) const;

    uint64_t node_to_index(node_index node) const;
    Kmer node_to_kmer(node_index node) const;
    node_index to_node(const Kmer &kmer) const;

    size_t k_;
    Mode mode_;
    mtg::kmer::KmerExtractor2Bit seq_encoder_;

    bit_vector_smart kmers_;

    bool complete_ = false;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_BITMAP_HPP__
