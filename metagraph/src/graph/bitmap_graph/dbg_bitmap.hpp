#ifndef __DBG_BITMAP_HPP__
#define __DBG_BITMAP_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"
#include "bit_vector.hpp"


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
    explicit DBGBitmap(size_t k, bool canonical_mode = false);

    // Initialize graph from builder
    explicit DBGBitmap(DBGBitmapConstructor *builder);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const;
    void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const;

    size_t outdegree(node_index node) const;

    size_t indegree(node_index node) const;

    node_index kmer_to_node(const std::string &kmer) const;

    std::string get_node_sequence(node_index node) const;

    inline size_t get_k() const { return k_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const;


    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    bool operator==(const DBGBitmap &other) const { return equals(other, false); }
    bool operator!=(const DBGBitmap &other) const { return !(*this == other); }

    virtual bool operator==(const DeBruijnGraph &other) const override final;

    bool equals(const DBGBitmap &other, bool verbose = false) const;

    inline bool is_complete() const { return complete_; }

    friend std::ostream& operator<<(std::ostream &out, const DBGBitmap &graph);

    const std::string &alphabet;

    using Chunk = bit_vector_sd;

  private:
    using Kmer = KmerExtractor2Bit::Kmer64;

    void add_sequence(const std::string &,
                      bit_vector_dyn *) {
        throw std::runtime_error("Not implemented");
    }

    Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                   bool canonical = false) const;

    uint64_t node_to_index(node_index node) const;
    Kmer node_to_kmer(node_index node) const;
    node_index to_node(const Kmer &kmer) const;

    size_t k_;
    bool canonical_mode_;
    KmerExtractor2Bit seq_encoder_;

    bit_vector_sd kmers_;

    bool complete_ = false;

    static constexpr auto kExtension = ".bitmapdbg";
};

#endif // __DBG_BITMAP_HPP__
