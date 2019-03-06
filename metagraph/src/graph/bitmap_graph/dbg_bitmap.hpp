#ifndef __DBG_BITMAP_HPP__
#define __DBG_BITMAP_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"
#include "bit_vector.hpp"


class DBGSDConstructor;

/**
 * Node-centric de Bruijn graph
 *
 * In canonical mode, for each k-mer in the graph, its
 * reverse complement is stored in the graph as well.
 */
class DBGSD : public DeBruijnGraph {
    friend DBGSDConstructor;

  public:
    // Initialize complete graph
    explicit DBGSD(size_t k, bool canonical_mode = false);

    // Initialize graph from builder
    explicit DBGSD(DBGSDConstructor *builder);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Map k-mers from sequence to nodes of the graph similarly to map_to_nodes
    // Guarantees that the k-mers from sequence are called in their natural order
    void map_kmers_sequentially(std::string::const_iterator,
                                std::string::const_iterator,
                                const std::function<void(node_index)> &,
                                const std::function<bool()> &) const {
        // TODO: Complete map_sequence_sequentially for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const {
        // TODO: Complete call_outgoing_kmers for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const {
        // TODO: Complete call_incoming_kmers for DBGSD.
        throw std::runtime_error("Not implemented");
    }

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

    node_index kmer_to_node(const std::string &kmer) const;
    std::string node_to_kmer(node_index i) const;

    std::string get_node_sequence(node_index node) const {
        return node_to_kmer(node);
    }

    char get_kmer_last_char(node_index) const {
        // TODO: Complete get_kmer_last_char for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    inline size_t get_k() const { return k_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }
    bool is_single_outgoing(node_index) const {
        // TODO: Complete is_single_outgoing for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    inline uint64_t num_nodes() const { return kmers_.num_set_bits() - 1; }
    inline uint64_t capacity() const { return kmers_.size() - 1; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    template <class... T>
    using Call = typename std::function<void(T...)>;

    // traverse all nodes in graph
    void call_kmers(Call<node_index, const std::string&> callback) const;

    // call paths (or simple paths if |split_to_contigs| is true) that cover
    // exactly all kmers in graph
    void call_sequences(Call<const std::string&> callback,
                        bool split_to_contigs = false) const;

    bool operator==(const DBGSD &other) const { return equals(other, false); }
    bool operator!=(const DBGSD &other) const { return !(*this == other); }

    bool equals(const DBGSD &other, bool verbose = false) const;

    friend std::ostream& operator<<(std::ostream &out, const DBGSD &graph);

    const std::string &alphabet;

    using Chunk = bit_vector_sd;

  private:
    KmerExtractor2Bit seq_encoder_;
    using Kmer = KmerExtractor2Bit::Kmer64;

    // Not supported
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL);

    Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                   bool to_canonical = false) const;

    // get the node index if |kmer| is in the graph and npos otherwise
    node_index get_node(const Kmer &kmer) const;

    // translate index |node| from [1...1+|A|^k] to a Kmer in A^k
    Kmer to_kmer(node_index node) const;
    // translate |kmer| from A^k to a node index in [1...1+|A|^k]
    node_index to_index(const Kmer &kmer) const;

    typedef uint8_t TAlphabet;

    void call_paths(Call<const std::vector<node_index>,
                    const std::vector<TAlphabet>&> callback,
                    bool split_to_contigs) const;

    size_t k_;
    bool canonical_mode_;

    bit_vector_sd kmers_;

    static constexpr auto kExtension = ".bitmapdbg";
};

#endif // __DBG_BITMAP_HPP__
