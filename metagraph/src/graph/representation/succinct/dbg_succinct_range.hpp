#ifndef __DBG_SUCCINCT_RANGE_HPP__
#define __DBG_SUCCINCT_RANGE_HPP__

#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "dbg_succinct.hpp"

#include <tsl/ordered_set.h>

namespace mtg {
namespace graph {

/**
 * DBGSuccinctRange is a wrapper around a DBGSuccinct object which partially
 * implements a variable-order de Bruijn graph. It extents the node space of its
 * graph with virtual nodes which represent ranges of BOSS nodes that share a common
 * suffix. Each range corresponds to two virtual nodes: one representing a set of
 * dummy source k-mers (e.g., $$$AGT, having an offset of 3) and the other
 * representing dummy sink k-mers (e.g., AGT$$$). Forward traversal on source
 * k-mers and backwards traversal on sink k-mers increase the matched suffix length.
 * To circularize the graph, the root (e.g., $$$$$) is represented by two nodes
 * (source and sink), where the only child of the sink root node is the source root node.
 * Only source ranges are explicitly stored, sink ranges are implicit.
 */
class DBGSuccinctRange : public DeBruijnGraph {
  public:
    using node_index = DBGSuccinct::node_index;
    using EdgeRange = std::tuple<boss::BOSS::edge_index /* first edge */,
                                 boss::BOSS::edge_index /* last edge */,
                                 size_t /* offset - 1 */>;

    struct EdgeRangeHash {
        uint64_t operator()(const EdgeRange &range) const {
            // combine hashes of offset and the range end
            // https://stackoverflow.com/questions/8513911/how-to-create-a-good-hash-combine-with-64-bit-output-inspired-by-boosthash-co
            auto [first, last, offset] = range;
            constexpr uint64_t kMul = 0x9ddfea08eb382d69ULL;
            uint64_t a = (offset ^ last) * kMul;
            a ^= (a >> 47);
            uint64_t b = (last ^ a) * kMul;
            b ^= (b >> 47);
            return b * kMul;
        }
    };

    using EdgeDescriptor = std::pair<EdgeRange, bool /* is_sink_k-mer */>;
    using EdgeRangeStorage = tsl::ordered_set<EdgeRange, EdgeRangeHash,
                                              std::equal_to<EdgeRange>,
                                              std::allocator<EdgeRange>,
                                              std::vector<EdgeRange>,
                                              uint64_t>;

    DBGSuccinctRange(std::shared_ptr<const DBGSuccinct> graph)
          : dbg_succ_ptr_(graph),
            dbg_succ_(*dbg_succ_ptr_),
            offset_(dbg_succ_.max_index() + 1) {}

    DBGSuccinctRange(std::shared_ptr<const DeBruijnGraph> graph)
          : graph_ptr_(graph),
            dbg_succ_ptr_(std::dynamic_pointer_cast<const DBGSuccinct>(graph)),
            dbg_succ_(*dbg_succ_ptr_),
            offset_(dbg_succ_.max_index() + 1) {}

    DBGSuccinctRange(const DBGSuccinct &graph)
          : dbg_succ_(graph), offset_(dbg_succ_.max_index() + 1) {}

    virtual ~DBGSuccinctRange() override {}

    virtual size_t get_k() const override { return dbg_succ_.get_k(); }
    virtual bool is_canonical_mode() const override { return dbg_succ_.is_canonical_mode(); }

    /**
     * Traverse the outgoing edge. If node is a source k-mer, increase the suffix
     * match length by 1.
     */
    virtual node_index traverse(node_index node, char next_char) const override;

    /**
     * Traverse the incoming edge. If node is a sink k-mer, increase the suffix
     * match length by 1.
     */
    virtual node_index traverse_back(node_index node, char prev_char) const override;

    /**
     * Traverse graph mapping sequence to the graph nodes
     * and run callback for each node until the termination condition is satisfied.
     * Guarantees that nodes are called in the same order as the input sequence.
     * In canonical mode, non-canonical k-mers are not mapped to canonical ones.
     * If a k-mer in the sequence is not found, the longest matching prefix of
     * the k-mer is found and its corresponding source k-mer is called.
     * If the last k-mer is not mapped, each suffix of the sequence of
     * length >= get_dbg_succ().get_boss().get_indexed_suffix_length() is also
     * mapped and called.
     */
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override;

    /**
     * Traverse graph mapping sequence to the graph nodes
     * and run callback for each node until the termination condition is satisfied.
     * If a k-mer in the sequence is not found, the longest matching prefix of
     * the k-mer is found and its corresponding source k-mer is called.
     * If the last k-mer is not mapped, each suffix of the sequence of
     * length >= get_dbg_succ().get_boss().get_indexed_suffix_length() is also
     * mapped and called.
     */
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Given a starting node and a sequence of edge labels, traverse the graph
    // forward. The traversal is terminated once terminate() returns true or
    // when the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones.
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const override;

    typedef std::function<void(const std::string&, const std::vector<node_index>&)> CallPath;

    /**
     * Call contigs, a set of sequences covering each node in the graph exactly once.
     * @param num_threads number of threads to use for graph traversal
     * @param kmers_in_single_form if true, output each k-mer only in one of its forms
     * (canonical/non-canonical). That is, skip a k-mer if its reverse-complement has been
     * extracted.
     */
    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override {
        dbg_succ_.call_sequences(callback, num_threads, kmers_in_single_form);
    }
    /**
     * Call all unitigs except short tips, where tips are
     * the unitigs with InDegree(first) + OutDegree(last) < 2.
     * If |kmers_in_single_form| is true, output each k-mer only in one of its
     * forms (canonical/non-canonical). That is, skip a k-mer if its
     * reverse-complement has been extracted.
     */
    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override {
        dbg_succ_.call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override {
        dbg_succ_.call_kmers(callback);
    }

    virtual size_t outdegree(node_index) const override;
    virtual size_t indegree(node_index) const override;
    virtual bool has_single_outgoing(node_index) const override;
    virtual bool has_multiple_outgoing(node_index) const override;
    virtual bool has_no_incoming(node_index) const override;
    virtual bool has_single_incoming(node_index) const override;

    // Given a k-mer, return its corresponding node in the graph.
    // Strings representing source and sink k-mers (e.g., $$$ATG, ATG$$$) are
    // mapped to their corresponding range nodes.
    virtual node_index kmer_to_node(std::string_view kmer) const override;

    using OutgoingEdgeCallback = std::function<void(node_index /* target_kmer */,
                                                    char /* last_target_char */)>;
    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override;

    using IncomingEdgeCallback = std::function<void(node_index /* source_kmer */,
                                                    char /* first_source_char */)>;
    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;
    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override {
        return dbg_succ_.find(sequence, discovery_fraction);
    }

    virtual const std::string& alphabet() const override { return dbg_succ_.alphabet(); }

    // Call the source root node (e.g., $$$$$)
    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override;

    virtual std::string get_node_sequence(node_index node) const override;

    virtual void add_sequence(std::string_view, const std::function<void(node_index)> &) override {
        throw std::runtime_error("Not implemented");
    }

    virtual bool load(const std::string&) override {
        throw std::runtime_error("Not implemented");
    }

    virtual void serialize(const std::string&) const override {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string file_extension() const override {
        throw std::runtime_error("Not implemented");
    }

    /**
     * Returns the number of nodes in the graph, which is equal to the number of
     * nodes in the underlying graph plus the number of cached source and sink nodes.
     */
    virtual uint64_t num_nodes() const override {
        // Only ranges corresponding to sources are explicitly stored
        return dbg_succ_.num_nodes() + edge_pairs_.size() * 2;
    }

    virtual uint64_t max_index() const override {
        // Define a very loose upper bound
        return (1llu << 63) - 1;
    };

    const DBGSuccinct& get_dbg_succ() const { return dbg_succ_; }

    // Get the number of non-matching characters of the node.
    size_t get_offset(node_index node) const;

    // Map a sink node to its corresponding source node and vice-versa. Normal
    // nodes are returned as-is.
    node_index toggle_node_sink_source(node_index node) const;

    // Return true iff it is a sink node.
    bool is_sink(node_index node) const;

    // Compute the reverse complement of seq and compute its corresponding nodes.
    // Source nodes are mapped to their corresponding sink nodes and vice-versa.
    void reverse_complement(std::string &seq, std::vector<node_index> &path) const;

    // Given a node, call all corresponding nodes in the underlying graph.
    void call_nodes_in_range(node_index node,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate = []() { return false; }) const;

    void reset_range_storage() {
        std::lock_guard<std::mutex> lock(edge_pair_mutex_);
        edge_pairs_ = EdgeRangeStorage();
    }

  private:
    const std::shared_ptr<const DeBruijnGraph> graph_ptr_;
    const std::shared_ptr<const DBGSuccinct> dbg_succ_ptr_;
    const DBGSuccinct &dbg_succ_;

    mutable EdgeRangeStorage edge_pairs_;
    mutable std::mutex edge_pair_mutex_;
    size_t offset_;

    node_index kmer_to_node(const boss::BOSS::TAlphabet *begin,
                            const boss::BOSS::TAlphabet *end,
                            bool require_exact_length = false) const;

    EdgeDescriptor fetch_edge_range(node_index node) const;

    void call_left_tightened_ranges(node_index node,
                                    boss::BOSS::edge_index first,
                                    boss::BOSS::edge_index last,
                                    size_t offset,
                                    const std::function<void(node_index, boss::BOSS::TAlphabet)> &callback,
                                    boss::BOSS::TAlphabet s = std::numeric_limits<boss::BOSS::TAlphabet>::max()) const;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_RANGE_HPP__
