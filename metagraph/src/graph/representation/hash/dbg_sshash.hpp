#ifndef __DBG_SSHASH_HPP__
#define __DBG_SSHASH_HPP__

#include <iostream>
#include <tsl/ordered_set.h>
#include "common/utils/string_utils.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace sshash{
class dictionary;
}
namespace mtg {
namespace graph {

class DBGSSHash : public DeBruijnGraph {
  public:
    explicit DBGSSHash(size_t k);
    DBGSSHash(std::string const& input_filename, size_t k);

    ~DBGSSHash();

// SequenceGraph overrides
    void add_sequence(
            std::string_view sequence,
            const std::function<void(node_index)> &on_insertion = [](node_index) {}) override;

    void map_to_nodes(
            std::string_view sequence,
            const std::function<void(node_index)> &callback,
            const std::function<bool()> &terminate = []() { return false; }) const override;

    void map_to_nodes_sequentially(
            std::string_view sequence,
            const std::function<void(node_index)> &callback,
            const std::function<bool()> &terminate = []() { return false; }) const override;

    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const override;

    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const override;

    uint64_t num_nodes() const override;

    bool load(std::istream &in);
    bool load(const std::string &filename) override;

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const override;

    static constexpr auto kExtension = ".sshashdbg";
    std::string file_extension() const override;

    std::string get_node_sequence(node_index node) const override;

// DeBruijnGraph overrides
    size_t get_k() const override ;

    // TODO: add the support for the canonical mode
    Mode get_mode() const override;
    
    node_index traverse(node_index node, char next_char) const override;
    node_index traverse_back(node_index node, char prev_char) const override;

    void call_kmers(const std::function<void(node_index, const std::string &)> &callback) const override;

    size_t outdegree(node_index) const override;
    bool has_single_outgoing(node_index) const override;
    bool has_multiple_outgoing(node_index) const override;

    size_t indegree(node_index) const override;
    bool has_no_incoming(node_index) const override;
    bool has_single_incoming(node_index) const override;

    node_index kmer_to_node(std::string_view kmer) const override;
    
    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const override;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const override;


    bool operator==(const DeBruijnGraph &other) const override;

    const std::string &alphabet() const override;

  private:
    static const std::string alphabet_;
    std::unique_ptr<sshash::dictionary> dict_;
    size_t k_;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_SSHASH_HPP__
