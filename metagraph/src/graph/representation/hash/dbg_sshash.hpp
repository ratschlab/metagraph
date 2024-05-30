#ifndef __DBG_SSHASH_HPP__
#define __DBG_SSHASH_HPP__

#include <iostream>

#include <tsl/ordered_set.h>
#include <dictionary.hpp>

#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg::graph {

class DBGSSHash : public DeBruijnGraph {
    using kmer_t = sshash::dna_uint_kmer_t<uint64_t>;

  public:
    explicit DBGSSHash(size_t k, Mode mode = BASIC);
    DBGSSHash(std::string const& input_filename, size_t k, Mode mode = BASIC);

    // SequenceGraph overrides
    void add_sequence(
            std::string_view sequence,
            const std::function<void(node_index)>& on_insertion = [](node_index) {}) override;

    void map_to_nodes(
            std::string_view sequence,
            const std::function<void(node_index)>& callback,
            const std::function<bool()>& terminate = []() { return false; }) const override;

    void map_to_nodes_sequentially(
            std::string_view sequence,
            const std::function<void(node_index)>& callback,
            const std::function<bool()>& terminate = []() { return false; }) const override;

    void map_to_nodes_with_rc(
            std::string_view sequence,
            const std::function<void(node_index, bool)>& callback,
            const std::function<bool()>& terminate = []() { return false; }) const;

    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)>& callback) const override;

    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)>& callback) const override;

    uint64_t num_nodes() const override;

    bool load(std::istream& in);
    bool load(const std::string& filename) override;

    void serialize(std::ostream& out) const;
    void serialize(const std::string& filename) const override;

    static constexpr auto kExtension = ".sshashdbg";
    std::string file_extension() const override;

    std::string get_node_sequence(node_index node) const override;

    // DeBruijnGraph overrides
    size_t get_k() const override final { return k_; }

    Mode get_mode() const override final { return mode_; }

    node_index traverse(node_index node, char next_char) const override;
    node_index traverse_back(node_index node, char prev_char) const override;

    void call_nodes(const std::function<void(node_index)>& callback,
                    const std::function<bool()> &terminate = [](){ return false; }) const override;

    size_t outdegree(node_index) const override;

    size_t indegree(node_index) const override;

    node_index kmer_to_node(std::string_view kmer) const override;
    std::pair<node_index, bool> kmer_to_node_with_rc(std::string_view kmer) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback& callback) const override;

    void call_outgoing_kmers_with_rc(
            node_index node,
            const std::function<void(node_index, char, bool)>& callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback& callback) const override;

    void call_incoming_kmers_with_rc(
            node_index node,
            const std::function<void(node_index, char, bool)>& callback) const;


    bool operator==(const DeBruijnGraph& other) const override;

    const std::string& alphabet() const override final { return alphabet_; }

    const sshash::dictionary<kmer_t>& data() const { return dict_; }

    node_index reverse_complement(node_index node) const;

  private:
    static const std::string alphabet_;
    sshash::dictionary<kmer_t> dict_;
    size_t k_;
    size_t num_nodes_;
    Mode mode_;
};

} // namespace mtg::graph

#endif // __DBG_SSHASH_HPP__
