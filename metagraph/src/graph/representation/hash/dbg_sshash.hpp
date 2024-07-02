#ifndef __DBG_SSHASH_HPP__
#define __DBG_SSHASH_HPP__

#include <iostream>
#include <variant>

#include <dictionary.hpp>
#include <sdsl/uint256_t.hpp>

#include "graph/representation/base/sequence_graph.hpp"

namespace mtg::graph {

class DBGSSHash : public DeBruijnGraph {
  public:
    using KmerInt64 = uint64_t;
    using KmerInt128 = __uint128_t;
    using KmerInt256 = sdsl::uint256_t;

#if _PROTEIN_GRAPH
    template <typename KmerInt>
    using kmer_t = sshash::aa_uint_kmer_t<KmerInt>;
#else
    template <typename KmerInt>
    using kmer_t = sshash::dna_uint_kmer_t<KmerInt>;
#endif

    using dict_t = std::variant<
                        sshash::dictionary<kmer_t<KmerInt64>>,
                        sshash::dictionary<kmer_t<KmerInt128>>,
                        sshash::dictionary<kmer_t<KmerInt256>>>;

    explicit DBGSSHash(size_t k, Mode mode = BASIC);
    DBGSSHash(const std::string &input_filename, size_t k, Mode mode = BASIC, size_t num_chars = 0);

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

    template <bool with_rc = true>
    void map_to_nodes_with_rc(
            std::string_view sequence,
            const std::function<void(node_index, bool)>& callback,
            const std::function<bool()>& terminate = []() { return false; }) const;

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

    template <bool with_rc = true>
    std::pair<node_index, bool> kmer_to_node_with_rc(std::string_view kmer) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback& callback) const override;

    template <bool with_rc = true>
    void call_outgoing_kmers_with_rc(
            node_index node,
            const std::function<void(node_index, char, bool)>& callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback& callback) const override;

    template <bool with_rc = true>
    void call_incoming_kmers_with_rc(
            node_index node,
            const std::function<void(node_index, char, bool)>& callback) const;


    const std::string& alphabet() const override final { return alphabet_; }

    const dict_t& data() const { return dict_; }

    node_index reverse_complement(node_index node) const;

  private:
    static const std::string alphabet_;
    dict_t dict_;
    size_t k_;
    size_t num_nodes_;
    Mode mode_;

    size_t dict_size() const;
};

} // namespace mtg::graph

#endif // __DBG_SSHASH_HPP__
