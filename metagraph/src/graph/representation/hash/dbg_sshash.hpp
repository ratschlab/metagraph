#ifndef __DBG_SSHASH_HPP__
#define __DBG_SSHASH_HPP__

#include <iostream>
#include <limits>
#include <memory>
#include <variant>
#include <filesystem>

#include <dictionary.hpp>
#include <offsets.hpp>
#include <kmer.hpp>
#include <sdsl/uint256_t.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

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

    using dict_t = std::variant<sshash::dictionary<kmer_t<KmerInt64>, sshash::decoded_offsets>,
                                sshash::dictionary<kmer_t<KmerInt128>, sshash::decoded_offsets>,
                                sshash::dictionary<kmer_t<KmerInt256>, sshash::decoded_offsets>>;

    explicit DBGSSHash(size_t k, Mode mode = BASIC);
    // Construction uses scratch files; an empty swap_dir selects the system temp directory.
    DBGSSHash(const std::string& input_filename,
              size_t k,
              Mode mode = BASIC,
              size_t num_chars = 0,
              const std::filesystem::path &swap_dir = {});

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
    uint64_t max_index() const override;

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

    void call_nodes(
            const std::function<void(node_index)>& callback,
            const std::function<bool()>& terminate = []() { return false; },
            size_t num_threads = 1,
            size_t batch_size = 1'000'000) const override;

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

    std::vector<seq_index> get_sequence_ids(node_index node) const override;


    const std::string& alphabet() const override final { return alphabet_; }

    const dict_t& data() const { return dict_; }

    node_index reverse_complement(node_index node) const;
    node_index get_base_node(node_index node) const;

    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)>& callback) const override;

    template <bool with_rc = true>
    void
    adjacent_outgoing_nodes_with_rc(node_index node,
                                    const std::function<void(node_index, bool)>& callback) const;

    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)>& callback) const override;

    template <bool with_rc = true>
    void
    adjacent_incoming_nodes_with_rc(node_index node,
                                    const std::function<void(node_index, bool)>& callback) const;

    sdsl::bit_vector get_invalid_char_mask(std::string_view sequence) const;

    // Call each node whose prefix matches at least minimum query characters,
    // reporting its longest matching prefix. Invalid characters end a match.
    // Return the total number of matching nodes, even when it exceeds maximum
    // and suppresses all callbacks. The prefix index is lazy
    // and is not serialized. Canonical DNA includes both orientations; protein
    // requires BASIC mode. A zero minimum is treated as one.
    size_t call_nodes_with_prefix_matching(
            std::string_view query,
            size_t minimum,
            const std::function<void(node_index, size_t)>& callback,
            size_t maximum = std::numeric_limits<size_t>::max()) const;

    // Call nodes whose suffix matches the longest prefix of query (at most k
    // characters). Reject invalid input and suppress the entire longest-match
    // bucket if it exceeds maximum. Return its size, capped at maximum + 1.
    // The suffix order is lazy and is not serialized. A zero minimum is one.
    size_t call_nodes_with_suffix_matching_longest_prefix(
            std::string_view query,
            std::function<void(node_index, uint64_t)> callback,
            size_t minimum = 1,
            size_t maximum = std::numeric_limits<size_t>::max()) const;

    static constexpr node_index sshash_to_graph_index(uint64_t idx) {
        return idx + 1;
    }

    static constexpr uint64_t graph_index_to_sshash(node_index idx) {
        return idx - 1;
    }

  private:
    struct SeedIndexCache;
    // Copies share only immutable node ordering; a successful load replaces it.
    std::shared_ptr<SeedIndexCache> seed_index_cache_;

    static const std::string alphabet_;
    dict_t dict_;
    size_t k_;
    size_t num_nodes_;
    Mode mode_;
    bit_vector_smart succ_is_next_;
    bit_vector_smart pred_is_prev_;

    size_t dict_size() const;

    char get_last_char(node_index node) const;
    char get_first_char(node_index node) const;
};

} // namespace mtg::graph

#endif // __DBG_SSHASH_HPP__
