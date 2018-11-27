#ifndef __DBG_SD_HPP__
#define __DBG_SD_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "utils.hpp"
#include "bit_vector.hpp"


class DBGSD : public DeBruijnGraph {
    using KmerExtractor = utils::KmerExtractor2Bit;
    using Kmer = KmerExtractor::Kmer;
    using KmerIndex = bit_vector_sd;
  public:
    explicit DBGSD(size_t k, bool canonical_only = false);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Check whether graph contains fraction of nodes from the sequence
    bool find(const std::string &sequence,
              double kmer_discovery_fraction = 1) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    node_index kmer_to_node(const std::string &kmer) const;
    std::string node_to_kmer(node_index i) const;

    size_t get_k() const { return k_; }
    bool is_canonical_mode() const { return canonical_only_; }

    uint64_t num_nodes() const { return kmers_.size() - 1; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

  private:
    Vector<Kmer> sequence_to_kmers(const std::string &sequence) const;
    inline node_index get_index(const Kmer &kmer) const;
    inline Kmer get_kmer(node_index node) const;

    size_t k_;
    bool canonical_only_;

    KmerIndex kmers_;
    KmerExtractor seq_encoder_;

    uint64_t offset_;

    static constexpr auto kExtension = ".sddbg";
};

#endif // __DBG_SD_HPP__
