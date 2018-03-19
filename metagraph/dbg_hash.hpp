#ifndef __DBG_HASH_HPP__
#define __DBG_HASH_HPP__

#include <fstream>
#include <unordered_map>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "dbg_bloom_annotator.hpp"


class DBGHash : public annotate::DeBruijnGraphWrapper {
  public:
    DBGHash(const size_t k) : k_(k) {}

    size_t get_k() const { return k_; }

    edge_index first_edge() const { return 0; }
    edge_index last_edge() const { return indices_.size() - 1; }

    std::string encode_sequence(const std::string &sequence) const;

    std::string get_node_kmer(edge_index i) const;

    char get_edge_label(edge_index i) const { return kmers_[i]->back(); }

    // Check if the source k-mer for this edge has the only outgoing edge
    bool has_the_only_outgoing_edge(edge_index i) const;

    bool has_the_only_incoming_edge(edge_index i) const;

    bool is_dummy_edge(const std::string &kmer) const;
    bool is_dummy_label(char edge_label) const;

    edge_index next_edge(edge_index i, char edge_label) const;

    edge_index prev_edge(edge_index i) const;

    void add_sequence(const std::string &sequence, bool rooted = false);

    size_t get_num_edges() const;

    void serialize(const std::string &filename) const;

    void load(const std::string &filename);

    std::string transform_sequence(const std::string &sequence, bool rooted = false) const;

  private:
    size_t k_;
    std::unordered_map<std::string, size_t> indices_;
    std::vector<const std::string*> kmers_;
};


#endif // __DBG_HASH_HPP__
