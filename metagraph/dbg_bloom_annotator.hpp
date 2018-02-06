#ifndef __DBG_BLOOM_ANNOTATOR_HPP__
#define __DBG_BLOOM_ANNOTATOR_HPP__

#include "hashers.hpp"

namespace annotate {

class DeBruijnGraphWrapper {
  public:
    typedef uint64_t edge_index;

    virtual ~DeBruijnGraphWrapper() {}

    virtual size_t get_k() const = 0;

    virtual edge_index first_edge() const = 0;
    virtual edge_index last_edge() const = 0;

    // Transform sequence to the same kind as the de bruijn graph stores
    virtual std::string encode_sequence(const std::string &sequence) const {
        return sequence;
    }

    virtual std::string get_node_kmer(edge_index i) const = 0;
    virtual char get_edge_label(edge_index i) const = 0;

    // Check if the source k-mer for this edge has the only outgoing edge
    virtual bool has_the_only_outgoing_edge(edge_index i) const = 0;
    virtual bool has_the_only_incoming_edge(edge_index i) const = 0;

    virtual bool is_dummy_edge(const std::string &kmer) const = 0;

    virtual edge_index next_edge(edge_index i, char edge_label) const = 0;
    virtual edge_index prev_edge(edge_index i) const = 0;
};


class PreciseAnnotator {
  public:
    PreciseAnnotator(const DeBruijnGraphWrapper &graph) : graph_(graph) {}

    void add_sequence(const std::string &sequence, size_t column);

    void add_column(const std::string &sequence);

    std::vector<size_t> annotation_from_kmer(const std::string &kmer) const;

  private:
    annotate::HashAnnotation<annotate::ExactFilter> annotation_exact;
    const DeBruijnGraphWrapper &graph_;
};


class BloomAnnotator {
  public:
    BloomAnnotator(size_t num_hash_functions,
                   const DeBruijnGraphWrapper &graph,
                   double bloom_size_factor,
                   bool verbose = false);

    void add_sequence(const std::string &sequence, size_t column);

    void add_column(const std::string &sequence);

    std::vector<size_t> annotation_from_kmer(const std::string &kmer) const;

    std::vector<size_t> get_annotation(DeBruijnGraphWrapper::edge_index i) const;

    std::vector<size_t> get_annotation_corrected(DeBruijnGraphWrapper::edge_index i,
                                                 size_t path_cutoff = 50) const;

    void test_fp_all(const PreciseAnnotator &annotation_exact, size_t step = 1) const;

    void serialize(std::ostream &out) const;

    void serialize(const std::string &filename) const;

    static std::vector<size_t> unpack(const std::vector<size_t> &packed);

  private:
    std::string kmer_from_index(DeBruijnGraphWrapper::edge_index index) const;

    std::vector<uint8_t> test_fp(DeBruijnGraphWrapper::edge_index i,
                                 const PreciseAnnotator &annotation_exact) const;

    annotate::HashAnnotation<annotate::BloomFilter> annotation;

    const DeBruijnGraphWrapper &graph_;
    double bloom_size_factor_;

    //TODO: get rid of this if not using degree Bloom filter
    std::vector<size_t> sizes_v;

    size_t total_traversed_;
    bool verbose_;
};

} // namespace annotate

#endif // __DBG_BLOOM_ANNOTATOR_HPP__
