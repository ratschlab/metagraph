#ifndef __DBG_BLOOM_ANNOTATOR_HPP__
#define __DBG_BLOOM_ANNOTATOR_HPP__

#include <libmaus2/util/NumberSerialisation.hpp>

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

    virtual std::string transform_sequence(const std::string &sequence, bool rooted = false) const {
        return rooted ? sequence : sequence;
    }

    virtual std::string get_node_kmer(edge_index i) const = 0;
    virtual char get_edge_label(edge_index i) const = 0;

    // Check if the source k-mer for this edge has the only outgoing edge
    virtual bool has_the_only_outgoing_edge(edge_index i) const = 0;
    virtual bool has_the_only_incoming_edge(edge_index i) const = 0;

    virtual bool is_dummy_label(char edge_label) const = 0;
    virtual bool is_dummy_edge(const std::string &kmer) const = 0;

    virtual edge_index next_edge(edge_index i, char edge_label) const = 0;
    virtual edge_index prev_edge(edge_index i) const = 0;

};


class PreciseAnnotator {
  public:
    PreciseAnnotator(const DeBruijnGraphWrapper &graph) : graph_(graph) {}

    void add_sequence(const std::string &sequence, size_t column = -1llu, bool rooted = false);

    void add_column(const std::string &sequence, bool rooted = false);

    std::vector<uint64_t> annotation_from_kmer(const std::string &kmer) const;

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    void load(std::istream &in);
    void load(const std::string &filename);

    void export_rows(std::ostream &out) const;
    void export_rows(const std::string &filename) const;

  private:
    ExactHashAnnotation annotation_exact;
    const DeBruijnGraphWrapper &graph_;
};


class BloomAnnotator {
  public:
    // Computes optimal `bloom_size_factor` and `num_hash_functions` automatically
    BloomAnnotator(const DeBruijnGraphWrapper &graph,
                   double bloom_fpp,
                   bool verbose = false);

    // If not provided, computes optimal `num_hash_functions` automatically
    BloomAnnotator(const DeBruijnGraphWrapper &graph,
                   double bloom_size_factor,
                   size_t num_hash_functions,
                   bool verbose = false);

    void add_sequence(const std::string &sequence, size_t column, size_t num_elements = 0);

    void add_column(const std::string &sequence, size_t num_elements = 0);

    std::vector<uint64_t> get_annotation(DeBruijnGraphWrapper::edge_index i) const;

    std::vector<uint64_t> get_annotation_corrected(DeBruijnGraphWrapper::edge_index i,
                                                   bool check_both_directions = false,
                                                   size_t path_cutoff = 50) const;

    template <typename Annotator>
    void test_fp_all(const Annotator &annotation_exact, size_t num = 0, bool check_both_directions = false) const;

    void serialize(std::ostream &out) const;

    void serialize(const std::string &filename) const;

    static std::vector<size_t> unpack(const std::vector<uint64_t> &packed);

    size_t num_hash_functions() const;

    double size_factor() const;

    double approx_false_positive_rate() const;

    size_t get_size(size_t i) const;

  private:
    std::vector<uint64_t> annotation_from_kmer(const std::string &kmer) const;

    std::string kmer_from_index(DeBruijnGraphWrapper::edge_index index) const;

    template <typename Annotator>
    std::vector<size_t> test_fp(DeBruijnGraphWrapper::edge_index i,
                                 const Annotator &annotation_exact,
                                 bool check_both_directions = false) const;

    const DeBruijnGraphWrapper &graph_;
    double bloom_size_factor_;
    double bloom_fpp_;
    BloomHashAnnotation annotation;

    //TODO: get rid of this if not using degree Bloom filter
    std::vector<size_t> sizes_v;

    size_t total_traversed_;
    bool verbose_;
};

} // namespace annotate

#endif // __DBG_BLOOM_ANNOTATOR_HPP__
