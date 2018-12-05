#ifndef __DBG_SD_HPP__
#define __DBG_SD_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "utils.hpp"
#include "kmer_extractor.hpp"
#include "bit_vector.hpp"


template <typename KMER>
class SDChunkConstructor;
class DBGSDConstructor;

class DBGSD : public DeBruijnGraph {
    using KmerExtractor = KmerExtractor2Bit;
    using KmerIndex = bit_vector_sd;
    using Kmer = KmerExtractor::Kmer64;
    using TAlphabet = KmerExtractor::TAlphabet;
    friend DBGSDConstructor;
    friend SDChunkConstructor<Kmer>;

  public:
    using Chunk = bit_vector_sd;

    explicit DBGSD(size_t k, bool canonical_mode = false);

    explicit DBGSD(DBGSDConstructor *builder);

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

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    node_index kmer_to_node(const std::string &kmer) const;
    std::string node_to_kmer(node_index i) const;

    inline size_t get_k() const { return k_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }

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
    void call_edges(Call<node_index, const std::vector<TAlphabet>&> callback) const;
    // call paths (or simple paths if |split_to_contigs| is true) that cover
    // exactly all edges in graph
    void call_paths(Call<const std::vector<node_index>,
                         const std::vector<TAlphabet>&> callback,
                    bool split_to_contigs = false) const;
    void call_sequences(Call<const std::string&> callback,
                        bool split_to_contigs = false) const;

    inline bool operator==(const DBGSD &other) const {
        return equals_internally(other, false);
    }

    inline bool operator!=(const DBGSD &other) const {
        return !(*this == other);
    }

    bool equals_internally(const DBGSD &other, bool verbose = false) const;

    friend std::ostream& operator<<(std::ostream &out, const DBGSD &graph);

    //TODO: get rid of this
    template <typename KMER>
    static typename KMER::KMerWordType kmer_to_index(const KMER &kmer);

    static size_t capacity(size_t k, size_t kLogSigma);

    const KmerIndex& data() const { return kmers_; }

    const TAlphabet alph_size;
    const std::string &alphabet;

  private:

    Vector<Kmer> sequence_to_kmers(const std::string &sequence) const;

    node_index get_index(const Kmer &kmer) const;
    Kmer get_kmer(node_index node) const;

    static size_t infer_k(size_t kmer_index_size, size_t kLogSigma);

    // Traverse the outgoing edge
    node_index outgoing(node_index node, TAlphabet next_char) const;
    // Traverse the incoming edge
    node_index incoming(node_index node, TAlphabet prev_char) const;

    std::vector<node_index> outgoing(node_index node) const;

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;

    static constexpr auto kExtension = ".sddbg";
};

#endif // __DBG_SD_HPP__
