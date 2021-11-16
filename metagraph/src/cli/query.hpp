#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <variant>
#include <optional>
#include <json/json.h>
#include <sdsl/int_vector.hpp>
#include "common/vector.hpp"
#include "graph/alignment/aligner_alignment.hpp"

class ThreadPool;

namespace mtg {

namespace seq_io {
    class FastaParser;
}

namespace graph {
    class AnnotatedDBG;
    namespace align {
        class DBGAlignerConfig;
    }
}


namespace cli {

class Config;

using StringGenerator = std::function<void(std::function<void(const std::string &)>)>;

/**
 * Construct a query graph and augment it with neighboring paths around it
 * if `config` is specified.
 * @param anno_graph annotated de Bruijn graph (the index)
 * @param call_sequences generate sequences to be queried against anno_graph
 * @param num_threads number of threads to use
 * @param canonical if true, the returned query graph is a canonical graph
 * @param config a pointer to a Config to determine parameters of the hull
 */
std::unique_ptr<graph::AnnotatedDBG>
construct_query_graph(const graph::AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      size_t num_threads,
                      const Config *config = nullptr);


/**
 * Simple struct to wrap a query sequence.
 */
struct QuerySequence {
    /**
     * Simple struct to wrap alignment results for a query sequence
     */
    struct Alignment {
        graph::align::DBGAlignerConfig::score_t score;  // Alignment score
        std::string cigar;                              // Cigar string

        Alignment(graph::align::DBGAlignerConfig::score_t score, std::string&& cigar)
            : score(score), cigar(std::move(cigar)) {};
    };

    size_t id;                          // Sequence ID
    std::string name;                   // Sequence name
    std::string sequence;               // Sequence string representation
    std::optional<Alignment> alignment; // Optional wrapper for alignment struct

    /**
     * Simple constructor for QuerySequence which moves sequence into struct.
     * Sets alignment to std::nullopt
     */
    QuerySequence(size_t id, const std::string &name, const std::string &sequence)
        : id(id), name(name), sequence(std::move(sequence)) {
        alignment = std::nullopt;
    }
};


/**
 * Search result for an individual sequence. Stores the search result in one of several possible
 * vector types (kmer counts, label counts, quantiles, coordinates, signatures) represented
 * as a typesafe union std::variant.
 *
 * Note that this is an intermediate storage class and does not (1) store the actual sequence
 * nor (2) perform any of the search itself.
 */
class SeqSearchResult {

public:

    typedef std::string label;
    typedef std::vector<label> label_vec;
    typedef std::vector<std::pair<label, size_t>> label_count_vec;
    typedef std::vector<std::pair<label, sdsl::bit_vector>> label_sig_vec;
    typedef std::vector<std::pair<label, std::vector<size_t>>> label_quantile_vec;
    typedef std::vector<std::pair<label, std::vector<SmallVector<uint64_t>>>> label_coord_vec;

    typedef std::variant<
            label_vec,
            label_count_vec,
            label_sig_vec,
            label_quantile_vec,
            label_coord_vec
            > result_type;

    /**
     * Construct an instance of SeqSearchResult by providing QuerySequence instance describing
     * sequence and result from that query on an annotated DBG. Will move both into this instance.
     *
     * @param sequence  rvalue reference to QuerySequence (WILL BE MOVED)
     * @param result    rvalue reference to result_type (WILL BE MOVED)
     */
    SeqSearchResult(QuerySequence&& sequence, result_type&& result)
            : sequence(std::move(sequence)),
            result(std::move(result)) {};

    /** Const reference getters */
    const QuerySequence& get_sequence() const { return sequence; }
    const result_type& get_result() const { return result; }

    /**
     * Returns a Json object representing the individual query result for the
     * represented sequence.
     *
     * @param counts_kmers      should counts be labeled kmer (t) or label (f) counts?
     * @param anno_graph        reference to annotated dbg for kmer presence mask scoring
     * @param expand_coords     do not collapse continuous ranges of coords (query_coords)
     * @return  Json::Value instance representing sequence result
     */
    Json::Value to_json(bool counts_kmers,
                        bool expand_coords,
                        const graph::AnnotatedDBG &anno_graph) const;

    /**
     * Returns a string representing the individual query result for the represented sequence.
     * Follows the format '<seq_name>\t<label>:<info>\t<label>:<info>\t ...'
     *
     * Done this way to allow for some output related config (as opposed to ostream override)
     *
     * @param delimiter             the delimiter between labels for that sequence
     * @param suppress_unlabeled    do not print seq_name if sequence is unlabeled
     * @param expand_coords         do not collapse continuous ranges of coords (query_coords)
     * @param anno_graph            reference to annotated dbg for kmer presence mask scoring
     */
    std::string to_string(std::string delimiter,
                          bool suppress_unlabeled,
                          bool expand_coords,
                          const graph::AnnotatedDBG &anno_graph) const;

private:

    QuerySequence sequence;     // query sequence this result represents
    result_type result;         // result vector of labels and additional info

};


class QueryExecutor {

public:

    QueryExecutor(const Config &config,
                  const graph::AnnotatedDBG &anno_graph,
                  std::unique_ptr<graph::align::DBGAlignerConfig>&& aligner_config,
                  ThreadPool &thread_pool) :
                  config_(config), anno_graph_(anno_graph),
                  aligner_config_(std::move(aligner_config)),
                  thread_pool_(thread_pool) {};

    /**
     * Query sequences from a FASTA file on the stored QueryExecutor::anno_graph.
     *
     * Provide a callback which receives a reference to a SeqSearchResult instance which the
     * callback function should use as desired (typically generate a string/JSON output).
     *
     * @param file_path     path to FASTA file
     * @param callback      callback function
     */
    void query_fasta(const std::string &file_path,
                     const std::function<void(const SeqSearchResult &)> &callback);

    static SeqSearchResult execute_query(QuerySequence&& sequence,
                                         bool count_labels,
                                         bool print_signature,
                                         size_t num_top_labels,
                                         double discovery_fraction,
                                         double presence_fraction,
                                         const graph::AnnotatedDBG &anno_graph,
                                         bool with_kmer_counts = false,
                                         const std::vector<double> &count_quantiles = {},
                                         bool query_coords = false);

private:

    const Config &config_;
    const graph::AnnotatedDBG &anno_graph_;
    std::unique_ptr<graph::align::DBGAlignerConfig> aligner_config_;
    ThreadPool &thread_pool_;

    void batched_query_fasta(mtg::seq_io::FastaParser &fasta_parser,
                             const std::function<void(const SeqSearchResult &)> &callback);

};


int query_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __QUERY_GRAPH_HPP__
