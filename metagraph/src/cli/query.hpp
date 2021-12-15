#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <vector>

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


class QueryExecutor {
  public:
    QueryExecutor(const Config &config,
                  const graph::AnnotatedDBG &anno_graph,
                  std::unique_ptr<graph::align::DBGAlignerConfig>&& aligner_config,
                  ThreadPool &thread_pool);

    void query_fasta(const std::string &file_path,
                     const std::function<void(const std::string &)> &callback);

    static std::string execute_query(const std::string &seq_name,
                                     const std::string &sequence,
                                     bool count_labels,
                                     bool print_signature,
                                     bool suppress_unlabeled,
                                     size_t num_top_labels,
                                     double discovery_fraction,
                                     double presence_fraction,
                                     std::string anno_labels_delimiter,
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
                             const std::function<void(const std::string &)> &callback);
};


int query_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __QUERY_GRAPH_HPP__
