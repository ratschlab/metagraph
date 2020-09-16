#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>

class ThreadPool;

namespace mtg {

namespace seq_io {
    class FastaParser;
}

namespace graph {
    class AnnotatedDBG;
    namespace align {
        class IDBGAligner;
        class DBGAlignerConfig;
    }
}


namespace cli {

class Config;

using StringGenerator = std::function<void(std::function<void(const std::string &)>)>;

/**
 * Construct a query graph
 * @param anno_graph the input annotated de Bruijn graph
 * @param call_sequences generate sequences to be queried against anno_graph
 * @param num_threads number of threads to use
 * @param canonical if true, the returned query graph is a canonical graph
 * @param config a pointer to a Config to determine construction parameters
 */
std::unique_ptr<graph::AnnotatedDBG>
construct_query_graph(const graph::AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      size_t num_threads,
                      bool canonical = false,
                      const Config *config = nullptr);


class QueryExecutor {
  public:
    QueryExecutor(const Config &config,
                  const graph::AnnotatedDBG &anno_graph,
                  const graph::align::DBGAlignerConfig *aligner_config,
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
                                     std::string anno_labels_delimiter,
                                     const graph::AnnotatedDBG &anno_graph);

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
