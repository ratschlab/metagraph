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

std::unique_ptr<graph::AnnotatedDBG>
construct_query_graph(const graph::AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads,
                      bool canonical = false,
                      size_t sub_k = std::numeric_limits<size_t>::max(),
                      size_t max_fork_count = 0,
                      size_t max_traversal_distance = 0,
                      double max_traversed_nodes_per_seq_char = 0.0);


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
