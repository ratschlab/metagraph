#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>

class AnnotatedDBG;
class IDBGAligner;
class Config;
class ThreadPool;

namespace mtg {
namespace seq_io {
class FastaParser;
} // namespace seq_io
} // namespace mtg

using StringGenerator = std::function<void(std::function<void(const std::string &)>)>;

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads);


class QueryExecutor {
  public:
    QueryExecutor(const Config &config,
                  const AnnotatedDBG &anno_graph,
                  const IDBGAligner *aligner,
                  ThreadPool &thread_pool)
      : config_(config),
        anno_graph_(anno_graph),
        aligner_(aligner),
        thread_pool_(thread_pool) {}

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
                                     const AnnotatedDBG &anno_graph);

  private:
    const Config &config_;
    const AnnotatedDBG &anno_graph_;
    const IDBGAligner *aligner_;
    ThreadPool &thread_pool_;

    void batched_query_fasta(mtg::seq_io::FastaParser &fasta_parser,
                             const std::function<void(const std::string &)> &callback);
};


int query_graph(Config *config);

#endif // __QUERY_GRAPH_HPP__
