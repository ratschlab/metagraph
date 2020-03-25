#ifndef __QUERY_GRAPH_HPP__
#define __QUERY_GRAPH_HPP__

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>

#include "common/threads/threading.hpp"
#include "seq_io/sequence_io.hpp"

class AnnotatedDBG;
class IDBGAligner;
class Config;

std::string execute_query(const std::string &seq_name,
                          const std::string &sequence,
                          bool count_labels,
                          bool print_signature,
                          bool suppress_unlabeled,
                          size_t num_top_labels,
                          double discovery_fraction,
                          std::string anno_labels_delimiter,
                          const AnnotatedDBG &anno_graph);


using StringGenerator = std::function<void(std::function<void(const std::string &)>)>;

std::unique_ptr<AnnotatedDBG> construct_query_graph(const AnnotatedDBG &anno_graph,
                                                    StringGenerator call_sequences,
                                                    double discovery_fraction,
                                                    size_t num_threads);

int query_graph(Config *config);

class QueryExecutor {
  public:
    void process_fasta_file(const std::string &file_path, std::ostream &out);

    QueryExecutor(const Config *config,
                  const AnnotatedDBG *anno_graph,
                  const IDBGAligner *aligner,
                  ThreadPool *thread_pool)
        : config_(config),
          anno_graph_(anno_graph),
          aligner_(aligner),
          thread_pool_(thread_pool) {}

  private:
    const Config *config_;
    const AnnotatedDBG *anno_graph_;
    const IDBGAligner *aligner_;
    ThreadPool *thread_pool_;

    void process_fasta_file_fast(FastaParser &fasta_parser,
                                 const std::string &file_name,
                                 std::ostream &out);

    void forward_query(size_t id,
                       const std::string &name,
                       const std::string &seq,
                       const AnnotatedDBG *graph_to_query,
                       std::ostream &out,
                       std::mutex &stream_mutex);
};


#endif // __QUERY_GRAPH_HPP__
