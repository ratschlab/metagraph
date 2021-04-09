#include "tax_classify.hpp"

#include "annotation/taxonomy/tax_classifier.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

uint64_t QUERY_SEQ_BATCH_SIZE = 10000;

void execute_fasta_seq(const std::string sequence,
                       const std::string seq_label,
                       const mtg::graph::DeBruijnGraph &graph,
                       const mtg::annot::TaxClassifier &tax_classifier,
                       std::mutex &result_mutex,
                       std::vector<std::pair<std::string, uint64_t> > *results) {
    uint64_t taxid = tax_classifier.assign_class(graph, sequence);
    std::lock_guard<std::mutex> lock(result_mutex);
    results->push_back(std::pair<std::string, uint64_t>{seq_label, taxid});
}

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
                        const mtg::annot::TaxClassifier &tax_classifier,
                        std::mutex &result_mutex,
                        std::vector<std::pair<std::string, uint64_t> > *results) {
    logger->trace("Parsing query sequences from file {}.", file);

    seq_io::FastaParser fasta_parser(file);

    std::vector<std::pair<std::string, std::string> > seq_batch;

    uint64_t cnt_queries_executed = 0;
    for (const seq_io::kseq_t &kseq : fasta_parser) {
        seq_batch.push_back({std::string(kseq.seq.s), std::string(kseq.name.s)});

        if (seq_batch.size() != QUERY_SEQ_BATCH_SIZE) {
            continue;
        }
        thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
            for (std::pair<std::string, std::string> &seq : sequences) {
                execute_fasta_seq(seq.first, seq.second, graph, tax_classifier, result_mutex, results);
            }
        }, std::move(seq_batch));

        cnt_queries_executed ++;
        if (cnt_queries_executed % 100000 == 0) {
            logger->trace("Started to executed the first {} queries from file {}.", cnt_queries_executed, file);
        }
        seq_batch.clear();
    }

    thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string>> sequences){
        for (std::pair<std::string, std::string> &seq : sequences) {
            execute_fasta_seq(seq.first, seq.second, graph, tax_classifier, result_mutex, results);
        }
    }, std::move(seq_batch));
}

int taxonomic_classification(Config *config) {
    assert(config);

    const std::vector<std::string> &files = config->fnames;

    Timer timer;
    logger->trace("Loading TaxonomyDB...");
    mtg::annot::TaxClassifier tax_classifier(config->taxonomic_tree,
                                             config->lca_coverage_fraction,
                                             config->discovery_fraction);
    logger->trace("Finished loading TaxonomyDB in {}s.", timer.elapsed());

    timer.reset();
    logger->trace("Graph loading...");
    auto graph = load_critical_dbg(config->infbase);
    logger->trace("Finished graph loading in {}s.", timer.elapsed());

    timer.reset();
    logger->trace("Processing the classification...");
    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);
    std::vector<std::pair<std::string, uint64_t> > results;
    std::mutex result_mutex;
    for (const std::string &file : files) {
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
                           tax_classifier,
                           result_mutex,
                           &results);
    }
    thread_pool.join();
    for (const std::pair<std::string, uint64_t> &result : results) {
        std::string output = fmt::format("Sequence '{}' was classified with Tax ID '{}'\n",
                                         result.first, result.second);
        std::cout << output << std::endl;
    }
    logger->trace("Finished all the queries in {}s.", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
