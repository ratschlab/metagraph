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

std::vector<std::pair<std::string, uint64_t> > pair_label_taxid;
std::mutex tax_mutex;

void append_new_result(const std::string &seq_label, const uint64_t taxid) {
    std::scoped_lock<std::mutex> guard(tax_mutex);
    pair_label_taxid.emplace_back(seq_label, taxid);
}

void print_all_results(const std::function<void(const std::string &, const uint64_t &)> &callback) {
    for (const std::pair<std::string, uint64_t> &label_taxid : pair_label_taxid) {
        callback(label_taxid.first, label_taxid.second);
    }
}

void run_sequence_batch(const std::vector<std::pair<std::string, std::string> > &seq_batch,
                        const mtg::annot::TaxClassifier &tax_classifier,
                        const mtg::graph::DeBruijnGraph &graph,
                        ThreadPool &thread_pool) {
    thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
        for (std::pair<std::string, std::string> &seq : sequences) {
            append_new_result(seq.second, tax_classifier.assign_class(graph, seq.first));
        }
    }, std::move(seq_batch));
}

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
                        const mtg::annot::TaxClassifier &tax_classifier) {
    logger->trace("Parsing query sequences from file {}.", file);

    seq_io::FastaParser fasta_parser(file);

    std::vector<std::pair<std::string, std::string> > seq_batch;

    uint64_t cnt_queries_started = 0;
    for (const seq_io::kseq_t &kseq : fasta_parser) {
        seq_batch.push_back({std::string(kseq.seq.s), std::string(kseq.name.s)});

        if (seq_batch.size() != QUERY_SEQ_BATCH_SIZE) {
            continue;
        }

        run_sequence_batch(seq_batch, tax_classifier, graph, thread_pool);

        cnt_queries_started ++;
        if (cnt_queries_started % 100000 == 0) {
            logger->trace("Processing the first {} queries from file {}.", cnt_queries_started, file);
        }
        seq_batch.clear();
    }
    run_sequence_batch(seq_batch, tax_classifier, graph, thread_pool);
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
    for (const std::string &file : files) {
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
                           tax_classifier);
    }
    thread_pool.join();

    std::cout << "time_spent_map_to_nodes = " << tax_classifier.time_spent_map_to_nodes << "\n";
    std::cout << "algorithm               = " << tax_classifier.time_spent_assign_class - tax_classifier.time_spent_map_to_nodes << "\n";
    std::cout << "time_spent_assign_class = " << tax_classifier.time_spent_assign_class << "\n";

    print_all_results(
            [](const std::string name_seq, const uint64_t &taxid) {
                std::string result = fmt::format(
                      "Sequence '{}' was classified with Tax ID '{}'\n",
                      name_seq, taxid);
                std::cout << result << std::endl;
            });

    logger->trace("Finished all the queries in {}s.", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
