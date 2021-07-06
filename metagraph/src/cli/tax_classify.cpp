#include "tax_classify.hpp"

#include "annotation/taxonomy/tax_classifier.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "seq_io/sequence_io.hpp"

#include "annotation/taxonomy/taxonomic_db.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

const uint64_t QUERY_SEQ_BATCH_SIZE = 1000;

void append_new_result(const std::string &seq_label,
                       const uint64_t taxid,
                       std::vector<std::pair<std::string, uint64_t> > *pair_label_taxid,
                       std::mutex *tax_mutex) {
    std::scoped_lock<std::mutex> guard(*tax_mutex);
    (*pair_label_taxid).emplace_back(seq_label, taxid);
}

void print_all_results(const std::vector<std::pair<std::string, uint64_t> > &pair_label_taxid,
                       const std::function<void(const std::string &, const uint64_t &)> &callback) {
    for (const std::pair<std::string, uint64_t> &label_taxid : pair_label_taxid) {
        callback(label_taxid.first, label_taxid.second);
    }
}

void execute_fasta_file(const string &file,
                        std::function<void(const std::vector<std::pair<std::string, std::string> > &)> &callback) {
    logger->trace("Parsing query sequences from file {}.", file);

    seq_io::FastaParser fasta_parser(file);
    std::vector<std::pair<std::string, std::string> > seq_batch;

    uint64_t cnt_queries_started = 0;
    for (const seq_io::kseq_t &kseq : fasta_parser) {
        seq_batch.push_back({std::string(kseq.seq.s), std::string(kseq.name.s)});

        if (seq_batch.size() != QUERY_SEQ_BATCH_SIZE) {
            continue;
        }
        callback(seq_batch);

        cnt_queries_started ++;
        if (cnt_queries_started % 100000 == 0) {
            logger->trace("Processing the first {} queries from file {}.", cnt_queries_started, file);
        }
        seq_batch.clear();
    }
    callback(seq_batch);
}

int taxonomic_classification(Config *config) {
    assert(config);

    const std::vector<std::string> &files = config->fnames;

    Timer timer;
    logger->trace("Graph loading...");
    auto graph = load_critical_dbg(config->infbase);
    logger->trace("Finished graph loading after {}s.", timer.elapsed());

    timer.reset();
    logger->trace("Processing the classification...");
    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    std::function<void(const std::vector<std::pair<std::string, std::string> > &)> callback;

    std::vector<std::pair<std::string, uint64_t> > pair_label_taxid;
    std::mutex tax_mutex;

    std::unique_ptr<mtg::annot::TaxClassifier> tax_classifier;
    std::unique_ptr<graph::AnnotatedDBG> anno_graph;
    std::unique_ptr<annot::TaxonomyDB> taxonomy;

    if (config->taxonomic_db != "") {
        // Given a taxdb filepath, use the tax_class version with taxonomic database.
        timer.reset();
        logger->trace("Loading TaxonomyDB...");
        tax_classifier = std::make_unique<annot::TaxClassifier>(config->taxonomic_db,
                                                                config->lca_coverage_fraction,
                                                                config->discovery_fraction);

        logger->trace("Finished loading TaxonomyDB in {}s.", timer.elapsed());
        callback = [&](const std::vector<std::pair<std::string, std::string> > &seq_batch){
                thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
                    for (std::pair<std::string, std::string> &seq : sequences) {
                        append_new_result(seq.second, tax_classifier->assign_class(*graph, seq.first),
                                          &pair_label_taxid, &tax_mutex);
                    }
                }, std::move(seq_batch));
        };
    } else {
        // Use tax_class without any precomputed database.
        if (config->infbase_annotators.size() == 0) {
            logger->error("The annotation matrix is missing from the command line, please use '-a' flag for the annotation matrix filepath.");
            std::exit(1);
        }
        timer.reset();
        logger->trace("Graph and Annotation loading...");
        graph = load_critical_dbg(config->infbase);
        anno_graph = initialize_annotated_dbg(graph, *config);
        logger->trace("Finished graph annotation loading after {}s.", timer.elapsed());

        timer.reset();
        logger->trace("Constructing TaxonomyDB...");
        // Create a 'TaxonomyDB' object here. The constructor will initialize all the taxonomic tree statistics that are further used by taxonomy.classify().
        taxonomy = std::make_unique<annot::TaxonomyDB>(config->taxonomic_tree, "", tsl::hopscotch_set<std::string>{});
        logger->trace("Finished TaxonomyDB construction after {}s.", timer.elapsed());

        if (config->top_label_fraction > 0) {
            // Use tax_class version that is returning the LCA of the top labels among the kmers.
            // This version is fast, but less precise.
            callback = [&](const std::vector<std::pair<std::string, std::string> > &seq_batch){
                thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
                    for (std::pair<std::string, std::string> &seq : sequences) {
                        append_new_result(seq.second, taxonomy->assign_class_toplabels(*anno_graph, seq.first, config->top_label_fraction),
                                          &pair_label_taxid, &tax_mutex);
                    }
                }, std::move(seq_batch));
            };
        } else {
            // Use tax_class version that computest the LCA taxid for each kmer.
            // The prediction result will be identical to the one using taxdb, but the computations will be slower.
            callback = [&](const std::vector<std::pair<std::string, std::string> > &seq_batch){
                thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
                    for (std::pair<std::string, std::string> &seq : sequences) {
                        append_new_result(seq.second, taxonomy->assign_class_getrows(*anno_graph, seq.first,
                                                                                     config->lca_coverage_fraction,
                                                                                     config->discovery_fraction),
                                          &pair_label_taxid, &tax_mutex);
                    }
                }, std::move(seq_batch));
            };
        }
    }

    for (const std::string &file : files) {
        execute_fasta_file(file, callback);
    }
    thread_pool.join();

    print_all_results(pair_label_taxid, [](const std::string name_seq, const uint64_t &taxid) {
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
