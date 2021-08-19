#include "tax_class.hpp"

#include "annotation/taxonomy/tax_classifier.hpp"
#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"

#include "graph/annotated_dbg.hpp"

#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "load/load_annotation.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

const uint32_t QUERY_SEQ_BATCH_SIZE = 100000;

void append_new_result(const std::string &seq_label,
                       const uint32_t taxid,
                       std::vector<std::pair<std::string, uint32_t> > *pair_label_taxid,
                       std::mutex *tax_mutex) {
    std::scoped_lock<std::mutex> guard(*tax_mutex);
    (*pair_label_taxid).emplace_back(seq_label, taxid);
}

void print_all_results(const std::vector<std::pair<std::string, uint32_t> > &pair_label_taxid,
                       const std::function<void(const std::string &, const uint32_t &)> &callback) {
    for (const std::pair<std::string, uint32_t> &label_taxid : pair_label_taxid) {
        callback(label_taxid.first, label_taxid.second);
    }
}

void execute_fasta_file(const string &file,
                        std::function<void(const std::vector<std::pair<std::string, std::string> > &)> &callback) {
    logger->trace("Parsing query sequences from file {}.", file);

    seq_io::FastaParser fasta_parser(file);
    std::vector<std::pair<std::string, std::string> > seq_batch;

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        seq_batch.push_back({std::string(kseq.seq.s), std::string(kseq.name.s)});

        if (seq_batch.size() != QUERY_SEQ_BATCH_SIZE) {
            continue;
        }
        callback(seq_batch);

        logger->trace("Processing an another bucket of {} queries from file {}.", QUERY_SEQ_BATCH_SIZE, file);
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
    logger->trace("Finished graph loading after {} sec.", timer.elapsed());

    timer.reset();
    logger->trace("Processing the classification...");
    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    std::function<void(const std::vector<std::pair<std::string, std::string> > &)> callback;

    std::vector<std::pair<std::string, uint32_t> > pair_label_taxid;
    std::mutex tax_mutex;

    std::unique_ptr<mtg::annot::TaxonomyClsAnno> taxonomy;
    std::unique_ptr<graph::AnnotatedDBG> anno_graph;


    timer.reset();
    logger->trace("Graph and Annotation loading...");
    // graph = load_critical_dbg(config->infbase);
    anno_graph = initialize_annotated_dbg(graph, *config);
    logger->trace("Finished graph annotation loading after {} sec.", timer.elapsed());


    std::shared_ptr<mtg::graph::DBGBitmap> graph2 = load_critical_dbg_dd(config->infbase2);
//    config->infbase = config->infbase2;
//    config->infbase_annotators = config->infbase_annotators2;

    uint64_t max_index_anno2 = graph2->max_index();

    std::unique_ptr<annot::RowCompressed<std::string>> anno_graph2 = initialize_annotation_dd(config->infbase_annotators2.at(0), max_index_anno2);
    bool loaded = anno_graph2->load(config->infbase_annotators2.at(0));
    if (!loaded) {
        logger->error("Cannot load annotations2, file corrupted");
        exit(1);
    }

    if (config->taxonomic_db != "") {
        throw std::runtime_error("internal error: taxonomic classification with taxDB is not implemented.");
    } else {
        // Use tax_class without any precomputed database.
        if (config->infbase_annotators.size() == 0) {
            logger->error("The annotation matrix is missing from the command line, "
                          "please use '-a' flag for the annotation matrix filepath.");
            std::exit(1);
        }
        std::cerr << "\n\nbefore load anno_matrix2\n" << std::endl;

        
        // std::shared_ptr<mtg::graph::DBGSuccinct> graph2 = std::make_shared<mtg::graph::DBGSuccinct>(graph2_dbg);


        std::cerr << "\n\nafter load anno_matrix2\n" << std::endl;

        timer.reset();
        logger->trace("Constructing TaxonomyClsAnno...");
        taxonomy = std::make_unique<annot::TaxonomyClsAnno>(*anno_graph, config->lca_coverage_fraction,
                                                            config->discovery_fraction, config->taxonomic_tree,
                                                            config->label_taxid_map);
        logger->trace("Finished TaxonomyDB construction after {} sec.", timer.elapsed());

        if (config->top_label_fraction > 0) {
            // Use tax_class version that is returning the LCA of the top labels among the kmers.
            // This version is fast, but less precise.
            callback = [&](const std::vector<std::pair<std::string, std::string> > &seq_batch){
                thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
                    for (std::pair<std::string, std::string> &seq : sequences) {
                        append_new_result(seq.second, taxonomy->assign_class_toplabels(
                                seq.first, config->top_label_fraction), &pair_label_taxid, &tax_mutex);
                    }
                }, std::move(seq_batch));
            };
        } else {
            // Use tax_class version that computest the LCA taxid for each kmer.
            // The prediction result will be identical to the one using taxdb, but the computations will be slower.
            callback = [&](const std::vector<std::pair<std::string, std::string> > &seq_batch){
                thread_pool.enqueue([&](std::vector<std::pair<std::string, std::string> > sequences){
                    for (std::pair<std::string, std::string> &seq : sequences) {
                        append_new_result(seq.second, taxonomy->assign_class(seq.first, *graph2, *anno_graph2, seq.second), &pair_label_taxid, &tax_mutex);
                    }
                }, std::move(seq_batch));
            };
        }
    }

    for (const std::string &file : files) {
        execute_fasta_file(file, callback);
    }
    thread_pool.join();
    uint64_t num_hits = 0;

    print_all_results(pair_label_taxid, [&](const std::string name_seq, const uint32_t &taxid) {
        std::string result = fmt::format(
                "Sequence '{}' was classified with Tax ID '{}'\n",
                name_seq, taxid);
        std::cout << result << std::endl;
        if (utils::split_string(name_seq, "|")[1] == to_string(taxid)) {
            num_hits += 1;
        }
    });

    logger->trace("Finished all the queries in {} sec.", timer.elapsed());

    std::cerr << "num hits = " << num_hits << "\n total results =" << pair_label_taxid.size() << "\n";
    std::cerr << "hit rate = " << (double)num_hits / pair_label_taxid.size() << "\n";

    return 0;
}

} // namespace cli
} // namespace mtg
