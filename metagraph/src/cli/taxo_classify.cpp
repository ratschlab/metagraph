#include "taxo_classify.hpp"

#include "common/logger.hpp"
#include "config/config.hpp"
#include "common/unix_tools.hpp"
#include "load/load_graph.hpp"
#include "annotation/taxonomy/taxo_classifier.hpp"
#include "common/threads/threading.hpp"
#include "seq_io/sequence_io.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
                        mtg::annot::TaxoClassifier &taxo_classifier,
                        const std::function<void(const std::string &, const uint64_t &)> &callback,
                        const Config &config) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file);

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        thread_pool.enqueue([&](const auto&... args){
            uint64_t taxid = taxo_classifier.assign_class(graph,
                                                          std::string(kseq.seq.s),
                                                          config.lca_coverage_threshold);
            callback(std::string(kseq.name.s), taxid);
        });
    }
}

int taxonomic_classification(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    Timer timer;
    logger->trace("Loading TaxonomyDB ...");
    mtg::annot::TaxoClassifier taxo_classifier(config->taxonomic_tree);
    logger->trace("Finished loading TaxonomyDB in {}", timer.elapsed());

    timer.reset();
    logger->trace("Graph loading ...");
    auto graph = load_critical_dbg(config->infbase);
    logger->trace("Finished graph loading in {}", timer.elapsed());

    timer.reset();
    logger->trace("Processing the classification");
    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);
    for (const auto &file: files) {
        Timer curr_timer;
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
                           taxo_classifier,
                           [](const std::string name_seq, const uint64_t &taxid) {
                               std::string result = fmt::format("Sequence '{}' was classified with Tax ID '{}'\n",
                                                                name_seq,
                                                                taxid);
                               std::cout << result << std::endl;
                           },
                           *config);
        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }
    logger->trace("Finished processing all the classification in {}", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
