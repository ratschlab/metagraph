#include "classify.hpp"

#include "common/logger.hpp"
#include "config/config.hpp"
#include "common/unix_tools.hpp"
#include "load/load_graph.hpp"
#include "annotation/taxonomic/classifier.hpp"
#include "common/threads/threading.hpp"
#include "seq_io/sequence_io.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

void execute_fasta_file(const string &file,
                        ThreadPool &thread_pool,
                        const mtg::graph::DeBruijnGraph &graph,
                        mtg::annot::Classifier &taxo_classifier,
                        const std::function<void(const std::string &)> &callback) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file);
    std::cout << " finished with fastaParser construction\n";

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        std::cout << "fasta for ->" << kseq.name.s << "\n";
        std::cout << std::string(kseq.seq.s) << endl;
        thread_pool.enqueue([&](const auto&... args){
            // TODO add this lca_coverage_threshold as optional param
            std::cout << "inside thread pool?" << endl;
//          const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);
//          graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), canonical);

            std::string label = taxo_classifier.assign_class(graph,
                                                             std::string(kseq.seq.s),
                                                             0.66);
            std::cout << "label = " << label << std::endl;
            callback(label);
        });
    }
    std::cout << "last line in execute_fasta_file" << std::endl;
}

int classify(Config *config) {
    assert(config);

    const auto &files = config->fnames;
//    assert(config->infbase_annotators.size() == 1);

    Timer timer;
    logger->trace("Loading taxonomy data ...", timer.elapsed());
    mtg::annot::Classifier taxo_classifier(config->taxonomic_tree );
    std::cout << "after Classifier constructor" << std::endl;
//    config->infbase_annotators.pop_back();

    logger->trace("Graph loading ...");
    auto graph = load_critical_dbg(config->infbase);
//    auto anno_graph = initialize_annotated_dbg(graph, *config);


    logger->trace("Processing the classification", timer.elapsed());

//    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);
    ThreadPool thread_pool(1u);
    for (const auto &file: files) {
        Timer curr_timer;
        execute_fasta_file(file,
                           thread_pool,
                           *graph,
                           taxo_classifier,
                           [](const std::string &result) { std::cout << result << std::endl; });
        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}

} // namespace cli
} // namespace mtg
