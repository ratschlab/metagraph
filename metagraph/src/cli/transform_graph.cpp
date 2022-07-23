#include "transform_graph.hpp"

#include <filesystem>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_rc.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;


int transform_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    if (config->initialize_bloom)
        std::filesystem::remove(utils::make_suffix(config->outfbase, ".bloom"));

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    auto dbg_succ = std::dynamic_pointer_cast<graph::DBGSuccinct>(graph);

    if (!dbg_succ.get())
        throw std::runtime_error("Only implemented for DBGSuccinct");

    if (config->adjrc) {
        graph::NodeRC(*dbg_succ, true).serialize(config->outfbase + dbg_succ->file_extension());
        return 0;
    }

    if (config->initialize_bloom) {
        assert(config->bloom_fpp > 0.0 && config->bloom_fpp <= 1.0);
        assert(config->bloom_bpk >= 0.0);
        assert(config->bloom_fpp < 1.0 || config->bloom_bpk > 0.0);

        logger->trace("Construct Bloom filter for nodes...");

        timer.reset();

        if (config->bloom_fpp < 1.0) {
            dbg_succ->initialize_bloom_filter_from_fpr(
                config->bloom_fpp,
                config->bloom_max_num_hash_functions
            );
        } else {
            dbg_succ->initialize_bloom_filter(
                config->bloom_bpk,
                config->bloom_max_num_hash_functions
            );
        }

        logger->trace("Bloom filter constructed in {} sec", timer.elapsed());

        assert(dbg_succ->get_bloom_filter());

        auto fname = utils::make_suffix(config->outfbase, dbg_succ->bloom_filter_file_extension());
        std::ofstream bloom_outstream(fname, std::ios::binary);
        if (!bloom_outstream)
            throw std::ios_base::failure("Can't write to file " + fname);

        dbg_succ->get_bloom_filter()->serialize(bloom_outstream);

        return 0;
    }

    if (config->clear_dummy) {
        logger->trace("Traverse the tree of source dummy edges and remove redundant ones...");
        timer.reset();

        // remove redundant dummy edges and mark all other dummy edges
        dbg_succ->mask_dummy_kmers(get_num_threads(), true);

        logger->trace("The tree of source dummy edges traversed in {} sec", timer.elapsed());
        timer.reset();
    }

    if (config->node_suffix_length != dbg_succ->get_boss().get_indexed_suffix_length()) {
        size_t suffix_length = std::min((size_t)config->node_suffix_length,
                                        dbg_succ->get_boss().get_k());

        if (suffix_length * log2(dbg_succ->get_boss().alph_size - 1) > 63) {
            logger->error("Node ranges for k-mer suffixes longer than {} cannot be indexed",
                          static_cast<int>(63 / log2(dbg_succ->get_boss().alph_size - 1)));
            exit(1);
        }

        logger->trace("Index all node ranges for suffixes of length {} in {:.2f} MB",
                      suffix_length,
                      std::pow(dbg_succ->get_boss().alph_size - 1, suffix_length)
                            * 2. * sizeof(uint64_t) * 1e-6);
        logger->trace("Compressed node ranges to {:.2f} MB",
                      dbg_succ->get_boss().get_suffix_ranges_index_size() / 8e6);
        timer.reset();

        dbg_succ->get_boss().index_suffix_ranges(suffix_length, get_num_threads());

        logger->trace("Indexing of node ranges took {} sec", timer.elapsed());
        timer.reset();
    }

    if (config->to_adj_list) {
        logger->trace("Converting graph to adjacency list...");

        auto *boss = &dbg_succ->get_boss();
        timer.reset();

        std::ofstream outstream(config->outfbase + ".adjlist");
        boss->print_adj_list(outstream);

        logger->trace("Conversion done in {} sec", timer.elapsed());

        return 0;
    }

    if (config->graph_mode == graph::DeBruijnGraph::PRIMARY
            && dbg_succ->get_mode() == graph::DeBruijnGraph::BASIC) {
        logger->info("Changing graph mode from basic to primary");
        logger->warn("FYI: This doesn't rebuild the graph. Apply with caution"
                     " and only to graphs constructed from primary contigs!");
        // keep the graph state (representation) unchanged
        config->state = dbg_succ->get_state();
        graph::boss::BOSS* boss = dbg_succ->release_boss();
        dbg_succ.reset(new graph::DBGSuccinct(boss, graph::DeBruijnGraph::PRIMARY));
        logger->info("Graph mode changed to primary");
    }

    if (config->state != dbg_succ->get_state()) {
        logger->trace("Converting graph to state {}", Config::state_to_string(config->state));
        timer.reset();

        dbg_succ->switch_state(config->state);

        logger->trace("Conversion done in {} sec", timer.elapsed());
    }

    logger->trace("Serializing transformed graph...");
    timer.reset();
    dbg_succ->serialize(config->outfbase);
    logger->trace("Serialization done in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
