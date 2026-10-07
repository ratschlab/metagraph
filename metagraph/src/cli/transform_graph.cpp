#include "transform_graph.hpp"

#include <filesystem>
#include <fstream>

#include <unistd.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

namespace {

/**
 * transform --mask-dummy (DESIGN-pattern-search.md §4): give an existing succinct graph the
 * dummy-edge mask `build --mask-dummy` would have written, as the file <graph>.edgemask that
 * DBGSuccinct::load reads beside it. Only that file is written: the mask prunes nothing, so
 * node ids and the annotation stay valid, and a graph of hundreds of GB is not rewritten for a
 * file a small fraction of its size. Without a mask every dummy edge passes for a k-mer, and
 * the pattern search cannot count (mask_required).
 */
int write_dummy_mask(const std::string &graph_path, const Config &config) {
    using graph::DBGSuccinct;

    if (!utils::ends_with(graph_path, DBGSuccinct::kExtension)) {
        logger->error("--mask-dummy: '{}' is not a succinct graph ({}): only a succinct graph "
                      "has a dummy-edge mask", graph_path, DBGSuccinct::kExtension);
        return 1;
    }
    const std::string prefix = utils::remove_suffix(graph_path, DBGSuccinct::kExtension);
    const std::string mask_path = prefix + DBGSuccinct::kDummyMaskExtension;
    const bool replacing = std::filesystem::exists(mask_path);
    if (replacing && !config.force) {
        // a mask written by build or by an earlier run is not overwritten unseen
        logger->error("--mask-dummy: {} exists: the graph has its dummy-edge mask already; "
                      "pass --force to build it again and replace it", mask_path);
        return 1;
    }

    Timer timer;
    logger->info("Loading the graph {} without its mask...", graph_path);
    // without the mask: a replaced one is not read (it may be the broken file being replaced)
    DBGSuccinct graph(2);
    bool loaded = false;
    try {
        loaded = graph.load_without_mask(graph_path);
    } catch (const std::exception &e) {
        logger->error("Cannot load graph from '{}': {}", graph_path, e.what());
        return 1;
    }
    if (!loaded) {
        logger->error("Cannot load graph from '{}': {}", graph_path,
                      utils::file_read_failure_detail(graph_path));
        return 1;
    }
    logger->info("Graph loaded in {:.3f} s: k = {}, mode {}, {} edges", timer.elapsed(),
                 graph.get_k(), Config::graphmode_to_string(graph.get_mode()),
                 graph.max_index());

    DummyMaskCounts counts;
    try {
        counts = mask_dummy_edges(&graph, get_num_threads());
    } catch (const std::exception &e) {
        // (memory, above all, on a large graph): nothing was written
        logger->error("--mask-dummy: the mask could not be built: {}", e.what());
        return 1;
    }
    logger->info("Dummy edges marked in {:.3f} s with {} threads: {} edges, {} source dummies "
                 "(the main dummy edge included), {} sink dummies, {} k-mers",
                 counts.seconds, get_num_threads(), counts.edges, counts.source_dummy,
                 counts.sink_dummy, counts.kmers);

    // written to a temporary file in the same directory and renamed over the target, so that a
    // loader never reads a partial mask: DBGSuccinct::load refuses the graph with a mask it
    // cannot read, and an interrupted run would otherwise leave such a file beside the graph
    timer.reset();
    const std::string tmp_path = mask_path + ".tmp." + std::to_string(getpid());
    try {
        {
            std::ofstream out = utils::open_new_ofstream(tmp_path);
            if (!out.good())
                throw std::ios_base::failure("cannot open " + tmp_path + " for writing");
            graph.get_mask()->serialize(out);
            out.close();
            if (!out)
                throw std::ios_base::failure("cannot write " + tmp_path);
        }
        std::filesystem::rename(tmp_path, mask_path);
    } catch (const std::exception &e) {
        std::error_code ignored;
        std::filesystem::remove(tmp_path, ignored);
        logger->error("--mask-dummy: the mask was not written: {}", e.what());
        return 1;
    }
    logger->info("{} {} ({} bytes) in {:.3f} s; {} is unchanged", replacing ? "Replaced" : "Wrote",
                 mask_path, std::filesystem::file_size(mask_path), timer.elapsed(), graph_path);

    // what the new file changes for the graph's other readers, which load it from now on
    logger->info("Every loader of {} reads the mask from now on: the graph states its k-mers "
                 "as nodes ({} instead of {} edges), and an index manifest "
                 "(--index-manifest) must list the new file", graph_path, counts.kmers,
                 counts.edges);
    const std::string bloom_path = prefix + DBGSuccinct::kBloomFilterExtension;
    if (std::filesystem::exists(bloom_path)) {
        // DBGSuccinct::load reads the Bloom filter only together with the mask
        logger->warn("{} exists beside the graph: DBGSuccinct::load reads a Bloom filter only "
                     "when the mask is present, so it is loaded with the graph from now on",
                     bloom_path);
    }
    return 0;
}

} // namespace


int transform_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);

    if (config->mark_dummy_kmers)
        return write_dummy_mask(files.at(0), *config);

    assert(config->outfbase.size());

    if (config->initialize_bloom)
        std::filesystem::remove(utils::make_suffix(config->outfbase, ".bloom"));

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    auto dbg_succ = std::dynamic_pointer_cast<graph::DBGSuccinct>(graph);

    if (!dbg_succ.get()) {
        logger->warn("Transformations only implemented for DBGSuccinct, serializing graph and exiting");
        graph->serialize(config->outfbase);
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
        std::ofstream bloom_out = utils::open_new_ofstream(fname);
        if (!bloom_out)
            throw std::ios_base::failure("Can't write to file " + fname);

        dbg_succ->get_bloom_filter()->serialize(bloom_out);

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
        timer.reset();
        dbg_succ->get_boss().index_suffix_ranges(suffix_length, get_num_threads());
        logger->trace("Indexing of node ranges took {} sec", timer.elapsed());
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
