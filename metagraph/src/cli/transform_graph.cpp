#include "transform_graph.hpp"

#include <filesystem>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "clean.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;


// Build an SSHash graph from the contigs of `graph` and save it with the row map
// that `transform_anno --remap-graph` needs to move annotations to it.
static int convert_to_sshash(Config *config, std::shared_ptr<graph::DeBruijnGraph> graph) {
    using graph::DBGSSHash;

    if (graph->get_mode() == graph::DeBruijnGraph::CANONICAL) {
        logger->error("Only graphs in basic or primary mode can be converted to SSHash");
        exit(1);
    }
    const size_t k = graph->get_k();
    const auto mode = graph->get_mode();
    // without --mode, config->graph_mode is basic
    if (config->graph_mode != graph::DeBruijnGraph::BASIC && config->graph_mode != mode) {
        logger->error("The SSHash graph keeps the mode of the input graph."
                      " Change the mode of the input graph first");
        exit(1);
    }

    auto tmp_dir = utils::create_temp_dir(
            config->tmp_dir.empty() ? std::filesystem::path(config->outfbase).remove_filename()
                                    : config->tmp_dir,
            "convert");

    // Extract the contigs as `transform --to-fasta --row-map` does
    Config extract(*config);
    extract.to_fasta = true;
    extract.row_map = true;
    extract.outfbase = tmp_dir / "contigs";
    // For a basic succinct graph, the node degrees read during the extraction give
    // the succ/pred bitmaps, which SSHash otherwise computes with k-mer lookups.
    std::pair<sdsl::bit_vector, sdsl::bit_vector> succ_pred;
    const bool shortcut = mode == graph::DeBruijnGraph::BASIC
                            && dynamic_cast<const graph::DBGSuccinct*>(graph.get());
    clean_graph(&extract, std::move(graph), shortcut ? &succ_pred : nullptr);

    Timer timer;
    try {
        DBGSSHash sshash(extract.outfbase + ".fasta.gz", k, mode,
                         config->num_chars, config->memory_available,
                         shortcut ? &succ_pred : nullptr);
        logger->trace("SSHash graph built in {} sec", timer.elapsed());
        sshash.serialize(config->outfbase);

        const std::string row_map_fname
                = utils::remove_suffix(config->outfbase, DBGSSHash::kExtension) + ".rowmap";
        std::error_code ec;
        std::filesystem::rename(extract.outfbase + ".rowmap", row_map_fname, ec);
        // the temporary directory can be on another file system
        if (ec) {
            std::filesystem::copy_file(extract.outfbase + ".rowmap", row_map_fname,
                                       std::filesystem::copy_options::overwrite_existing);
        }
        logger->trace("Row map saved to {}", row_map_fname);
    } catch (const std::exception &e) {
        logger->error("{}", e.what());
        utils::remove_temp_dir(tmp_dir);
        exit(1);
    }
    utils::remove_temp_dir(tmp_dir);

    return 0;
}


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

    if (config->graph_type == Config::GraphType::SSHASH)
        return convert_to_sshash(config, std::move(graph));

    auto dbg_succ = std::dynamic_pointer_cast<graph::DBGSuccinct>(graph);

    if (!dbg_succ.get()) {
        logger->warn("Transformations require a succinct graph. Saving the graph without transformations and exiting.");
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
