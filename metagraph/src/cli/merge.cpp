#include "merge.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_merge.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "stats.hpp"

using mtg::common::logger;
using mtg::common::get_verbose;


int merge_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    BOSS *graph = NULL;

    Timer timer;

    std::vector<std::shared_ptr<DBGSuccinct>> dbg_graphs;
    std::vector<const BOSS*> graphs;

    config->canonical = true;

    for (const auto &file : files) {
        logger->info("Opening file '{}'", file);

        dbg_graphs.emplace_back(load_critical_graph_from_file<DBGSuccinct>(file));

        graphs.push_back(&dbg_graphs.back()->get_boss());

        if (get_verbose())
            print_boss_stats(*graphs.back());

        config->canonical &= dbg_graphs.back()->is_canonical_mode();
    }

    logger->info("Graphs are loaded in {} sec", timer.elapsed());

    if (config->dynamic) {
        logger->info("Start merging traversal");
        timer.reset();

        graph = dbg_graphs.at(0)->release_boss();

        if (graph->get_state() != BOSS::State::DYN) {
            logger->trace("Switching state of succinct graph to dynamic...");

            graph->switch_state(BOSS::State::DYN);

            logger->trace("Switching done in {} sec", timer.elapsed());
        }

        for (size_t i = 1; i < graphs.size(); ++i) {
            graph->merge(dbg_graphs.at(i)->get_boss());

            logger->info("Graph '{}' merged in {} sec", files[i], timer.elapsed());

            dbg_graphs.at(i).reset();
        }
    } else if (get_num_threads() > 1 || config->parts_total > 1) {
        logger->info("Start merging blocks");
        timer.reset();

        auto *chunk = merge::merge_blocks_to_chunk(
            graphs,
            config->part_idx,
            config->parts_total,
            get_num_threads(),
            config->num_bins_per_thread,
            get_verbose()
        );
        if (!chunk) {
            logger->error("ERROR when building chunk {}", config->part_idx);
            exit(1);
        }
        logger->info("Blocks merged in {} sec", timer.elapsed());

        if (config->parts_total > 1) {
            chunk->serialize(config->outfbase
                              + "." + std::to_string(config->part_idx)
                              + "_" + std::to_string(config->parts_total));
        } else {
            graph = new BOSS(graphs[0]->get_k());
            chunk->initialize_boss(graph);
        }
        delete chunk;
    } else {
        logger->info("Start merging graphs");
        timer.reset();

        graph = merge::merge(graphs, get_verbose());
    }
    dbg_graphs.clear();

    assert(graph);

    logger->info("Graphs merged in {} sec", timer.elapsed());

    // graph output
    DBGSuccinct(graph, config->canonical).serialize(config->outfbase);

    return 0;
}
