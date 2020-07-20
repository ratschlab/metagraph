#include "augment.hpp"

#include <string>
#include <vector>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "graph/annotated_dbg.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotation.hpp"
#include "config/config.hpp"
#include "parse_sequences.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;


int augment_graph(Config *config) {
    assert(config);
    assert(config->infbase_annotators.size() <= 1);

    const auto &files = config->fnames;

    Timer timer;

    // load graph
    auto graph = load_critical_dbg(config->infbase);

    auto node_weights = graph->load_extension<graph::NodeWeights>(config->infbase);
    // TODO: fix extension of DBGSuccinct with k-mer counts
    //       DBGSuccinct with mask of dummy edges initialized uses
    //       contiguous indexes that are not compatible with node weights,
    //       which are indexed by the rows of the BOSS table.
    //       This can be fixed by using the same indexes in all cases
    //       (non-contiguous indexing)
    if (!node_weights->is_compatible(*graph)) {
        logger->error("Node weights are not compatible with graph '{}' "
                      "and will not be updated", config->infbase);
        node_weights.reset();
    }

    logger->trace("De Bruijn graph with k-mer length k={} was loaded in {} sec",
                  graph->get_k(), timer.elapsed());
    timer.reset();

    if (dynamic_cast<graph::DBGSuccinct*>(graph.get())) {
        auto &succinct_graph = dynamic_cast<graph::DBGSuccinct&>(*graph);

        if (succinct_graph.get_state() != graph::BOSS::State::DYN) {
            logger->trace("Switching state of succinct graph to dynamic...");

            succinct_graph.switch_state(graph::BOSS::State::DYN);

            logger->trace("State switching done in {} sec", timer.elapsed());
        }
    }

    std::unique_ptr<bit_vector_dyn> inserted_nodes;
    if (config->infbase_annotators.size() || node_weights)
        inserted_nodes.reset(new bit_vector_dyn(graph->max_index() + 1, 0));

    timer.reset();

    logger->trace("Start graph augmentation");

    if (graph->is_canonical_mode())
        config->forward_and_reverse = false;

    config->canonical = graph->is_canonical_mode();

    std::function<void(uint64_t)> on_node_insert = [](uint64_t) {};
    if (inserted_nodes)
        on_node_insert = [&](uint64_t new_node) { inserted_nodes->insert_bit(new_node, 1); };

    for (const auto &file : files) {
        parse_sequences(file, *config,
            [&](std::string_view seq) {
                graph->add_sequence(seq, on_node_insert);
            },
            [&](std::string_view seq, uint32_t /*count*/) {
                graph->add_sequence(seq, on_node_insert);
            }
        );
        logger->trace("Extracted all sequences from file '{}' in {} sec",
                      file, timer.elapsed());
    }
    assert(!inserted_nodes || inserted_nodes->size() == graph->max_index() + 1);

    logger->trace("Graph augmentation done in {} sec", timer.elapsed());
    timer.reset();

    if (node_weights) {
        node_weights->insert_nodes(*inserted_nodes);

        assert(node_weights->is_compatible(*graph));

        if (graph->is_canonical_mode())
            config->forward_and_reverse = true;

        for (const auto &file : files) {
            parse_sequences(file, *config,
                [&graph,&node_weights](std::string_view seq) {
                    graph->map_to_nodes_sequentially(seq,
                        [&](auto node) { node_weights->add_weight(node, 1); }
                    );
                },
                [&graph,&node_weights](std::string_view seq, uint32_t count) {
                    graph->map_to_nodes_sequentially(seq,
                        [&](auto node) { node_weights->add_weight(node, count); }
                    );
                }
            );
            logger->trace("Extracted all sequences from file '{}' in {} sec",
                          file, timer.elapsed());
        }
    }

    logger->trace("Node weights updated in {} sec", timer.elapsed());

    assert(config->outfbase.size());

    // serialize graph
    timer.reset();

    graph->serialize(config->outfbase);
    graph->serialize_extensions(config->outfbase);
    graph.reset();

    logger->trace("Serialized in {} sec", timer.elapsed());

    timer.reset();

    if (!config->infbase_annotators.size())
        return 0;

    auto annotation = initialize_annotation(config->infbase_annotators.at(0), *config);

    if (!annotation->load(config->infbase_annotators.at(0))) {
        logger->error("Cannot load graph annotation from '{}'",
                      config->infbase_annotators.at(0));
        exit(1);
    } else {
        logger->trace("Annotation was loaded in {} sec", timer.elapsed());
    }

    timer.reset();

    assert(inserted_nodes);

    if (annotation->num_objects() + 1
            != inserted_nodes->size() - inserted_nodes->num_set_bits()) {
        logger->error("Graph and annotation are incompatible");
        exit(1);
    }

    logger->trace("Insert empty rows in annotation matrix...");

    // transform indexes of the inserved k-mers to the annotation format
    std::vector<uint64_t> inserted_rows;
    inserted_nodes->call_ones([&](auto i) {
        inserted_rows.push_back(graph::AnnotatedDBG::graph_to_anno_index(i));
    });
    annotation->insert_rows(inserted_rows);

    logger->trace("Rows inserted in {} sec", timer.elapsed());

    annotation->serialize(config->outfbase);

    return 0;
}

} // namespace cli
} // namespace mtg
