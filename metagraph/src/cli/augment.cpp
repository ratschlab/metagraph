#include "augment.hpp"

#include <string>
#include <vector>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "graph/annotated_dbg.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotation.hpp"
#include "config/config.hpp"
#include "sequence_reader.hpp"

using mg::common::logger;
using utils::get_verbose;


int augment_graph(Config *config) {
    assert(config);
    assert(config->infbase_annotators.size() <= 1);

    const auto &files = config->fnames;

    Timer timer;

    // load graph
    auto graph = load_critical_dbg(config->infbase);

    auto node_weights = graph->load_extension<NodeWeights>(config->infbase);
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

    if (dynamic_cast<DBGSuccinct*>(graph.get())) {
        auto &succinct_graph = dynamic_cast<DBGSuccinct&>(*graph);

        if (succinct_graph.get_state() != BOSS::State::DYN) {
            logger->trace("Switching state of succinct graph to dynamic...");

            succinct_graph.switch_state(BOSS::State::DYN);

            logger->trace("State switching done in {} sec", timer.elapsed());
        }
    }

    std::unique_ptr<bit_vector_dyn> inserted_edges;
    if (config->infbase_annotators.size() || node_weights)
        inserted_edges.reset(new bit_vector_dyn(graph->max_index() + 1, 0));

    timer.reset();

    logger->trace("Start graph augmentation");

    if (graph->is_canonical_mode())
        config->forward_and_reverse = false;

    config->canonical = graph->is_canonical_mode();

    parse_sequences(files, *config, timer,
        [&graph,&inserted_edges](std::string&& seq) {
            graph->add_sequence(seq, inserted_edges.get());
        },
        [&graph,&inserted_edges](std::string&& kmer, uint32_t /*count*/) {
            graph->add_sequence(kmer, inserted_edges.get());
        },
        [&graph,&inserted_edges](const auto &loop) {
            loop([&graph,&inserted_edges](const char *seq) {
                graph->add_sequence(seq, inserted_edges.get());
            });
        }
    );

    logger->trace("Graph augmentation done in {} sec", timer.elapsed());
    timer.reset();

    if (node_weights) {
        node_weights->insert_nodes(*inserted_edges);

        assert(node_weights->is_compatible(*graph));

        if (graph->is_canonical_mode())
            config->forward_and_reverse = true;

        parse_sequences(files, *config, timer,
            [&graph,&node_weights](std::string&& seq) {
                graph->map_to_nodes_sequentially(seq.begin(), seq.end(),
                    [&](auto node) { node_weights->add_weight(node, 1); }
                );
            },
            [&graph,&node_weights](std::string&& kmer, uint32_t count) {
                node_weights->add_weight(graph->kmer_to_node(kmer), count);
            },
            [&graph,&node_weights](const auto &loop) {
                loop([&graph,&node_weights](const char *seq) {
                    std::string seq_str(seq);
                    graph->map_to_nodes_sequentially(seq_str.begin(), seq_str.end(),
                        [&](auto node) { node_weights->add_weight(node, 1); }
                    );
                });
            }
        );
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

    assert(inserted_edges);

    if (annotation->num_objects() + 1
            != inserted_edges->size() - inserted_edges->num_set_bits()) {
        logger->error("Graph and annotation are incompatible");
        exit(1);
    }

    logger->trace("Insert empty rows in annotation matrix...");

    // transform indexes of the inserved k-mers to the annotation format
    std::vector<uint64_t> inserted_rows;
    inserted_edges->call_ones([&](auto i) {
        inserted_rows.push_back(AnnotatedDBG::graph_to_anno_index(i));
    });
    annotation->insert_rows(inserted_rows);

    logger->trace("Rows inserted in {} sec", timer.elapsed());

    annotation->serialize(config->outfbase);

    return 0;
}
