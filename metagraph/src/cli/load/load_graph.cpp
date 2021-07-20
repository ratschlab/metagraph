#include "load_graph.hpp"

#include <string>
#include <cassert>

#include "common/logger.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/dbg_succinct_cached.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;

using mtg::common::logger;


Config::GraphType parse_graph_type(const std::string &filename) {
    if (utils::ends_with(filename, DBGSuccinct::kExtension)) {
        return Config::GraphType::SUCCINCT;

    } else if (utils::ends_with(filename, DBGHashOrdered::kExtension)) {
        return Config::GraphType::HASH;

    } else if (utils::ends_with(filename, DBGHashString::kExtension)) {
        return Config::GraphType::HASH_STR;

    } else if (utils::ends_with(filename, DBGHashFast::kExtension)) {
        return Config::GraphType::HASH_FAST;

    } else if (utils::ends_with(filename, graph::DBGBitmap::kExtension)) {
        return Config::GraphType::BITMAP;

    } else {
        return Config::GraphType::INVALID;
    }
}

std::shared_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename) {
    auto graph_type = parse_graph_type(filename);
    switch (graph_type) {
        case Config::GraphType::SUCCINCT:
            return load_critical_graph_from_file<DBGSuccinct>(filename);

        case Config::GraphType::HASH:
            return load_critical_graph_from_file<DBGHashOrdered>(filename);

        case Config::GraphType::HASH_PACKED:
            return load_critical_graph_from_file<DBGHashOrdered>(filename);

        case Config::GraphType::HASH_STR:
            return load_critical_graph_from_file<DBGHashString>(filename);

        case Config::GraphType::HASH_FAST:
            return load_critical_graph_from_file<DBGHashFast>(filename);

        case Config::GraphType::BITMAP:
            return load_critical_graph_from_file<graph::DBGBitmap>(filename);

        case Config::GraphType::INVALID:
            logger->error("Cannot load graph from file '{}', needs a valid file extension",
                          filename);
            exit(1);
    }
    assert(false);
    exit(1);
}

std::shared_ptr<DeBruijnGraph> primary_to_canonical(std::shared_ptr<DeBruijnGraph> graph,
                                                    size_t cache_size) {
    if (graph->get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only primary mode graphs can be wrapped into canonical mode.");
        exit(1);
    }

    logger->trace("Primary graph wrapped into canonical");
    return std::make_shared<CanonicalDBG>(graph, cache_size);
}

std::shared_ptr<DeBruijnGraph> make_cached_graph(std::shared_ptr<DeBruijnGraph> graph,
                                                 const Config &config,
                                                 size_t cache_size) {
    // if alignment in both directions is not required, then there's no need to cache
    if (graph->get_mode() != DeBruijnGraph::CANONICAL && config.align_one_strand)
        return graph;

    auto base_graph = graph;
    auto canonical = std::dynamic_pointer_cast<CanonicalDBG>(graph);

    if (canonical) {
        base_graph = canonical->get_graph_ptr();
        assert(base_graph && "CanonicalDBG should be built from a non-const graph");
    }

    if (auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(base_graph)) {
        graph = make_cached_dbgsuccinct(base_graph, cache_size);
    } else {
        // graphs other than DBGSuccinct can't be cached
        return graph;
    }

    if (graph->get_mode() == DeBruijnGraph::PRIMARY)
        graph = primary_to_canonical(graph, cache_size);

    return graph;
}

} // namespace cli
} // namespace mtg
