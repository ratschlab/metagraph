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

template <class Graph>
std::shared_ptr<CanonicalDBG> primary_to_canonical(std::shared_ptr<Graph> graph,
                                                   size_t cache_size) {
    if (graph->get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only primary mode graphs can be wrapped into canonical mode.");
        exit(1);
    }

    logger->trace("Primary graph wrapped into canonical");
    return std::make_shared<CanonicalDBG>(graph, cache_size);
}

template std::shared_ptr<CanonicalDBG> primary_to_canonical(std::shared_ptr<const DeBruijnGraph>, size_t);
template std::shared_ptr<CanonicalDBG> primary_to_canonical(std::shared_ptr<DeBruijnGraph>, size_t);

std::shared_ptr<const DeBruijnGraph> make_cached_graph(const DeBruijnGraph &graph,
                                                       const Config &config,
                                                       size_t cache_size) {
    std::shared_ptr<const DeBruijnGraph> ret_graph(std::shared_ptr<const DeBruijnGraph>{}, &graph);

    // if alignment in both directions is not required, then there's no need to cache
    if (ret_graph->get_mode() != DeBruijnGraph::CANONICAL && config.align_only_forwards)
        return ret_graph;

    auto base_graph = ret_graph;

    if (auto canonical = std::dynamic_pointer_cast<const CanonicalDBG>(ret_graph))
        base_graph = canonical->get_graph_ptr();

    if (auto dbg_succ = std::dynamic_pointer_cast<const DBGSuccinct>(base_graph)) {
        ret_graph = dbg_succ->get_cached_view(cache_size);
    } else {
        // graphs other than DBGSuccinct can't be cached
        return ret_graph;
    }

    if (ret_graph->get_mode() == DeBruijnGraph::PRIMARY)
        ret_graph = primary_to_canonical(ret_graph, cache_size);

    return ret_graph;
}

} // namespace cli
} // namespace mtg
