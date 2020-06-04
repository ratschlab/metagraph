#include "load_graph.hpp"

#include <string>
#include <cassert>

#include "common/logger.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

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

    } else if (utils::ends_with(filename, bitmap_graph::DBGBitmap::kExtension)) {
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
            return load_critical_graph_from_file<bitmap_graph::DBGBitmap>(filename);

        case Config::GraphType::INVALID:
            logger->error("Cannot load graph from file '{}', needs a valid file extension",
                          filename);
            exit(1);
    }
    assert(false);
    exit(1);
}

} // namespace cli
} // namespace mtg
