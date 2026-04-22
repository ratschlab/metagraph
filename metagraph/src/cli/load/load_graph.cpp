#include "load_graph.hpp"

#include <string>
#include <vector>

#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using mtg::common::logger;

using GraphLoader = std::shared_ptr<DeBruijnGraph> (*)(const std::string&);

template <class Graph>
std::shared_ptr<DeBruijnGraph> load_as(const std::string &filename) {
    return load_critical_graph_from_file<Graph>(filename);
}

#define GRAPH_LOADER_ENTRY(GRAPH) { GRAPH::kExtension, &load_as<GRAPH> }
static const std::vector<std::pair<std::string, GraphLoader>> kGraphLoaders = {
    GRAPH_LOADER_ENTRY(DBGSuccinct),
    GRAPH_LOADER_ENTRY(DBGHashOrdered),
    GRAPH_LOADER_ENTRY(DBGHashString),
    GRAPH_LOADER_ENTRY(DBGHashFast),
    GRAPH_LOADER_ENTRY(DBGBitmap),
    GRAPH_LOADER_ENTRY(DBGSSHash),
};

std::shared_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename) {
    for (const auto &[extension, load] : kGraphLoaders) {
        if (utils::ends_with(filename, extension))
            return load(filename);
    }
    auto exts = utils::get_firsts<std::vector<std::string>>(kGraphLoaders);
    logger->error("Cannot load graph from '{}': unrecognized file extension (expected one of {})",
                  filename, fmt::join(exts, "/"));
    exit(1);
}

} // namespace cli
} // namespace mtg
