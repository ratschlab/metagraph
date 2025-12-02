#ifndef __METAGRAPH_HPP__
#define __METAGRAPH_HPP__

/**
 * When using MetaGraph via ExternalProject_Add or as a compiled library,
 * this header should be included FIRST in your source files. It provides essential
 * compile-time configuration definitions that affect how other MetaGraph headers
 * are interpreted (e.g., _DNA_GRAPH, _PROTEIN_GRAPH, _USE_FOLLY, etc.).
 */

// Compile-time configuration definitions (must be included first)
#include "common/compile_definitions.hpp"

// Core graph and annotation structures
#include "graph/annotated_dbg.hpp"

// Command-line configuration and runtime parameters
#include "cli/config/config.hpp"

// Loading functionality
#include "cli/load/load_graph.hpp"
#include "cli/load/load_annotated_graph.hpp"

// Query and alignment operations
#include "cli/query.hpp"
#include "cli/align.hpp"

#endif // __METAGRAPH_HPP__
