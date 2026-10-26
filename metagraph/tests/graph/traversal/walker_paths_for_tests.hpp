#ifndef __TESTS_GRAPH_TRAVERSAL_WALKER_PATHS_FOR_TESTS_HPP__
#define __TESTS_GRAPH_TRAVERSAL_WALKER_PATHS_FOR_TESTS_HPP__

/**
 * A walked path of the segment DAG as the tests read it (PathResult stores only its leaf;
 * walk_path_leaf_first visits its first-parent chain).
 */

#include <algorithm>
#include <string>
#include <vector>

#include "graph/traversal/walker.hpp"


namespace mtg {
namespace graph {
namespace traversal {

// The flank of |path| in natural orientation, from the segments
inline std::string spell_path(const ArmResult &arm, const PathResult &path) {
    // walked leaf -> root through first parents: natural orientation is root -> leaf on
    // the right arm (so the pieces are reversed at the end) and leaf -> root on the left
    // (whose segment sequences are stored natural, i.e. reversed walking order)
    std::vector<const std::string*> pieces;
    walk_path_leaf_first(arm, path, [&](size_t s) {
        pieces.push_back(&arm.segments[s].sequence);
        return true;
    });
    if (arm.arm == Arm::RIGHT)
        std::reverse(pieces.begin(), pieces.end());
    std::string out;
    for (const std::string *p : pieces) {
        out += *p;
    }
    return out;
}

// The segments of |path|, root -> leaf, through first parents at joins: the chain that
// PathResult does not store. O(path depth in segments) per call.
inline std::vector<size_t> path_segments(const ArmResult &arm, const PathResult &path) {
    std::vector<size_t> chain;
    walk_path_leaf_first(arm, path, [&](size_t s) {
        chain.push_back(s);
        return true;
    });
    std::reverse(chain.begin(), chain.end());
    return chain;
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TESTS_GRAPH_TRAVERSAL_WALKER_PATHS_FOR_TESTS_HPP__
