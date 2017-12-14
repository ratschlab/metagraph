#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>


class DBG_succ;

namespace merge {

    /*
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    DBG_succ* merge(const std::vector<DBG_succ*> &Gv,
                    std::vector<uint64_t> kv,
                    std::vector<uint64_t> nv);

    /**
    * Heavily borrowing from the graph sequence traversal, this function gets a graph pointer |mergeable| and merges its
    * nodes into the target graph object |target|. The edges of |mergeable| are fully traversed and nodes are added to
    * G_t if not existing yet. This function is well suited to merge small graphs into large ones.
    */
    void merge(DBG_succ *target, const DBG_succ &mergeable);

}

#endif // __MERGE_HPP__
