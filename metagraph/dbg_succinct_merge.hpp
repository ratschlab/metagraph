#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>
#include <string>


class DBG_succ;
class Config;

namespace merge {

    DBG_succ* build_chunk(const std::vector<const DBG_succ*> &Gv, Config *config);

    DBG_succ* merge_chunks(const std::string &filenamebase, size_t num_chunks);

    /*
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    DBG_succ* merge(const std::vector<const DBG_succ*> &Gv,
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
