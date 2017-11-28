#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>

#include "dbg_succinct.hpp"


namespace merge {

    typedef uint64_t TAlphabet;

    /*
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    void merge(DBG_succ *Gt,
               std::vector<DBG_succ*> Gv,
               std::vector<uint64_t> kv,
               std::vector<uint64_t> nv);

    /**
    * Heavily borrowing from the graph sequence traversal, this function gets a graph pointer G_m and merges its
    * nodes into the target graph object G_t. The edges of G_m are fully traversed and nodes are added to
    * G_t if not existing yet. This function is well suited to merge small graphs into large ones.
    */
    void merge(DBG_succ *G_t, DBG_succ *G_m);

    /*
     * Helper function to determine the bin boundaries, given
     * a number of threads.
     */
    std::vector<std::pair<uint64_t, uint64_t>> get_bins(
        DBG_succ* G, uint64_t bins
    );

    std::vector<std::pair<uint64_t, uint64_t>> get_bins_relative(
        DBG_succ* G_from,
        DBG_succ* G_to,
        std::vector<std::pair<uint64_t, uint64_t> > ref_bins,
        uint64_t first_pos, uint64_t last_pos
    );

}

#endif // __MERGE_HPP__
