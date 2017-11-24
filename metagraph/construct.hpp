#ifndef __CONSTRUCT_HPP__
#define __CONSTRUCT_HPP__

#include <cstdint>
#include <string>
#include <parallel/algorithm>

#include "dbg_succinct_libmaus.hpp"


namespace construct {

    // add a full sequence to the graph
    void add_seq(DBG_succ *G, kstring_t &seq, bool append = true);
    void add_seq_fast(DBG_succ *G, kstring_t &seq, const std::string &name, bool add_bridge = true,
                      unsigned int parallel = 1, std::string suffix = "", bool add_anno = false);
    void construct_succ(DBG_succ *G, unsigned int parallel = 1, bool add_anno = false);

    /** This function takes a character c and appends it to the end of the graph sequence
     * given that the corresponding note is not part of the graph yet.
     */
    uint64_t append_pos(DBG_succ *G, uint64_t c, uint64_t *ckmer = NULL, uint64_t i = 0);

    /** This function takes a pointer to a graph structure and concatenates the arrays W, last
     * and F to this graph's arrays. In almost all cases this will not produce a valid graph and
     * should only be used as a helper in the parallel merge procedure.
     */
    void append_graph(DBG_succ *G_t, DBG_succ *G_s);

    /**
     * This function takes a pointer to a graph structure and concatenates the arrays W, last
     * and F to this graph's static containers last_stat and W_stat. In almost all cases
     * this will not produce a valid graph and should only be used as a helper in the
     * parallel merge procedure.
     */
    void append_graph_static(DBG_succ *G_t, DBG_succ *G_s);

    uint64_t remove_edges(DBG_succ *G, std::set<uint64_t> &edges, uint64_t ref_point = 0);

    std::deque<std::string> generate_suffices(DBG_succ *G, unsigned int nsplits);

    void add_sink(DBG_succ* G, unsigned int parallel = 1, std::string suffix = "", bool add_anno = false);

    // translate from ascii into talphabet
    //TODO: move this to dbg_succinct_libmaus.hpp
    const char nt_lookup[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5
    };

} // namespace construct

#endif // __CONSTRUCT_HPP__
