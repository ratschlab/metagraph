#ifndef __CONSTRUCT_HPP__
#define __CONSTRUCT_HPP__

#include <cstdint>
#include <string>
#include "kseq.h"

#include "dbg_succinct_libmaus.hpp"

namespace construct {
    
    typedef uint64_t TAlphabet;

    // add a full sequence to the graph
    void add_seq (DBG_succ* G ,kstring_t &seq);
    void add_seq_alt (DBG_succ* G, kstring_t &seq, bool bridge=true, unsigned int parallel=1, std::string suffix="");
    void construct_succ(DBG_succ* G, unsigned int parallel=1);

    /** This function takes a character c and appends it to the end of the graph sequence
     * given that the corresponding note is not part of the graph yet.
     */
    void append_pos(DBG_succ* G, TAlphabet c);

    /** This function takes a pointer to a graph structure and concatenates the arrays W, last 
     * and F to this graph's arrays. In almost all cases this will not produce a valid graph and 
     * should only be used as a helper in the parallel merge procedure.
     */
    void append_graph(DBG_succ *G_t, DBG_succ* G_s);

    /** 
     * This function takes a pointer to a graph structure and concatenates the arrays W, last 
     * and F to this graph's static containers last_stat and W_stat. In almost all cases 
     * this will not produce a valid graph and should only be used as a helper in the 
     * parallel merge procedure.
     */
    void append_graph_static(DBG_succ *G_t, DBG_succ* G_s);

}
#endif
