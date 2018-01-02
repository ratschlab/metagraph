#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>
#include <string>

#include "dbg_succinct.hpp"


namespace merge {

    /**
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    DBG_succ* merge(const std::vector<const DBG_succ*> &graphs);


    class graph_chunk {
      public:
        virtual ~graph_chunk() {}

        virtual void push_back(TAlphabet W, TAlphabet F, bool last) = 0;

        virtual TAlphabet get_W_back() const = 0;

        virtual void alter_W_back(TAlphabet W) = 0;

        virtual void alter_last_back(bool last) = 0;

        virtual uint64_t size() const = 0;

        virtual void extend(const graph_chunk &other) = 0;

        virtual void initialize_graph(DBG_succ *graph) = 0;

        virtual bool load(const std::string &filename_base) = 0;

        virtual void serialize(const std::string &filename_base) const = 0;
    };

    class dynamic_graph_chunk;

    class vector_graph_chunk;

    graph_chunk* build_chunk(const std::vector<const DBG_succ*> &graphs,
                             size_t chunk_idx,
                             size_t num_chunks,
                             size_t num_threads,
                             size_t num_bins_per_thread);

    /**
     * Merge graph chunks from the vector passed and release the chunks afterwards.
     * If the null pointers are passed,
     * load chunks from files "filenamebase.<chunk_idx>_<num_chunks>".
     */
    DBG_succ* merge_chunks(size_t k,
                           const std::vector<graph_chunk*> &graph_chunks,
                           const std::string &filenamebase = "");

} // namespace merge

#endif // __MERGE_HPP__
