#ifndef __BOSS_MERGE_HPP__
#define __BOSS_MERGE_HPP__

#include <cstdint>
#include <vector>
#include <string>

#include "boss.hpp"
#include "boss_chunk.hpp"


namespace mtg {
namespace graph {
namespace boss {

    /**
     * Given a list of boss tables, this function
     * merges all of them into a new one.
     */
    BOSS* merge(const std::vector<const BOSS*> &graphs,
                bool verbose = false);

    BOSS::Chunk merge_blocks_to_chunk(const std::vector<const BOSS*> &graphs,
                                      size_t chunk_idx,
                                      size_t num_chunks,
                                      size_t num_threads,
                                      size_t num_bins_per_thread,
                                      bool verbose = false);

} // namespace boss
} // namespace graph
} // namespace mtg

#endif // __BOSS_MERGE_HPP__
