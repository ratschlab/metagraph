#ifndef __BRWT_BUILDERS_HPP__
#define __BRWT_BUILDERS_HPP__

#include <vector>
#include <unordered_map>
#include <memory>
#include <filesystem>

#include "common/vectors/bit_vector.hpp"
#include "common/threads/threading.hpp"
#include "brwt.hpp"


namespace mtg {
namespace annot {
namespace binmat {

class BRWTBottomUpBuilder {
  public:
    typedef std::vector<const bit_vector *> VectorPtrs;
    typedef std::vector<std::vector<BRWT::Column>> Partition;
    typedef std::function<Partition(const VectorPtrs &)> Partitioner;

    static Partitioner get_basic_partitioner(size_t arity = 2);

    // Build the Multi-BRWT compressed representation of a binary matrix
    static BRWT build(std::vector<std::unique_ptr<bit_vector>>&& columns,
                      Partitioner partitioner = get_basic_partitioner(),
                      size_t num_nodes_parallel = 1,
                      size_t num_threads = 1);

    using CallColumn
        = std::function<void(uint64_t, std::unique_ptr<bit_vector>&&)>;

    static BRWT build(const std::function<void(const CallColumn &)> &get_columns,
                      const std::vector<std::vector<uint64_t>> &linkage,
                      const std::filesystem::path &tmp_dir,
                      size_t num_nodes_parallel = 1,
                      size_t num_threads = 1);

    // Merge multiple binary matrices compressed with Multi-BRWT
    static BRWT merge(std::vector<BRWT>&& matrices,
                      Partitioner partitioner = get_basic_partitioner(),
                      size_t num_nodes_parallel = 1,
                      size_t num_threads = 1);

  private:
    // Concatenate multiple Multi-BRWT submatrices
    static BRWT concatenate(std::vector<BRWT>&& submatrices,
                            sdsl::bit_vector *buffer,
                            ThreadPool &thread_pool);
    // Concatenate multiple Multi-BRWT submatrices
    static BRWT concatenate_sparse(std::vector<BRWT>&& submatrices,
                                   sdsl::bit_vector *buffer,
                                   ThreadPool &thread_pool);
};


class BRWTOptimizer {
  public:
    virtual ~BRWTOptimizer() {}

    // remove some internal nodes to make the tree
    // smaller and increase the arity
    static void relax(BRWT *brwt_matrix,
                      uint64_t max_arity = -1,
                      size_t num_threads = 1);
  private:
    // check if removing this node is going to reduce the size
    static bool should_prune(const BRWT &node);
    // remove the node and reassign all its children to its parent
    static void reassign(size_t node_rank, BRWT *parent, size_t num_threads);
    // estimate delta between the transformed tree and the current one
    static double pruning_delta(const BRWT &node);
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __BRWT_BUILDERS_HPP__
