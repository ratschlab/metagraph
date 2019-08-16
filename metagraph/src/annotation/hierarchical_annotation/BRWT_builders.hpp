#ifndef __BRWT_BUILDERS_HPP__
#define __BRWT_BUILDERS_HPP__

#include <vector>
#include <unordered_map>
#include <memory>

#include "BRWT.hpp"
#include "bit_vector.hpp"


class BRWTBuilder {
    friend class BRWTOptimizer;

  public:
    virtual ~BRWTBuilder() {}

    using Row = typename BRWT::Row;
    using Column = typename BRWT::Column;
    using Child = typename BRWT::Child;

  protected:
    struct NodeBRWT {
        // indexes of columns in the initial matrix
        // compressed in subtree with the root in this node
        std::vector<Column> column_arrangement;
        // sizes of groups in affine partition of |column_arrangement|
        std::vector<size_t> group_sizes;
        // compressed submatrices of the current matrix without zero rows
        // with columns distributed accoring to partition
        std::vector<std::unique_ptr<BinaryMatrix>> child_nodes;
    };

    static BRWT initialize(NodeBRWT&& node, bit_vector&& nonzero_rows);
};


class BRWTBottomUpBuilder : public BRWTBuilder {
  public:
    typedef std::vector<std::unique_ptr<bit_vector>> VectorsPtr;
    typedef std::vector<Column> Group;
    typedef std::vector<Group> Partition;
    typedef std::function<Partition(const VectorsPtr &)> Partitioner;

    static Partitioner get_basic_partitioner(size_t arity);

    static BRWT build(VectorsPtr&& columns,
                      Partitioner partitioner = get_basic_partitioner(2),
                      size_t num_threads = 1);

  private:
    static std::pair<NodeBRWT, std::unique_ptr<bit_vector>>
    merge(std::vector<NodeBRWT> &&nodes,
          std::vector<std::unique_ptr<bit_vector>> &&nonzero_rows);
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
    using NodeBRWT = BRWTBuilder::NodeBRWT;

    static void add_submatrix(std::unique_ptr<BinaryMatrix>&& node,
                              NodeBRWT *parent,
                              uint64_t max_delta_arity,
                              size_t num_threads);

    // removes a node and reassigns all its children to its parent
    static void reassign(std::unique_ptr<BRWT>&& node,
                         NodeBRWT *parent,
                         size_t num_threads);
    // estimate delta between the transformed tree and the current one
    static double pruning_delta(const BRWT &node);
};


#endif // __BRWT_BUILDERS_HPP__
