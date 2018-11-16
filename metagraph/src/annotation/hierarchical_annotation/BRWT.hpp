#ifndef __BRWT_HPP__
#define __BRWT_HPP__

#include <vector>
#include <unordered_map>
#include <memory>

#include "binary_matrix.hpp"
#include "bit_vector.hpp"
#include "utils.hpp"


// Compress sparse binary matrix (binary relations)
// using the BRWT structure described in:
// Barbay, Jérémy, Francisco Claude, and Gonzalo Navarro.
// "Compact binary relation representations with rich functionality."
// Information and Computation 232 (2013): 19-37.
class BRWT : public BinaryMatrix {
    friend class BRWTBuilder;
    friend class BRWTOptimizer;

    typedef uint32_t Child;

  public:
    BRWT() = default;

    BRWT(const BRWT &other) = default;
    BRWT& operator=(const BRWT &other) = default;

    BRWT(BRWT&& other) = default;
    BRWT& operator=(BRWT&& other) = default;

    uint64_t num_columns() const { return assignments_.size(); }
    uint64_t num_rows() const { return nonzero_rows_.size(); }

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

    // internal stats functions
    double avg_arity() const;
    uint64_t num_nodes() const;
    double shrinking_rate() const;
    uint64_t total_column_size() const;
    uint64_t total_num_set_bits() const;

  private:
    // breadth-first traversal
    void BFT(std::function<void(const BRWT &node)> callback) const;

    // assigns columns to the child nodes
    utils::RangePartition assignments_;
    bit_vector_rrr<> nonzero_rows_;
    std::vector<std::unique_ptr<BinaryMatrix>> child_nodes_;
};

#endif // __BRWT_HPP__
