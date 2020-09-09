#pragma once

#include <fstream>

#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <sdsl/rrr_vector.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vector.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg {
namespace annot {
namespace binmat {

/**
 * Row-major representation of annotations where rows are stored as differences vs a
 * predecessor row. The predecessor of a row is determined by a path in an external graph
 * structure, a #graph::DBGSuccinct, which has the property that rows on a path are likely
 * very similar.
 * The row-diff binary matrix is defined by three data structures:
 *   1. #diffs_ stores the concatenated diffs row by row. Since each row-diff has variable
 *      length, the boundary_ vector marks the end of each row-diff.
 *   2. boundary_ marks the end of a row diff in #diffs_
 *   3. terminal_ Rows marked as terminal are stored in full.
 * Retrieving data from RowDiff requires the associated #graph_. In order to get the i-th
 * row, we start traversing the node corresponding to i in #graph_ and accumulate the
 * diffs until we hit a terminal node, which is stored in full.
 */
class RowDiff : public BinaryMatrix {
  public:
    RowDiff() {}

    RowDiff(const uint64_t num_columns,
            const uint64_t num_relations,
            const graph::DBGSuccinct *graph,
            const sdsl::enc_vector<> &diffs,
            const sdsl::bit_vector &boundary,
            const sdsl::bit_vector &terminal)
            : num_columns_(num_columns),
              num_relations_(num_relations),
              graph_(graph),
              diffs_(diffs),
              boundary_(boundary),
              terminal_(terminal) {
        sdsl::util::init_support(sboundary_, &boundary_);
    }

    uint64_t num_columns() const override { return num_columns_; }

    /**
     * Returns the number of set bits in the matrix.
     */
    uint64_t num_relations() const override { return num_relations_; }

    uint64_t num_rows() const override { return terminal_.size(); }

    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

    bool get(Row row, Column column) const override;

    /**
     * Returns the given column.
     * WARNING: this operation is slow. Prefer using other formats for column operations.
     */
    std::vector<Row> get_column(Column column) const override;


    SetBitPositions get_row(Row row) const override;

    bool load(std::istream &f) override;
    void serialize(std::ostream &f) const override;

    void serialize(const std::string &name) const;
    bool load(const std::string &name);

    const sdsl::rrr_vector<> &terminal() const { return terminal_; }
    const sdsl::rrr_vector<> &boundary() const { return boundary_; }
    const sdsl::enc_vector<> &diffs() const { return diffs_; }

    Vector<uint64_t> get_diff(uint64_t node_id) const;

  private:
    static void merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2);

    uint64_t num_columns_ = 0;

    uint64_t num_relations_ = 0;

    const graph::DBGSuccinct *graph_;

    sdsl::enc_vector<> diffs_;
    sdsl::rrr_vector<> boundary_;
    sdsl::rrr_vector<> terminal_;
    sdsl::rrr_vector<>::select_1_type sboundary_;
};
} // namespace binmat
} // namespace annot
} // namespace mtg
