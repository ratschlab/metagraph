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

class RowDiff : public BinaryMatrix {
  public:
    RowDiff() {}

    RowDiff(const uint64_t num_columns,
            const graph::DBGSuccinct *graph,
            const Vector<uint64_t> &diffs,
            const sdsl::bit_vector &boundary,
            const sdsl::bit_vector &terminal)
        : num_columns_(num_columns),
          graph_(graph),
          diffs_(diffs),
          boundary_(boundary),
          terminal_(terminal) {
        sdsl::util::init_support(sboundary_, &boundary_);
    }

    RowDiff(const uint64_t num_columns,
            const graph::DBGSuccinct *graph,
            const sdsl::enc_vector<> &diffs,
            const sdsl::bit_vector &boundary,
            const sdsl::bit_vector &terminal)
            : num_columns_(num_columns),
              graph_(graph),
              diffs_(diffs),
              boundary_(boundary),
              terminal_(terminal) {
        sdsl::util::init_support(sboundary_, &boundary_);
    }

    uint64_t num_columns() const override { return num_columns_; }

    uint64_t num_rows() const override { return terminal_.size(); }

    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

    bool get(Row row, Column column) const override;

    /**
     * Returns the given column.
     * WARNING: this operation is slow. Prefer using other formats for column operations.
     */
    std::vector<Row> get_column(Column column) const override;

    /**
     * Returns the number of set bits in the matrix.
     * TODO: this is very slow, consider storing/loading the value
     */
    uint64_t num_relations() const override;

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

    const graph::DBGSuccinct *graph_;

    sdsl::enc_vector<> diffs_;
    sdsl::rrr_vector<> boundary_;
    sdsl::rrr_vector<> terminal_;
    sdsl::rrr_vector<>::select_1_type sboundary_;
};
} // namespace binmat
} // namespace annot
} // namespace mtg
