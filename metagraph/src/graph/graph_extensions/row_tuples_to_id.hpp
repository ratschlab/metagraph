#ifndef __NODE_WEIGHTS_HPP__
#define __NODE_WEIGHTS_HPP__

#include <string>
#include <vector>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg::graph {

// Rename to CoordToAccession?
class RowTuplesToId : public SequenceGraph::GraphExtension {
  public:
    using Column = mtg::annot::matrix::BinaryMatrix::Column;
    using Tuple = mtg::annot::matrix::MultiIntMatrix::Tuple;
    using RowTuples = mtg::annot::matrix::MultiIntMatrix::RowTuples;

    RowTuplesToId() {}
    // Constructor from vector of (sequence_header, num_kmers) pairs per annotation column.
    // The column names are stored in col_names.
    RowTuplesToId(const std::vector<std::vector<std::pair<std::string, uint64_t>>> &accessions,
                  const std::vector<std::string> &col_names);
    // Merge multiple serialized RowTuplesToId objects into one and remove them once merged.
    RowTuplesToId(const std::vector<std::string> &fnames);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>>
    rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuple) const;

    uint64_t num_columns() const { return seq_delims_.size(); }

    static constexpr auto kRowTuplesExtension = ".seqs";

  private:
    std::vector<std::vector<std::string>> seq_id_labels_;
    std::vector<bit_vector_sd> seq_delims_;
};

} // namespace mtg::graph

#endif // __NODE_WEIGHTS_HPP__
