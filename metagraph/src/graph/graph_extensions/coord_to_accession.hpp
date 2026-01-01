#ifndef __COORD_TO_ACCESSION_HPP__
#define __COORD_TO_ACCESSION_HPP__

#include <string>
#include <vector>

#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vectors/bit_vector_sd.hpp"


namespace mtg::graph {

// A mapping from k-mer coordinates in columns to sequence accessions (headers).
// Maps k-mer coordinates within annotation columns to sequence accessions (sequence headers).
// When querying with flag `--accessions`, this mapping is used to transform coordinates and labels
// in query results to show sequence-based (e.g., "<seq1>:0-0-3\t<seq3>:0-1-5") instead of file-
// based coordinates (e.g., "<file.fa>:0-0-3:0-23-27"). This allows querying coordinate annotations
// but displaying results as if sequences were annotated by their headers instead of source files.
class CoordToAccession : public SequenceGraph::GraphExtension {
  public:
    using Column = mtg::annot::matrix::BinaryMatrix::Column;
    using Tuple = mtg::annot::matrix::MultiIntMatrix::Tuple;
    using RowTuples = mtg::annot::matrix::MultiIntMatrix::RowTuples;

    CoordToAccession() {}
    // Constructor from vector of (sequence_header, num_kmers) pairs per annotation column.
    // The column names are stored in col_names.
    CoordToAccession(const std::vector<std::vector<std::pair<std::string, uint64_t>>> &accessions,
                     const std::vector<std::string> &col_names);
    // Merge multiple serialized CoordToAccession objects into one and remove them once merged.
    CoordToAccession(const std::vector<std::string> &fnames, size_t num_threads);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>>
    rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuple) const;

    uint64_t num_columns() const { return seq_delims_.size(); }

    static constexpr auto kExtension = ".seqs";

  private:
    std::vector<std::vector<std::string>> seq_id_labels_;
    std::vector<bit_vector_sd> seq_delims_;
};

} // namespace mtg::graph

#endif // __COORD_TO_ACCESSION_HPP__
