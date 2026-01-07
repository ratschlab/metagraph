#ifndef __COORD_TO_HEADER_HPP__
#define __COORD_TO_HEADER_HPP__

#include <string>
#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vectors/bit_vector_sd.hpp"


namespace mtg {
namespace annot {

/**
 * An index mapping k-mer coordinates in annotation columns to sequence headers.
 * This mapping can be used to transform coordinates and labels in query results to show
 * sequence-based (e.g., "<seq1>:0-0-3\t<seq3>:0-1-5") instead of file-based coordinates
 * (e.g., "<file.fa>:0-0-3:0-23-27"). This allows showing results of querying coordinate
 * annotations as if sequences were annotated by their headers instead of source files.
 */
class CoordToHeader {
  public:
    using Column = mtg::annot::matrix::BinaryMatrix::Column;
    using Tuple = mtg::annot::matrix::MultiIntMatrix::Tuple;
    using RowTuples = mtg::annot::matrix::MultiIntMatrix::RowTuples;

    CoordToHeader() {}
    // Constructor from sequence names and number of k-mers in them, per annotation column.
    CoordToHeader(std::vector<std::vector<std::string>> &&headers,
                  std::vector<std::vector<uint64_t>> &&num_kmers);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    // Transforms global coords to sequence-based coords and outputs them as a list of tuples
    // (header, num_matches, coords). `coords` contains all coordinates from `rows_tuple`
    // corresponding to that sequence/header and mapped to the sequence-based positions and
    // `num_matches` is the number of rows with matches for that sequence/header (the number of
    // non-zero tuples in `coords`). Also applies filtering (`num_matches` must be >= `min_count`).
    std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>>
    rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuple, size_t min_count) const;

    uint64_t num_columns() const { return seq_delims_.size(); }
    uint64_t num_sequences(Column column) const { return seq_id_labels_[column].size(); }
    // Get number of k-mers/coordinates in a specific column
    uint64_t num_kmers(Column column) const { return seq_delims_[column].size(); }
    const std::vector<std::string>& get_headers(Column column) const {
        return seq_id_labels_[column];
    }

    static constexpr auto kExtension = ".seqs";

  private:
    std::vector<std::vector<std::string>> seq_id_labels_;
    std::vector<bit_vector_sd> seq_delims_;
};

} // namespace annot
} // namespace mtg

#endif // __COORD_TO_HEADER_HPP__
