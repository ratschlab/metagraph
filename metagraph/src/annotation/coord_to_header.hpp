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

    /**
     * Transforms global coordinates to local (sequence-based) coordinates in-place.
     *
     * For each coordinate in `rows_tuples`, this function:
     *   1. Determines which sequence (header) the coordinate belongs to within its column
     *   2. Converts the global coordinate to a local coordinate within that sequence
     *   3. Encodes both the local coordinate and header index into a single value
     *
     * Encoding scheme: For a coordinate in column `j`, the transformed value is:
     *   `encoded_coord = local_coord * num_headers[j] + header_id`
     *
     * Decoding: To extract the header index and local coordinate from an encoded value:
     *   `header_id = encoded_coord % num_headers[j]`
     *   `local_coord = encoded_coord / num_headers[j]`
     */
    void map_to_local_coords(std::vector<RowTuples> *rows_tuples) const;

    uint64_t num_columns() const { return coord_offsets_.size(); }
    uint64_t num_sequences(Column column) const { return headers_[column].size(); }
    // Get number of k-mers/coordinates in a specific column
    uint64_t num_kmers(Column column) const { return coord_offsets_[column].size(); }
    const std::vector<std::string>& get_headers(Column column) const { return headers_[column]; }
    size_t num_headers(Column column) const { return headers_[column].size(); }

    static constexpr auto kExtension = ".seqs";

  private:
    std::vector<std::vector<std::string>> headers_;
    std::vector<bit_vector_sd> coord_offsets_;
};

} // namespace annot
} // namespace mtg

#endif // __COORD_TO_HEADER_HPP__
