#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <memory>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace annot {

template <typename RowType = matrix::BinaryMatrix::SetBitPositions>
class StreamRows;

// TODO: implement this as an annotation matrix
// StaticBinRelAnnotator<VectorRowBinMat>
// TODO: remove?
template <typename Label = std::string>
class RowCompressed : public MultiLabelAnnotation<Label> {
  public:
    using Index = typename MultiLabelAnnotation<Label>::Index;
    using VLabels = typename MultiLabelAnnotation<Label>::VLabels;

    RowCompressed(uint64_t num_rows = 0, bool sparse = false);

    template <typename RowType>
    RowCompressed(Vector<RowType>&& annotation_rows,
                  const std::vector<Label> &column_labels);

    void reinitialize(uint64_t num_rows);

    void add_labels(const std::vector<Index> &indices, const VLabels &labels);
    void add_labels_fast(const std::vector<Index> &indices, const VLabels &labels);

    void insert_rows(const std::vector<Index> &rows);

    uint64_t num_objects() const;
    uint64_t num_relations() const;

    /************************* Serialization *************************/

    void serialize(const std::string &filename) const;
    bool load(const std::string &filename) { return merge_load({ filename }); }
    bool merge_load(const std::vector<std::string> &filenames);

    static void serialize(const std::string &filename,
                          const LabelEncoder<Label> &label_encoder,
                          const std::function<void(matrix::BinaryMatrix::RowCallback)> &call_rows);

    static void read_shape(const std::string &filename,
                           uint64_t *num_objects,
                           uint64_t *num_relations);

    static LabelEncoder<Label> read_label_encoder(const std::string &filename);

    static StreamRows<matrix::BinaryMatrix::SetBitPositions>
    get_row_streamer(const std::string &filename);

    /*****************************************************************/

    const matrix::BinaryMatrixRowDynamic& get_matrix() const { return *matrix_; };

    std::string file_extension() const { return kExtension; }

    static constexpr auto kExtension = ".row.annodbg";

  private:
    using MultiLabelAnnotation<Label>::label_encoder_;

    std::unique_ptr<matrix::BinaryMatrixRowDynamic> matrix_;

    static LabelEncoder<Label> read_label_encoder(std::istream &instream);
};


// Row streamer -- read rows from a serialized row major binary matrix
template <typename RowType>
class StreamRows {
  public:
    StreamRows(const std::string &filename, size_t offset);

    //TODO: implement constructor from stream once
    //      it's implemented for sdsl::int_vector_buffer<>.
    //      Then, use StreamRows to simplify load functions.
    // StreamRows(std::istream &instream);

    // return nullptr after all rows have been called
    RowType* next_row();

  private:
    RowType row_;
    sdsl::int_vector_buffer<> inbuf_;
    uint64_t i_ = 0;
};

// Write matrix to the end
void append_row_major(const std::string &filename,
                      const std::function<void(matrix::BinaryMatrix::RowCallback)> &call_rows,
                      uint64_t num_cols);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
