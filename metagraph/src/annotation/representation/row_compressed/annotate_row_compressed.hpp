#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <memory>

#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace annot {

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

    static matrix::StreamRows<matrix::BinaryMatrix::SetBitPositions>
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

} // namespace annot
} // namespace mtg

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
