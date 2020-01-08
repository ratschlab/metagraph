#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <memory>

#include <Eigen/Sparse>

#include "annotation/representation/base/annotation.hpp"


namespace annotate {

// TODO: implement this as an annotation matrix
// StaticBinRelAnnotator<VectorRowBinMat>
template <typename Label = std::string>
class RowCompressed : public MultiLabelEncoded<Label> {
  public:
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;
    using SetBitPositions = typename MultiLabelEncoded<Label>::SetBitPositions;

    RowCompressed(uint64_t num_rows = 0, bool sparse = false);

    typedef std::function<void(Index, SetBitPositions&&)> CallRow;
    RowCompressed(uint64_t num_rows,
                  const std::vector<Label> &labels,
                  std::function<void(CallRow)> call_rows);

    void reinitialize(uint64_t num_rows);

    void set(Index i, const VLabels &labels);

    void add_labels(const std::vector<Index> &indices, const VLabels &labels);
    void add_labels_fast(const std::vector<Index> &indices, const VLabels &labels);

    void insert_rows(const std::vector<Index> &rows);

    bool has_label(Index i, const Label &label) const;
    bool has_labels(Index i, const VLabels &labels) const;

    SetBitPositions get_label_codes(Index i) const { return matrix_->get_row(i); }

    uint64_t num_objects() const;
    uint64_t num_relations() const;

    /************************* Serialization *************************/

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    static void load_shape(const std::string &filename,
                           uint64_t *num_objects,
                           uint64_t *num_relations);

    static std::unique_ptr<LabelEncoder<Label>>
    load_label_encoder(const std::string &filename);

    template <typename RowType = BinaryMatrix::SetBitPositions>
    class StreamRows {
      public:
        explicit StreamRows(std::string filename);

        // return null after all rows have been called
        RowType* next_row() { return sr_->next_row(); }

      private:
        std::unique_ptr<::StreamRows<RowType>> sr_;
    };

    // TODO: move to BinaryMatrixRowDynamic
    static void write_rows(std::string filename,
                           const LabelEncoder<Label> &label_encoder,
                           const std::function<void(BinaryMatrix::RowCallback)> &call_rows);

    /*****************************************************************/

    const BinaryMatrixRowDynamic& get_matrix() const { return *matrix_; };

    std::string file_extension() const { return kExtension; }

    static constexpr auto kExtension = ".row.annodbg";

  private:
    using MultiLabelEncoded<Label>::label_encoder_;

    std::unique_ptr<BinaryMatrixRowDynamic> matrix_;

    static std::unique_ptr<LabelEncoder<Label>>
    load_label_encoder(std::istream &instream);
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
