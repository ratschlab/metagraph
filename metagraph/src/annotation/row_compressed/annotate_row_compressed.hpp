#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <memory>

#include <Eigen/Sparse>

#include "annotate.hpp"
#include "binary_matrix.hpp"
#include "vector_row_binmat.hpp"


namespace annotate {

const char kRowAnnotatorExtension[] = ".row.annodbg";

template <typename Label>
class ColumnCompressed;


template <typename Label = std::string>
class RowCompressed : public MultiLabelEncoded<uint64_t, Label> {
    friend ColumnCompressed<Label>;

    template <class A, typename L>
    friend std::unique_ptr<A> convert(RowCompressed<L>&&);

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    RowCompressed(uint64_t num_rows, bool sparse = false);

    void set_labels(Index i, const VLabels &labels);
    VLabels get_labels(Index i) const;

    void add_label(Index i, const Label &label);
    void add_labels(Index i, const VLabels &labels);
    void add_labels(const std::vector<Index> &indices, const VLabels &labels);
    void add_labels_fast(const std::vector<Index> &indices, const VLabels &labels);

    bool has_label(Index i, const Label &label) const;
    bool has_labels(Index i, const VLabels &labels) const;

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    void insert_rows(const std::vector<Index> &rows);

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    VLabels get_labels(const std::vector<Index> &indices,
                       double presence_ratio) const;

    uint64_t num_objects() const;
    size_t num_labels() const;
    uint64_t num_relations() const;

    static LabelEncoder<Label>* load_label_encoder(const std::string &filename);
    static void stream_counts(std::string filename,
                              uint64_t &num_objects,
                              uint64_t &num_relations,
                              bool sparse = false);
    class StreamRows {
      public:
        StreamRows(std::string filename, bool sparse);
        // TODO: decode and return VLabels?
        // or, make private, and make rowcompressed merge a friend method of RowCompressed class
        std::unique_ptr<std::vector<VectorRowBinMat::Row> > next_row() { return sr_->next_row(); };
      private:
        std::unique_ptr<VectorRowBinMat::StreamRows> sr_;
    };
    static uint64_t write_rows(std::string filename,
                           const LabelEncoder<Label> &label_encoder,
                           const std::function<void (const std::function<void (const std::vector<uint64_t> &)>&)> &callback,
                           //BinaryMatrix::GetRow &callback,
                           bool sparse);

  private:
    void reinitialize(uint64_t num_rows);

    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    std::unique_ptr<BinaryMatrixRowDynamic> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    static constexpr auto kExtension = kRowAnnotatorExtension;
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
