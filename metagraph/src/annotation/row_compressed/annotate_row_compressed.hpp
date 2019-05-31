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
    template <class A>
    friend std::unique_ptr<A> convert(const std::string&);
    template <class A, typename L>
    friend void merge(const std::vector<const MultiLabelEncoded<uint64_t, L>*>&,
                      const std::vector<std::string>&,
                      const std::string&);
    template <typename L>
    friend void merge_rows(const std::vector<const LabelEncoder<L> *>&,
                           std::function<const std::vector<uint64_t>(uint64_t)>,
                           uint64_t,
                           const std::string&);

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    RowCompressed(uint64_t num_rows = 0, bool sparse = false);

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

    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const;

    uint64_t num_objects() const;
    size_t num_labels() const;
    uint64_t num_relations() const;

  private:
    void reinitialize(uint64_t num_rows);

    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    std::unique_ptr<BinaryMatrixRowDynamic> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    static std::unique_ptr<LabelEncoder<Label>>
    load_label_encoder(const std::string &filename);
    static std::unique_ptr<LabelEncoder<Label>>
    load_label_encoder(std::istream &instream);

    static void stream_counts(const std::string &filename,
                              uint64_t *num_objects,
                              uint64_t *num_relations);

    class StreamRows {
      public:
        explicit StreamRows(std::string filename);
        // return null after all rows have been called
        std::vector<VectorRowBinMat::Column>* next_row() { return sr_->next_row(); };
      private:
        std::unique_ptr<VectorRowBinMat::StreamRows> sr_;
    };

    static void write_rows(std::string filename,
                           const LabelEncoder<Label> &label_encoder,
                           const std::function<void(BinaryMatrix::RowCallback&)> &call_rows);

    std::vector<uint64_t> get_label_indexes(Index i) const {
        return matrix_->get_row(i);
    }

    static constexpr auto kExtension = kRowAnnotatorExtension;
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
