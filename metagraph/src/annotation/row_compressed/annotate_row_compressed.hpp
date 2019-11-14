#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <memory>

#include <Eigen/Sparse>

#include "annotate.hpp"
#include "binary_matrix.hpp"


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
    friend void merge(std::vector<std::unique_ptr<MultiLabelEncoded<uint64_t, L>>>&&,
                      const std::vector<std::string>&,
                      const std::string&);

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;
    using IterateRows = typename MultiLabelEncoded<uint64_t, Label>::IterateRows;
    using SetBitPositions = typename MultiLabelEncoded<uint64_t, Label>::SetBitPositions;

    RowCompressed(uint64_t num_rows = 0, bool sparse = false);

    typedef std::function<void(Index, SetBitPositions&&)> CallRow;
    RowCompressed(uint64_t num_rows,
                  const std::vector<Label> &labels,
                  std::function<void(CallRow)> call_rows);

    void set_labels(Index i, const VLabels &labels);

    void add_label(Index i, const Label &label);
    void add_labels(Index i, const VLabels &labels);
    void add_labels(const std::vector<Index> &indices, const VLabels &labels);
    void add_labels_fast(const std::vector<Index> &indices, const VLabels &labels);

    bool has_label(Index i, const Label &label) const;
    bool has_labels(Index i, const VLabels &labels) const;

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    void insert_rows(const std::vector<Index> &rows);

    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const;

    uint64_t num_objects() const;
    size_t num_labels() const;
    uint64_t num_relations() const;

    std::string file_extension() const { return kExtension; }

    static void write_rows(std::string filename,
                           const LabelEncoder<Label> &label_encoder,
                           const std::function<void(BinaryMatrix::RowCallback)> &call_rows);

  private:
    void reinitialize(uint64_t num_rows);

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

    template <typename RowType = BinaryMatrix::SetBitPositions>
    class StreamRows {
      public:
        explicit StreamRows(std::string filename);

        // return null after all rows have been called
        RowType* next_row() { return sr_->next_row(); }

      private:
        std::unique_ptr<::StreamRows<RowType>> sr_;
    };

    SetBitPositions get_label_codes(Index i) const {
        return matrix_->get_row(i);
    }

    static constexpr auto kExtension = kRowAnnotatorExtension;
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
