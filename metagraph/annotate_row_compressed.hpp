#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include <Eigen/Sparse>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Label>
class ColumnCompressed;


class RowMajorSparseBinaryMatrix {
  public:
    virtual ~RowMajorSparseBinaryMatrix() {}

    virtual void set_bit(size_t i, size_t j) = 0;
    virtual bool is_set_bit(size_t i, size_t j) const = 0;

    virtual size_t select(size_t i, size_t k) const = 0;

    virtual size_t size() const = 0;
    virtual size_t size(size_t i) const = 0;

    virtual void clear(size_t i) = 0;

    virtual void reinitialize(size_t num_rows) = 0;

    virtual void insert_rows(const std::vector<uint64_t> &rows) = 0;
};


template <typename Label = std::string>
class RowCompressed : public MultiLabelEncoded<uint64_t, Label> {
    friend ColumnCompressed<Label>;

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    RowCompressed(uint64_t num_rows, bool sparse = false);

    void set_labels(Index i, const VLabels &labels);
    VLabels get_labels(Index i) const;

    void add_label(Index i, const Label &label);
    void add_labels(Index i, const VLabels &labels);
    void add_labels(const std::vector<Index> &indices, const VLabels &labels);

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
    double sparsity() const;

  private:
    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    std::unique_ptr<RowMajorSparseBinaryMatrix> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    static constexpr auto kExtension = ".row.annodbg";
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
