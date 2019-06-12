#ifndef __ANNOTATE_STATIC_HPP__
#define __ANNOTATE_STATIC_HPP__

#include <memory>
#include <vector>

#include "annotate.hpp"


namespace annotate {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelEncoded<uint64_t, Label> {
  public:
    typedef BinaryMatrixType binary_matrix_type;
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    StaticBinRelAnnotator() : matrix_(new BinaryMatrixType()) {}
    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder);

    bool has_label(Index i, const Label &label) const override;
    bool has_labels(Index i, const VLabels &labels) const override;

    VLabels get_labels(Index i) const override;
    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    VLabels get_labels(const std::vector<Index> &indices,
                       double presence_ratio) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;

    uint64_t num_objects() const override;
    size_t num_labels() const override;
    uint64_t num_relations() const override;

    void set_labels(Index, const VLabels &) override { except_dyn(); }
    void add_label(Index, const Label &) override { except_dyn(); }
    void add_labels(Index, const VLabels &) override { except_dyn(); }
    void add_labels(const std::vector<Index> &, const VLabels &) override { except_dyn(); }
    void insert_rows(const std::vector<Index> &) override { except_dyn(); }

    const BinaryMatrixType& data() const { return *matrix_; }

    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const override;

    std::string file_extension() const override { return kExtension; }

  private:
    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const override;

    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    std::vector<uint64_t> get_label_indexes(Index i) const override {
        return matrix_->get_row(i);
    }

    static const std::string kExtension;
};

} // namespace annotate

#endif // __ANNOTATE_STATIC_HPP__
