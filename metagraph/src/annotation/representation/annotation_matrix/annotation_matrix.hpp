#ifndef __ANNOTATION_MATRIX_HPP__
#define __ANNOTATION_MATRIX_HPP__

#include <memory>
#include <vector>

#include "annotation/binary_matrix/column_sparse//column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/base/annotation.hpp"

namespace mtg {
namespace annot {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelAnnotation<Label> {
  public:
    typedef BinaryMatrixType binary_matrix_type;
    using Index = typename MultiLabelAnnotation<Label>::Index;
    using VLabels = typename MultiLabelAnnotation<Label>::VLabels;

    StaticBinRelAnnotator(const StaticBinRelAnnotator&) = delete;
    StaticBinRelAnnotator(StaticBinRelAnnotator&&) = default;

    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder);

    template <typename... Args>
    StaticBinRelAnnotator(const LabelEncoder<Label> &label_encoder = {}, Args&&... args)
        : StaticBinRelAnnotator(std::make_unique<BinaryMatrixType>(std::forward<Args>(args)...),
                                label_encoder) {}

    void serialize(const std::string &filename) const override;
    bool load(const std::string &filename) override;
    static LabelEncoder<Label> read_label_encoder(const std::string &filename);
    // Dump columns to separate files in human-readable format
    bool dump_columns(const std::string &prefix, uint64_t num_threads = 1) const;

    uint64_t num_objects() const override { return matrix_->num_rows(); }
    uint64_t num_relations() const override { return matrix_->num_relations(); }

    void add_labels(const std::vector<Index> &, const VLabels &) override { except_dyn(); }
    void insert_rows(const std::vector<Index> &) override { except_dyn(); }

    const BinaryMatrixType& get_matrix() const override { return *matrix_; };

    std::unique_ptr<BinaryMatrixType> release_matrix()  { return std::move(matrix_); };

    std::string file_extension() const override { return kExtension; }

    static const std::string kExtension;

  private:
    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    using MultiLabelAnnotation<Label>::label_encoder_;
};

using ColumnCallback = std::function<void(uint64_t index,
                                          const std::string &label,
                                          std::unique_ptr<bit_vector>&& column)>;
bool merge_load_row_diff(const std::vector<std::string> &filenames,
                         const ColumnCallback &callback,
                         size_t num_threads);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_MATRIX_HPP__
