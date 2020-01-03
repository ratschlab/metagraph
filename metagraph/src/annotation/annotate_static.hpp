#ifndef __ANNOTATE_STATIC_HPP__
#define __ANNOTATE_STATIC_HPP__

#include <memory>
#include <vector>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"


namespace annotate {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelEncoded<Label> {
  public:
    typedef BinaryMatrixType binary_matrix_type;
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;
    using SetBitPositions = typename MultiLabelEncoded<Label>::SetBitPositions;

    explicit StaticBinRelAnnotator(size_t row_cache_size = 0);

    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder,
                          size_t row_cache_size = 0);

    bool has_label(Index i, const Label &label) const override;
    bool has_labels(Index i, const VLabels &labels) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;
    bool dump_columns(const std::string &prefix, uint64_t num_threads = 1) const;

    uint64_t num_objects() const override;
    uint64_t num_relations() const override;

    void set(Index, const VLabels &) override { except_dyn(); }
    void add_labels(const std::vector<Index> &, const VLabels &) override { except_dyn(); }
    void insert_rows(const std::vector<Index> &) override { except_dyn(); }

    const BinaryMatrixType& data() const { return *matrix_; }

    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const override;

    void reset_row_cache(size_t size);

    std::string file_extension() const override;

    const BinaryMatrix& get_matrix() const override { return *matrix_; };

  private:
    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    LabelEncoder<Label> &label_encoder_ { MultiLabelEncoded<Label>::label_encoder_ };

    SetBitPositions get_label_codes(Index i) const override;
    std::vector<SetBitPositions>
    get_label_codes(const std::vector<Index> &indices) const override;

    typedef caches::fixed_sized_cache<Index,
                                      SetBitPositions,
                                      caches::LRUCachePolicy<Index>> RowCacheType;
    mutable std::unique_ptr<RowCacheType> cached_rows_;

    static const std::string kExtension;
};

} // namespace annotate

#endif // __ANNOTATE_STATIC_HPP__
