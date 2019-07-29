#ifndef __ANNOTATE_STATIC_HPP__
#define __ANNOTATE_STATIC_HPP__

#include <memory>
#include <vector>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"


namespace annotate {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelEncoded<uint64_t, Label> {
  public:
    typedef BinaryMatrixType binary_matrix_type;
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    StaticBinRelAnnotator(size_t row_cache_size = 0);
    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder,
                          size_t row_cache_size = 0);

    bool has_label(Index i, const Label &label) const override;
    bool has_labels(Index i, const VLabels &labels) const override;

    VLabels get_labels(Index i) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;
    void dump_columns(const std::string &prefix,
                      bool binary = false,
                      uint64_t num_threads = 1) const;

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

    void reset_row_cache(size_t size);

    std::string file_extension() const override;

  private:
    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    std::vector<uint64_t> get_label_codes(Index i) const override;

    typedef caches::fixed_sized_cache<Index,
                                      std::vector<uint64_t>,
                                      caches::LRUCachePolicy<Index>> RowCacheType;
    mutable std::unique_ptr<RowCacheType> cached_rows_;

    static const std::string kExtension;
};

} // namespace annotate

#endif // __ANNOTATE_STATIC_HPP__
