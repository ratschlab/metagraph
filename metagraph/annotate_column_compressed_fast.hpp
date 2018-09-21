#ifndef __ANNOTATE_COLUMN_COMPRESSED_FAST_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_FAST_HPP__

#include "annotate_column_compressed.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Label = std::string>
class FastColumnCompressed : public MultiLabelAnnotation<uint64_t, Label> {
  public:
    using Index = typename MultiLabelAnnotation<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelAnnotation<uint64_t, Label>::VLabels;

    FastColumnCompressed(uint64_t num_rows,
                         size_t num_columns_cached = 1,
                         bool verbose = false);

    // Initialize from ColumnCompressed annotator
    FastColumnCompressed(ColumnCompressed<Label>&& annotator,
                         size_t num_columns_cached = 1,
                         bool verbose = false,
                         bool build_index = true);

    FastColumnCompressed(const FastColumnCompressed&) = delete;
    FastColumnCompressed& operator=(const FastColumnCompressed&) = delete;

    ~FastColumnCompressed();

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

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1)) const;

    size_t num_labels() const;
    double sparsity() const;

    // chooses an optimal value for |num_aux_cols| by default
    void rebuild_index(size_t num_aux_cols = -1);

  private:
    bool get_entry(Index i, size_t j) const;
    bool get_index_entry(Index i, size_t t) const;

    std::vector<size_t> get_row(Index i) const;
    std::vector<size_t> filter_row(Index i, const std::vector<bool> &filter) const;
    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    void flush() const;
    void flush(size_t j, const std::vector<bool> &annotation_curr);
    void flush_index(size_t t, const std::vector<bool> &annotation_curr);

    std::vector<bool>& decompress(size_t j);
    std::vector<bool>& decompress_index(size_t t);

    void update_index();

    uint64_t num_rows_;

    bool to_update_ = false;
    bool to_update_index_ = false;

    std::vector<std::unique_ptr<bit_vector_small>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              std::vector<bool>*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    std::vector<size_t> column_to_index_;
    std::vector<std::vector<uint64_t>> index_to_columns_;

    std::vector<std::unique_ptr<bit_vector_small>> index_;

    caches::fixed_sized_cache<size_t,
                              std::vector<bool>*,
                              caches::LRUCachePolicy<size_t>> cached_index_;

    LabelEncoder<Label> label_encoder_;

    bool verbose_;

    static constexpr auto kExtension = ".column.annodbg";
    static constexpr auto kIndexExtension = ".colindex.annodbg";
};

} // namespace annotate

#endif // __ANNOTATE_COLUMN_COMPRESSED_FAST_HPP__
