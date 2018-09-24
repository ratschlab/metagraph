#ifndef __ANNOTATE_COLUMN_COMPRESSED_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Label>
class FastColumnCompressed;

template <typename Label>
class RowCompressed;


template <typename Label = std::string>
class ColumnCompressed : public MultiLabelEncoded<uint64_t, Label> {
    friend FastColumnCompressed<Label>;

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    ColumnCompressed(uint64_t num_rows,
                     size_t num_columns_cached = 1,
                     bool verbose = false);

    ColumnCompressed(const ColumnCompressed&) = delete;
    ColumnCompressed& operator=(const ColumnCompressed&) = delete;

    ~ColumnCompressed();

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

    // For each pair (first, second) in the dictionary, renames
    // column |first| with |second| and merges the columns with matching names.
    void rename_columns(const std::map<std::string, std::string> &dict);

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    VLabels get_labels(const std::vector<Index> &indices,
                       double presence_ratio) const;

    size_t num_labels() const;
    double sparsity() const;

    void convert_to_row_annotator(RowCompressed<Label> *annotator,
                                  size_t num_threads = 1) const;

  private:
    void set(Index i, size_t j, bool value);
    bool is_set(Index i, size_t j) const;
    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    void add_labels(uint64_t begin, uint64_t end,
                    RowCompressed<Label> *annotator) const;
    void release();
    void flush() const;
    void flush(size_t j, const std::vector<bool> &annotation_curr);
    std::vector<bool>& decompress(size_t j);

    uint64_t num_rows_;

    std::vector<std::unique_ptr<bit_vector_small>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              std::vector<bool>*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    LabelEncoder<Label> &label_encoder_
        = MultiLabelEncoded<uint64_t, Label>::label_encoder_;

    bool verbose_;

    static constexpr auto kExtension = ".column.annodbg";
};

} // namespace annotate

#endif // __ANNOTATE_COLUMN_COMPRESSED_HPP__
