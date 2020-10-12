#ifndef __ANNOTATE_COLUMN_COMPRESSED_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_HPP__

#include <mutex>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "common/vectors/bit_vector.hpp"
#include "common/vector.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"


namespace mtg {
namespace annot {

/**
 * Multithreading:
 *  The non-const methods must be called sequentially.
 *  Then, any subset of the public const methods can be called concurrently.
 */
template <typename Label = std::string>
class ColumnCompressed : public MultiLabelEncoded<Label> {
  public:
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;

    ColumnCompressed(uint64_t num_rows = 0,
                     size_t num_columns_cached = 1);

    ColumnCompressed(const ColumnCompressed&) = delete;
    ColumnCompressed& operator=(const ColumnCompressed&) = delete;

    ~ColumnCompressed();

    void set(Index i, const VLabels &labels) override;

    void add_labels(const std::vector<Index> &indices,
                    const VLabels &labels) override;
    void add_label_counts(const std::vector<Index> &indices,
                          const VLabels &labels,
                          const std::vector<uint32_t> &counts) override;

    bool has_label(Index i, const Label &label) const override;
    bool has_labels(Index i, const VLabels &labels) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;
    using ColumnCallback = std::function<void(uint64_t offset,
                                              const Label &,
                                              std::unique_ptr<bit_vector>&&)>;
    static bool merge_load(const std::vector<std::string> &filenames,
                           const ColumnCallback &callback,
                           size_t num_threads = 1);
    // Dump columns to separate files in human-readable format
    bool dump_columns(const std::string &prefix, size_t num_threads = 1) const;

    void insert_rows(const std::vector<Index> &rows) override;

    // For each pair (first, second) in the dictionary, renames
    // column |first| with |second| and merges the columns with matching names.
    void rename_labels(const tsl::hopscotch_map<Label, Label> &dict) override;

    uint64_t num_objects() const override;
    size_t num_labels() const override { return bitmatrix_.size(); }
    uint64_t num_relations() const override;

    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const override;

    /**
     * Return all labels for which counts are greater than or equal to |min_count|.
     * Stop counting if count is greater than |count_cap|.
     */
    std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
    count_labels(const std::vector<std::pair<Index, size_t>> &index_counts,
                 size_t min_count = 1,
                 size_t count_cap = std::numeric_limits<size_t>::max()) const override;

    const bitmap& get_column(const Label &label) const;

    /**
     * Returns a view of the current annotation matrix.
     * Warning: The returned object doesn't own its data and will become invalid when the
     * current object is destroyed.
     */
    const binmat::ColumnMajor& get_matrix() const override;

    /**
     * Returns the current annotation matrix. The data is moved into the return value,
     * which leaves the current object empty.
     */
    binmat::ColumnMajor release_matrix();

    std::string file_extension() const override { return kExtension; }

    static constexpr auto kExtension = ".column.annodbg";

  private:
    void set(Index i, size_t j, bool value);
    void flush() const;
    void flush(size_t j, const bitmap_builder &column_builder);
    bitmap_builder& decompress_builder(size_t j);
    bitmap_dyn& decompress_bitmap(size_t j);
    const bitmap& get_column(size_t j) const;

    uint64_t num_rows_;

    std::vector<std::unique_ptr<bit_vector>> bitmatrix_;
    mutable binmat::ColumnMajor annotation_matrix_view_;

    mutable std::mutex bitmap_conversion_mu_;
    mutable bool flushed_ = true;

    caches::fixed_sized_cache<size_t,
                              bitmap_builder*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    std::vector<sdsl::int_vector<>> relation_counts_;

    using MultiLabelEncoded<Label>::label_encoder_;
};

} // namespace annot
} // namespace mtg

#endif // __ANNOTATE_COLUMN_COMPRESSED_HPP__
