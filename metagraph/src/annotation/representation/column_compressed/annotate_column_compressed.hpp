#ifndef __ANNOTATE_COLUMN_COMPRESSED_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_HPP__

#include <mutex>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "common/vectors/bit_vector.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector.hpp"
#include "common/sorted_vector.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"


namespace mtg {
namespace annot {

/**
 * Multithreading:
 *  The non-const methods must be called sequentially (except add_label_counts).
 *  Then, any subset of the public const methods can be called concurrently.
 */
template <typename Label = std::string>
class ColumnCompressed : public MultiLabelAnnotation<Label> {
  public:
    typedef matrix::ColumnMajor binary_matrix_type;
    using Index = typename MultiLabelAnnotation<Label>::Index;
    using VLabels = typename MultiLabelAnnotation<Label>::VLabels;

    // |swap_dir| specifies a location on disk for tempprary buffers during the
    // column construction. If empty, no swap is used and all is stored in RAM.
    // |buffer_size_bytes| sets the buffer size per column (1 GB by default)
    ColumnCompressed(uint64_t num_rows = 0,
                     size_t num_columns_cached = 1,
                     const std::string &swap_dir = "",
                     uint64_t buffer_size_bytes = 1e7,
                     uint8_t count_width = 8,
                     size_t max_chunks_open = 2000);

    ColumnCompressed(sdsl::bit_vector&& column,
                     const std::string &column_label,
                     size_t num_columns_cached = 1,
                     const std::string &swap_dir = "",
                     uint64_t buffer_size_bytes = 1e7,
                     uint8_t count_width = 8,
                     size_t max_chunks_open = 2000);

    ColumnCompressed(std::vector<std::unique_ptr<bit_vector>>&& columns,
                     const LabelEncoder<Label> &label_encoder,
                     size_t num_columns_cached = 1,
                     const std::string &swap_dir = "",
                     uint64_t buffer_size_bytes = 1e7,
                     uint8_t count_width = 8,
                     size_t max_chunks_open = 2000);

    ColumnCompressed(const ColumnCompressed&) = delete;
    ColumnCompressed& operator=(const ColumnCompressed&) = delete;

    ~ColumnCompressed();

    void add_labels(const std::vector<Index> &indices,
                    const VLabels &labels) override;
    // for each label and index 'indices[i]' add count 'counts[i]'
    // thread-safe
    void add_label_counts(const std::vector<Index> &indices,
                          const VLabels &labels,
                          const std::vector<uint64_t> &counts) override;
    // for each label and index 'i' add numeric attribute 'coord'
    void add_label_coord(Index i, const VLabels &labels, uint64_t coord) override;
    void add_label_coords(const std::vector<std::pair<Index, uint64_t>> &coords,
                          const VLabels &labels) override;

    void serialize(const std::string &filename) const override;
    bool load(const std::string &filename) override;
    // the order of the columns may be changed when merging multiple annotators
    bool merge_load(const std::vector<std::string> &filenames);
    using ColumnCallback = std::function<void(uint64_t offset,
                                              const Label &,
                                              std::unique_ptr<bit_vector>&&)>;
    static bool merge_load(const std::vector<std::string> &filenames,
                           const ColumnCallback &callback,
                           size_t num_threads = 1);
    static size_t read_num_labels(const std::string &filename);
    static LabelEncoder<Label> read_label_encoder(const std::string &filename);

    using ValuesCallback = std::function<void(uint64_t offset,
                                              const Label &,
                                              sdsl::int_vector<>&&)>;
    static void load_column_values(const std::vector<std::string> &filenames,
                                   const ValuesCallback &callback,
                                   size_t num_threads = 1);

    using ColumnsValuesCallback = std::function<void(uint64_t offset,
                                                     const Label &,
                                                     std::unique_ptr<bit_vector>&&,
                                                     sdsl::int_vector<>&&)>;
    static void load_columns_and_values(const std::vector<std::string> &filenames,
                                        const ColumnsValuesCallback &callback,
                                        size_t num_threads = 1);

    using ColumnsDelimsValuesCallback = std::function<void(uint64_t offset,
                                                      const Label &,
                                                      std::unique_ptr<bit_vector>&&,
                                                      bit_vector_smart&&,
                                                      sdsl::int_vector<>&&)>;
    static void load_columns_delims_and_values(const std::vector<std::string> &filenames,
                                               const ColumnsDelimsValuesCallback &callback,
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

    const bitmap& get_column(const Label &label) const;

    /**
     * Returns a view of the current annotation matrix.
     * Warning: The returned object doesn't own its data and will become invalid when the
     * current object is destroyed.
     */
    const matrix::ColumnMajor& get_matrix() const override;

    /**
     * Returns the current annotation matrix. The data is moved into the return value,
     * which leaves the current object empty.
     */
    std::unique_ptr<matrix::ColumnMajor> release_matrix();

    std::string file_extension() const override { return kExtension; }

    static constexpr auto kExtension = ".column.annodbg";
    static constexpr auto kCountExtension = ".column.annodbg.counts";
    static constexpr auto kCoordExtension = ".column.annodbg.coords";

  private:
    void set(Index i, size_t j, bool value);
    void flush() const;
    void flush(size_t j, bitmap_builder *column_builder);
    bitmap_builder& decompress_builder(size_t j);
    bitmap_dyn& decompress_bitmap(size_t j);
    const bitmap& get_column(size_t j) const;
    void serialize_counts(const std::string &filename) const;
    void serialize_coordinates(const std::string &filename) const;

    uint64_t num_rows_;

    const std::string swap_dir_;
    const uint64_t buffer_size_bytes_;

    matrix::ColumnMajor matrix_;
    std::vector<std::unique_ptr<bit_vector>> &bitmatrix_ { matrix_.data() };

    mutable std::mutex bitmap_conversion_mu_;
    mutable bool flushed_ = true;

    caches::fixed_sized_cache<size_t,
                              bitmap_builder*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    mutable std::mutex counts_mu_;
    uint8_t count_width_;
    uint64_t max_count_;
    std::vector<sdsl::int_vector<>> relation_counts_;
    // depending on parameters, coords are stored in RAM or in chunks dumped to disk
    std::vector<common::SortedVector<std::pair<Index, uint64_t>>> coords_;
    std::vector<uint64_t> max_coord_; // stores max coord stored in |coords_|
    size_t max_chunks_open_;

    using MultiLabelAnnotation<Label>::label_encoder_;
};

} // namespace annot
} // namespace mtg

#endif // __ANNOTATE_COLUMN_COMPRESSED_HPP__
