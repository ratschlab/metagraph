#ifndef __ANNOTATE_COLUMN_COMPRESSED_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_HPP__

#include <mutex>

#include <cache.hpp>
#include <lru_cache_policy.hpp>
#include <progress_bar.hpp>

#include "common/vectors/bit_vector.hpp"
#include "common/vector.hpp"
#include "annotation/annotate.hpp"


namespace annotate {

const char kColumnAnnotatorExtension[] = ".column.annodbg";

template <typename Label>
class RowCompressed;


/**
 * Multithreading:
 *  The non-const methods must be called sequentially.
 *  Then, any subset of the public const methods can be called concurrently.
 */
template <typename Label = std::string>
class ColumnCompressed : public MultiLabelEncoded<Label> {
    template <class A, typename L>
    friend std::unique_ptr<A> convert(ColumnCompressed<L>&&);

    template <class A, typename L, class P>
    friend std::unique_ptr<A> convert_to_BRWT(ColumnCompressed<L>&&, P, size_t, size_t);

  public:
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;
    using IterateRows = typename MultiLabelEncoded<Label>::IterateRows;
    using SetBitPositions = typename MultiLabelEncoded<Label>::SetBitPositions;

    ColumnCompressed(uint64_t num_rows = 0,
                     size_t num_columns_cached = 1,
                     bool verbose = false);

    ColumnCompressed(const ColumnCompressed&) = delete;
    ColumnCompressed& operator=(const ColumnCompressed&) = delete;

    ~ColumnCompressed();

    void set(Index i, const VLabels &labels) override;

    void add_labels(const std::vector<Index> &indices,
                    const VLabels &labels) override;

    bool has_label(Index i, const Label &label) const override;

    // For each index i in indices, check if i has the label.
    void call_relations(const std::vector<Index> &indices,
                        const Label &label,
                        std::function<void(Index, bool)> callback,
                        std::function<bool()> terminate = []() { return false; }) const;

    bool has_labels(Index i, const VLabels &labels) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;

    void insert_rows(const std::vector<Index> &rows) override;

    // For each pair (first, second) in the dictionary, renames
    // column |first| with |second| and merges the columns with matching names.
    void rename_labels(const std::unordered_map<Label, Label> &dict) override;

    uint64_t num_objects() const override;
    uint64_t num_relations() const override;
    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const override;

    /**
     * Return all labels for which counts are greater than or equal to |min_count|.
     * Stop counting if count is greater than |count_cap|.
     */
    std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
    count_labels(const std::unordered_map<Index, size_t> &index_counts,
                 size_t min_count = 1,
                 size_t count_cap = std::numeric_limits<size_t>::max()) const override;

    void convert_to_row_annotator(const std::string &outfbase) const;
    void convert_to_row_annotator(RowCompressed<Label> *annotator,
                                  size_t num_threads = 1) const;

    bool dump_columns(const std::string &prefix, size_t num_threads = 1) const;

    const auto& data() const { return bitmatrix_; };

    std::unique_ptr<IterateRows> iterator() const override;

    const bitmap& get_column(const Label &label) const;

    std::string file_extension() const override { return kExtension; }

  private:
    void set(Index i, size_t j, bool value);

    void add_labels(uint64_t begin, uint64_t end,
                    RowCompressed<Label> *annotator,
                    ProgressBar *progress_bar) const;
    void release();
    void flush() const;
    void flush(size_t j, const bitmap_builder &column_builder);
    bitmap_builder& decompress_builder(size_t j);
    bitmap_dyn& decompress_bitmap(size_t j);
    const bitmap& get_column(size_t j) const;

    SetBitPositions get_label_codes(Index i) const override;
    std::vector<SetBitPositions>
    get_label_codes(const std::vector<Index> &indices) const override;

    uint64_t num_rows_;

    std::vector<std::unique_ptr<bit_vector>> bitmatrix_;

    mutable std::mutex bitmap_conversion_mu_;
    mutable bool flushed_ = true;

    caches::fixed_sized_cache<size_t,
                              bitmap_builder*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    LabelEncoder<Label> &label_encoder_ { MultiLabelEncoded<Label>::label_encoder_ };

    bool verbose_;

    static constexpr auto kExtension = kColumnAnnotatorExtension;
};

} // namespace annotate

#endif // __ANNOTATE_COLUMN_COMPRESSED_HPP__
