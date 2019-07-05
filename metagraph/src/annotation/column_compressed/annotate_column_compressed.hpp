#ifndef __ANNOTATE_COLUMN_COMPRESSED_HPP__
#define __ANNOTATE_COLUMN_COMPRESSED_HPP__


#include <cache.hpp>
#include <lru_cache_policy.hpp>
#include <progress_bar.hpp>

#include "annotate.hpp"
#include "bit_vector.hpp"


namespace annotate {

const char kColumnAnnotatorExtension[] = ".column.annodbg";

template <typename Label>
class RowCompressed;


template <typename Label = std::string>
class ColumnCompressed : public MultiLabelEncoded<uint64_t, Label> {
    template <class A, typename L>
    friend std::unique_ptr<A> convert(ColumnCompressed<L>&&);

    template <class A, typename L, class P>
    friend std::unique_ptr<A> convert_to_BRWT(ColumnCompressed<L>&&, P, size_t);

  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    ColumnCompressed(uint64_t num_rows = 0,
                     size_t num_columns_cached = 1,
                     bool verbose = false);

    ColumnCompressed(const ColumnCompressed&) = delete;
    ColumnCompressed& operator=(const ColumnCompressed&) = delete;

    ~ColumnCompressed();

    using MultiLabelEncoded<uint64_t, Label>::set;
    void set_labels(Index i, const VLabels &labels) override;
    VLabels get_labels(Index i) const override;

    void add_label(Index i, const Label &label) override;
    void add_labels(Index i, const VLabels &labels) override;
    void add_labels(const std::vector<Index> &indices,
                    const VLabels &labels) override;

    bool has_label(Index i, const Label &label) const override;

    // For each index i in indices, check if i has the label. Return
    // true if the finished callback evaluates true during execution.
    bool call_relations(const std::vector<Index> &indices,
                        const Label &label,
                        std::function<void(Index)> object_callback,
                        std::function<bool()> terminate = []() { return false; }) const;

    bool has_labels(Index i, const VLabels &labels) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;

    void insert_rows(const std::vector<Index> &rows) override;

    // For each pair (first, second) in the dictionary, renames
    // column |first| with |second| and merges the columns with matching names.
    void rename_labels(const std::unordered_map<Label, Label> &dict) override;

    uint64_t num_objects() const override;
    size_t num_labels() const override;
    uint64_t num_relations() const override;
    void call_objects(const Label &label,
                      std::function<void(Index)> callback) const override;

    void convert_to_row_annotator(const std::string &outfbase) const;
    void convert_to_row_annotator(RowCompressed<Label> *annotator,
                                  size_t num_threads = 1) const;

    void dump_columns(const std::string &prefix) const;

    const auto& data() const { return bitmatrix_; };

    std::unique_ptr<IterateRows> iterator() const override;

    const bitmap& get_column(const Label &label) const;

    std::string file_extension() const override { return kExtension; }

  private:
    void set(Index i, size_t j, bool value);
    bool is_set(Index i, size_t j) const;

    void add_labels(uint64_t begin, uint64_t end,
                    RowCompressed<Label> *annotator,
                    ProgressBar *progress_bar) const;
    void release();
    void flush() const;
    void flush(size_t j, const bitmap &annotation_curr);
    bitmap_dyn& decompress(size_t j);
    const bitmap& get_column(size_t j) const;

    std::vector<uint64_t> get_label_codes(Index i) const override;

    uint64_t num_rows_;

    std::vector<std::unique_ptr<bit_vector>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              bitmap_dyn*,
                              caches::LRUCachePolicy<size_t>> cached_columns_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    bool verbose_;

    static constexpr auto kExtension = kColumnAnnotatorExtension;
};

} // namespace annotate

#endif // __ANNOTATE_COLUMN_COMPRESSED_HPP__
