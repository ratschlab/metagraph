#ifndef __ANNOTATE_COLOR_COMPRESSED_HPP__
#define __ANNOTATE_COLOR_COMPRESSED_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Color>
class FastColorCompressed;

template <typename Color>
class RowCompressed;


template <typename Color = std::string>
class ColorCompressed : public MultiColorAnnotation<uint64_t, Color> {
    friend FastColorCompressed<Color>;

  public:
    using Index = typename MultiColorAnnotation<uint64_t, Color>::Index;
    using Coloring = typename MultiColorAnnotation<uint64_t, Color>::Coloring;

    ColorCompressed(uint64_t num_rows,
                    size_t num_columns_cached = 1,
                    bool verbose = false);

    ColorCompressed(const ColorCompressed&) = delete;
    ColorCompressed& operator=(const ColorCompressed&) = delete;

    ~ColorCompressed();

    void set_coloring(Index i, const Coloring &coloring);
    Coloring get_coloring(Index i) const;

    void add_color(Index i, const Color &color);
    void add_colors(Index i, const Coloring &coloring);
    void add_colors(const std::vector<Index> &indices, const Coloring &coloring);

    bool has_color(Index i, const Color &color) const;
    bool has_colors(Index i, const Coloring &coloring) const;

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    void insert_rows(const std::vector<Index> &rows);

    // For each pair (first, second) in the dictionary, renames
    // column |first| with |second| and merges the columns with matching names.
    void rename_columns(const std::map<std::string, std::string> &dict);

    // Get colors that occur at least in |discovery_ratio| colorings.
    // If |discovery_ratio| = 0, return the union of colorings.
    Coloring aggregate_colors(const std::vector<Index> &indices,
                              double discovery_ratio = 1) const;

    // Count all colors collected from extracted colorings
    // and return top |num_top| with the counts computed.
    std::vector<std::pair<Color, size_t>>
    get_most_frequent_colors(const std::vector<Index> &indices,
                             size_t num_top = static_cast<size_t>(-1)) const;

    size_t num_colors() const;
    double sparsity() const;

    void convert_to_row_annotator(RowCompressed<Color> *annotator,
                                  size_t num_threads = 1) const;

  private:
    std::vector<uint64_t> count_colors(const std::vector<Index> &indices) const;

    void add_labels(uint64_t begin, uint64_t end,
                    RowCompressed<Color> *annotator) const;
    void release();
    void flush() const;
    void flush(size_t j, const std::vector<bool> &annotation_curr);
    std::vector<bool>& decompress(size_t j);

    uint64_t num_rows_;

    std::vector<std::unique_ptr<bit_vector_small>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              std::vector<bool>*,
                              caches::LRUCachePolicy<size_t>> cached_colors_;

    LabelEncoder<Color> color_encoder_;

    bool verbose_;

    static constexpr auto kExtension = ".color.annodbg";
};

} // namespace annotate

#endif // __ANNOTATE_COLOR_COMPRESSED_HPP__
