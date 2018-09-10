#ifndef __ANNOTATE_COLOR_COMPRESSED_HPP__
#define __ANNOTATE_COLOR_COMPRESSED_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Color, class Encoder>
class K2Compressed;

template <typename Color, class Encoder>
class FastColorCompressed;

template <typename Color, class Encoder>
class RowCompressed;


template <typename Color = std::string, class Encoder = StringEncoder>
class ColorCompressed : public MultiColorAnnotation<uint64_t, Color> {
    friend K2Compressed<Color, Encoder>;
    friend FastColorCompressed<Color, Encoder>;

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

    void convert_to_row_annotator(RowCompressed<Color, Encoder> *annotator,
                                  size_t num_threads = 1) const;

  private:
    std::vector<uint64_t> count_colors(const std::vector<Index> &indices) const;

    static void add_labels(const std::vector<sdsl::select_support_sd<>> *select_columns,
                           const std::vector<sdsl::rank_support_sd<>> *rank_columns,
                           uint64_t begin, uint64_t end,
                           RowCompressed<Color, Encoder> *annotator);
    void release();
    void flush() const;
    void flush(size_t j, sdsl::bit_vector *annotation_curr);
    sdsl::bit_vector& decompress(size_t j);

    uint64_t num_rows_;

    std::vector<std::unique_ptr<sdsl::sd_vector<>>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              sdsl::bit_vector*,
                              caches::LRUCachePolicy<size_t>> cached_colors_;

    std::unique_ptr<ColorEncoder<Color>> color_encoder_;

    bool verbose_;

    static constexpr auto kExtension = ".color.annodbg";
};

} // namespace annotate

#endif // __ANNOTATE_COLOR_COMPRESSED_HPP__
