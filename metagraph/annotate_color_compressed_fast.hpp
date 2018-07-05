#ifndef __ANNOTATE_COLOR_COMPRESSED_FAST_HPP__
#define __ANNOTATE_COLOR_COMPRESSED_FAST_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "annotate_color_compressed.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Color = std::string, class Encoder = StringEncoder>
class FastColorCompressed : public MultiColorAnnotation<uint64_t, Color> {
  public:
    using Index = typename MultiColorAnnotation<uint64_t, Color>::Index;
    using Coloring = typename MultiColorAnnotation<uint64_t, Color>::Coloring;

    FastColorCompressed(uint64_t num_rows,
                        size_t num_columns_cached = 1,
                        bool verbose = false);

    // Initialize from ColorCompressed annotator
    FastColorCompressed(ColorCompressed<Color, Encoder>&& annotator,
                        size_t num_columns_cached = 1,
                        bool verbose = false,
                        bool build_index = true);

    FastColorCompressed(const FastColorCompressed&) = delete;
    FastColorCompressed& operator=(const FastColorCompressed&) = delete;

    ~FastColorCompressed();

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

    // chooses an optimal value for |num_aux_cols| by default
    void rebuild_index(size_t num_aux_cols = -1);

  private:
    bool get_entry(Index i, size_t j) const;
    bool get_index_entry(Index i, size_t t) const;

    std::vector<size_t> get_row(Index i) const;
    std::vector<size_t> filter_row(Index i, const std::vector<bool> &filter) const;
    std::vector<uint64_t> count_colors(const std::vector<Index> &indices) const;

    void flush();
    void flush(size_t j, sdsl::bit_vector *annotation_curr);
    void flush_index(size_t t, sdsl::bit_vector *annotation_curr);

    sdsl::bit_vector& decompress(size_t j);
    sdsl::bit_vector& decompress_index(size_t t);

    void update_index();

    uint64_t num_rows_;

    bool to_update_ = false;
    bool to_update_index_ = false;

    std::vector<std::unique_ptr<sdsl::sd_vector<>>> bitmatrix_;

    caches::fixed_sized_cache<size_t,
                              sdsl::bit_vector*,
                              caches::LRUCachePolicy<size_t>> cached_colors_;

    std::vector<size_t> column_to_index_;
    std::vector<std::vector<uint64_t>> index_to_columns_;

    std::vector<std::unique_ptr<sdsl::sd_vector<>>> index_;

    caches::fixed_sized_cache<size_t,
                              sdsl::bit_vector*,
                              caches::LRUCachePolicy<size_t>> cached_index_;

    std::unique_ptr<ColorEncoder<Color>> color_encoder_;

    bool verbose_;

    static const std::string kExtension;
    static const std::string kIndexExtension;
};

} // namespace annotate

#endif // __ANNOTATE_COLOR_COMPRESSED_FAST_HPP__
