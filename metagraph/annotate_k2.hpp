#ifndef __ANNOTATE_K2_HPP__
#define __ANNOTATE_K2_HPP__

#include <sdsl/k2_tree.hpp>

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Color, class Encoder>
class ColorCompressed;


template <typename Color = std::string, class Encoder = StringEncoder>
class K2Compressed : public MultiColorAnnotation<uint64_t, Color> {
  public:
    using Index = typename MultiColorAnnotation<uint64_t, Color>::Index;
    using Coloring = typename MultiColorAnnotation<uint64_t, Color>::Coloring;

    K2Compressed();
    K2Compressed(const ColorCompressed<Color, Encoder> &color_compressed);

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

  private:
    sdsl::k2_tree<2> k2_tree_;
    std::unique_ptr<ColorEncoder<Color>> color_encoder_;
};

} // namespace annotate

#endif // __ANNOTATE_K2_HPP__
