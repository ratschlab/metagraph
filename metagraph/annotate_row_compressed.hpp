#ifndef __ANNOTATE_ROW_COMPRESSED_HPP__
#define __ANNOTATE_ROW_COMPRESSED_HPP__

#include "annotate.hpp"
#include "dbg_succinct.hpp"


namespace annotate {

template <typename Color = std::string, class Encoder = StringEncoder>
class RowCompressed : public MultiColorAnnotation<uint64_t, Color> {
  public:
    using Index = typename MultiColorAnnotation<uint64_t, Color>::Index;
    using Coloring = typename MultiColorAnnotation<uint64_t, Color>::Coloring;

    RowCompressed(uint64_t num_rows);

    void set_coloring(Index i, const Coloring &coloring);
    Coloring get_coloring(Index i) const;

    void add_color(Index i, const Color &color);
    void add_colors(Index i, const Coloring &coloring);
    void add_colors(const std::vector<Index> &indices, const Coloring &coloring);

    bool has_color(Index i, const Color &color) const;
    bool has_colors(Index i, const Coloring &coloring) const;

    bool load(const std::string &filename);
    void serialize(const std::string &filename) const;

    // Get colors that occur at least in |discovery_ratio| colorings.
    // If |discovery_ratio| = 0, return the union of colorings.
    Coloring aggregate_colors(const std::vector<Index> &indices,
                              double discovery_ratio = 1) const;

    // Count all colors collected from extracted colorings
    // and return top |num_top| with the counts computed.
    std::vector<std::pair<Color, size_t>>
    get_most_frequent_colors(const std::vector<Index> &indices,
                             size_t num_top = static_cast<size_t>(-1)) const;

  private:
    std::vector<std::vector<uint32_t>> encoded_colorings_;
    std::unique_ptr<ColorEncoder<Color>> color_encoder_;
};

} // namespace annotate

#endif // __ANNOTATE_ROW_COMPRESSED_HPP__
