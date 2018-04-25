#include "annotate_k2.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "annotate_color_compressed.hpp"


namespace annotate {

sdsl::k2_tree<2> build_k2_tree(const std::vector<sdsl::sd_vector<>*> &bitmatrix) {
    if (!bitmatrix.size())
        return sdsl::k2_tree<2>();

    //TODO: optimize this
    std::vector<std::tuple<size_t, size_t>> result;
    for (size_t j = 0; j < bitmatrix.size(); ++j) {
        auto slct = sdsl::select_support_sd<>(bitmatrix[j]);
        auto rank = sdsl::rank_support_sd<>(bitmatrix[j]);
        auto num_set_bits = rank(bitmatrix[j]->size());

        for (uint64_t i = 1; i <= num_set_bits; ++i) {
            result.push_back(std::make_pair(slct(i), j));
        }
    }

    return sdsl::k2_tree<2>(result, std::max(bitmatrix.size(),
                                             bitmatrix[0]->size()));
}

template <typename Color, class Encoder>
K2Compressed<Color, Encoder>::K2Compressed() : color_encoder_(new Encoder()) {}

template <typename Color, class Encoder>
K2Compressed<Color, Encoder>
::K2Compressed(const ColorCompressed<Color, Encoder> &annotator)
      : color_encoder_(
            new Encoder(dynamic_cast<Encoder&>(*annotator.color_encoder_))
        ) {
    const_cast<ColorCompressed<Color, Encoder>&>(annotator).flush();

    k2_tree_ = build_k2_tree(annotator.bitmatrix_);
}

template <typename Color, class Encoder>
void K2Compressed<Color, Encoder>::set_coloring(Index i, const Coloring &coloring) {
    //TODO
    std::ignore = i;
    std::ignore = coloring;
    throw std::runtime_error("To be implemented");
}

template <typename Color, class Encoder>
typename K2Compressed<Color, Encoder>::Coloring
K2Compressed<Color, Encoder>::get_coloring(Index i) const {
    Coloring coloring;

    auto neigh = k2_tree_.neigh(i);

    for (auto adjacent : neigh) {
        coloring.push_back(color_encoder_->decode(adjacent));
    }

    return coloring;
}

template <typename Color, class Encoder>
void K2Compressed<Color, Encoder>::add_color(Index i, const Color &color) {
    //TODO
    std::ignore = i;
    std::ignore = color;
    throw std::runtime_error("To be implemented");
}

template <typename Color, class Encoder>
void K2Compressed<Color, Encoder>::add_colors(Index i, const Coloring &coloring) {
    for (const auto &color : coloring) {
        add_color(i, color);
    }
}

template <typename Color, class Encoder>
void K2Compressed<Color, Encoder>::add_colors(const std::vector<Index> &indices,
                                              const Coloring &coloring) {
    for (Index i : indices) {
        add_colors(i, coloring);
    }
}

template <typename Color, class Encoder>
bool K2Compressed<Color, Encoder>::has_color(Index i, const Color &color) const {
    try {
        return k2_tree_.adj(i, color_encoder_->encode(color));
    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
bool K2Compressed<Color, Encoder>::has_colors(Index i, const Coloring &coloring) const {
    try {
        for (const auto &color : coloring) {
            if (!has_color(i, color))
                return false;
        }
        return true;
    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
void K2Compressed<Color, Encoder>::serialize(const std::string &filename) const {
    std::ofstream outstream(filename + ".k2.annodbg");
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    color_encoder_->serialize(outstream);
    k2_tree_.serialize(outstream);
}

template <typename Color, class Encoder>
bool K2Compressed<Color, Encoder>::merge_load(const std::vector<std::string> &filenames) {
    std::ifstream instream(filenames.at(0) + ".k2.annodbg");
    if (!instream.good())
        return false;

    try {
        if (!color_encoder_->load(instream))
            return false;

        k2_tree_.load(instream);
        return true;
    } catch (...) {
        return false;
    }
}

// Get colors that occur at least in |discovery_ratio| colorings.
// If |discovery_ratio| = 0, return the union of colorings.
template <typename Color, class Encoder>
typename K2Compressed<Color, Encoder>::Coloring
K2Compressed<Color, Encoder>::aggregate_colors(const std::vector<Index> &indices,
                                               double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    // const size_t max_colors_missing = indices.size() - min_colors_discovered;

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        auto neigh = k2_tree_.neigh(i);
        for (auto adjacent : neigh) {
            encoded_counter[adjacent]++;
        }
    }

    Coloring filtered_colors;

    for (auto it = encoded_counter.begin(); it != encoded_counter.end(); ++it) {
        if (it->second >= min_colors_discovered)
            filtered_colors.push_back(color_encoder_->decode(it->first));
    }

    return filtered_colors;
}

// Count all colors collected from extracted colorings
// and return top |num_top| with the counts computed.
template <typename Color, class Encoder>
std::vector<std::pair<Color, size_t>>
K2Compressed<Color, Encoder>::get_most_frequent_colors(const std::vector<Index> &indices,
                                                       size_t num_top) const {
    std::vector<size_t> encoded_counter(color_encoder_->size(), 0);

    for (Index i : indices) {
        auto neigh = k2_tree_.neigh(i);
        for (auto adjacent : neigh) {
            encoded_counter[adjacent]++;
        }
    }

    std::vector<std::pair<size_t, size_t>> counts;
    counts.reserve(encoded_counter.size());
    for (size_t j = 0; j < encoded_counter.size(); ++j) {
        if (encoded_counter[j])
            counts.emplace_back(j, encoded_counter[j]);
    }

    // sort in decreasing order
    std::sort(counts.begin(), counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    counts.resize(std::min(counts.size(), num_top));

    std::vector<std::pair<Color, size_t>> top_counts;
    for (const auto &encoded_pair : counts) {
        top_counts.emplace_back(color_encoder_->decode(encoded_pair.first),
                                encoded_pair.second);
    }

    return top_counts;
}

template <typename Color, class Encoder>
size_t K2Compressed<Color, Encoder>::num_colors() const {
    return color_encoder_->size();
}

template <typename Color, class Encoder>
double K2Compressed<Color, Encoder>::sparsity() const {
    throw std::runtime_error("To be implemented");
    return 0;
}

template class K2Compressed<std::string, StringEncoder>;

} // namespace annotate
