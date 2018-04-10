#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::StringSerialisation;


namespace annotate {

template <typename Color, class Encoder>
RowCompressed<Color, Encoder>::RowCompressed(uint64_t num_rows)
      : encoded_colorings_(num_rows), color_encoder_(new Encoder()) {}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::set_coloring(Index i, const Coloring &coloring) {
    assert(i < encoded_colorings_.size());

    encoded_colorings_[i].clear();
    for (const auto &color : coloring) {
        encoded_colorings_[i].push_back(color_encoder_->encode(color, true));
    }
}

template <typename Color, class Encoder>
typename RowCompressed<Color, Encoder>::Coloring
RowCompressed<Color, Encoder>::get_coloring(Index i) const {
    Coloring coloring;
    for (auto code : encoded_colorings_[i]) {
        coloring.push_back(color_encoder_->decode(code));
    }
    return coloring;
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_color(Index i, const Color &color) {
    auto code = color_encoder_->encode(color, true);
    if (std::find(encoded_colorings_[i].begin(),
                  encoded_colorings_[i].end(), code)
            == encoded_colorings_[i].end())
        encoded_colorings_[i].push_back(code);
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_colors(Index i, const Coloring &coloring) {
    for (const auto &color : coloring) {
        add_color(i, color);
    }
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_colors(const std::vector<Index> &indices,
                                               const Coloring &coloring) {
    for (Index i : indices) {
        add_colors(i, coloring);
    }
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::has_color(Index i, const Color &color) const {
    try {
        return std::find(encoded_colorings_[i].begin(),
                         encoded_colorings_[i].end(),
                         color_encoder_->encode(color))
                != encoded_colorings_[i].end();
    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::has_colors(Index i, const Coloring &coloring) const {
    std::set<size_t> querying_codes;
    try {
        for (const auto &color : coloring) {
            querying_codes.insert(color_encoder_->encode(color));
        }
    } catch (...) {
        return false;
    }
    std::set<size_t> encoded_coloring(encoded_colorings_[i].begin(),
                                      encoded_colorings_[i].end());
    return std::includes(encoded_coloring.begin(), encoded_coloring.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::serialize(const std::string &filename) const {
    std::ofstream outstream(filename + ".row.annodbg");

    NumberSerialisation::serialiseNumber(outstream, encoded_colorings_.size());

    color_encoder_->serialize(outstream);

    std::vector<uint32_t> full_vector;
    for (const auto &indices : encoded_colorings_) {
        for (auto j : indices) {
            full_vector.push_back(j + 1);
        }
        full_vector.push_back(0);
    }
    serialize_number_vector(outstream,
                            full_vector,
                            std::log2(color_encoder_->size() + 1) + 1);
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::load(const std::string &filename) {
    std::ifstream instream(filename + ".row.annodbg");
    if (!instream.good())
        return false;

    try {
        size_t num_rows = NumberSerialisation::deserialiseNumber(instream);
        encoded_colorings_.clear();
        encoded_colorings_.resize(num_rows);

        if (!color_encoder_->load(instream))
            return false;

        auto full_vector = load_number_vector<uint32_t>(instream);
        for (size_t j = 0, i = 0; i < full_vector.size(); ++i) {
            if (full_vector[i]) {
                encoded_colorings_[j].push_back(full_vector[i] - 1);
            } else {
                encoded_colorings_[j++].shrink_to_fit();
            }
        }
        return true;
    } catch (...) {
        return false;
    }
}

// Get colors that occur at least in |discovery_ratio| colorings.
// If |discovery_ratio| = 0, return the union of colorings.
template <typename Color, class Encoder>
typename RowCompressed<Color, Encoder>::Coloring
RowCompressed<Color, Encoder>::aggregate_colors(const std::vector<Index> &indices,
                                                double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    // const size_t max_colors_missing = indices.size() - min_colors_discovered;

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        for (auto code : encoded_colorings_[i]) {
            encoded_counter[code]++;
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
RowCompressed<Color, Encoder>::get_most_frequent_colors(const std::vector<Index> &indices,
                                                        size_t num_top) const {
    std::vector<size_t> encoded_counter(color_encoder_->size(), 0);

    for (Index i : indices) {
        for (auto code : encoded_colorings_[i]) {
            encoded_counter[code]++;
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

template class RowCompressed<std::string, StringEncoder>;

} // namespace annotate
