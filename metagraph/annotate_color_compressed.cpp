#include "annotate_color_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::StringSerialisation;


namespace annotate {

template <typename Color, class Encoder>
ColorCompressed<Color, Encoder>::ColorCompressed(uint64_t num_rows,
                                                 size_t num_columns_cached)
      : num_rows_(num_rows),
        cached_colors_(
            num_columns_cached,
            caches::LRUCachePolicy<size_t>(),
            [this](size_t j, sdsl::bit_vector *col_uncompressed) {
                this->flush(j, col_uncompressed);
                delete col_uncompressed;
            }
        ),
        color_encoder_(new Encoder()) {}

template <typename Color, class Encoder>
ColorCompressed<Color, Encoder>::~ColorCompressed() {
    cached_colors_.Clear();

    for (auto *column : bitmatrix_) {
        if (column)
            delete column;
    }
}

template <typename Color, class Encoder>
void
ColorCompressed<Color, Encoder>
::set_coloring(Index i, const Coloring &coloring) {
    assert(i < num_rows_);

    for (size_t j = 0; j < color_encoder_->size(); ++j) {
        uncompress(j)[i] = 0;
    }
    for (const auto &color : coloring) {
        uncompress(color_encoder_->encode(color, true))[i] = 1;
    }
}

template <typename Color, class Encoder>
typename ColorCompressed<Color, Encoder>::Coloring
ColorCompressed<Color, Encoder>::get_coloring(Index i) const {
    assert(i < num_rows_);

    const_cast<ColorCompressed*>(this)->flush();

    Coloring coloring;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j]);

        if ((*bitmatrix_[j])[i])
            coloring.push_back(color_encoder_->decode(j));
    }
    return coloring;
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::add_color(Index i, const Color &color) {
    assert(i < num_rows_);

    uncompress(color_encoder_->encode(color, true))[i] = 1;
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::add_colors(Index i, const Coloring &coloring) {
    for (const auto &color : coloring) {
        add_color(i, color);
    }
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::add_colors(const std::vector<Index> &indices,
                                                 const Coloring &coloring) {
    for (const auto &color : coloring) {
        for (Index i : indices) {
            add_color(i, color);
        }
    }
}

template <typename Color, class Encoder>
bool
ColorCompressed<Color, Encoder>
::has_color(Index i, const Color &color) const {
    try {
        auto j = color_encoder_->encode(color);
        if (cached_colors_.Cached(j)) {
            return (*cached_colors_.Get(j))[i];
        } else {
            return (*bitmatrix_[j])[i];
        }
    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
bool
ColorCompressed<Color, Encoder>
::has_colors(Index i, const Coloring &coloring) const {
    for (const auto &color : coloring) {
        if (!has_color(i, color))
            return false;
    }
    return true;
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::serialize(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();

    std::ofstream outstream(filename + ".color.annodbg");

    NumberSerialisation::serialiseNumber(outstream, num_rows_);

    color_encoder_->serialize(outstream);

    for (auto *column : bitmatrix_) {
        assert(column);
        column->serialize(outstream);
    }
}


template <typename Color, class Encoder>
bool ColorCompressed<Color, Encoder>
::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_colors_.Clear();

    for (auto *column : bitmatrix_) {
        if (column)
            delete column;
    }

    color_encoder_.reset(new Encoder());

    try {
        for (auto filename : filenames) {
            // load annotation from file
            std::ifstream instream(filename + ".color.annodbg");
            if (!instream.good())
                return false;

            if (filename == filenames.at(0)) {
                num_rows_ = NumberSerialisation::deserialiseNumber(instream);
            } else if (num_rows_ != NumberSerialisation::deserialiseNumber(instream)) {
                return false;
            }

            Encoder color_encoder_load;
            if (!color_encoder_load.load(instream))
                return false;

            // extend the color dictionary
            for (size_t c = 0; c < color_encoder_load.size(); ++c) {
                color_encoder_->encode(color_encoder_load.decode(c), true);
            }

            // update the existing and add some new columns
            bitmatrix_.resize(color_encoder_->size(), NULL);
            for (size_t c = 0; c < color_encoder_load.size(); ++c) {
                size_t col = color_encoder_->encode(color_encoder_load.decode(c));

                auto *new_column = new sdsl::sd_vector<>();
                new_column->load(instream);

                if (bitmatrix_.at(col)) {
                    auto &column = uncompress(col);

                    auto slct = sdsl::select_support_sd<>(new_column);
                    auto rank = sdsl::rank_support_sd<>(new_column);
                    auto num_set_bits = rank(new_column->size());

                    for (uint64_t i = 1; i <= num_set_bits; ++i) {
                        assert(slct(i) < column.size());

                        column[slct(i)] = 1;
                    }
                    delete new_column;
                } else {
                    bitmatrix_.at(col) = new_column;
                }
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
typename ColorCompressed<Color, Encoder>::Coloring
ColorCompressed<Color, Encoder>
::aggregate_colors(const std::vector<Index> &indices,
                   double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const_cast<ColorCompressed*>(this)->flush();

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    const size_t max_colors_missing = indices.size() - min_colors_discovered;

    Coloring filtered_colors;

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        uint64_t discovered = 0;
        uint64_t missing = 0;

        for (Index i : indices) {
            if ((*bitmatrix_[j])[i]) {
                if (++discovered >= min_colors_discovered) {
                    filtered_colors.push_back(color_encoder_->decode(j));
                    break;
                }
            } else if (++missing > max_colors_missing) {
                break;
            }
        }
    }

    return filtered_colors;
}

// Count all colors collected from extracted colorings
// and return top |num_top| with the counts computed.
template <typename Color, class Encoder>
std::vector<std::pair<Color, size_t>>
ColorCompressed<Color, Encoder>
::get_most_frequent_colors(const std::vector<Index> &indices,
                           size_t num_top) const {
    const_cast<ColorCompressed*>(this)->flush();

    std::vector<uint64_t> counter(bitmatrix_.size(), 0);
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        for (Index i : indices) {
            if ((*bitmatrix_[j])[i])
                counter[j]++;
        }
    }

    std::vector<std::pair<size_t, size_t>> counts;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        if (counter[j])
            counts.emplace_back(j, counter[j]);
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
void ColorCompressed<Color, Encoder>::flush() {
    for (const auto &cached_vector : cached_colors_) {
        flush(cached_vector.first, cached_vector.second);
    }
    assert(bitmatrix_.size() == color_encoder_->size());
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::flush(size_t j, sdsl::bit_vector *vector) {
    assert(vector);
    assert(cached_colors_.Cached(j));

    while (bitmatrix_.size() <= j) {
        bitmatrix_.push_back(NULL);
    }

    if (bitmatrix_[j])
        delete bitmatrix_[j];

    bitmatrix_[j] = new sdsl::sd_vector<>(*vector);
}

template <typename Color, class Encoder>
sdsl::bit_vector& ColorCompressed<Color, Encoder>::uncompress(size_t j) {
    assert(j < color_encoder_->size());

    try {
        return *cached_colors_.Get(j);
    } catch (...) {
        sdsl::bit_vector *bit_vector = new sdsl::bit_vector(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j]) {
            // inflate vector

            auto slct = sdsl::select_support_sd<>(bitmatrix_[j]);
            auto rank = sdsl::rank_support_sd<>(bitmatrix_[j]);
            auto num_set_bits = rank(bitmatrix_[j]->size());

            for (uint64_t i = 1; i <= num_set_bits; ++i) {
                assert(slct(i) < bit_vector->size());

                (*bit_vector)[slct(i)] = 1;
            }
        }

        cached_colors_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template class ColorCompressed<std::string, StringEncoder>;

} // namespace annotate
