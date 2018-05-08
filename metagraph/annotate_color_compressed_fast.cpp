#include "annotate_color_compressed_fast.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"

using libmaus2::util::NumberSerialisation;
using utils::remove_suffix;


namespace annotate {

template <typename Color, class Encoder>
const std::string FastColorCompressed<Color, Encoder>::kExtension = ".color.annodbg";

template <typename Color, class Encoder>
FastColorCompressed<Color, Encoder>::FastColorCompressed(uint64_t num_rows,
                                                         size_t num_columns_cached,
                                                         bool verbose)
      : num_rows_(num_rows),
        cached_colors_(
            num_columns_cached,
            caches::LRUCachePolicy<size_t>(),
            [this](size_t j, sdsl::bit_vector *col_uncompressed) {
                this->flush(j, col_uncompressed);
                delete col_uncompressed;
                index_.clear();
            }
        ),
        color_encoder_(new Encoder()),
        verbose_(verbose) {}

template <typename Color, class Encoder>
FastColorCompressed<Color, Encoder>
::FastColorCompressed(ColorCompressed<Color, Encoder>&& annotator,
                      size_t num_columns_cached,
                      bool verbose)
      : num_rows_(annotator.num_rows_),
        cached_colors_(
            num_columns_cached,
            caches::LRUCachePolicy<size_t>(),
            [this](size_t j, sdsl::bit_vector *col_uncompressed) {
                this->flush(j, col_uncompressed);
                delete col_uncompressed;
                index_.clear();
            }
        ),
        verbose_(verbose) {
    annotator.cached_colors_.Clear();
    annotator.flush();

    bitmatrix_ = std::move(annotator.bitmatrix_);

    color_encoder_ = std::move(annotator.color_encoder_);
    annotator.color_encoder_.reset(new Encoder());

    update_index();
}

template <typename Color, class Encoder>
FastColorCompressed<Color, Encoder>::~FastColorCompressed() {
    cached_colors_.Clear();

    for (auto *column : bitmatrix_) {
        assert(column);
        delete column;
    }
}

template <typename Color, class Encoder>
void
FastColorCompressed<Color, Encoder>
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
typename FastColorCompressed<Color, Encoder>::Coloring
FastColorCompressed<Color, Encoder>::get_coloring(Index i) const {
    assert(i < num_rows_);

    Coloring coloring;
    for (size_t j : get_row(i)) {
        coloring.push_back(color_encoder_->decode(j));
    }
    return coloring;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::add_color(Index i, const Color &color) {
    assert(i < num_rows_);

    uncompress(color_encoder_->encode(color, true))[i] = 1;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::add_colors(Index i, const Coloring &coloring) {
    for (const auto &color : coloring) {
        add_color(i, color);
    }
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::add_colors(const std::vector<Index> &indices,
                                                 const Coloring &coloring) {
    for (const auto &color : coloring) {
        for (Index i : indices) {
            add_color(i, color);
        }
    }
}

template <typename Color, class Encoder>
bool
FastColorCompressed<Color, Encoder>
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
FastColorCompressed<Color, Encoder>
::has_colors(Index i, const Coloring &coloring) const {
    for (const auto &color : coloring) {
        if (!has_color(i, color))
            return false;
    }
    return true;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::serialize(const std::string &filename) const {
    const_cast<FastColorCompressed*>(this)->flush();

    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    NumberSerialisation::serialiseNumber(outstream, num_rows_);

    color_encoder_->serialize(outstream);

    for (auto *column : bitmatrix_) {
        assert(column);
        column->serialize(outstream);
    }

    if (!index_.size())
        const_cast<FastColorCompressed*>(this)->update_index();

    NumberSerialisation::serialiseNumber(outstream, index_.size());
    for (auto &aux_column : index_) {
        aux_column.serialize(outstream);
    }
}


template <typename Color, class Encoder>
bool FastColorCompressed<Color, Encoder>
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
            if (verbose_) {
                std::cout << "Loading annotations from file "
                          << remove_suffix(filename, kExtension) + kExtension
                          << "..." << std::endl;
            }

            std::ifstream instream(remove_suffix(filename, kExtension) + kExtension);
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

                if (verbose_) {
                    auto rank = sdsl::rank_support_sd<>(new_column);
                    auto num_set_bits = rank(new_column->size());

                    std::cout << "Color <" << color_encoder_load.decode(c)
                                           << ">, "
                              << "density: "
                              << static_cast<double>(num_set_bits) / new_column->size()
                              << ", set bits: " << num_set_bits << std::endl;
                }

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

            if (filenames.size() == 1) {
                try {
                    index_.resize(NumberSerialisation::deserialiseNumber(instream));

                    for (auto &aux_column : index_) {
                        aux_column.load(instream);
                    }
                } catch (...) {}
            }
        }

        if (verbose_) {
            std::cout << "Annotation loading finished ("
                      << bitmatrix_.size() << " columns)" << std::endl;
        }

        if (!index_.size())
            update_index();

        if (verbose_) {
            for (size_t t = 0; t < index_.size(); ++t) {
                auto rank = sdsl::rank_support_sd<>(&index_[t]);
                auto num_set_bits = rank(index_[t].size());

                std::cout << "Aux column " << t << ", "
                          << "density: "
                          << static_cast<double>(num_set_bits) / index_[t].size()
                          << ", set bits: " << num_set_bits << std::endl;
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
typename FastColorCompressed<Color, Encoder>::Coloring
FastColorCompressed<Color, Encoder>
::aggregate_colors(const std::vector<Index> &indices,
                   double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    const size_t max_colors_missing = indices.size() - min_colors_discovered;

    Coloring filtered_colors;

    std::vector<bool> color_filter(bitmatrix_.size(), true);
    std::vector<size_t> discovered(bitmatrix_.size(), 0);
    uint64_t iteration = 0;

    for (Index i : indices) {
        for (size_t j : filter_row(i, color_filter)) {
            assert((*bitmatrix_[j])[i]);

            if (++discovered[j] >= min_colors_discovered) {
                filtered_colors.push_back(color_encoder_->decode(j));
                color_filter[j] = false;
            }
        }

        if (++iteration <= max_colors_missing)
            continue;

        for (size_t j = 0; j < color_filter.size(); ++j) {
            if (color_filter[j]
                    && iteration - discovered[j] > max_colors_missing) {
                color_filter[j] = false;
            }
        }
    }

    return filtered_colors;
}

// Count all colors collected from extracted colorings
// and return top |num_top| with the counts computed.
template <typename Color, class Encoder>
std::vector<std::pair<Color, size_t>>
FastColorCompressed<Color, Encoder>
::get_most_frequent_colors(const std::vector<Index> &indices,
                           size_t num_top) const {
    auto counter = count_colors(indices);

    std::vector<std::pair<size_t, size_t>> counts;
    for (size_t j = 0; j < counter.size(); ++j) {
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
size_t FastColorCompressed<Color, Encoder>::num_colors() const {
    return color_encoder_->size();
}

template <typename Color, class Encoder>
double FastColorCompressed<Color, Encoder>::sparsity() const {
    uint64_t num_set_bits = 0;

    const_cast<FastColorCompressed*>(this)->flush();

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        auto rank = sdsl::rank_support_sd<>(bitmatrix_[j]);
        num_set_bits += rank(bitmatrix_[j]->size());
    }

    return 1 - static_cast<double>(num_set_bits) / num_colors() / num_rows_;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::update_index(size_t num_aux_cols) {
    const_cast<FastColorCompressed*>(this)->flush();

    index_.clear();

    if (num_aux_cols == static_cast<size_t>(-1)) {
        // In the worst case, on average, we access
        // |num_aux_cols| auxiliary columns, plus
        // (|density| * |num_columns|) * (|num_columns| / |num_aux_cols|)
        // real columns per query.
        num_aux_cols = std::sqrt((1 - sparsity()) * bitmatrix_.size()
                                                  * bitmatrix_.size());
    }

    if (verbose_)
        std::cout << "Updating " << num_aux_cols << " index columns" << std::endl;

    size_t block_width = (bitmatrix_.size() + num_aux_cols - 1) / num_aux_cols;

    if (verbose_)
        std::cout << "Processed index columns..." << std::endl;

    for (size_t t = 0; t < num_aux_cols; ++t) {
        sdsl::bit_vector aux_col(num_rows_, 0);

        if (verbose_)
            std::cout << "Block " << t << ": " << std::flush;

        for (size_t j = t * block_width; j < std::min((t + 1) * block_width,
                                                      bitmatrix_.size()); ++j) {

            if (verbose_)
                std::cout << j << ".." << std::flush;

            sdsl::select_support_sd<> slct(bitmatrix_[j]);
            sdsl::rank_support_sd<> rank(bitmatrix_[j]);
            uint64_t num_set_bits = rank(bitmatrix_[j]->size());

            for (uint64_t i = 1; i <= num_set_bits; ++i) {
                assert(slct(i) < aux_col.size());

                aux_col[slct(i)] = 1;
            }
        }

        if (verbose_)
            std::cout << std::endl;

        index_.emplace_back(aux_col);
    }
}

template <typename Color, class Encoder>
std::vector<size_t> FastColorCompressed<Color, Encoder>::get_row(Index i) const {
    const_cast<FastColorCompressed*>(this)->flush();

    std::vector<size_t> set_bits;
    set_bits.reserve(bitmatrix_.size());

    if (!index_.size()) {
        // query all columns if auxiliary index is not initialized
        for (size_t j = 0; j < bitmatrix_.size(); ++j) {
            if ((*bitmatrix_[j])[i])
                set_bits.push_back(j);
        }

        return set_bits;
    }

    // use the initialized auxiliary index columns
    size_t block_width = (bitmatrix_.size() + index_.size() - 1) / index_.size();

    for (size_t t = 0; t < index_.size(); ++t) {
        if (!index_[t][i])
            continue;

        for (size_t j = t * block_width;
                    j < (t + 1) * block_width && j < bitmatrix_.size(); ++j) {
            if ((*bitmatrix_[j])[i])
                set_bits.push_back(j);
        }
    }

    return set_bits;
}

template <typename Color, class Encoder>
std::vector<size_t>
FastColorCompressed<Color, Encoder>
::filter_row(Index i, const std::vector<bool> &filter) const {
    const_cast<FastColorCompressed*>(this)->flush();

    assert(filter.size() == bitmatrix_.size());

    std::vector<size_t> set_bits;
    set_bits.reserve(bitmatrix_.size());

    if (!index_.size()) {
        // query all columns if auxiliary index is not initialized
        for (size_t j = 0; j < bitmatrix_.size(); ++j) {
            if (filter[j] && (*bitmatrix_[j])[i])
                set_bits.push_back(j);
        }

        return set_bits;
    }

    // use the initialized auxiliary index columns
    size_t block_width = (bitmatrix_.size() + index_.size() - 1) / index_.size();

    std::vector<size_t> bits_to_check;
    bits_to_check.reserve(bitmatrix_.size());

    for (size_t t = 0; t < index_.size(); ++t) {
        bits_to_check.clear();

        size_t block_begin = t * block_width;
        size_t block_end = std::min((t + 1) * block_width, bitmatrix_.size());

        for (size_t j = block_begin; j < block_end; ++j) {
            if (filter[j])
                bits_to_check.push_back(j);
        }

        if (bits_to_check.size() > 1 && !index_[t][i])
            continue;

        for (size_t j : bits_to_check) {
            if ((*bitmatrix_[j])[i])
                set_bits.push_back(j);
        }
    }

    return set_bits;
}

template <typename Color, class Encoder>
std::vector<uint64_t>
FastColorCompressed<Color, Encoder>
::count_colors(const std::vector<Index> &indices) const {
    std::vector<uint64_t> counter(bitmatrix_.size(), 0);

    for (Index i : indices) {
        for (size_t j : get_row(i)) {
            counter[j]++;
        }
    }

    return counter;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::flush() {
    for (const auto &cached_vector : cached_colors_) {
        flush(cached_vector.first, cached_vector.second);
    }
    assert(bitmatrix_.size() == color_encoder_->size());
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::flush(size_t j, sdsl::bit_vector *vector) {
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
sdsl::bit_vector& FastColorCompressed<Color, Encoder>::uncompress(size_t j) {
    assert(j < color_encoder_->size());

    index_.clear();

    try {
        return *cached_colors_.Get(j);
    } catch (...) {
        sdsl::bit_vector *bit_vector = new sdsl::bit_vector(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j]) {
            // inflate vector

            sdsl::select_support_sd<> slct(bitmatrix_[j]);
            sdsl::rank_support_sd<> rank(bitmatrix_[j]);
            uint64_t num_set_bits = rank(bitmatrix_[j]->size());

            for (uint64_t i = 1; i <= num_set_bits; ++i) {
                assert(slct(i) < bit_vector->size());

                (*bit_vector)[slct(i)] = 1;
            }
        }

        cached_colors_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template class FastColorCompressed<std::string, StringEncoder>;

} // namespace annotate
