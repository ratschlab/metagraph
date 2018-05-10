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
const std::string FastColorCompressed<Color, Encoder>::kIndexExtension = ".colindex.annodbg";

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
            }
        ),
        cached_index_(
            std::max(static_cast<size_t>(1), num_columns_cached / 2),
            caches::LRUCachePolicy<size_t>(),
            [this](size_t t, sdsl::bit_vector *col_uncompressed) {
                this->flush_index(t, col_uncompressed);
                delete col_uncompressed;
            }
        ),
        color_encoder_(new Encoder()),
        verbose_(verbose) {}

template <typename Color, class Encoder>
FastColorCompressed<Color, Encoder>
::FastColorCompressed(ColorCompressed<Color, Encoder>&& annotator,
                      size_t num_columns_cached,
                      bool verbose,
                      bool build_index)
      : num_rows_(annotator.num_rows_),
        cached_colors_(
            num_columns_cached,
            caches::LRUCachePolicy<size_t>(),
            [this](size_t j, sdsl::bit_vector *col_uncompressed) {
                this->flush(j, col_uncompressed);
                delete col_uncompressed;
            }
        ),
        cached_index_(
            std::max(static_cast<size_t>(1), num_columns_cached / 2),
            caches::LRUCachePolicy<size_t>(),
            [this](size_t t, sdsl::bit_vector *col_uncompressed) {
                this->flush_index(t, col_uncompressed);
                delete col_uncompressed;
            }
        ),
        verbose_(verbose) {
    annotator.cached_colors_.Clear();

    bitmatrix_ = std::move(annotator.bitmatrix_);

    color_encoder_ = std::move(annotator.color_encoder_);
    annotator.color_encoder_.reset(new Encoder());

    if (build_index)
        rebuild_index();
}

template <typename Color, class Encoder>
FastColorCompressed<Color, Encoder>::~FastColorCompressed() {
    cached_colors_.Clear();
    cached_index_.Clear();
}

template <typename Color, class Encoder>
void
FastColorCompressed<Color, Encoder>
::set_coloring(Index i, const Coloring &coloring) {
    assert(i < num_rows_);

    // add new colors
    for (const auto &color : coloring) {
        color_encoder_->encode(color, true);
    }

    // coloring as a row
    std::vector<bool> row(num_colors(), 0);
    for (const auto &color : coloring) {
        row[color_encoder_->encode(color, true)] = 1;
    }

    // reset coloring
    for (size_t j = 0; j < row.size(); ++j) {
        if (get_entry(i, j) != row[j])
            decompress(j)[i] = row[j];
    }

    update_index();
    for (size_t t = 0; t < index_to_columns_.size(); ++t) {
        bool non_empty_block = false;

        for (size_t j : index_to_columns_[t]) {
            if (row[j])
                non_empty_block = true;
        }

        if (get_index_entry(i, t) != non_empty_block)
            decompress_index(t)[i] = non_empty_block;
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

    size_t j = color_encoder_->encode(color, true);
    if (get_entry(i, j))
        return;

    decompress(j)[i] = 1;

    update_index();
    if (index_to_columns_.size()) {
        size_t t = column_to_index_[j];
        if (!get_index_entry(i, t))
            decompress_index(t)[i] = 1;
    }
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
        return get_entry(i, color_encoder_->encode(color));
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

    for (const auto &column : bitmatrix_) {
        column->serialize(outstream);
    }

    // serialize the column index
    outstream = std::ofstream(remove_suffix(filename, kExtension) + kIndexExtension);

    NumberSerialisation::serialiseNumber(outstream, index_.size());
    for (size_t t = 0; t < index_.size(); ++t) {
        index_[t]->serialize(outstream);
        serialize_number_vector(outstream, index_to_columns_[t]);
    }
}


template <typename Color, class Encoder>
bool FastColorCompressed<Color, Encoder>
::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_colors_.Clear();
    cached_index_.Clear();

    bitmatrix_.clear();
    index_.clear();
    column_to_index_.clear();
    index_to_columns_.clear();

    color_encoder_.reset(new Encoder());

    to_update_ = false;
    to_update_index_ = false;

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
            bitmatrix_.resize(color_encoder_->size());
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

                if (bitmatrix_.at(col).get()) {
                    utils::decompress_sd_vector(*new_column, &decompress(col));
                    delete new_column;
                } else {
                    bitmatrix_[col].reset(new_column);
                }
            }
        }

        if (verbose_) {
            std::cout << "Annotation loading finished ("
                      << bitmatrix_.size() << " columns)" << std::endl;
        }

        if (filenames.size() == 1) {
            if (verbose_) {
                std::cout << "Loading index from file "
                          << remove_suffix(filenames[0], kExtension) + kIndexExtension
                          << "..." << std::endl;
            }

            try {
                std::ifstream instream(remove_suffix(filenames[0], kExtension) + kIndexExtension);

                size_t index_size = NumberSerialisation::deserialiseNumber(instream);

                column_to_index_.resize(num_colors());

                for (size_t t = 0; t < index_size; ++t) {
                    index_.emplace_back(new sdsl::sd_vector<>());
                    index_.back()->load(instream);

                    index_to_columns_.emplace_back(load_number_vector<size_t>(instream));

                    for (size_t j : index_to_columns_.back()) {
                        column_to_index_[j] = t;
                    }
                }
            } catch (...) {
                std::cerr << "Warning: Can't load index from disk."
                          << " Initializing an index using the default strategy..." << std::endl;
                rebuild_index();
                std::cerr << "The default index has been initialized." << std::endl;
            }
        } else {
            std::cerr << "Warning: Annotation was loaded from multiple files."
                      << " Rebuild index using the default strategy..." << std::endl;
            rebuild_index();
            std::cerr << "The default index has been initialized." << std::endl;
        }

        if (verbose_) {
            std::cout << "Column index loading finished ("
                      << index_.size() << " index columns)" << std::endl;

            for (size_t t = 0; t < index_.size(); ++t) {
                auto rank = sdsl::rank_support_sd<>(index_[t].get());
                auto num_set_bits = rank(index_[t]->size());

                std::cout << "Index column " << t << ", "
                          << "density: "
                          << static_cast<double>(num_set_bits) / index_[t]->size()
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

    std::vector<bool> color_filter(num_colors(), true);
    std::vector<size_t> discovered(num_colors(), 0);
    uint64_t iteration = 0;

    for (Index i : indices) {
        for (size_t j : filter_row(i, color_filter)) {
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
        auto rank = sdsl::rank_support_sd<>(bitmatrix_[j].get());
        num_set_bits += rank(bitmatrix_[j]->size());
    }

    double density = static_cast<double>(num_set_bits) / num_colors() / num_rows_;
    return 1 - density;
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::rebuild_index(size_t num_aux_cols) {
    const_cast<FastColorCompressed*>(this)->flush();

    index_.clear();
    column_to_index_.clear();
    index_to_columns_.clear();

    if (num_aux_cols == static_cast<size_t>(-1)) {
        // In the worst case, on average, we access
        // |num_aux_cols| auxiliary columns, plus
        // (|density| * |num_columns|) * (|num_columns| / |num_aux_cols|)
        // real columns per query.
        num_aux_cols = std::sqrt((1 - sparsity()) * num_colors() * num_colors());
    }

    if (!num_aux_cols)
        return;

    index_.resize(num_aux_cols);
    column_to_index_.resize(bitmatrix_.size());
    index_to_columns_.resize(num_aux_cols);

    if (verbose_)
        std::cout << "Updating " << num_aux_cols << " index columns" << std::endl;

    size_t block_width = (num_colors() + num_aux_cols - 1) / num_aux_cols;

    if (verbose_)
        std::cout << "Processed index columns..." << std::endl;

    for (size_t t = 0; t < num_aux_cols; ++t) {
        sdsl::bit_vector &index_col = decompress_index(t);

        if (verbose_)
            std::cout << "Block " << t << ": " << std::flush;

        for (size_t j = t * block_width; j < std::min((t + 1) * block_width,
                                                      num_colors()); ++j) {
            column_to_index_[j] = t;
            index_to_columns_[t].push_back(j);

            if (verbose_)
                std::cout << j << ".." << std::flush;

            sdsl::select_support_sd<> slct(bitmatrix_[j].get());
            sdsl::rank_support_sd<> rank(bitmatrix_[j].get());
            uint64_t num_set_bits = rank(bitmatrix_[j]->size());

            for (uint64_t i = 1; i <= num_set_bits; ++i) {
                index_col[slct(i)] = 1;
            }
        }

        if (verbose_)
            std::cout << std::endl;
    }
}

template <typename Color, class Encoder>
std::vector<size_t> FastColorCompressed<Color, Encoder>::get_row(Index i) const {
    const_cast<FastColorCompressed*>(this)->flush();

    std::vector<size_t> set_bits;
    set_bits.reserve(bitmatrix_.size());

    if (!index_to_columns_.size()) {
        // query all columns if auxiliary index is not initialized
        for (size_t j = 0; j < bitmatrix_.size(); ++j) {
            if ((*bitmatrix_[j])[i])
                set_bits.push_back(j);
        }

        return set_bits;
    }

    // use the initialized auxiliary index columns
    for (size_t t = 0; t < index_to_columns_.size(); ++t) {
        if (!index_to_columns_[t].size() || !(*index_[t])[i])
            continue;

        for (size_t j : index_to_columns_[t]) {
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
    std::vector<size_t> bits_to_check;
    bits_to_check.reserve(bitmatrix_.size());

    for (size_t t = 0; t < index_.size(); ++t) {
        bits_to_check.clear();

        for (size_t j : index_to_columns_[t]) {
            if (filter[j])
                bits_to_check.push_back(j);
        }

        if (bits_to_check.size() > 1 && !(*index_[t])[i])
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
bool FastColorCompressed<Color, Encoder>::get_entry(Index i, size_t j) const {
    if (cached_colors_.Cached(j)) {
        return (*cached_colors_.Get(j))[i];
    } else if (j < bitmatrix_.size()) {
        assert(bitmatrix_[j].get());
        return (*bitmatrix_[j])[i];
    } else {
        return false;
    }
}

template <typename Color, class Encoder>
bool FastColorCompressed<Color, Encoder>::get_index_entry(Index i, size_t t) const {
    if (cached_index_.Cached(t)) {
        return (*cached_index_.Get(t))[i];
    } else {
        assert(t < index_.size() && index_[t].get());
        return (*index_[t])[i];
    }
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::flush() {
    if (to_update_) {
        for (const auto &cached_vector : cached_colors_) {
            flush(cached_vector.first, cached_vector.second);
        }
        to_update_ = false;
    }

    if (to_update_index_) {
        for (const auto &cached_index : cached_index_) {
            flush_index(cached_index.first, cached_index.second);
        }
        to_update_index_ = false;
    }
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::flush(size_t j, sdsl::bit_vector *vector) {
    assert(vector);
    assert(cached_colors_.Cached(j));

    if (!to_update_)
        return;

    while (bitmatrix_.size() <= j) {
        bitmatrix_.emplace_back();
    }

    bitmatrix_[j].reset();
    bitmatrix_[j].reset(new sdsl::sd_vector<>(*vector));
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::flush_index(size_t t, sdsl::bit_vector *vector) {
    assert(vector);
    assert(cached_index_.Cached(t));
    assert(t < index_.size());

    if (!to_update_index_)
        return;

    index_[t].reset();
    index_[t].reset(new sdsl::sd_vector<>(*vector));
}

template <typename Color, class Encoder>
sdsl::bit_vector& FastColorCompressed<Color, Encoder>::decompress(size_t j) {
    assert(j < num_colors());

    to_update_ = true;

    try {
        return *cached_colors_.Get(j);
    } catch (...) {
        sdsl::bit_vector *bit_vector = new sdsl::bit_vector(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j].get()) {
            utils::decompress_sd_vector(*bitmatrix_[j], bit_vector);
            bitmatrix_[j].reset();
        }

        cached_colors_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template <typename Color, class Encoder>
sdsl::bit_vector& FastColorCompressed<Color, Encoder>::decompress_index(size_t t) {
    assert(t < index_.size());

    to_update_index_ = true;

    try {
        return *cached_index_.Get(t);
    } catch (...) {
        sdsl::bit_vector *bit_vector = new sdsl::bit_vector(num_rows_, 0);

        if (t < index_.size() && index_[t].get()) {
            utils::decompress_sd_vector(*index_[t], bit_vector);
            index_[t].reset();
        }

        cached_index_.Put(t, bit_vector);
        return *bit_vector;
    }
}

template <typename Color, class Encoder>
void FastColorCompressed<Color, Encoder>::update_index() {
    if (!index_to_columns_.size())
        return;

    for (size_t j = column_to_index_.size(); j < num_colors(); ++j) {
        // distribute all new columns uniformly
        size_t t = j % index_to_columns_.size();
        column_to_index_.push_back(t);
        index_to_columns_[t].push_back(j);
    }
}

template class FastColorCompressed<std::string, StringEncoder>;

} // namespace annotate
