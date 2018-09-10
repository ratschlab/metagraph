#include "annotate_color_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"
#include "annotate_row_compressed.hpp"

using utils::remove_suffix;

size_t kNumRowsInBlock = 5'000'000;


namespace annotate {

template <typename Color, class Encoder>
ColorCompressed<Color, Encoder>::ColorCompressed(uint64_t num_rows,
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
        color_encoder_(new Encoder()),
        verbose_(verbose) {
    assert(num_columns_cached > 0);
}

template <typename Color, class Encoder>
ColorCompressed<Color, Encoder>::~ColorCompressed() {
    cached_colors_.Clear();
}

template <typename Color, class Encoder>
void
ColorCompressed<Color, Encoder>
::set_coloring(Index i, const Coloring &coloring) {
    assert(i < num_rows_);

    for (size_t j = 0; j < color_encoder_->size(); ++j) {
        decompress(j)[i] = 0;
    }
    for (const auto &color : coloring) {
        decompress(color_encoder_->encode(color, true))[i] = 1;
    }
}

template <typename Color, class Encoder>
typename ColorCompressed<Color, Encoder>::Coloring
ColorCompressed<Color, Encoder>::get_coloring(Index i) const {
    assert(i < num_rows_);

    flush();

    Coloring coloring;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j].get());

        if ((*bitmatrix_[j])[i])
            coloring.push_back(color_encoder_->decode(j));
    }
    return coloring;
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::add_color(Index i, const Color &color) {
    assert(i < num_rows_);

    decompress(color_encoder_->encode(color, true))[i] = 1;
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
    flush();

    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    serialize_number(outstream, num_rows_);

    color_encoder_->serialize(outstream);

    for (const auto &column : bitmatrix_) {
        assert(column.get());
        column->serialize(outstream);
    }
}


template <typename Color, class Encoder>
bool ColorCompressed<Color, Encoder>
::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_colors_.Clear();
    bitmatrix_.clear();

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
                num_rows_ = load_number(instream);
            } else if (num_rows_ != load_number(instream)) {
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
                    auto &column = decompress(col);

                    auto slct = sdsl::select_support_sd<>(new_column);
                    auto rank = sdsl::rank_support_sd<>(new_column);
                    auto num_set_bits = rank(new_column->size());

                    for (uint64_t i = 1; i <= num_set_bits; ++i) {
                        assert(slct(i) < column.size());

                        column[slct(i)] = 1;
                    }
                    delete new_column;
                } else {
                    bitmatrix_.at(col).reset(new_column);
                }
            }
        }

        if (verbose_) {
            std::cout << "Annotation loading finished ("
                      << bitmatrix_.size() << " columns)" << std::endl;
        }

        return true;

    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::insert_rows(const std::vector<Index> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));

    for (size_t j = 0; j < color_encoder_->size(); ++j) {
        auto &column = decompress(j);
        sdsl::bit_vector old = column;
        column = sdsl::bit_vector(old.size() + rows.size(), 0);

        uint64_t i = 0;
        uint64_t num_inserted = 0;

        for (auto next_inserted_row : rows) {
            while (i + num_inserted < next_inserted_row) {
                assert(i < old.size() && "Invalid indexes of inserted rows");
                assert(i + num_inserted < column.size() && "Invalid indexes of inserted rows");

                column[i + num_inserted] = old[i];
                i++;
            }
            // insert 0, not colored edge
            num_inserted++;
        }
        while (i < old.size()) {
            assert(i + num_inserted < column.size() && "Invalid indexes of inserted rows");

            column[i + num_inserted] = old[i];
            i++;
        }
    }
    num_rows_ += rows.size();
}

// For each pair (first, second) in the dictionary, renames
// column |first| with |second| and merges the columns with matching names.
template <typename Color, class Encoder>
void
ColorCompressed<Color, Encoder>
::rename_columns(const std::map<std::string, std::string> &dict) {
    cached_colors_.Clear();

    std::vector<std::string> old_index_to_label(color_encoder_->size());
    for (size_t i = 0; i < old_index_to_label.size(); ++i) {
        old_index_to_label[i] = color_encoder_->decode(i);
    }
    for (const auto &pair : dict) {
        old_index_to_label[color_encoder_->encode(pair.first)] = pair.second;
    }

    std::vector<std::string> index_to_label;
    std::unordered_map<std::string, std::set<size_t>> old_columns;
    for (size_t i = 0; i < old_index_to_label.size(); ++i) {
        if (!old_columns.count(old_index_to_label[i]))
            index_to_label.push_back(old_index_to_label[i]);

        old_columns[old_index_to_label[i]].insert(i);
    }

    std::vector<std::unique_ptr<sdsl::sd_vector<>>> old_bitmatrix;
    old_bitmatrix.swap(bitmatrix_);

    color_encoder_.reset(new Encoder());

    for (const auto &label : index_to_label) {
        color_encoder_->encode(label, true);

        const auto &cols = old_columns[label];

        assert(cols.size());
        if (cols.size() == 1) {
            bitmatrix_.emplace_back(std::move(old_bitmatrix[*cols.begin()]));
        } else {
            std::unique_ptr<sdsl::bit_vector> bit_vector(new sdsl::bit_vector(num_rows_, 0));
            for (size_t c : cols) {
                utils::decompress_sd_vector(*old_bitmatrix[c], bit_vector.get());
                old_bitmatrix[c].reset();
            }
            bitmatrix_.emplace_back(new sdsl::sd_vector<>(*bit_vector));
        }
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

    flush();

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
size_t ColorCompressed<Color, Encoder>::num_colors() const {
    return color_encoder_->size();
}

template <typename Color, class Encoder>
double ColorCompressed<Color, Encoder>::sparsity() const {
    uint64_t num_set_bits = 0;

    flush();

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        auto rank = sdsl::rank_support_sd<>(bitmatrix_[j].get());
        num_set_bits += rank(bitmatrix_[j]->size());
    }

    return 1 - static_cast<double>(num_set_bits) / num_colors() / num_rows_;
}

template <typename Color, class Encoder>
std::vector<uint64_t>
ColorCompressed<Color, Encoder>
::count_colors(const std::vector<Index> &indices) const {
    flush();

    std::vector<uint64_t> counter(bitmatrix_.size(), 0);

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        for (Index i : indices) {
            if ((*bitmatrix_[j])[i])
                counter[j]++;
        }
    }

    return counter;
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::flush() const {
    for (const auto &cached_vector : cached_colors_) {
        const_cast<ColorCompressed*>(this)->flush(
            cached_vector.first, cached_vector.second
        );
    }
    assert(bitmatrix_.size() == color_encoder_->size());
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::flush(size_t j, sdsl::bit_vector *vector) {
    assert(vector);
    assert(cached_colors_.Cached(j));

    while (bitmatrix_.size() <= j) {
        bitmatrix_.emplace_back();
    }

    bitmatrix_[j].reset();
    bitmatrix_[j].reset(new sdsl::sd_vector<>(*vector));
}

template <typename Color, class Encoder>
sdsl::bit_vector& ColorCompressed<Color, Encoder>::decompress(size_t j) {
    assert(j < color_encoder_->size());

    try {
        return *cached_colors_.Get(j);
    } catch (...) {
        sdsl::bit_vector *bit_vector = new sdsl::bit_vector(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j].get()) {
            utils::decompress_sd_vector(*bitmatrix_[j], bit_vector);
        }

        cached_colors_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::
convert_to_row_annotator(RowCompressed<Color, Encoder> *annotator,
                         size_t num_threads) const {
    flush();

    annotator->matrix_->reinitialize(num_rows_);
    annotator->color_encoder_.reset(new Encoder());
    for (size_t c = 0; c < color_encoder_->size(); ++c) {
        annotator->color_encoder_->encode(color_encoder_->decode(c), true);
    }

    std::vector<sdsl::select_support_sd<>> select_columns;
    std::vector<sdsl::rank_support_sd<>> rank_columns;

    for (const auto &v : bitmatrix_) {
        select_columns.emplace_back(v.get());
        rank_columns.emplace_back(v.get());
    }

    if (num_threads <= 1) {
        add_labels(&select_columns, &rank_columns, 0, num_rows_, annotator);
        return;
    }

    auto thread_pool = std::make_unique<utils::ThreadPool>(num_threads);
    for (uint64_t i = 0; i < num_rows_; i += kNumRowsInBlock) {
        thread_pool->enqueue(add_labels,
            &select_columns, &rank_columns,
            i, std::min(i + kNumRowsInBlock, num_rows_),
            annotator
        );
    }
    thread_pool.reset();
}

template <typename Color, class Encoder>
void ColorCompressed<Color, Encoder>::
add_labels(const std::vector<sdsl::select_support_sd<>> *select_columns,
           const std::vector<sdsl::rank_support_sd<>> *rank_columns,
           uint64_t begin, uint64_t end,
           RowCompressed<Color, Encoder> *annotator) {
    assert(begin <= end);
    assert(end <= annotator->matrix_->size());

    for (size_t j = 0; j < select_columns->size(); ++j) {
        // sdsl's rank computes the result excluding the query index
        uint64_t first = rank_columns->at(j)(begin) + 1;
        uint64_t last = rank_columns->at(j)(end);

        for (uint64_t r = first; r <= last; ++r) {
            annotator->matrix_->set_bit(select_columns->at(j)(r), j);
        }
    }
}

template class ColorCompressed<std::string, StringEncoder>;

} // namespace annotate
