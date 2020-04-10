#include "annotate_column_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "common/serialization.hpp"
#include "common/utils/string_utils.hpp"
#include "common/algorithms.hpp"
#include "common/vectors/bitmap_builder.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/threads/threading.hpp"

using utils::remove_suffix;

size_t kNumElementsReservedInBitmapBuilder = 10'000'000;


namespace annotate {

template <typename Label>
ColumnCompressed<Label>::ColumnCompressed(uint64_t num_rows,
                                          size_t num_columns_cached,
                                          bool verbose)
      : num_rows_(num_rows),
        cached_columns_(num_columns_cached,
                        caches::LRUCachePolicy<size_t>(),
                        [this](size_t j, bitmap_builder *column_builder) {
                            assert(column_builder);
                            this->flush(j, *column_builder);
                            delete column_builder;
                        }),
        verbose_(verbose) {
    assert(num_columns_cached > 0);
}

template <typename Label>
ColumnCompressed<Label>::~ColumnCompressed() {
    // Note: this is needed to make sure that everything is flushed to bitmatrix_
    //       BEFORE bitmatrix_ is destroyed
    cached_columns_.Clear();
}

template <typename Label>
void ColumnCompressed<Label>::set(Index i, const VLabels &labels) {
    assert(i < num_rows_);
    // add new labels
    for (const auto &label : labels) {
        label_encoder_.insert_and_encode(label);
    }
    // labels as a row
    std::vector<bool> row(label_encoder_.size(), 0);
    for (const auto &label : labels) {
        row[label_encoder_.encode(label)] = 1;
    }
    // set bits
    for (size_t j = 0; j < row.size(); ++j) {
        set(i, j, row[j]);
    }
}

template <typename Label>
void ColumnCompressed<Label>::add_labels(const std::vector<Index> &indices,
                                         const VLabels &labels) {
    for (const auto &label : labels) {
        const auto j = label_encoder_.insert_and_encode(label);
        decompress_builder(j).add_ones(indices.data(),
                                       indices.data() + indices.size());
    }
}

template <typename Label>
bool ColumnCompressed<Label>::has_label(Index i, const Label &label) const {
    try {
        return get_column(label)[i];
    } catch (...) {
        return false;
    }
}

template <typename Label>
bool ColumnCompressed<Label>::has_labels(Index i, const VLabels &labels) const {
    for (const auto &label : labels) {
        if (!has_label(i, label))
            return false;
    }
    return true;
}

template <typename Label>
void ColumnCompressed<Label>::serialize(const std::string &filename) const {
    flush();

    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension,
                            std::ios::binary);
    if (!outstream.good())
        throw std::ofstream::failure("Bad stream");

    serialize_number(outstream, num_rows_);

    label_encoder_.serialize(outstream);

    for (const auto &column : bitmatrix_) {
        assert(column.get());
        column->serialize(outstream);
    }
}

template <typename Label>
bool ColumnCompressed<Label>::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_columns_.Clear();
    bitmatrix_.clear();

    label_encoder_.clear();

    std::atomic<bool> error_occurred = false;

    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < filenames.size(); ++i) {
        if (error_occurred)
            continue;

        try {
            auto filename = remove_suffix(filenames[i], kExtension) + kExtension;

            if (verbose_) {
                std::cout << "Loading annotations from file " + filename + "\n" << std::flush;
            }

            std::ifstream instream(filename, std::ios::binary);
            if (!instream.good())
                throw std::ifstream::failure("can't open stream");

            const auto num_rows = load_number(instream);

            LabelEncoder<Label> label_encoder_load;
            if (!label_encoder_load.load(instream))
                throw std::ifstream::failure("can't load label encoder");

            if (!label_encoder_load.size())
                std::cerr << "No labels in " << filename << "\n" << std::flush;

            // update the existing and add some new columns
            for (size_t c = 0; c < label_encoder_load.size(); ++c) {
                std::unique_ptr<bit_vector> new_column { new bit_vector_smart() };

                auto pos = instream.tellg();

                if (!new_column->load(instream)) {
                    instream.seekg(pos, instream.beg);

                    new_column = std::make_unique<bit_vector_sd>();
                    if (!new_column->load(instream))
                        throw std::ifstream::failure("can't load next column");
                }

                if (new_column->size() != num_rows)
                    throw std::ifstream::failure("inconsistent column size");

                if (verbose_) {
                    auto num_set_bits = new_column->num_set_bits();

                    std::cout << "Column <" + label_encoder_load.decode(c) + ">"
                                + ", density: "
                                + std::to_string(static_cast<double>(num_set_bits) / new_column->size())
                                + ", set bits: " + std::to_string(num_set_bits) + "\n" << std::flush;
                }

                #pragma omp critical
                {
                    size_t col = label_encoder_.insert_and_encode(label_encoder_load.decode(c));

                    // set |num_rows_| with the first column inserted
                    if (!col)
                        num_rows_ = num_rows;

                    if (num_rows == num_rows_) {
                        assert(col <= bitmatrix_.size());

                        if (col == bitmatrix_.size()) {
                            bitmatrix_.emplace_back(std::move(new_column));
                        } else if (!bitmatrix_.at(col).get()) {
                            bitmatrix_.at(col) = std::move(new_column);
                        } else {
                            decompress_bitmap(col) |= *new_column;
                        }
                    } else {
                        error_occurred = true;
                    }
                }
            }
        } catch (...) {
            error_occurred = true;
        }
    }

    if (error_occurred)
        return false;

    if (verbose_) {
        std::cout << "Annotation loading finished ("
                  << bitmatrix_.size() << " columns)" << std::endl;
    }

    return true;
}

template <typename Label>
void ColumnCompressed<Label>::insert_rows(const std::vector<Index> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));
    for (size_t j = 0; j < label_encoder_.size(); ++j) {
        decompress_bitmap(j).insert_zeros(rows);
    }
    num_rows_ += rows.size();
}

template <typename Label>
void ColumnCompressed<Label>
::call_objects(const Label &label,
               std::function<void(Index)> callback) const {
    try {
        get_column(label).call_ones(callback);
    } catch (...) {
        return;
    }
}

template <typename Label>
std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
ColumnCompressed<Label>
::count_labels(const std::vector<std::pair<Index, size_t>> &index_counts,
               size_t min_count,
               size_t count_cap) const {

    assert(count_cap >= min_count);

    if (!count_cap)
        return {};

    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<std::pair<uint64_t, size_t>> label_counts;
    label_counts.reserve(num_labels());

    for (size_t j = 0; j < num_labels(); ++j) {
        size_t total_checked = 0;
        size_t total_matched = 0;

        const auto &column = get_column(j);

        for (auto [i, count] : index_counts) {
            total_checked += count;
            total_matched += count * column[i];

            if (total_matched >= count_cap
                    || total_matched + (total_sum_count - total_checked) < min_count)
                break;
        }

        if (total_matched >= min_count)
            label_counts.emplace_back(j, std::min(total_matched, count_cap));
    }

    return label_counts;
}

// For each pair (first, second) in the dictionary, renames
// column |first| with |second| and merges the columns with matching names.
template <typename Label>
void ColumnCompressed<Label>
::rename_labels(const tsl::hopscotch_map<Label, Label> &dict) {
    std::vector<Label> index_to_label(label_encoder_.size());
    // old labels
    for (size_t i = 0; i < index_to_label.size(); ++i) {
        index_to_label[i] = label_encoder_.decode(i);
    }
    // new labels
    for (const auto &pair : dict) {
        try {
            index_to_label[label_encoder_.encode(pair.first)] = pair.second;
        } catch (const std::out_of_range &) {
            std::cerr << "Warning: label '" << pair.first << "' not"
                      << " found in annotation. Skipping instruction"
                      << " '" << pair.first << " -> " << pair.second << "'."
                      << std::endl;
        }
    }

    std::vector<Label> new_index_to_label;
    tsl::hopscotch_map<Label, std::set<size_t>> old_columns;
    for (size_t i = 0; i < index_to_label.size(); ++i) {
        if (!old_columns.count(index_to_label[i]))
            new_index_to_label.push_back(index_to_label[i]);

        old_columns[index_to_label[i]].insert(i);
    }

    cached_columns_.Clear();
    std::vector<std::unique_ptr<bit_vector>> old_bitmatrix;
    old_bitmatrix.swap(bitmatrix_);

    label_encoder_.clear();

    for (const auto &label : new_index_to_label) {
        label_encoder_.insert_and_encode(label);

        const auto &cols = old_columns[label];

        assert(cols.size());
        if (cols.size() == 1) {
            bitmatrix_.emplace_back(std::move(old_bitmatrix[*cols.begin()]));
        } else {
            sdsl::bit_vector bit_vector(num_rows_, 0);
            for (size_t c : cols) {
                old_bitmatrix[c]->add_to(&bit_vector);
                old_bitmatrix[c].reset();
            }
            bitmatrix_.emplace_back(new bit_vector_smart(bit_vector));
        }
    }
}

template <typename Label>
uint64_t ColumnCompressed<Label>::num_objects() const {
    return num_rows_;
}

template <typename Label>
uint64_t ColumnCompressed<Label>::num_relations() const {
    uint64_t num_rels = 0;
    for (size_t i = 0; i < num_labels(); ++i) {
        num_rels += get_column(i).num_set_bits();
    }
    return num_rels;
}

template <typename Label>
void ColumnCompressed<Label>::set(Index i, size_t j, bool value) {
    assert(i < num_rows_);

    // We update the value if:
    //  * we are inserting a new column
    //  * found uncompressed, hence it's easier to blindly update
    //  * the column exists but compressed -- update only if value differs
    if (j >= bitmatrix_.size()
            || cached_columns_.Cached(j)
            || (*bitmatrix_[j])[i] != value) {  // only compressed -- check value
        decompress_bitmap(j).set(i, value);
    }
}

template <typename Label>
void ColumnCompressed<Label>::flush() const {
    std::lock_guard<std::mutex> lock(bitmap_conversion_mu_);

    if (!flushed_) {
        for (const auto &cached_vector : cached_columns_) {
            const_cast<ColumnCompressed*>(this)->flush(
                cached_vector.first, *cached_vector.second
            );
        }
        flushed_ = true;
    }
    assert(bitmatrix_.size() == label_encoder_.size());
}

template <typename Label>
void ColumnCompressed<Label>::flush(size_t j, const bitmap_builder &builder) {
    assert(j < bitmatrix_.size());

    // Note: asserting that j is cached cannot be done here when this function
    //       is invovled as part of the OnEraseCallback, since the mutex locking
    //       in the caches library would cause the check to be done after it has
    //       been erased.

    bitmatrix_[j].reset();

    auto initialization_data = builder.get_initialization_data();
    bitmatrix_[j].reset(new bit_vector_smart(initialization_data.call_ones,
                                             initialization_data.size,
                                             initialization_data.num_set_bits));

    assert(initialization_data.size == bitmatrix_[j]->size());
    assert(initialization_data.num_set_bits == bitmatrix_[j]->num_set_bits());
    assert(bitmatrix_[j]->num_set_bits() == bitmatrix_[j]->rank1(bitmatrix_[j]->size() - 1));
}

/**
 * Get a reference to an uncompressed column in cache.
 * If not cached, decompress the compressed one and cache it.
 * The uncompressed column can be used for any modifications.
 */
template <typename Label>
bitmap_dyn& ColumnCompressed<Label>::decompress_bitmap(size_t j) {
    // get the bitmap builder (or bitmap) from cache
    auto &builder = decompress_builder(j);
    flushed_ = false;

    assert(j < bitmatrix_.size());

    if (bitmap_dyn *uncompressed = dynamic_cast<bitmap_dyn*>(&builder))
        return *uncompressed;

    // if the column is new and we have only its builder, build the column
    auto initialization_data = builder.get_initialization_data();

    sdsl::bit_vector column_data(num_rows_, 0);
    initialization_data.call_ones([&](uint64_t i) { column_data[i] = 1; });

    auto *vector = new bitmap_vector(std::move(column_data));
    cached_columns_.Put(j, vector);

    return *vector;
}

/**
 * Get a column builder (or uncompressed column) for inserting new set bits.
 *
 *  - If a new column is being inserted, initialize its fast builder.
 *  - If the column (or its builder) is cached, return a reference to it.
 *  - Otherwise, if the column exists but not cached, decompress it
 *    and add to cache. Then, return a reference to it.
 */
template <typename Label>
bitmap_builder& ColumnCompressed<Label>::decompress_builder(size_t j) {
    assert(j < label_encoder_.size());

    flushed_ = false;

    try {
        // check the  the cached bitmap builder
        return *cached_columns_.Get(j);
    } catch (...) {
        assert(j <= bitmatrix_.size());

        bitmap_builder *vector;

        if (j == bitmatrix_.size()) {
            // the column is new, create an efficient builder for it
            bitmatrix_.emplace_back();
            // Work with the full uncompressed bitmap if it takes less space
            // than the buffer in its builder.
            if (num_rows_ < kNumElementsReservedInBitmapBuilder * 64) {
                vector = new bitmap_vector(num_rows_, 0);
            } else {
                // For large bitmaps, use the efficient builder, using only
                // space proportional to the number of bits set in the bitmap.
                vector = new bitmap_builder_set(num_rows_, get_num_threads(),
                                                kNumElementsReservedInBitmapBuilder);
            }
        } else {
            // otherwise, decompress the existing column and initialize a bitmap
            vector = new bitmap_vector(bitmatrix_[j]->template convert_to<sdsl::bit_vector>());
            bitmatrix_[j].reset();
        }

        cached_columns_.Put(j, vector);
        return *cached_columns_.Get(j);
    }
}

template <typename Label>
const bitmap& ColumnCompressed<Label>::get_column(size_t j) const {
    if (!cached_columns_.Cached(j)) {
        assert(j < bitmatrix_.size() && bitmatrix_[j].get());
        return (*bitmatrix_[j]);
    }

    // lock the mutex in case the bitmap conversion happens
    std::lock_guard<std::mutex> lock(bitmap_conversion_mu_);

    const auto &builder = *cached_columns_.Get(j);

    if (const bitmap *uncompressed = dynamic_cast<const bitmap*>(&builder))
        return *uncompressed;

    // if the column is new and we have only its builder, build the column
    auto initialization_data = builder.get_initialization_data();

    sdsl::bit_vector column_data(num_rows_, 0);
    initialization_data.call_ones([&](uint64_t i) { column_data[i] = 1; });

    auto *vector = new bitmap_vector(std::move(column_data));
    const_cast<ColumnCompressed*>(this)->cached_columns_.Put(j, vector);

    return *vector;
}

template <typename Label>
const bitmap& ColumnCompressed<Label>::get_column(const Label &label) const {
    return get_column(label_encoder_.encode(label));
}

template <typename Label>
const ColumnMajor& ColumnCompressed<Label>::get_matrix() const {
    flush();
    return annotation_matrix_view_;
}

template <typename Label>
bool ColumnCompressed<Label>
::dump_columns(const std::string &prefix, size_t num_threads) const {
    bool success = true;

    #pragma omp parallel for num_threads(num_threads)
    for (uint64_t j = 0; j < num_labels(); ++j) {
        std::ofstream outstream(
            remove_suffix(prefix, kExtension)
                + "." + std::to_string(j)
                + ".text.annodbg"
        );

        if (!outstream.good()) {
            std::cerr << "ERROR: dumping column " << j << " failed" << std::endl;
            success = false;
            continue;
        }

        const auto &column = get_column(j);

        outstream << num_objects() << " " << column.num_set_bits() << "\n";

        column.call_ones([&](const auto &pos) {
            outstream << pos << "\n";
        });
    }

    return success;
}

template class ColumnCompressed<std::string>;

} // namespace annotate
