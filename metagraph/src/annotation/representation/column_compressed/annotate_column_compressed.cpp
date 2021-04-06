#include "annotate_column_compressed.hpp"

#include <string>
#include <numeric>
#include <algorithm>
#include <stdexcept>

#include "common/serialization.hpp"
#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "common/vectors/bitmap_builder.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/threads/threading.hpp"


namespace mtg {
namespace annot {

using utils::remove_suffix;
using mtg::common::logger;

const uint8_t kCountBits = 8;
const uint32_t kMaxCount = sdsl::bits::lo_set[kCountBits];


template <typename Label>
ColumnCompressed<Label>::ColumnCompressed(uint64_t num_rows,
                                          size_t num_columns_cached,
                                          const std::string &swap_dir,
                                          uint64_t buffer_size_bytes)
      : num_rows_(num_rows),
        swap_dir_(swap_dir),
        buffer_size_bytes_(buffer_size_bytes),
        cached_columns_(std::max(num_columns_cached, (size_t)1),
                        caches::LRUCachePolicy<size_t>(),
                        [this](size_t j, bitmap_builder *column_builder) {
                            assert(column_builder);
                            this->flush(j, column_builder);
                            delete column_builder;
                        }) {}

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

// for each label and index 'indices[i]' add count 'counts[i]'
template <typename Label>
void ColumnCompressed<Label>::add_label_counts(const std::vector<Index> &indices,
                                               const VLabels &labels,
                                               const std::vector<uint32_t> &counts) {
    assert(indices.size() == counts.size());

    const auto &columns = get_matrix().data();

    if (relation_counts_.size() != columns.size())
        relation_counts_.resize(columns.size());

    for (const auto &label : labels) {
        const auto j = label_encoder_.insert_and_encode(label);

        if (!relation_counts_[j].size()) {
            relation_counts_[j] = sdsl::int_vector<>(columns[j]->num_set_bits(), 0, kCountBits);

        } else if (relation_counts_[j].size() != columns[j]->num_set_bits()) {
            logger->error("Binary relation matrix was changed while adding relation counts");
            exit(1);
        }

        for (size_t i = 0; i < indices.size(); ++i) {
            if (uint64_t rank = columns[j]->conditional_rank1(indices[i])) {
                uint32_t count = std::min(counts[i], kMaxCount);
                sdsl::int_vector_reference<sdsl::int_vector<>> ref = relation_counts_[j][rank - 1];
                ref = std::min((uint32_t)ref, kMaxCount - count) + count;

            } else {
                logger->warn("Trying to add count {} for non-annotated object {}."
                             " The count was ignored.", counts[i], indices[i]);
            }
        }
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
    if (!outstream.good()) {
        logger->trace("Could not open {}", remove_suffix(filename, kExtension) + kExtension);
        throw std::ofstream::failure("Bad stream");
    }

    serialize_number(outstream, num_rows_);

    label_encoder_.serialize(outstream);

    for (const auto &column : bitmatrix_) {
        assert(column.get());
        column->serialize(outstream);
    }

    if (!relation_counts_.size())
        return;

    outstream.close();

    outstream.open(remove_suffix(filename, kExtension) + kExtension + ".counts",
                   std::ios::binary);
    if (!outstream.good()) {
        logger->trace("Could not open {}",
                      remove_suffix(filename, kExtension) + kExtension + ".counts");
        throw std::ofstream::failure("Bad stream");
    }

    uint64_t num_counts = 0;
    uint64_t sum_counts = 0;

    for (size_t j = 0; j < relation_counts_.size(); ++j) {
        if (!relation_counts_[j].size()) {
            sdsl::int_vector<>(bitmatrix_[j]->num_set_bits(), 0, kCountBits).serialize(outstream);

        } else {
            assert(relation_counts_[j].size() == bitmatrix_[j]->num_set_bits()
                && "Binary relation matrix must not be changed while adding relation counts");

            for (uint64_t v : relation_counts_[j]) {
                num_counts++;
                sum_counts += v;
            }

            relation_counts_[j].serialize(outstream);
        }
    }

    for (size_t j = relation_counts_.size(); j < bitmatrix_.size(); ++j) {
        sdsl::int_vector<>(bitmatrix_[j]->num_set_bits(), 0, kCountBits).serialize(outstream);
    }

    logger->info("Num relation counts: {}", num_counts);
    logger->info("Total relation count: {}", sum_counts);
}

template <typename Label>
bool ColumnCompressed<Label>::load(const std::string &filename) {
    // release the columns stored
    cached_columns_.Clear();
    bitmatrix_.clear();

    label_encoder_.clear();

    const std::string &f = remove_suffix(filename, kExtension) + kExtension;
    logger->trace("Loading annotations from file {}", f);

    try {
        std::ifstream instream(f, std::ios::binary);
        if (!instream.good())
            throw std::ifstream::failure("can't open file");

        num_rows_ = load_number(instream);

        if (!label_encoder_.load(instream))
            throw std::ifstream::failure("can't load label encoder");

        if (!label_encoder_.size())
            logger->warn("No columns in {}", f);

        for (size_t c = 0; c < label_encoder_.size(); ++c) {
            auto column = std::make_unique<bit_vector_smart>();

            if (!column->load(instream))
                throw std::ifstream::failure("can't load next column");

            if (column->size() != num_rows_)
                throw std::ifstream::failure("inconsistent column size");

            uint64_t num_set_bits = column->num_set_bits();
            logger->trace("Column: {}, Density: {}, Set bits: {}",
                          label_encoder_.decode(c),
                          static_cast<double>(num_set_bits) / column->size(),
                          num_set_bits);

            bitmatrix_.emplace_back(std::move(column));
        }

    } catch (const std::exception &e) {
        logger->error("Caught exception when loading columns from {}: {}", f, e.what());
        return false;
    } catch (...) {
        logger->error("Unknown exception when loading columns from {}", f);
        return false;
    }

    logger->trace("Annotation loading finished ({} columns)", bitmatrix_.size());
    return true;
}

template <typename Label>
bool ColumnCompressed<Label>::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_columns_.Clear();
    bitmatrix_.clear();

    label_encoder_.clear();

    bool no_errors = true;

    bool merge_successful = merge_load(filenames,
        [&](uint64_t, const Label &label, std::unique_ptr<bit_vector>&& column) {
            uint64_t num_set_bits = column->num_set_bits();
            logger->trace("Column: {}, Density: {}, Set bits: {}", label,
                          static_cast<double>(num_set_bits) / column->size(),
                          num_set_bits);

            #pragma omp critical
            {
                // set |num_rows_| with the first column inserted
                if (!bitmatrix_.size())
                    num_rows_ = column->size();

                if (column->size() == num_rows_) {
                    size_t col = label_encoder_.insert_and_encode(label);
                    assert(col <= bitmatrix_.size());

                    if (col == bitmatrix_.size()) {
                        bitmatrix_.emplace_back(std::move(column));
                    } else {
                        assert(bitmatrix_.at(col).get());
                        decompress_bitmap(col) |= *column;
                    }
                } else {
                    no_errors = false;
                }
            }
        },
        filenames.size() > 1u ? get_num_threads() : 0
    );

    if (merge_successful && no_errors) {
        logger->trace("Annotation loading finished ({} columns)", bitmatrix_.size());
        return true;
    } else {
        return false;
    }
}

template <typename Label>
bool ColumnCompressed<Label>::merge_load(const std::vector<std::string> &filenames,
                                         const ColumnCallback &callback,
                                         size_t num_threads) {
    std::atomic<bool> error_occurred = false;

    std::vector<uint64_t> offsets(filenames.size(), 0);

    // load labels
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 1; i < filenames.size(); ++i) {
        auto filename = remove_suffix(filenames[i - 1], kExtension) + kExtension;

        std::ifstream instream(filename, std::ios::binary);
        if (!instream.good()) {
            logger->error("Can't read from {}", filename);
            error_occurred = true;
        }
        std::ignore = load_number(instream);

        LabelEncoder<Label> label_encoder;
        if (!label_encoder.load(instream)) {
            logger->error("Can't load label encoder from {}", filename);
            error_occurred = true;
        }

        offsets[i] = label_encoder.size();
    }

    if (error_occurred)
        return false;

    // compute global offsets (partial sums)
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());

    // load annotations
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < filenames.size(); ++i) {
        const auto &filename = remove_suffix(filenames[i], kExtension) + kExtension;
        logger->trace("Loading annotations from file {}", filename);
        try {
            std::ifstream instream(filename, std::ios::binary);
            if (!instream.good())
                throw std::ifstream::failure("can't open file");

            const auto num_rows = load_number(instream);

            LabelEncoder<Label> label_encoder_load;
            if (!label_encoder_load.load(instream))
                throw std::ifstream::failure("can't load label encoder");

            if (!label_encoder_load.size())
                logger->warn("No labels in {}", filename);

            // update the existing and add some new columns
            for (size_t c = 0; c < label_encoder_load.size(); ++c) {
                auto new_column = std::make_unique<bit_vector_smart>();

                if (!new_column->load(instream))
                    throw std::ifstream::failure("can't load next column");

                if (new_column->size() != num_rows)
                    throw std::ifstream::failure("inconsistent column size");

                callback(offsets[i] + c,
                         label_encoder_load.decode(c),
                         std::move(new_column));
            }
        } catch (const std::exception &e) {
            logger->error("Caught exception when loading columns from {}: {}", filename, e.what());
            error_occurred = true;
        } catch (...) {
            logger->error("Unknown exception when loading columns from {}", filename);
            error_occurred = true;
        }
    }

    return !error_occurred;
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
            logger->warn("Label '{}' not found in annotation."
                         " Skipping instruction '{} -> {}'.",
                         pair.first, pair.first, pair.second);
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
        // the columns are flushed automatically on erasing from cache
        const_cast<ColumnCompressed*>(this)->cached_columns_.Clear();
        flushed_ = true;
    }
    assert(bitmatrix_.size() == label_encoder_.size());
}

template <typename Label>
void ColumnCompressed<Label>::flush(size_t j, bitmap_builder *builder) {
    assert(builder);
    assert(j < bitmatrix_.size());

    // Note: asserting that j is cached cannot be done here when this function
    //       is invovled as part of the OnEraseCallback, since the mutex locking
    //       in the caches library would cause the check to be done after it has
    //       been erased.

    bitmatrix_[j].reset();

    auto initialization_data = builder->get_initialization_data();
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
            if (num_rows_ < buffer_size_bytes_ * 8) {
                vector = new bitmap_vector(num_rows_, 0);

            // For large bitmaps, use a more space efficient builder that
            // requires only 64 bits per each 1-bit in the bitmap.
            } else if (swap_dir_.size()) {
                // use a fixed size buffer and disk swap
                vector = new bitmap_builder_set_disk(num_rows_, get_num_threads(),
                                                     buffer_size_bytes_ / 8, swap_dir_);
            } else {
                // all in RAM, the buffer is automatically resized for large bitmaps
                vector = new bitmap_builder_set(num_rows_, get_num_threads(),
                                                buffer_size_bytes_ / 8);
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

    auto &builder = *cached_columns_.Get(j);

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
const binmat::ColumnMajor& ColumnCompressed<Label>::get_matrix() const {
    flush();
    annotation_matrix_view_.update_pointer(bitmatrix_);
    return annotation_matrix_view_;
}

template <typename Label>
binmat::ColumnMajor ColumnCompressed<Label>::release_matrix() {
    flush();
    label_encoder_.clear();
    return binmat::ColumnMajor(std::move(bitmatrix_));
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
            logger->error("Dumping column {} failed", j);
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

} // namespace annot
} // namespace mtg
