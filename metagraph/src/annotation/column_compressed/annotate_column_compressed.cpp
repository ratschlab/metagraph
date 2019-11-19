#include "annotate_column_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include <fmt/format.h>

#include "serialization.hpp"
#include "string_utils.hpp"
#include "algorithms.hpp"
#include "bitmap_mergers.hpp"
#include "threading.hpp"
#include "annotate_row_compressed.hpp"

using utils::remove_suffix;

size_t kNumRowsInBlock = 50'000;


namespace annotate {

template <typename Label>
ColumnCompressed<Label>::ColumnCompressed(uint64_t num_rows,
                                          size_t num_columns_cached,
                                          bool verbose)
      : num_rows_(num_rows),
        cached_columns_(num_columns_cached,
                        caches::LRUCachePolicy<size_t>(),
                        [this](size_t j, bitmap *col_uncompressed) {
                            assert(col_uncompressed);
                            this->flush(j, *col_uncompressed);
                            delete col_uncompressed;
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
void ColumnCompressed<Label>::set_labels(Index i, const VLabels &labels) {
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
typename ColumnCompressed<Label>::SetBitPositions
ColumnCompressed<Label>::get_label_codes(Index i) const {
    assert(i < num_rows_);

    SetBitPositions label_indices;
    for (size_t j = 0; j < num_labels(); ++j) {
        if (is_set(i, j))
            label_indices.push_back(j);
    }
    return label_indices;
}

template <typename Label>
std::vector<typename ColumnCompressed<Label>::SetBitPositions>
ColumnCompressed<Label>::get_label_codes(const std::vector<Index> &indices) const {
    std::vector<SetBitPositions> rows(indices.size());

    for (size_t j = 0; j < num_labels(); ++j) {
        for (size_t i = 0; i < indices.size(); ++i) {
            assert(indices[i] < num_rows_);

            if (is_set(indices[i], j))
                rows[i].push_back(j);
        }
    }

    return rows;
}

template <typename Label>
void ColumnCompressed<Label>::add_label(Index i, const Label &label) {
    assert(i < num_rows_);

    set(i, label_encoder_.insert_and_encode(label), 1);
}

template <typename Label>
void ColumnCompressed<Label>::add_labels(Index i, const VLabels &labels) {
    for (const auto &label : labels) {
        add_label(i, label);
    }
}

template <typename Label>
void ColumnCompressed<Label>::add_labels(const std::vector<Index> &indices,
                                         const VLabels &labels) {
    for (const auto &label : labels) {
        for (Index i : indices) {
            add_label(i, label);
        }
    }
}

template <typename Label>
bool ColumnCompressed<Label>::has_label(Index i, const Label &label) const {
    try {
        return is_set(i, label_encoder_.encode(label));
    } catch (...) {
        return false;
    }
}

template <typename Label>
void ColumnCompressed<Label>
::call_relations(const std::vector<Index> &indices,
                 const Label &label,
                 std::function<void(Index, bool)> callback,
                 std::function<bool()> terminate) const {
    try {
        size_t label_code = label_encoder_.encode(label);

        for (Index i : indices) {
            if (terminate())
                return;

            callback(i, is_set(i, label_code));
        }

    } catch (...) {

        for (Index i : indices) {
            if (terminate())
                return;

            callback(i, false);
        }
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
                            decompress(col) |= *new_column;
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
        decompress(j).insert_zeros(rows);
    }
    num_rows_ += rows.size();
}

template <typename Label>
void ColumnCompressed<Label>
::call_objects(const Label &label,
               std::function<void(Index)> callback) const {
    size_t col;
    try {
        col = label_encoder_.encode(label);
    } catch (...) {
        return;
    }

    get_column(col).call_ones(callback);
}

template <typename Label>
std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
ColumnCompressed<Label>
::count_labels(const std::unordered_map<Index, size_t> &index_counts,
               size_t min_count,
               size_t count_cap) const {

    min_count = std::max(min_count, size_t(1));

    assert(count_cap >= min_count);

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<uint64_t> indices(index_counts.size());
    std::transform(index_counts.begin(), index_counts.end(), indices.begin(),
                   [](const auto &pair) { return pair.first; });

    std::vector<std::pair<uint64_t, size_t>> label_counts;
    label_counts.reserve(num_labels());

    // TODO: get rid of label encoder from here
    for (const auto &label : this->get_all_labels()) {
        size_t total_checked = 0;
        size_t total_matched = 0;
        call_relations(
            indices,
            label,
            [&](auto i, bool matched) {
                auto count = index_counts.at(i);
                total_checked += count;
                total_matched += count * matched;
            },
            [&]() { return total_matched >= count_cap
                        || total_matched + (total_sum_count - total_checked) < min_count; }
        );

        if (total_matched >= min_count)
            label_counts.emplace_back(label_encoder_.encode(label),
                                      std::min(total_matched, count_cap));
    }

    return label_counts;
}

// For each pair (first, second) in the dictionary, renames
// column |first| with |second| and merges the columns with matching names.
template <typename Label>
void ColumnCompressed<Label>
::rename_labels(const std::unordered_map<Label, Label> &dict) {
    std::vector<Label> index_to_label(label_encoder_.size());
    // old labels
    for (size_t i = 0; i < index_to_label.size(); ++i) {
        index_to_label[i] = label_encoder_.decode(i);
    }
    // new labels
    for (const auto &pair : dict) {
        try {
            index_to_label[label_encoder_.encode(pair.first)] = pair.second;
        } catch (const std::runtime_error&) {
            std::cerr << "Warning: label '" << pair.first << "' not"
                      << " found in annotation. Skipping instruction"
                      << " '" << pair.first << " -> " << pair.second << "'."
                      << std::endl;
        }
    }

    std::vector<Label> new_index_to_label;
    std::unordered_map<Label, std::set<size_t>> old_columns;
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
size_t ColumnCompressed<Label>::num_labels() const {
    assert(bitmatrix_.size() == label_encoder_.size());
    return label_encoder_.size();
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
    if (j >= bitmatrix_.size() || is_set(i, j) != value)
        decompress(j).set(i, value);
}

template <typename Label>
bool ColumnCompressed<Label>::is_set(Index i, size_t j) const {
    return get_column(j)[i];
}

template <typename Label>
void ColumnCompressed<Label>::flush() const {
    for (const auto &cached_vector : cached_columns_) {
        const_cast<ColumnCompressed*>(this)->flush(
            cached_vector.first, *cached_vector.second
        );
    }
    assert(bitmatrix_.size() == label_encoder_.size());
}

template <typename Label>
void ColumnCompressed<Label>::flush(size_t j, const bitmap &vector) {
    assert(j < bitmatrix_.size());

    // Note: asserting that j is cached cannot be done here when this function
    //       is invovled as part of the OnEraseCallback, since the mutex locking
    //       in the caches library would cause the check to be done after it has
    //       been erased.

    bitmatrix_[j].reset();
    bitmatrix_[j].reset(new bit_vector_smart(
        [&](const auto &callback) { vector.call_ones(callback); },
        vector.size(),
        vector.num_set_bits()
    ));

    assert(vector.size() == bitmatrix_[j]->size());
    assert(vector.num_set_bits() == bitmatrix_[j]->num_set_bits());
    assert(bitmatrix_[j]->num_set_bits() == bitmatrix_[j]->rank1(bitmatrix_[j]->size() - 1));
}

template <typename Label>
bitmap_dyn& ColumnCompressed<Label>::decompress(size_t j) {
    assert(j < label_encoder_.size());

    try {
        return *cached_columns_.Get(j);
    } catch (...) {
        assert(j <= bitmatrix_.size());
        // initialize the bit_vector if it exists in bitmatrix,
        // or if num_rows_ is small
        if (j == bitmatrix_.size())
            bitmatrix_.emplace_back();

        auto *vector = new bitmap_adaptive(num_rows_, 0);

        if (bitmatrix_[j].get())
            *vector |= *bitmatrix_[j];

        cached_columns_.Put(j, vector);
        return *cached_columns_.Get(j);
    }
}

template <typename Label>
const bitmap& ColumnCompressed<Label>::get_column(size_t j) const {
    if (cached_columns_.Cached(j)) {
        return (*cached_columns_.Get(j));
    } else {
        assert(j < bitmatrix_.size() && bitmatrix_[j].get());
        return (*bitmatrix_[j]);
    }
}

template <typename Label>
const bitmap& ColumnCompressed<Label>::get_column(const Label &label) const {
    return get_column(label_encoder_.encode(label));
}

template <typename Label>
void ColumnCompressed<Label>
::convert_to_row_annotator(const std::string &outfbase) const {
    flush();

    ProgressBar progress_bar(num_rows_, "Serialized rows", std::cerr, !utils::get_verbose());

    RowCompressed<Label>::write_rows(
        outfbase,
        label_encoder_,
        [&](BinaryMatrix::RowCallback write_row) {

            #pragma omp parallel for ordered schedule(dynamic) num_threads(get_num_threads())
            for (uint64_t i = 0; i < num_rows_; i += kNumRowsInBlock) {

                uint64_t begin = i;
                uint64_t end = std::min(i + kNumRowsInBlock, num_rows_);

                std::vector<SetBitPositions> rows(end - begin);

                assert(begin <= end);
                assert(end <= num_rows_);

                // TODO: use RowsFromColumnsTransformer
                for (size_t j = 0; j < bitmatrix_.size(); ++j) {
                    bitmatrix_[j]->call_ones_in_range(begin, end,
                        [&](uint64_t idx) { rows[idx - begin].push_back(j); }
                    );
                }

                #pragma omp ordered
                {
                    for (const auto &row : rows) {
                        write_row(row);
                        ++progress_bar;
                    }
                }
            }
        }
    );
}

template <typename Label>
void ColumnCompressed<Label>
::convert_to_row_annotator(RowCompressed<Label> *annotator,
                           size_t num_threads) const {
    assert(annotator);

    flush();

    ProgressBar progress_bar(num_rows_, "Processed rows", std::cerr, !utils::get_verbose());

    annotator->reinitialize(num_rows_);
    annotator->label_encoder_ = label_encoder_;

    if (num_threads <= 1) {
        add_labels(0, num_rows_, annotator, &progress_bar);
        return;
    }

    ThreadPool thread_pool(num_threads);
    for (uint64_t i = 0; i < num_rows_; i += kNumRowsInBlock) {
        thread_pool.enqueue(
            [this](auto... args) { this->add_labels(args...); },
            i, std::min(i + kNumRowsInBlock, num_rows_),
            annotator,
            &progress_bar
        );
    }
    thread_pool.join();
}

template <typename Label>
void ColumnCompressed<Label>::add_labels(uint64_t begin, uint64_t end,
                                         RowCompressed<Label> *annotator,
                                         ProgressBar *progress_bar) const {
    assert(begin <= end);
    assert(end <= annotator->matrix_->num_rows());

    // TODO: use RowsFromColumnsTransformer
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        bitmatrix_[j]->call_ones_in_range(begin, end,
            [&](uint64_t idx) { annotator->matrix_->set(idx, j); }
        );
    }
    if (progress_bar)
        *progress_bar += end - begin;
}

template <typename Label>
bool ColumnCompressed<Label>
::dump_columns(const std::string &prefix, size_t num_threads) const {
    bool success = true;

    #pragma omp parallel for num_threads(num_threads)
    for (uint64_t i = 0; i < bitmatrix_.size(); ++i) {
        const auto &column = get_column(i);
        const std::string outfile = remove_suffix(prefix, kExtension)
            + fmt::format(".{}.text.annodbg", i);

        try {
            VectorFileOutStream outstream(outfile);
            outstream.write_value(num_objects());
            outstream.write_value(column.num_set_bits());
            column.call_ones([&](const auto &pos) { outstream.write_value(pos); });
        } catch (...) {
            std::cerr << "Failed to write to " << outfile << std::endl;
            success = false;
        }
    }

    return success;
}

template <class Annotator>
class IterateTransposed : public Annotator::IterateRows {
  public:
    IterateTransposed(std::unique_ptr<utils::RowsFromColumnsTransformer> transformer)
          : row_iterator_(std::move(transformer)) {}

    typename Annotator::SetBitPositions next_row() override final {
        return row_iterator_.next_row<typename Annotator::SetBitPositions>();
    }

  private:
    utils::RowsFromColumnsIterator row_iterator_;
};

template <typename Label>
std::unique_ptr<typename ColumnCompressed<Label>::IterateRows>
ColumnCompressed<Label>::iterator() const {
    flush();

    return std::make_unique<IterateTransposed<ColumnCompressed<Label>>>(
        std::make_unique<utils::RowsFromColumnsTransformer>(bitmatrix_)
    );
};

template class ColumnCompressed<std::string>;

} // namespace annotate
