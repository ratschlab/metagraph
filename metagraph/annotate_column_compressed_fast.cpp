#include "annotate_column_compressed_fast.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"

using utils::remove_suffix;


namespace annotate {

template <typename Label>
FastColumnCompressed<Label>::FastColumnCompressed(uint64_t num_rows,
                                                  size_t num_columns_cached,
                                                  bool verbose)
      : num_rows_(num_rows),
        cached_columns_(
            num_columns_cached,
            caches::LRUCachePolicy<size_t>(),
            [this](size_t j, std::vector<bool> *col_uncompressed) {
                this->flush(j, *col_uncompressed);
                delete col_uncompressed;
            }
        ),
        cached_index_(
            std::max(static_cast<size_t>(1), num_columns_cached / 2),
            caches::LRUCachePolicy<size_t>(),
            [this](size_t t, std::vector<bool> *col_uncompressed) {
                this->flush_index(t, *col_uncompressed);
                delete col_uncompressed;
            }
        ),
        verbose_(verbose) {}

template <typename Label>
FastColumnCompressed<Label>::FastColumnCompressed(ColumnCompressed<Label>&& annotator,
                                                  size_t num_columns_cached,
                                                  bool verbose,
                                                  bool build_index)
      : FastColumnCompressed(annotator.num_rows_, num_columns_cached, verbose) {
    annotator.cached_columns_.Clear();

    bitmatrix_ = std::move(annotator.bitmatrix_);
    label_encoder_ = std::move(annotator.label_encoder_);

    if (build_index)
        rebuild_index();
}

template <typename Label>
FastColumnCompressed<Label>::~FastColumnCompressed() {
    cached_columns_.Clear();
    cached_index_.Clear();
}

template <typename Label>
void FastColumnCompressed<Label>::set_labels(Index i, const VLabels &labels) {
    assert(i < num_rows_);

    // add new labels
    for (const auto &label : labels) {
        label_encoder_.insert_and_encode(label);
    }

    // labels as a row
    std::vector<bool> row(num_labels(), 0);
    for (const auto &label : labels) {
        row[label_encoder_.insert_and_encode(label)] = 1;
    }

    // reset labels
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

template <typename Label>
typename FastColumnCompressed<Label>::VLabels
FastColumnCompressed<Label>::get_labels(Index i) const {
    assert(i < num_rows_);

    VLabels labels;
    for (size_t j : get_row(i)) {
        labels.push_back(label_encoder_.decode(j));
    }
    return labels;
}

template <typename Label>
void FastColumnCompressed<Label>::add_label(Index i, const Label &label) {
    assert(i < num_rows_);

    size_t j = label_encoder_.insert_and_encode(label);
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

template <typename Label>
void FastColumnCompressed<Label>::add_labels(Index i, const VLabels &labels) {
    for (const auto &label : labels) {
        add_label(i, label);
    }
}

template <typename Label>
void FastColumnCompressed<Label>::add_labels(const std::vector<Index> &indices,
                                             const VLabels &labels) {
    for (const auto &label : labels) {
        for (Index i : indices) {
            add_label(i, label);
        }
    }
}

template <typename Label>
bool FastColumnCompressed<Label>::has_label(Index i, const Label &label) const {
    try {
        return get_entry(i, label_encoder_.encode(label));
    } catch (...) {
        return false;
    }
}

template <typename Label>
bool FastColumnCompressed<Label>::has_labels(Index i, const VLabels &labels) const {
    for (const auto &label : labels) {
        if (!has_label(i, label))
            return false;
    }
    return true;
}

template <typename Label>
void FastColumnCompressed<Label>::serialize(const std::string &filename) const {
    flush();

    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    serialize_number(outstream, num_rows_);

    label_encoder_.serialize(outstream);

    for (const auto &column : bitmatrix_) {
        column->serialize(outstream);
    }

    // serialize the column index
    outstream = std::ofstream(remove_suffix(filename, kExtension) + kIndexExtension);

    serialize_number(outstream, index_.size());
    for (size_t t = 0; t < index_.size(); ++t) {
        index_[t]->serialize(outstream);
        serialize_number_vector(outstream, index_to_columns_[t]);
    }
}


template <typename Label>
bool FastColumnCompressed<Label>::merge_load(const std::vector<std::string> &filenames) {
    // release the columns stored
    cached_columns_.Clear();
    cached_index_.Clear();

    bitmatrix_.clear();
    index_.clear();
    column_to_index_.clear();
    index_to_columns_.clear();

    label_encoder_ = decltype(label_encoder_)();

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
                num_rows_ = load_number(instream);
            } else if (num_rows_ != load_number(instream)) {
                return false;
            }

            decltype(label_encoder_) label_encoder_load;
            if (!label_encoder_load.load(instream))
                return false;

            // extend the label dictionary
            for (size_t c = 0; c < label_encoder_load.size(); ++c) {
                label_encoder_.insert_and_encode(label_encoder_load.decode(c));
            }

            // update the existing and add some new columns
            bitmatrix_.resize(label_encoder_.size());
            for (size_t c = 0; c < label_encoder_load.size(); ++c) {
                size_t col = label_encoder_.encode(label_encoder_load.decode(c));

                auto *new_column = new bit_vector_small();
                new_column->load(instream);

                if (verbose_) {
                    auto num_set_bits = new_column->num_set_bits();

                    std::cout << "Label <" << label_encoder_load.decode(c)
                                           << ">, "
                              << "density: "
                              << static_cast<double>(num_set_bits) / new_column->size()
                              << ", set bits: " << num_set_bits << std::endl;
                }

                if (bitmatrix_.at(col).get()) {
                    new_column->add_to(&decompress(col));
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

                size_t index_size = load_number(instream);

                column_to_index_.resize(num_labels());

                for (size_t t = 0; t < index_size; ++t) {
                    index_.emplace_back(new bit_vector_small());
                    index_.back()->load(instream);

                    index_to_columns_.emplace_back(load_number_vector<uint64_t>(instream));

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
                auto num_set_bits = index_[t]->num_set_bits();

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

template <typename Label>
void FastColumnCompressed<Label>::insert_rows(const std::vector<Index> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));

    for (size_t j = 0; j < label_encoder_.size(); ++j) {
        auto &column = decompress(j);
        std::vector<bool> old(num_rows_ + rows.size(), 0);
        old.swap(column);

        uint64_t i = 0;
        uint64_t num_inserted = 0;

        for (auto next_inserted_row : rows) {
            while (i + num_inserted < next_inserted_row) {
                assert(i < old.size() && "Invalid indexes of inserted rows");
                assert(i + num_inserted < column.size() && "Invalid indexes of inserted rows");

                column[i + num_inserted] = old[i];
                i++;
            }
            // insert 0, not labeled edge
            num_inserted++;
        }
        while (i < old.size()) {
            assert(i + num_inserted < column.size() && "Invalid indexes of inserted rows");

            column[i + num_inserted] = old[i];
            i++;
        }
    }

    for (size_t t = 0; t < index_to_columns_.size(); ++t) {
        auto &column = decompress_index(t);
        std::vector<bool> old(num_rows_ + rows.size(), 0);
        old.swap(column);

        uint64_t i = 0;
        uint64_t num_inserted = 0;

        for (auto next_inserted_row : rows) {
            while (i + num_inserted < next_inserted_row) {
                assert(i < old.size() && "Invalid indexes of inserted rows");
                assert(i + num_inserted < column.size() && "Invalid indexes of inserted rows");

                column[i + num_inserted] = old[i];
                i++;
            }
            // insert 0, not labeled edge
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


// Get labels that occur at least in |presence_ratio| rows.
// If |presence_ratio| = 0, return all occurring labels.
template <typename Label>
typename FastColumnCompressed<Label>::VLabels
FastColumnCompressed<Label>::get_labels(const std::vector<Index> &indices,
                                        double presence_ratio) const {
    assert(presence_ratio >= 0 && presence_ratio <= 1);

    const size_t min_labels_discovered =
                        presence_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * presence_ratio);
    const size_t max_labels_missing = indices.size() - min_labels_discovered;

    VLabels filtered_labels;

    std::vector<bool> label_filter(num_labels(), true);
    std::vector<size_t> discovered(num_labels(), 0);
    uint64_t iteration = 0;

    for (Index i : indices) {
        for (size_t j : filter_row(i, label_filter)) {
            if (++discovered[j] >= min_labels_discovered) {
                filtered_labels.push_back(label_encoder_.decode(j));
                label_filter[j] = false;
            }
        }

        if (++iteration <= max_labels_missing)
            continue;

        for (size_t j = 0; j < label_filter.size(); ++j) {
            if (label_filter[j]
                    && iteration - discovered[j] > max_labels_missing) {
                label_filter[j] = false;
            }
        }
    }

    return filtered_labels;
}

// Count all labels collected from the given rows
// and return top |num_top| with the their counts.
template <typename Label>
std::vector<std::pair<Label, size_t>>
FastColumnCompressed<Label>::get_top_labels(const std::vector<Index> &indices,
                                            size_t num_top) const {
    auto counter = count_labels(indices);

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

    std::vector<std::pair<Label, size_t>> top_counts;
    for (const auto &encoded_pair : counts) {
        top_counts.emplace_back(label_encoder_.decode(encoded_pair.first),
                                encoded_pair.second);
    }

    return top_counts;
}

template <typename Label>
size_t FastColumnCompressed<Label>::num_labels() const {
    return label_encoder_.size();
}

template <typename Label>
double FastColumnCompressed<Label>::sparsity() const {
    uint64_t num_set_bits = 0;

    flush();

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        num_set_bits += bitmatrix_[j]->num_set_bits();
    }

    double density = static_cast<double>(num_set_bits) / num_labels() / num_rows_;
    return 1 - density;
}

template <typename Label>
void FastColumnCompressed<Label>::rebuild_index(size_t num_aux_cols) {
    flush();

    index_.clear();
    column_to_index_.clear();
    index_to_columns_.clear();

    if (num_aux_cols == static_cast<size_t>(-1)) {
        // In the worst case, on average, we access
        // |num_aux_cols| auxiliary columns, plus
        // (|density| * |num_columns|) * (|num_columns| / |num_aux_cols|)
        // real columns per query.
        num_aux_cols = std::sqrt((1 - sparsity()) * num_labels() * num_labels());
    }

    if (!num_aux_cols)
        return;

    index_.resize(num_aux_cols);
    column_to_index_.resize(bitmatrix_.size());
    index_to_columns_.resize(num_aux_cols);

    if (verbose_)
        std::cout << "Updating " << num_aux_cols << " index columns" << std::endl;

    size_t block_width = (num_labels() + num_aux_cols - 1) / num_aux_cols;

    if (verbose_)
        std::cout << "Processed index columns..." << std::endl;

    for (size_t t = 0; t < num_aux_cols; ++t) {
        auto &index_col = decompress_index(t);

        if (verbose_)
            std::cout << "Block " << t << ": " << std::flush;

        for (size_t j = t * block_width; j < std::min((t + 1) * block_width,
                                                      num_labels()); ++j) {
            column_to_index_[j] = t;
            index_to_columns_[t].push_back(j);

            if (verbose_)
                std::cout << j << ".." << std::flush;

            bitmatrix_[j]->add_to(&index_col);
        }

        if (verbose_)
            std::cout << std::endl;
    }
}

template <typename Label>
std::vector<size_t> FastColumnCompressed<Label>::get_row(Index i) const {
    flush();

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

template <typename Label>
std::vector<size_t>
FastColumnCompressed<Label>
::filter_row(Index i, const std::vector<bool> &filter) const {
    flush();

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

template <typename Label>
std::vector<uint64_t>
FastColumnCompressed<Label>
::count_labels(const std::vector<Index> &indices) const {
    std::vector<uint64_t> counter(bitmatrix_.size(), 0);

    for (Index i : indices) {
        for (size_t j : get_row(i)) {
            counter[j]++;
        }
    }

    return counter;
}

template <typename Label>
bool FastColumnCompressed<Label>::get_entry(Index i, size_t j) const {
    if (cached_columns_.Cached(j)) {
        return (*cached_columns_.Get(j))[i];
    } else if (j < bitmatrix_.size()) {
        assert(bitmatrix_[j].get());
        return (*bitmatrix_[j])[i];
    } else {
        return false;
    }
}

template <typename Label>
bool FastColumnCompressed<Label>::get_index_entry(Index i, size_t t) const {
    if (cached_index_.Cached(t)) {
        return (*cached_index_.Get(t))[i];
    } else {
        assert(t < index_.size() && index_[t].get());
        return (*index_[t])[i];
    }
}

template <typename Label>
void FastColumnCompressed<Label>::flush() const {
    if (to_update_) {
        for (const auto &cached_vector : cached_columns_) {
            const_cast<FastColumnCompressed*>(this)->flush(
                cached_vector.first, *cached_vector.second
            );
        }
        const_cast<FastColumnCompressed*>(this)->to_update_ = false;
    }

    if (to_update_index_) {
        for (const auto &cached_index : cached_index_) {
            const_cast<FastColumnCompressed*>(this)->flush_index(
                cached_index.first, *cached_index.second
            );
        }
        const_cast<FastColumnCompressed*>(this)->to_update_index_ = false;
    }
}

template <typename Label>
void FastColumnCompressed<Label>::flush(size_t j, const std::vector<bool> &vector) {
    assert(cached_columns_.Cached(j));

    if (!to_update_)
        return;

    while (bitmatrix_.size() <= j) {
        bitmatrix_.emplace_back();
    }

    bitmatrix_[j].reset();
    bitmatrix_[j].reset(new bit_vector_small(vector));
}

template <typename Label>
void FastColumnCompressed<Label>::flush_index(size_t t, const std::vector<bool> &vector) {
    assert(cached_index_.Cached(t));
    assert(t < index_.size());

    if (!to_update_index_)
        return;

    index_[t].reset();
    index_[t].reset(new bit_vector_small(vector));
}

template <typename Label>
std::vector<bool>& FastColumnCompressed<Label>::decompress(size_t j) {
    assert(j < num_labels());

    to_update_ = true;

    try {
        return *cached_columns_.Get(j);
    } catch (...) {
        auto *bit_vector = new std::vector<bool>(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j].get()) {
            bitmatrix_[j]->add_to(bit_vector);
            bitmatrix_[j].reset();
        }

        cached_columns_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template <typename Label>
std::vector<bool>& FastColumnCompressed<Label>::decompress_index(size_t t) {
    assert(t < index_.size());

    to_update_index_ = true;

    try {
        return *cached_index_.Get(t);
    } catch (...) {
        std::vector<bool> *bit_vector = new std::vector<bool>(num_rows_, 0);;

        if (t < index_.size() && index_[t].get()) {
            index_[t]->add_to(bit_vector);
            index_[t].reset();
        }

        cached_index_.Put(t, bit_vector);
        return *bit_vector;
    }
}

template <typename Label>
void FastColumnCompressed<Label>::update_index() {
    if (!index_to_columns_.size())
        return;

    for (size_t j = column_to_index_.size(); j < num_labels(); ++j) {
        // distribute all new columns uniformly
        size_t t = j % index_to_columns_.size();
        column_to_index_.push_back(t);
        index_to_columns_[t].push_back(j);
    }
}

template class FastColumnCompressed<std::string>;

} // namespace annotate
