#include "annotate_column_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"
#include "annotate_row_compressed.hpp"

using utils::remove_suffix;

size_t kNumRowsInBlock = 5'000'000;


namespace annotate {

template <typename Label>
ColumnCompressed<Label>::ColumnCompressed(uint64_t num_rows,
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
        verbose_(verbose) {
    assert(num_columns_cached > 0);
}

template <typename Label>
ColumnCompressed<Label>::~ColumnCompressed() {
    cached_columns_.Clear();
}

template <typename Label>
void ColumnCompressed<Label>::set_labels(Index i, const VLabels &labels) {
    assert(i < num_rows_);

    for (size_t j = 0; j < label_encoder_.size(); ++j) {
        decompress(j)[i] = 0;
    }
    for (const auto &label : labels) {
        decompress(label_encoder_.insert_and_encode(label))[i] = 1;
    }
}

template <typename Label>
typename ColumnCompressed<Label>::VLabels
ColumnCompressed<Label>::get_labels(Index i) const {
    assert(i < num_rows_);

    flush();

    VLabels labels;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j].get());

        if ((*bitmatrix_[j])[i])
            labels.push_back(label_encoder_.decode(j));
    }
    return labels;
}

template <typename Label>
void ColumnCompressed<Label>::add_label(Index i, const Label &label) {
    assert(i < num_rows_);

    decompress(label_encoder_.insert_and_encode(label))[i] = 1;
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
        auto j = label_encoder_.encode(label);
        if (cached_columns_.Cached(j)) {
            return (*cached_columns_.Get(j))[i];
        } else {
            return (*bitmatrix_[j])[i];
        }
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

    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

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

    label_encoder_ = decltype(label_encoder_)();

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

template <typename Label>
void ColumnCompressed<Label>::insert_rows(const std::vector<Index> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));

    for (size_t j = 0; j < label_encoder_.size(); ++j) {
        std::vector<bool> &column = decompress(j);
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

// For each pair (first, second) in the dictionary, renames
// column |first| with |second| and merges the columns with matching names.
template <typename Label>
void ColumnCompressed<Label>::rename_columns(const std::map<std::string, std::string> &dict) {
    cached_columns_.Clear();

    std::vector<std::string> old_index_to_label(label_encoder_.size());
    for (size_t i = 0; i < old_index_to_label.size(); ++i) {
        old_index_to_label[i] = label_encoder_.decode(i);
    }
    for (const auto &pair : dict) {
        old_index_to_label[label_encoder_.encode(pair.first)] = pair.second;
    }

    std::vector<std::string> index_to_label;
    std::unordered_map<std::string, std::set<size_t>> old_columns;
    for (size_t i = 0; i < old_index_to_label.size(); ++i) {
        if (!old_columns.count(old_index_to_label[i]))
            index_to_label.push_back(old_index_to_label[i]);

        old_columns[old_index_to_label[i]].insert(i);
    }

    std::vector<std::unique_ptr<bit_vector_small>> old_bitmatrix;
    old_bitmatrix.swap(bitmatrix_);

    label_encoder_ = decltype(label_encoder_)();

    for (const auto &label : index_to_label) {
        label_encoder_.insert_and_encode(label);

        const auto &cols = old_columns[label];

        assert(cols.size());
        if (cols.size() == 1) {
            bitmatrix_.emplace_back(std::move(old_bitmatrix[*cols.begin()]));
        } else {
            std::vector<bool> bit_vector(num_rows_, 0);
            for (size_t c : cols) {
                old_bitmatrix[c]->add_to(&bit_vector);
                old_bitmatrix[c].reset();
            }
            bitmatrix_.emplace_back(new bit_vector_small(bit_vector));
        }
    }
}

// Get labels that occur at least in |presence_ratio| rows.
// If |presence_ratio| = 0, return all occurring labels.
template <typename Label>
typename ColumnCompressed<Label>::VLabels
ColumnCompressed<Label>::get_labels(const std::vector<Index> &indices,
                                    double presence_ratio) const {
    assert(presence_ratio >= 0 && presence_ratio <= 1);

    flush();

    const size_t min_labels_discovered =
                        presence_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * presence_ratio);
    const size_t max_labels_missing = indices.size() - min_labels_discovered;

    VLabels filtered_labels;

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        uint64_t discovered = 0;
        uint64_t missing = 0;

        for (Index i : indices) {
            if ((*bitmatrix_[j])[i]) {
                if (++discovered >= min_labels_discovered) {
                    filtered_labels.push_back(label_encoder_.decode(j));
                    break;
                }
            } else if (++missing > max_labels_missing) {
                break;
            }
        }
    }

    return filtered_labels;
}

// Count all labels collected from the given rows
// and return top |num_top| with the their counts.
template <typename Label>
std::vector<std::pair<Label, size_t>>
ColumnCompressed<Label>::get_top_labels(const std::vector<Index> &indices,
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
size_t ColumnCompressed<Label>::num_labels() const {
    return label_encoder_.size();
}

template <typename Label>
double ColumnCompressed<Label>::sparsity() const {
    uint64_t num_set_bits = 0;

    flush();

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        num_set_bits += bitmatrix_[j]->num_set_bits();
    }

    return 1 - static_cast<double>(num_set_bits) / num_labels() / num_rows_;
}

template <typename Label>
std::vector<uint64_t>
ColumnCompressed<Label>::count_labels(const std::vector<Index> &indices) const {
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
void ColumnCompressed<Label>::flush(size_t j, const std::vector<bool> &vector) {
    assert(cached_columns_.Cached(j));

    while (bitmatrix_.size() <= j) {
        bitmatrix_.emplace_back();
    }

    bitmatrix_[j].reset();
    bitmatrix_[j].reset(new bit_vector_small(vector));
}

template <typename Label>
std::vector<bool>& ColumnCompressed<Label>::decompress(size_t j) {
    assert(j < label_encoder_.size());

    try {
        return *cached_columns_.Get(j);
    } catch (...) {
        std::vector<bool> *bit_vector = new std::vector<bool>(num_rows_, 0);

        if (j < bitmatrix_.size() && bitmatrix_[j].get()) {
            bitmatrix_[j]->add_to(bit_vector);
        }

        cached_columns_.Put(j, bit_vector);
        return *bit_vector;
    }
}

template <typename Label>
void ColumnCompressed<Label>::convert_to_row_annotator(RowCompressed<Label> *annotator,
                                                       size_t num_threads) const {
    flush();

    annotator->matrix_->reinitialize(num_rows_);
    annotator->label_encoder_ = label_encoder_;

    if (num_threads <= 1) {
        add_labels(0, num_rows_, annotator);
        return;
    }

    utils::ThreadPool thread_pool(num_threads);
    for (uint64_t i = 0; i < num_rows_; i += kNumRowsInBlock) {
        thread_pool.enqueue(
            [this](const auto&... args) { this->add_labels(args...); },
            i, std::min(i + kNumRowsInBlock, num_rows_),
            annotator
        );
    }
    thread_pool.join();
}

template <typename Label>
void ColumnCompressed<Label>::add_labels(uint64_t begin, uint64_t end,
                                         RowCompressed<Label> *annotator) const {
    assert(begin <= end);
    assert(end <= annotator->matrix_->size());

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        uint64_t first = begin == 0
                            ? 1
                            : bitmatrix_[j]->rank1(begin - 1) + 1;
        uint64_t last = bitmatrix_[j]->rank1(end - 1);

        for (uint64_t r = first; r <= last; ++r) {
            annotator->matrix_->set_bit(bitmatrix_[j]->select1(r), j);
        }
    }
}

template class ColumnCompressed<std::string>;

} // namespace annotate
