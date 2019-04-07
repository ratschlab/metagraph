#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"
#include "vector_row_binmat.hpp"
#include "eigen_spmat.hpp"
#include "binary_matrix.hpp"

using utils::remove_suffix;


namespace annotate {


template <typename Label>
RowCompressed<Label>::RowCompressed(uint64_t num_rows, bool sparse)  {
    if (sparse) {
        matrix_.reset(new EigenSpMat(num_rows));
    } else {
        matrix_.reset(new VectorRowBinMat(num_rows));
    }
}

template <typename Label>
void RowCompressed<Label>::reinitialize(uint64_t num_rows) {
    if (dynamic_cast<EigenSpMat*>(matrix_.get())) {
        matrix_.reset(new EigenSpMat(num_rows));
    } else {
        matrix_.reset(new VectorRowBinMat(num_rows));
    }

    label_encoder_.clear();
}

template <typename Label>
void RowCompressed<Label>::set_labels(Index i, const VLabels &labels) {
    assert(i < matrix_->num_rows());

    matrix_->clear_row(i);
    add_labels(i, labels);
}

template <typename Label>
typename RowCompressed<Label>::VLabels
RowCompressed<Label>::get_labels(Index i) const {
    VLabels labels;
    for (auto col : matrix_->get_row(i)) {
        labels.push_back(label_encoder_.decode(col));
    }
    return labels;
}

template <typename Label>
void RowCompressed<Label>::add_label(Index i, const Label &label) {
    matrix_->set(i, label_encoder_.insert_and_encode(label));
}

template <typename Label>
void RowCompressed<Label>::add_labels(Index i, const VLabels &labels) {
    for (const auto &label : labels) {
        add_label(i, label);
    }
}

template <typename Label>
void RowCompressed<Label>::add_labels(const std::vector<Index> &indices,
                                      const VLabels &labels) {
    for (Index i : indices) {
        add_labels(i, labels);
    }
}

template <typename Label>
void RowCompressed<Label>::add_labels_fast(const std::vector<Index> &indices,
                                           const VLabels &labels) {
    std::vector<uint64_t> col_ids;
    col_ids.reserve(labels.size());
    for (const auto &label : labels) {
        col_ids.push_back(label_encoder_.insert_and_encode(label));
    }
    std::sort(col_ids.begin(), col_ids.end());
    col_ids.erase(std::unique(col_ids.begin(), col_ids.end()), col_ids.end());

    auto unique_indices = indices;
    std::sort(unique_indices.begin(), unique_indices.end());
    unique_indices.erase(
        std::unique(unique_indices.begin(), unique_indices.end()),
        unique_indices.end()
    );

    if (dynamic_cast<VectorRowBinMat*>(matrix_.get())) {
        for (Index i : unique_indices) {
            for (auto j : col_ids) {
                dynamic_cast<VectorRowBinMat&>(*matrix_).force_set(i, j);
            }
        }
    } else {
        for (Index i : unique_indices) {
            for (auto j : col_ids) {
                matrix_->set(i, j);
            }
        }
    }
}

template <typename Label>
bool RowCompressed<Label>::has_label(Index i, const Label &label) const {
    try {
        return matrix_->get(i, label_encoder_.encode(label));
    } catch (...) {
        return false;
    }
}

template <typename Label>
bool RowCompressed<Label>::has_labels(Index i, const VLabels &labels) const {
    std::set<size_t> querying_codes;
    try {
        for (const auto &label : labels) {
            querying_codes.insert(label_encoder_.encode(label));
        }
    } catch (...) {
        return false;
    }
    std::set<size_t> encoded_labels;
    for (auto col : matrix_->get_row(i)) {
        encoded_labels.insert(col);
    }
    return std::includes(encoded_labels.begin(), encoded_labels.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <typename Label>
void RowCompressed<Label>::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension,
                            std::ios::binary);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }
    label_encoder_.serialize(outstream);
    matrix_->serialize(outstream);
}

template <typename Label>
bool RowCompressed<Label>::merge_load(const std::vector<std::string> &filenames) {
    std::ifstream instream(remove_suffix(filenames.at(0), kExtension) + kExtension,
                           std::ios::binary);
    if (!instream.good())
        return false;

    try {
        bool is_successfully_loaded = label_encoder_.load(instream)
                                            && matrix_->load(instream);
        if (filenames.size() == 1 || !is_successfully_loaded)
            return is_successfully_loaded;

        assert(filenames.size() > 1);

        if (!dynamic_cast<VectorRowBinMat*>(matrix_.get())) {
            std::cerr << "Error: loading from multiple row annotators is supported"
                      << " only for the VectorRowBinMat representation" << std::endl;
            exit(1);
        }

        auto &matrix = dynamic_cast<VectorRowBinMat&>(*matrix_);
        auto next_block = std::make_unique<VectorRowBinMat>(matrix_->num_rows());

        for (auto filename : filenames) {
            if (filename == filenames[0])
                continue;

            std::ifstream instream(remove_suffix(filename, kExtension) + kExtension,
                                   std::ios::binary);
            if (!instream.good())
                return false;

            LabelEncoder<Label> label_encoder_load;
            if (!label_encoder_load.load(instream)
                    || !next_block->load(instream)
                    || next_block->num_rows() != matrix_->num_rows())
                return false;

            // add new columns
            std::vector<uint64_t> new_column_positions(label_encoder_load.size());

            for (size_t c = 0; c < label_encoder_load.size(); ++c) {
                try {
                    new_column_positions[c]
                        = label_encoder_.encode(label_encoder_load.decode(c));
                    std::cerr << "Warning: the merged annotations have common"
                              << " labels. The merged result will include"
                              << " duplicates and be bigger than expected."
                              << " Detected label: "
                              << label_encoder_load.decode(c) << std::endl;
                } catch (...) {
                    // the label is new
                    new_column_positions[c]
                        = label_encoder_.insert_and_encode(label_encoder_load.decode(c));
                }
            }

            // set all bits from the next block
            for (uint64_t i = 0; i < next_block->num_rows(); ++i) {
                for (auto j : next_block->get_row(i)) {
                    matrix.force_set(i, new_column_positions[j]);
                }
            }
        }
    } catch (...) {
        return false;
    }
    return true;
}

template <typename Label>
void RowCompressed<Label>::insert_rows(const std::vector<Index> &rows) {
    matrix_->insert_rows(rows);
}

// Get labels that occur at least in |presence_ratio| rows.
// If |presence_ratio| = 0, return all occurring labels.
template <typename Label>
typename RowCompressed<Label>::VLabels
RowCompressed<Label>::get_labels(const std::vector<Index> &indices,
                                 double presence_ratio) const {
    assert(presence_ratio >= 0 && presence_ratio <= 1);

    const size_t min_labels_discovered =
                        presence_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * presence_ratio);
    // const size_t max_labels_missing = indices.size() - min_labels_discovered;

    auto counts = count_labels(indices);

    VLabels filtered_labels;

    for (size_t i = 0; i < counts.size(); ++i) {
        if (counts[i] && counts[i] >= min_labels_discovered)
            filtered_labels.push_back(label_encoder_.decode(i));
    }

    return filtered_labels;
}

// Count all labels collected from the given rows.
template <typename Label>
std::vector<uint64_t>
RowCompressed<Label>::count_labels(const std::vector<Index> &indices) const {
    std::vector<uint64_t> counter(num_labels(), 0);

    for (Index i : indices) {
        for (auto col : matrix_->get_row(i)) {
            counter[col]++;
        }
    }

    return counter;
}

template <typename Label>
uint64_t RowCompressed<Label>::num_objects() const {
    return matrix_->num_rows();
}

template <typename Label>
size_t RowCompressed<Label>::num_labels() const {
    assert(label_encoder_.size() == matrix_->num_columns());
    return label_encoder_.size();
}

template <typename Label>
uint64_t RowCompressed<Label>::num_relations() const {
    return matrix_->num_relations();
}

template <typename Label>
LabelEncoder<Label>* RowCompressed<Label>::load_label_encoder(const std::string &filename) {
    std::ifstream instream(remove_suffix(filename, kExtension) + kExtension,
                           std::ios::binary);

    auto *label_encoder = new LabelEncoder<Label>();
    if (!label_encoder->load(instream))
        throw std::ifstream::failure("Bad stream");

    return label_encoder;
}

template <typename Label>
bool RowCompressed<Label>::stream_counts(const std::string &filename,
                                         uint64_t &num_objects,
                                         uint64_t &num_relations,
                                         bool sparse) {
    throw std::runtime_error("Not implemented");
    //return stream_rows(filename, [](const std::vector<uint64_t> &row){}, sparse);
}

template <typename Label>
RowCompressed<Label>
::StreamRows::StreamRows(const std::string &filename, bool sparse) {
    std::ifstream instream(remove_suffix(filename, kExtension) + kExtension,
                           std::ios::binary);

    if (!instream.good())
        throw std::runtime_error("Bad stream");

    // skip header
    LabelEncoder<Label> label_encoder_load;
    if (!label_encoder_load.load(instream))
        throw std::runtime_error("Could not load label encoder");

    // rows
    if (sparse) {
        throw std::runtime_error("Not implemented");
    } else {
        sr_ = std::make_unique<VectorRowBinMat::StreamRows>(instream);
    }
}

template <typename Label>
void RowCompressed<Label>
::write_rows(const std::string &filename,
             const LabelEncoder<Label> &label_encoder,
             const uint64_t num_rows,
             const std::function<void (const std::function<void(void)>&, const std::function<void (const std::vector<uint64_t> &)>&, const std::function<void(void)>&)> &callback,
             //BinaryMatrix::GetRow &callback,
             bool sparse) {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension,
                            std::ios::binary);
    if (!outstream.good())
        throw std::ofstream::failure("Bad stream");

    label_encoder.serialize(outstream);

    if (sparse) {
        throw std::runtime_error("Not implemented");
    } else {
        uint64_t num_cols = label_encoder.size();
        VectorRowBinMat::write_rows(outstream, callback, num_rows, num_cols);
    }
}

template class RowCompressed<std::string>;

} // namespace annotate
