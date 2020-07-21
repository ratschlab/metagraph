#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>
#include <set>

#include "common/serialization.hpp"
#include "common/utils/string_utils.hpp"
#include "annotation/binary_matrix/row_vector/vector_row_binmat.hpp"
#include "annotation/binary_matrix/row_vector/eigen_spmat.hpp"


namespace mtg {
namespace annot {

using utils::remove_suffix;
using mtg::annot::binmat::BinaryMatrix;
using mtg::annot::binmat::StreamRows;
using mtg::annot::binmat::EigenSpMat;
using mtg::annot::binmat::VectorRowBinMat;


template <typename Label>
RowCompressed<Label>::RowCompressed(uint64_t num_rows, bool sparse)  {
    if (sparse) {
        matrix_.reset(new EigenSpMat(num_rows));
    } else {
        matrix_.reset(new VectorRowBinMat<>(num_rows));
    }
}

template <typename Label>
template <typename RowType>
RowCompressed<Label>::RowCompressed(std::vector<RowType>&& annotation_rows,
                                    const std::vector<Label> &labels)
      : matrix_(new VectorRowBinMat<RowType>(std::move(annotation_rows), labels.size())) {
    for (const auto &label : labels) {
        label_encoder_.insert_and_encode(label);
    }
}

template RowCompressed<std::string>::RowCompressed(std::vector<SmallVector<uint32_t>>&&, const std::vector<std::string> &);
template RowCompressed<std::string>::RowCompressed(std::vector<Vector<uint64_t>>&&, const std::vector<std::string> &);

template <typename Label>
void RowCompressed<Label>::reinitialize(uint64_t num_rows) {
    if (dynamic_cast<EigenSpMat*>(matrix_.get())) {
        matrix_.reset(new EigenSpMat(num_rows));
    } else {
        matrix_.reset(new VectorRowBinMat<>(num_rows));
    }

    label_encoder_.clear();
}

template <typename Label>
void RowCompressed<Label>::set(Index i, const VLabels &labels) {
    assert(i < matrix_->num_rows());

    matrix_->clear_row(i);

    for (const auto &label : labels) {
        matrix_->set(i, label_encoder_.insert_and_encode(label));
    }
}

template <typename Label>
void RowCompressed<Label>::add_labels(const std::vector<Index> &indices,
                                      const VLabels &labels) {
    std::vector<uint64_t> columns(labels.size());
    for (size_t t = 0; t < labels.size(); ++t) {
        columns[t] = label_encoder_.insert_and_encode(labels[t]);
    }

    for (Index i : indices) {
        for (size_t j : columns) {
            matrix_->set(i, j);
        }
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

    for (Index i : unique_indices) {
        for (auto j : col_ids) {
            matrix_->force_set(i, j);
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

        if (!dynamic_cast<VectorRowBinMat<>*>(matrix_.get())) {
            std::cerr << "Error: loading from multiple row annotators is supported"
                      << " only for the VectorRowBinMat representation" << std::endl;
            exit(1);
        }

        auto &matrix = dynamic_cast<VectorRowBinMat<>&>(*matrix_);
        auto next_block = std::make_unique<VectorRowBinMat<>>(matrix_->num_rows());

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
            for (uint64_t i = 0, n_rows = next_block->num_rows(); i < n_rows; ++i) {
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

template <typename Label>
uint64_t RowCompressed<Label>::num_objects() const {
    return matrix_->num_rows();
}

template <typename Label>
uint64_t RowCompressed<Label>::num_relations() const {
    return matrix_->num_relations();
}

template <typename Label>
std::unique_ptr<LabelEncoder<Label>>
RowCompressed<Label>::load_label_encoder(const std::string &filebase) {
    auto filename = remove_suffix(filebase, kExtension) + kExtension;

    std::ifstream instream(filename, std::ios::binary);
    if (!instream.good())
        throw std::ifstream::failure("Cannot read from file " + filename);

    try {
        return load_label_encoder(instream);
    } catch (...) {
        throw std::ifstream::failure("Cannot load label encoder from file " + filename);
    }
}

template <typename Label>
std::unique_ptr<LabelEncoder<Label>>
RowCompressed<Label>::load_label_encoder(std::istream &instream) {
    if (!instream.good())
        throw std::istream::failure("Bad stream");

    auto label_encoder = std::make_unique<LabelEncoder<Label>>();
    if (!label_encoder->load(instream))
        throw std::istream::failure("Cannot load label encoder");

    return label_encoder;
}

template <typename Label>
void RowCompressed<Label>::load_shape(const std::string &filename,
                                      uint64_t *num_objects_,
                                      uint64_t *num_relations_) {
    assert(num_objects_);
    assert(num_relations_);
    uint64_t &num_objects = *num_objects_;
    uint64_t &num_relations = *num_relations_;
    num_objects = 0;
    num_relations = 0;

    auto sr = get_row_streamer(filename);

    while (auto *row = sr.next_row()) {
        num_objects++;
        num_relations += row->size();
    }
}

template <typename Label>
StreamRows<BinaryMatrix::SetBitPositions>
RowCompressed<Label>::get_row_streamer(const std::string &filebase) {
    std::string filename = remove_suffix(filebase, kExtension) + kExtension;
    std::ifstream instream(filename, std::ios::binary);
    // skip header
    load_label_encoder(instream);
    // rows
    return StreamRows<BinaryMatrix::SetBitPositions>(filename, instream.tellg());
}

template <typename Label>
void RowCompressed<Label>
::serialize(const std::string &filebase,
            const LabelEncoder<Label> &label_encoder,
            const std::function<void(BinaryMatrix::RowCallback)> &call_rows) {
    auto filename = remove_suffix(filebase, kExtension) + kExtension;

    std::ofstream outstream(filename, std::ios::binary);
    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    label_encoder.serialize(outstream);
    outstream.close();

    binmat::append_row_major(filename, call_rows, label_encoder.size());
}

template class RowCompressed<std::string>;

} // namespace annot
} // namespace mtg
