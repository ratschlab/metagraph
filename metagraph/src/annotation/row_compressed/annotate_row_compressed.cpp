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
        matrix_.reset(new VectorRowBinMat<>(num_rows));
    }
}

template <typename Label>
RowCompressed<Label>::RowCompressed(uint64_t num_rows,
                                    const std::vector<Label> &labels,
                                    std::function<void(CallRow)> call_rows) {
    for (const auto &label : labels) {
        label_encoder_.insert_and_encode(label);
    }

    matrix_.reset(new VectorRowBinMat<>(num_rows, labels.size(), [&](auto call_row) {
        call_rows([&](Index i, const VLabels &row_labels) {
            VectorRowBinMat<>::row_type row(row_labels.size());
            for (size_t j = 0; j < row_labels.size(); ++j) {
                row[j] = label_encoder_.encode(row_labels[j]);
            }
            call_row(i, std::move(row));
        });
    }));
}

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
void RowCompressed<Label>::set_labels(Index i, const VLabels &labels) {
    assert(i < matrix_->num_rows());

    matrix_->clear_row(i);
    add_labels(i, labels);
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

template <typename Label>
void RowCompressed<Label>
::call_objects(const Label &label,
               std::function<void(Index)> callback) const {
    uint64_t col;
    try {
        col = label_encoder_.encode(label);
    } catch (...) {
        return;
    }

    for (Index index : matrix_->get_column(col)) {
        callback(index);
    }
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
void RowCompressed<Label>::stream_counts(const std::string &filename,
                                         uint64_t *num_objects_,
                                         uint64_t *num_relations_) {
    assert(num_objects_);
    assert(num_relations_);
    uint64_t &num_objects = *num_objects_;
    uint64_t &num_relations = *num_relations_;
    num_objects = 0;
    num_relations = 0;

    StreamRows<> sr(filename);

    while (auto *row = sr.next_row()) {
        num_objects++;
        num_relations += row->size();
    }
}

template <typename Label>
template <typename RowType>
RowCompressed<Label>::StreamRows<RowType>::StreamRows(std::string filename) {
    filename = remove_suffix(filename, kExtension) + kExtension;
    std::ifstream instream(filename, std::ios::binary);
    // skip header
    load_label_encoder(instream);
    // rows
    sr_ = std::make_unique<::StreamRows<RowType>>(filename, instream.tellg());
}

template <typename Label>
void RowCompressed<Label>
::write_rows(std::string filename,
             const LabelEncoder<Label> &label_encoder,
             const std::function<void(BinaryMatrix::RowCallback)> &call_rows) {
    filename = remove_suffix(filename, kExtension) + kExtension;

    std::ofstream outstream(filename, std::ios::binary);
    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    label_encoder.serialize(outstream);
    outstream.close();

    append_row_major(filename, call_rows, label_encoder.size());
}

template class RowCompressed<std::string>;

} // namespace annotate
