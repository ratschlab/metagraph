#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"
#include "utils.hpp"

using utils::remove_suffix;


namespace annotate {

class VectorVectorMatrix : public RowMajorSparseBinaryMatrix {
  public:
    VectorVectorMatrix(size_t num_rows) : vector_(num_rows) {}

    void set_bit(size_t i, size_t j) {
        assert(i < vector_.size());
        if (!is_set_bit(i, j))
            vector_[i].push_back(j);
    }
    bool is_set_bit(size_t i, size_t j) const {
        assert(i < vector_.size());
        return std::find(vector_[i].begin(), vector_[i].end(), j)
                    != vector_[i].end();
    }

    size_t select(size_t i, size_t k) const { return vector_[i][k]; }

    size_t size() const { return vector_.size(); }
    size_t size(size_t i) const { return vector_[i].size(); }

    void clear(size_t i) { vector_[i].clear(); }

    void reinitialize(size_t num_rows) {
        vector_ = decltype(vector_)();
        vector_.resize(num_rows);
    }

    void insert_rows(const std::vector<uint64_t> &rows) {
        assert(std::is_sorted(rows.begin(), rows.end()));

        utils::insert_default_values(rows, &vector_);
    }

  private:
    std::vector<SmallVector> vector_;
};

class EigenSparserMatrix : public RowMajorSparseBinaryMatrix {
    typedef Eigen::SparseMatrix<bool, Eigen::RowMajor> SparseMatrix;

  public:
    EigenSparserMatrix(size_t num_rows) { reinitialize(num_rows); }

    void set_bit(size_t i, size_t j) { mat_.coeffRef(i, j) = 1; }

    bool is_set_bit(size_t i, size_t j) const { return mat_.coeff(i, j); }

    size_t select(size_t i, size_t k) const {
        return mat_.innerIndexPtr()[mat_.outerIndexPtr()[i] + k];
    }

    size_t size() const { return mat_.rows(); }
    size_t size(size_t i) const { return mat_.innerVector(i).nonZeros(); }

    void clear(size_t i) {
        for (SparseMatrix::InnerIterator it(mat_, i); it; ++it) {
            it.valueRef() = 0;
        }
    }

    void reinitialize(size_t num_rows) {
        mat_ = Eigen::SparseMatrix<bool, Eigen::RowMajor>(num_rows, kMaxNumCols);
        mat_.reserve(num_rows);
    }

    void insert_rows(const std::vector<uint64_t>&) {
        throw std::runtime_error("Error: Not implemented");
    }

  private:
    const size_t kMaxNumCols = 10'000'000;

    SparseMatrix mat_;
};


template <typename Label>
RowCompressed<Label>::RowCompressed(uint64_t num_rows, bool sparse)  {
    if (sparse) {
        matrix_.reset(new EigenSparserMatrix(num_rows));
    } else {
        matrix_.reset(new VectorVectorMatrix(num_rows));
    }
}

template <typename Label>
void RowCompressed<Label>::set_labels(Index i, const VLabels &labels) {
    assert(i < matrix_->size());

    matrix_->clear(i);
    add_labels(i, labels);
}

template <typename Label>
typename RowCompressed<Label>::VLabels
RowCompressed<Label>::get_labels(Index i) const {
    VLabels labels;
    for (size_t k = 0; k < matrix_->size(i); ++k) {
        labels.push_back(label_encoder_.decode(matrix_->select(i, k)));
    }
    return labels;
}

template <typename Label>
void RowCompressed<Label>::add_label(Index i, const Label &label) {
    matrix_->set_bit(i, label_encoder_.insert_and_encode(label));
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
bool RowCompressed<Label>::has_label(Index i, const Label &label) const {
    try {
        return matrix_->is_set_bit(i, label_encoder_.encode(label));
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
    for (size_t k = 0; k < matrix_->size(i); ++k) {
        encoded_labels.insert(matrix_->select(i, k));
    }
    return std::includes(encoded_labels.begin(), encoded_labels.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <typename Label>
void RowCompressed<Label>::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    serialize_number(outstream, matrix_->size());

    label_encoder_.serialize(outstream);

    uint64_t dump_vector_size = 0;
    for (size_t i = 0; i < matrix_->size(); ++i) {
        dump_vector_size += matrix_->size(i) + 1;
    }
    sdsl::int_vector<> full_vector(dump_vector_size, 0,
                                   std::log2(label_encoder_.size() + 1) + 1);

    for (uint64_t i = 0, p = 0; i < matrix_->size(); ++i) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            assert(matrix_->select(i, k) + 1 < label_encoder_.size() + 1);

            full_vector[p++] = matrix_->select(i, k) + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

template <typename Label>
bool RowCompressed<Label>::merge_load(const std::vector<std::string> &filenames) {
    std::ifstream instream(remove_suffix(filenames.at(0), kExtension) + kExtension);
    if (!instream.good())
        return false;

    try {
        size_t num_rows = load_number(instream);
        matrix_->reinitialize(num_rows);

        if (!label_encoder_.load(instream))
            return false;

        sdsl::int_vector<> full_vector;
        full_vector.load(instream);

        for (size_t k = 0, i = 0; k < full_vector.size(); ++k) {
            if (full_vector[k]) {
                matrix_->set_bit(i, full_vector[k] - 1);
            } else {
                i++;
            }
        }

        return true;
    } catch (...) {
        return false;
    }
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

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            encoded_counter[matrix_->select(i, k)]++;
        }
    }

    VLabels filtered_labels;

    for (auto it = encoded_counter.begin(); it != encoded_counter.end(); ++it) {
        if (it->second >= min_labels_discovered)
            filtered_labels.push_back(label_encoder_.decode(it->first));
    }

    return filtered_labels;
}

// Count all labels collected from the given rows.
template <typename Label>
std::vector<uint64_t>
RowCompressed<Label>::count_labels(const std::vector<Index> &indices) const {
    std::vector<uint64_t> counter(num_labels(), 0);

    for (Index i : indices) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            counter[matrix_->select(i, k)]++;
        }
    }

    return counter;
}

template <typename Label>
uint64_t RowCompressed<Label>::num_objects() const {
    return matrix_->size();
}

template <typename Label>
size_t RowCompressed<Label>::num_labels() const {
    return label_encoder_.size();
}

template <typename Label>
double RowCompressed<Label>::sparsity() const {
    uint64_t num_set_bits = 0;

    for (uint64_t i = 0; i < matrix_->size(); ++i) {
        num_set_bits += matrix_->size(i);
    }

    return 1 - static_cast<double>(num_set_bits) / label_encoder_.size()
                                                 / matrix_->size();
}

template class RowCompressed<std::string>;

} // namespace annotate
