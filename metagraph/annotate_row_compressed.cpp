#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#if _USE_FOLLY
#include <folly/small_vector.h>
#endif

#include "serialization.hpp"
#include "utils.hpp"

using libmaus2::util::StringSerialisation;
using utils::remove_suffix;

#if _USE_FOLLY
typedef folly::small_vector<uint32_t, 2, uint32_t> SmallVector;
#else
typedef std::vector<uint32_t> SmallVector;
#endif


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

        std::vector<SmallVector> old;
        old.swap(vector_);
        // vector_.reserve(old.size() + rows.size());

        size_t i = 0;
        size_t num_inserted = 0;

        for (auto next_inserted_row : rows) {
            while (i + num_inserted < next_inserted_row) {
                assert(i < old.size() && "Invalid indexes of inserted rows");
                vector_.push_back(std::move(old[i++]));
            }
            vector_.push_back(SmallVector());
            num_inserted++;
        }
        while (i < old.size()) {
            vector_.push_back(std::move(old[i++]));
        }
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


template <typename Color, class Encoder>
RowCompressed<Color, Encoder>::RowCompressed(uint64_t num_rows, bool sparse)
      : color_encoder_(new Encoder()) {
    if (sparse) {
        matrix_.reset(new EigenSparserMatrix(num_rows));
    } else {
        matrix_.reset(new VectorVectorMatrix(num_rows));
    }
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::set_coloring(Index i, const Coloring &coloring) {
    assert(i < matrix_->size());

    matrix_->clear(i);
    for (const auto &color : coloring) {
        matrix_->set_bit(i, color_encoder_->encode(color, true));
    }
}

template <typename Color, class Encoder>
typename RowCompressed<Color, Encoder>::Coloring
RowCompressed<Color, Encoder>::get_coloring(Index i) const {
    Coloring coloring;
    for (size_t k = 0; k < matrix_->size(i); ++k) {
        coloring.push_back(color_encoder_->decode(matrix_->select(i, k)));
    }
    return coloring;
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_color(Index i, const Color &color) {
    matrix_->set_bit(i, color_encoder_->encode(color, true));
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_colors(Index i, const Coloring &coloring) {
    for (const auto &color : coloring) {
        add_color(i, color);
    }
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::add_colors(const std::vector<Index> &indices,
                                               const Coloring &coloring) {
    for (Index i : indices) {
        add_colors(i, coloring);
    }
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::has_color(Index i, const Color &color) const {
    try {
        return matrix_->is_set_bit(i, color_encoder_->encode(color));
    } catch (...) {
        return false;
    }
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::has_colors(Index i, const Coloring &coloring) const {
    std::set<size_t> querying_codes;
    try {
        for (const auto &color : coloring) {
            querying_codes.insert(color_encoder_->encode(color));
        }
    } catch (...) {
        return false;
    }
    std::set<size_t> encoded_coloring;
    for (size_t k = 0; k < matrix_->size(i); ++k) {
        encoded_coloring.insert(matrix_->select(i, k));
    }
    return std::includes(encoded_coloring.begin(), encoded_coloring.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }

    serialize_number(outstream, matrix_->size());

    color_encoder_->serialize(outstream);

    uint64_t dump_vector_size = 0;
    for (size_t i = 0; i < matrix_->size(); ++i) {
        dump_vector_size += matrix_->size(i) + 1;
    }
    sdsl::int_vector<> full_vector(dump_vector_size, 0,
                                   std::log2(color_encoder_->size() + 1) + 1);

    for (uint64_t i = 0, p = 0; i < matrix_->size(); ++i) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            assert(matrix_->select(i, k) + 1 < color_encoder_->size() + 1);

            full_vector[p++] = matrix_->select(i, k) + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

template <typename Color, class Encoder>
bool RowCompressed<Color, Encoder>::merge_load(const std::vector<std::string> &filenames) {
    std::ifstream instream(remove_suffix(filenames.at(0), kExtension) + kExtension);
    if (!instream.good())
        return false;

    try {
        size_t num_rows = load_number(instream);
        matrix_->reinitialize(num_rows);

        if (!color_encoder_->load(instream))
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

template <typename Color, class Encoder>
void RowCompressed<Color, Encoder>::insert_rows(const std::vector<Index> &rows) {
    matrix_->insert_rows(rows);
}

// Get colors that occur at least in |discovery_ratio| colorings.
// If |discovery_ratio| = 0, return the union of colorings.
template <typename Color, class Encoder>
typename RowCompressed<Color, Encoder>::Coloring
RowCompressed<Color, Encoder>::aggregate_colors(const std::vector<Index> &indices,
                                                double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    // const size_t max_colors_missing = indices.size() - min_colors_discovered;

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            encoded_counter[matrix_->select(i, k)]++;
        }
    }

    Coloring filtered_colors;

    for (auto it = encoded_counter.begin(); it != encoded_counter.end(); ++it) {
        if (it->second >= min_colors_discovered)
            filtered_colors.push_back(color_encoder_->decode(it->first));
    }

    return filtered_colors;
}

// Count all colors collected from extracted colorings
// and return top |num_top| with the counts computed.
template <typename Color, class Encoder>
std::vector<std::pair<Color, size_t>>
RowCompressed<Color, Encoder>::get_most_frequent_colors(const std::vector<Index> &indices,
                                                        size_t num_top) const {
    std::vector<size_t> encoded_counter(color_encoder_->size(), 0);

    for (Index i : indices) {
        for (size_t k = 0; k < matrix_->size(i); ++k) {
            encoded_counter[matrix_->select(i, k)]++;
        }
    }

    std::vector<std::pair<size_t, size_t>> counts;
    counts.reserve(encoded_counter.size());
    for (size_t j = 0; j < encoded_counter.size(); ++j) {
        if (encoded_counter[j])
            counts.emplace_back(j, encoded_counter[j]);
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
size_t RowCompressed<Color, Encoder>::num_colors() const {
    return color_encoder_->size();
}

template <typename Color, class Encoder>
double RowCompressed<Color, Encoder>::sparsity() const {
    uint64_t num_set_bits = 0;

    for (uint64_t i = 0; i < matrix_->size(); ++i) {
        num_set_bits += matrix_->size(i);
    }

    return 1 - static_cast<double>(num_set_bits) / color_encoder_->size()
                                                 / matrix_->size();
}

template class RowCompressed<std::string, StringEncoder>;

} // namespace annotate
