#include "annotate_static.hpp"

#include "utils.hpp"
#include "static_annotators_def.hpp"

using utils::remove_suffix;


namespace annotate {

template <>
const std::string RowFlatAnnotator::kExtension = kRowPackedExtension;
template <>
const std::string RainbowfishAnnotator::kExtension = kRainbowfishExtension;
template <>
const std::string BRWTCompressed<std::string>::kExtension = kBRWTExtension;
template <>
const std::string BinRelWT_sdslAnnotator::kExtension = kBinRelWT_sdslExtension;
template <>
const std::string BinRelWTAnnotator::kExtension = kBinRelWTExtension;

template <class BinaryMatrixType, typename Label>
std::string
StaticBinRelAnnotator<BinaryMatrixType, Label>
::file_extension() const { return kExtension; }

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::has_label(Index i, const Label &label) const {
    assert(i < num_objects());
    try {
        return matrix_->get(i, label_encoder_.encode(label));
    } catch (...) {
        return false;
    }
}

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::has_labels(Index i, const VLabels &labels) const {
    assert(i < num_objects());

    std::set<size_t> querying_codes;
    try {
        for (const auto &label : labels) {
            querying_codes.insert(label_encoder_.encode(label));
        }
    } catch (...) {
        return false;
    }
    std::set<size_t> encoded_labels;
    for (auto col : get_label_indices(i)) {
        encoded_labels.insert(col);
    }
    return std::includes(encoded_labels.begin(), encoded_labels.end(),
                         querying_codes.begin(), querying_codes.end());
}

template <class BinaryMatrixType, typename Label>
auto StaticBinRelAnnotator<BinaryMatrixType, Label>::get_labels(Index i) const
-> VLabels {
    assert(i < num_objects());

    VLabels labels;
    for (auto col : get_label_indices(i)) {
        labels.push_back(label_encoder_.decode(col));
    }
    return labels;
}

template <class BinaryMatrixType, typename Label>
void
StaticBinRelAnnotator<BinaryMatrixType, Label>
::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kExtension) + kExtension,
                            std::ios::binary);
    if (!outstream.good()) {
        throw std::ofstream::failure("Bad stream");
    }
    label_encoder_.serialize(outstream);
    matrix_->serialize(outstream);
}

template <class BinaryMatrixType, typename Label>
bool
StaticBinRelAnnotator<BinaryMatrixType, Label>
::merge_load(const std::vector<std::string> &filenames) {
    if (filenames.size() > 1)
        std::cerr << "Warning: Can't merge static annotators."
                     " Only the first will be loaded." << std::endl;

    std::ifstream instream(remove_suffix(filenames.at(0), kExtension) + kExtension,
                           std::ios::binary);
    if (!instream.good())
        return false;

    try {
        assert(matrix_.get());
        return label_encoder_.load(instream) && matrix_->load(instream);
    } catch (...) {
        return false;
    }
}

template <class BinaryMatrixType, typename Label>
uint64_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_objects() const {
    return matrix_->num_rows();;
}

template <class BinaryMatrixType, typename Label>
size_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_labels() const {
    assert(label_encoder_.size() == matrix_->num_columns());
    return label_encoder_.size();
}

template <class BinaryMatrixType, typename Label>
uint64_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_relations() const {
    return matrix_->num_relations();
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>
::call_objects(const Label &label,
               std::function<void(Index)> callback) const {
    uint64_t encoding;
    try {
        encoding = label_encoder_.encode(label);
    } catch (...) {
        return;
    }

    for (Index index : matrix_->get_column(encoding)) {
        callback(index);
    }
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>::except_dyn() {
    throw std::runtime_error("Dynamic actions are not supported"
                             " in static representation");
}

template <class BinaryMatrixType, typename Label>
std::vector<uint64_t> StaticBinRelAnnotator<BinaryMatrixType, Label>
::get_label_indices(Index i) const {
    if (cached_rows_.get()) {
        try {
            return cached_rows_->Get(i);
        } catch (...) {
            auto row = matrix_->get_row(i);
            cached_rows_->Put(i, row);
            return row;
        }
    } else {
        return matrix_->get_row(i);
    }
}

// For each index i in indices, check of i has the label. Return
// true if the finished callback evaluates true during execution.
template <class BinaryMatrixType, typename Label>
bool StaticBinRelAnnotator<BinaryMatrixType, Label>
::call_indices_until(const std::vector<Index> &indices,
                     const Label &label,
                     std::function<void(Index)> index_callback,
                     std::function<bool()> finished) const {
    auto column = matrix_->get_column(label_encoder_.encode(label));

    std::unordered_set<Index> column_set(column.begin(), column.end());
    for (Index i : indices) {
        if (column_set.find(i) != column_set.end()) {
            index_callback(i);

            if (finished())
                return true;
        }
    }

    return false;
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>
::reset_row_cache(size_t size) {
    cached_rows_.reset(size ? new RowCacheType(size) : nullptr);
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>::StaticBinRelAnnotator(size_t row_cache_size)
      : matrix_(new BinaryMatrixType()),
        cached_rows_({ row_cache_size ? new RowCacheType(row_cache_size) : nullptr }) {}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                        const LabelEncoder<Label> &label_encoder,
                        size_t row_cache_size)
      : cached_rows_({ row_cache_size ? new RowCacheType(row_cache_size) : nullptr }) {
    assert(matrix.get());
    matrix_ = std::move(matrix);
    label_encoder_ = label_encoder;
}

template class StaticBinRelAnnotator<RowConcatenated<>, std::string>;

template class StaticBinRelAnnotator<Rainbowfish, std::string>;

template class StaticBinRelAnnotator<BRWT, std::string>;

template class StaticBinRelAnnotator<BinRelWT_sdsl, std::string>;

template class StaticBinRelAnnotator<BinRelWT, std::string>;

} // namespace annotate
