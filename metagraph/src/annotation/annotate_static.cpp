#include "annotate_static.hpp"

#include "utils.hpp"
#include "static_annotators_def.hpp"

using utils::remove_suffix;


namespace annotate {

template <>
const std::string RowFlatAnnotator::kFileExtension = kRowPackedExtension;
template <>
const std::string RainbowfishAnnotator::kFileExtension = kRainbowfishExtension;
template <>
const std::string BRWTCompressed<std::string>::kFileExtension = kBRWTExtension;
template <>
const std::string BinRelWT_sdslAnnotator::kFileExtension = kBinRelWT_sdslExtension;
template <>
const std::string BinRelWTAnnotator::kFileExtension = kBinRelWTExtension;


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
    for (auto col : matrix_->get_row(i)) {
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
    for (auto col : matrix_->get_row(i)) {
        labels.push_back(label_encoder_.decode(col));
    }
    return labels;
}

// Get labels that occur at least in |presence_ratio| rows.
// If |presence_ratio| = 0, return all occurring labels.
template <class BinaryMatrixType, typename Label>
auto
StaticBinRelAnnotator<BinaryMatrixType, Label>
::get_labels(const std::vector<Index> &indices,
             double presence_ratio) const -> VLabels {
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

template <class BinaryMatrixType, typename Label>
void
StaticBinRelAnnotator<BinaryMatrixType, Label>
::serialize(const std::string &filename) const {
    std::ofstream outstream(remove_suffix(filename, kFileExtension)
                                                         + kFileExtension,
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

    std::ifstream instream(remove_suffix(filenames.at(0), kFileExtension)
                                                            + kFileExtension,
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

// Count all labels collected from the given rows.
template <class BinaryMatrixType, typename Label>
std::vector<uint64_t>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::count_labels(const std::vector<Index> &indices) const {
    std::vector<uint64_t> counter(num_labels(), 0);

    for (Index i : indices) {
        for (auto col : matrix_->get_row(i)) {
            counter[col]++;
        }
    }

    return counter;
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>::except_dyn() {
    throw std::runtime_error("Dynamic actions are not supported"
                             " in static representation");
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                        const LabelEncoder<Label> &label_encoder) {
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
