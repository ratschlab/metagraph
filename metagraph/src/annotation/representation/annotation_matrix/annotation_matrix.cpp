#include "annotation_matrix.hpp"

#include "common/utils/string_utils.hpp"
#include "common/serialization.hpp"
#include "static_annotators_def.hpp"


namespace mtg {
namespace anno {

using utils::remove_suffix;


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
    for (uint64_t col_id : get_matrix().get_row(i)) {
        encoded_labels.insert(col_id);
    }
    return std::includes(encoded_labels.begin(), encoded_labels.end(),
                         querying_codes.begin(), querying_codes.end());
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
uint64_t StaticBinRelAnnotator<BinaryMatrixType, Label>::num_relations() const {
    return matrix_->num_relations();
}

template <class BinaryMatrixType, typename Label>
void StaticBinRelAnnotator<BinaryMatrixType, Label>::except_dyn() {
    throw std::runtime_error("Dynamic actions are not supported"
                             " in static representation");
}

template <class BinaryMatrixType, typename Label>
bool StaticBinRelAnnotator<BinaryMatrixType, Label>
::dump_columns(const std::string &prefix, uint64_t num_threads) const {
    size_t m = this->num_labels();
    bool success = true;

    #pragma omp parallel for num_threads(num_threads)
    for (uint64_t i = 0; i < m; ++i) {
        std::ofstream outstream(
            remove_suffix(prefix, kExtension)
                + "." + std::to_string(i)
                + ".text.annodbg"
        );

        if (!outstream.good()) {
            std::cerr << "ERROR: dumping column " << i << " failed" << std::endl;
            success = false;
            continue;
        }

        outstream << num_objects() << " ";

        auto column = matrix_->get_column(i);

        outstream << column.size() << "\n";
        for (auto pos : column) {
            outstream << pos << "\n";
        }
    }

    return success;
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                        const LabelEncoder<Label> &label_encoder)
      : matrix_(std::move(matrix)) {
    assert(matrix_.get());
    label_encoder_ = label_encoder;
}

template class StaticBinRelAnnotator<binmat::RowConcatenated<>, std::string>;

template class StaticBinRelAnnotator<binmat::Rainbowfish, std::string>;

template class StaticBinRelAnnotator<binmat::BRWT, std::string>;

template class StaticBinRelAnnotator<binmat::BinRelWT_sdsl, std::string>;

template class StaticBinRelAnnotator<binmat::BinRelWT, std::string>;

template class StaticBinRelAnnotator<UniqueRowBinmat, std::string>;

template class StaticBinRelAnnotator<Rainbow<BRWT>, std::string>;

} // namespace anno
} // namespace mtg
