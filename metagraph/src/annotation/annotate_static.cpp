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
    for (auto col : get_label_codes(i)) {
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
    for (auto col : get_label_codes(i)) {
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
::get_label_codes(Index i) const {
    if (!cached_rows_.get())
        return matrix_->get_row(i);

    try {
        return cached_rows_->Get(i);
    } catch (...) {
        auto row = matrix_->get_row(i);
        cached_rows_->Put(i, row);
        return row;
    }
}

template <class BinaryMatrixType, typename Label>
bool StaticBinRelAnnotator<BinaryMatrixType, Label>
::dump_columns(const std::string &prefix, bool binary, uint64_t num_threads) const {
    size_t m = num_labels();
    bool success = true;

    if (binary) {
        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = 0; i < m; ++i) {
            std::ofstream outstream(
                remove_suffix(prefix, kExtension)
                    + "." + std::to_string(i)
                    + ".raw.annodbg",
                std::ios::binary
            );

            if (!outstream.good()) {
                std::cerr << "ERROR: dumping column " << i << " failed" << std::endl;
                success = false;
                continue;
            }

            serialize_number(outstream, num_objects());

            auto column = matrix_->get_column(i);

            serialize_number(outstream, column.size());
            for (auto pos : column) {
                serialize_number(outstream, pos);
            }
        }
    } else {
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

            outstream << column.size() << std::endl;
            for (auto pos : column) {
                outstream << pos << std::endl;
            }
        }
    }

    return success;
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(size_t row_cache_size) : matrix_(new BinaryMatrixType()) {
    MultiLabelEncoded<uint64_t, Label>::reset_row_cache(row_cache_size);
}

template <class BinaryMatrixType, typename Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>
::StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                        const LabelEncoder<Label> &label_encoder,
                        size_t row_cache_size) : matrix_(std::move(matrix)) {
    assert(matrix_.get());
    MultiLabelEncoded<uint64_t, Label>::reset_row_cache(row_cache_size);
    label_encoder_ = label_encoder;
}

template class StaticBinRelAnnotator<RowConcatenated<>, std::string>;

template class StaticBinRelAnnotator<Rainbowfish, std::string>;

template class StaticBinRelAnnotator<BRWT, std::string>;

template class StaticBinRelAnnotator<BinRelWT_sdsl, std::string>;

template class StaticBinRelAnnotator<BinRelWT, std::string>;

} // namespace annotate
