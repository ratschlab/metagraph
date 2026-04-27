#include "annotation_matrix.hpp"

#include <filesystem>

#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/serialization.hpp"
#include "static_annotators_def.hpp"

namespace mtg {
namespace annot {

using namespace mtg::annot::matrix;

using utils::remove_suffix;
using utils::make_suffix;
using mtg::common::logger;


template <class BinaryMatrixType, typename Label>
void
StaticBinRelAnnotator<BinaryMatrixType, Label>
::serialize(const std::string &filename) const {
    const auto &fname = make_suffix(filename, kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + fname);
    label_encoder_.serialize(out);
    matrix_->serialize(out);
}

template <class BinaryMatrixType, typename Label>
bool StaticBinRelAnnotator<BinaryMatrixType, Label>::load(const std::string &filename) {
    const auto &fname = make_suffix(filename, kExtension);
    std::unique_ptr<std::ifstream> in = utils::open_ifstream(fname, utils::with_mmap());
    if (!in->good()) {
        logger->error("Cannot open annotation file '{}': {}", fname,
                      utils::file_read_failure_detail(fname));
        return false;
    }

    try {
        assert(matrix_.get());
        if (!label_encoder_.load(*in)) {
            logger->error("Cannot load annotation label encoder from '{}': {}", fname,
                          utils::file_read_failure_detail(fname));
            return false;
        }
        if (!matrix_->load(*in)) {
            logger->error("Cannot load annotation matrix from '{}': {}", fname,
                          utils::file_read_failure_detail(fname));
            return false;
        }
        return true;
    } catch (const std::exception &e) {
        logger->error("Cannot load annotation from '{}': {} (caught: {})", fname,
                      utils::file_read_failure_detail(fname), e.what());
        return false;
    } catch (...) {
        logger->error("Cannot load annotation from '{}': {} (caught unknown exception)",
                      fname, utils::file_read_failure_detail(fname));
        return false;
    }
}

template <class BinaryMatrixType, typename Label>
LabelEncoder<Label>
StaticBinRelAnnotator<BinaryMatrixType, Label>::read_label_encoder(const std::string &filename) {
    const auto &fname = make_suffix(filename, kExtension);
    std::ifstream in(fname, std::ios::binary);
    LabelEncoder<Label> label_encoder;
    if (!label_encoder.load(in))
        throw std::ofstream::failure("Can't load label encoder from " + fname);
    return label_encoder;
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
        std::ofstream outstream(remove_suffix(prefix, kExtension)
                                    + fmt::format(".{}.text.annodbg", i));
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

// TODO: make row-diff annotations .column.annodbg format and remove this function
bool merge_load_row_diff(const std::vector<std::string> &filenames,
                         const ColumnCallback &callback,
                         size_t num_threads) {
    std::atomic<bool> error_occurred = false;

    std::vector<uint64_t> offsets(filenames.size(), 0);

    // load labels
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 1; i < filenames.size(); ++i) {
        try {
            offsets[i] = RowDiffColumnAnnotator::read_label_encoder(filenames[i - 1]).size();
        } catch (...) {
            common::logger->error("Can't load label encoder from {}", filenames[i - 1]);
            error_occurred = true;
        }
    }

    // compute global offsets (partial sums)
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());

    // load annotations
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (size_t i = 0; i < filenames.size(); ++i) {
        if (error_occurred)
            continue;

        try {
            common::logger->trace("Loading annotations from file {}", filenames[i]);
            LabelEncoder<std::string> label_encoder;
            std::unique_ptr<std::ifstream> in = utils::open_ifstream(filenames[i]);
            RowDiff<ColumnMajor> matrix;
            if (!label_encoder.load(*in) || !matrix.load(*in)) {
                common::logger->error("Can't load {}", filenames[i]);
                std::exit(1);
            }
            if (!label_encoder.size())
                common::logger->warn("No labels in {}", filenames[i]);

            std::vector<std::unique_ptr<bit_vector>> &cols = matrix.diffs().data();

            for (uint32_t j = 0; j < cols.size(); ++j) {
                callback(offsets[i] + j, label_encoder.decode(j), std::move(cols[j]));
            }

        } catch (...) {
            error_occurred = true;
        }
    }

    return !error_occurred;
}

template class StaticBinRelAnnotator<RowFlat<>, std::string>;

template class StaticBinRelAnnotator<Rainbowfish, std::string>;

template class StaticBinRelAnnotator<BRWT, std::string>;

template class StaticBinRelAnnotator<RowSparse, std::string>;

template class StaticBinRelAnnotator<BinRelWT, std::string>;

template class StaticBinRelAnnotator<UniqueRowBinmat, std::string>;

template class StaticBinRelAnnotator<Rainbow<BRWT>, std::string>;

template class StaticBinRelAnnotator<RowDiff<ColumnMajor>, std::string>;

template class StaticBinRelAnnotator<RowDiff<BRWT>, std::string>;

template class StaticBinRelAnnotator<RowDiff<RowFlat<>>, std::string>;

template class StaticBinRelAnnotator<RowDiff<RowSparse>, std::string>;

template class StaticBinRelAnnotator<RowDiff<RowDisk>, std::string>;

template class StaticBinRelAnnotator<IntRowDiff<IntRowDisk>, std::string>;

template class StaticBinRelAnnotator<TupleRowDiff<CoordRowDisk>, std::string>;

template class StaticBinRelAnnotator<CSCMatrix<BRWT, CountsVector>, std::string>;

template class StaticBinRelAnnotator<IntRowDiff<CSCMatrix<BRWT, CountsVector>>, std::string>;

template class StaticBinRelAnnotator<CSRMatrix, std::string>;

template class StaticBinRelAnnotator<TupleCSCMatrix<ColumnMajor>, std::string>;
template class StaticBinRelAnnotator<TupleCSCMatrix<BRWT>, std::string>;

template class StaticBinRelAnnotator<TupleRowDiff<TupleCSCMatrix<ColumnMajor>>, std::string>;
template class StaticBinRelAnnotator<TupleRowDiff<TupleCSCMatrix<BRWT>>, std::string>;

} // namespace annot
} // namespace mtg
