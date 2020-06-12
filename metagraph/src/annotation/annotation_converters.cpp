#include "annotation_converters.hpp"

#include <cassert>
#include <vector>
#include <functional>
#include <filesystem>

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "binary_matrix/row_vector/vector_row_binmat.hpp"
#include "binary_matrix/multi_brwt/brwt_builders.hpp"
#include "binary_matrix/multi_brwt/clustering.hpp"
#include "representation/annotation_matrix/static_annotators_def.hpp"
#include "representation/column_compressed/annotate_column_compressed.hpp"
#include "representation/row_compressed/annotate_row_compressed.hpp"

namespace annotate {

using namespace mtg;
using mtg::common::logger;

typedef LabelEncoder<std::string> LEncoder;

size_t kNumRowsInBlock = 50'000;


template <class RowCallback>
void call_rows(const BinaryMatrix &row_major_matrix,
               const RowCallback &callback,
               bool sort = false) {
    const auto num_rows = row_major_matrix.num_rows();

    for (size_t i = 0; i < num_rows; ++i) {
        auto row = row_major_matrix.get_row(i);
        if (sort)
            std::sort(row.begin(), row.end());
        callback(row);
    }
}


// RowCompressed -> other

template <>
std::unique_ptr<RowFlatAnnotator>
convert<RowFlatAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    ProgressBar progress_bar(num_rows, "Processing rows",
                             std::cerr, !common::get_verbose());

    if (const auto *mat = dynamic_cast<const VectorRowBinMat<>*>(&annotator.get_matrix()))
        const_cast<VectorRowBinMat<>*>(mat)->standardize_rows();

    auto matrix = std::make_unique<RowConcatenated<>>(
        [&](auto callback) {
            call_rows(annotator.get_matrix(),
                [&](const auto &row) {
                    assert(std::is_sorted(row.begin(), row.end()));
                    callback(row);
                    ++progress_bar;
                }
            );
        },
        num_columns,
        num_rows,
        num_set_bits
    );

    return std::make_unique<RowFlatAnnotator>(std::move(matrix),
                                              annotator.get_label_encoder());
}

template <>
std::unique_ptr<RainbowfishAnnotator>
convert<RainbowfishAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<Rainbowfish>([&](auto callback) {
        call_rows(annotator.get_matrix(), callback, true);
    }, num_columns);

    return std::make_unique<RainbowfishAnnotator>(std::move(matrix),
                                                  annotator.get_label_encoder());
}

template <>
std::unique_ptr<BinRelWT_sdslAnnotator>
convert<BinRelWT_sdslAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT_sdsl>(
        [&](auto callback) {
            call_rows(annotator.get_matrix(), callback);
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWT_sdslAnnotator>(std::move(matrix),
                                                    annotator.get_label_encoder());
}

template <>
std::unique_ptr<BinRelWTAnnotator>
convert<BinRelWTAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT>(
        [&](auto callback) {
            call_rows(annotator.get_matrix(), callback);
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWTAnnotator>(std::move(matrix),
                                               annotator.get_label_encoder());
}

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation> convert(const std::string &filename) {
    using MatrixType = typename StaticAnnotation::binary_matrix_type;
    using Label = typename StaticAnnotation::Label;

    auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);
    if (!label_encoder.get())
        throw std::iostream::failure("Cannot load label encoder");

    uint64_t num_rows;
    uint64_t num_relations;
    RowCompressed<Label>::load_shape(filename, &num_rows, &num_relations);

    constexpr size_t num_passes = std::is_same_v<MatrixType, Rainbowfish> ? 2u : 1u;
    ProgressBar progress_bar(num_rows * num_passes, "Processing rows",
                             std::cerr, !common::get_verbose());

    auto call_rows = [&](BinaryMatrix::RowCallback callback) {
        auto row_streamer = RowCompressed<Label>::get_row_streamer(filename);
        for (uint64_t r = 0; r < num_rows; ++r) {
            auto row = row_streamer.next_row();
            // TODO: remove sort?
            std::sort(row->begin(), row->end());
            callback(*row);
            ++progress_bar;
        }
    };

    std::unique_ptr<MatrixType> matrix;

    if constexpr(std::is_same_v<MatrixType, RowConcatenated<>>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder->size(), num_rows, num_relations);

    } else if constexpr(std::is_same_v<MatrixType, Rainbowfish>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder->size());

    } else if constexpr(std::is_same_v<MatrixType, BinRelWT>) {
        matrix = std::make_unique<MatrixType>(call_rows, num_relations, label_encoder->size());

    } else if constexpr(std::is_same_v<MatrixType, BinRelWT_sdsl>) {
        matrix = std::make_unique<MatrixType>(call_rows, num_relations, label_encoder->size());

    } else {
        static_assert(utils::dependent_false<StaticAnnotation>::value);
    }

    return std::make_unique<StaticAnnotation>(std::move(matrix), *label_encoder);
}

template std::unique_ptr<RowFlatAnnotator> convert(const std::string &filename);
template std::unique_ptr<RainbowfishAnnotator> convert(const std::string &filename);
template std::unique_ptr<BinRelWTAnnotator> convert(const std::string &filename);
template std::unique_ptr<BinRelWT_sdslAnnotator> convert(const std::string &filename);


// ColumnCompressed -> other

template <>
std::unique_ptr<RowFlatAnnotator>
convert<RowFlatAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowConcatenated<>>([&](auto callback) {
        utils::RowsFromColumnsTransformer(annotator.get_matrix().data()).call_rows(callback);
    }, num_columns, num_rows, num_set_bits);

    return std::make_unique<RowFlatAnnotator>(std::move(matrix),
                                              annotator.get_label_encoder());
}

template <>
std::unique_ptr<RainbowfishAnnotator>
convert<RainbowfishAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<Rainbowfish>([&](auto callback) {
        utils::RowsFromColumnsTransformer(annotator.get_matrix().data()).call_rows(callback);
    }, num_columns);

    return std::make_unique<RainbowfishAnnotator>(std::move(matrix),
                                                  annotator.get_label_encoder());
}

template <>
std::unique_ptr<UniqueRowAnnotator>
convert<UniqueRowAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<UniqueRowBinmat>([&](auto callback) {
        utils::RowsFromColumnsTransformer(annotator.get_matrix().data()).call_rows(callback);
    }, num_columns);

    return std::make_unique<UniqueRowAnnotator>(std::move(matrix),
                                                annotator.get_label_encoder());
}

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_BRWT(ColumnCompressed<Label>&& annotator,
                BRWTBottomUpBuilder::Partitioner partitioning,
                size_t num_parallel_nodes,
                size_t num_threads) {
    // we are going to take the columns from the annotator and thus
    // have to replace them with empty columns to keep the structure valid
    std::vector<std::unique_ptr<bit_vector>> columns(annotator.num_labels());
    columns.swap(const_cast<std::vector<std::unique_ptr<bit_vector>>&>(
        annotator.get_matrix().data()
    ));

    auto matrix = std::make_unique<BRWT>(
        BRWTBottomUpBuilder::build(std::move(columns),
                                   partitioning,
                                   num_parallel_nodes,
                                   num_threads)
    );

    return std::make_unique<StaticAnnotation>(std::move(matrix),
                                              annotator.get_label_encoder());
}

template <>
std::unique_ptr<MultiBRWTAnnotator>
convert_to_greedy_BRWT<MultiBRWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotation,
                                                        size_t num_parallel_nodes,
                                                        size_t num_threads,
                                                        uint64_t num_rows_subsampled) {
    return convert_to_BRWT<MultiBRWTAnnotator>(
        std::move(annotation),
        [num_threads,num_rows_subsampled](const auto &columns) {
            std::vector<sdsl::bit_vector> subvectors
                = random_submatrix(columns, num_rows_subsampled, num_threads);
            return greedy_matching(subvectors, num_threads);
        },
        num_parallel_nodes,
        num_threads
    );
}

template <>
std::unique_ptr<MultiBRWTAnnotator>
convert_to_simple_BRWT<MultiBRWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotation,
                                                        size_t grouping_arity,
                                                        size_t num_parallel_nodes,
                                                        size_t num_threads) {
    return convert_to_BRWT<MultiBRWTAnnotator>(
        std::move(annotation),
        BRWTBottomUpBuilder::get_basic_partitioner(grouping_arity),
        num_parallel_nodes,
        num_threads
    );
}

std::vector<std::vector<uint64_t>>
parse_linkage_matrix(const std::string &filename) {
    std::ifstream in(filename);

    std::vector<std::vector<uint64_t>> linkage;
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> parts = utils::split_string(line, " ");
        if (parts.empty())
            continue;

        try {
            if (parts.size() != 4)
                throw std::runtime_error("Invalid format");

            uint64_t first = std::stoi(parts.at(0));
            uint64_t second = std::stoi(parts.at(1));
            uint64_t merged = std::stoi(parts.at(3));

            if (first == second || first >= merged || second >= merged) {
                logger->error("Invalid format of the linkage matrix."
                              " Indexes of parent clusters must be larger than"
                              " indexes of the objects/clusters the include");
                exit(1);
            }

            while (linkage.size() <= merged) {
                linkage.push_back({});
            }

            linkage[merged].push_back(first);
            linkage[merged].push_back(second);

        } catch (const std::exception &e) {
            logger->error("Possibly invalid format of the linkage matrix."
                          " Each line must contsin exactly 4 values:"
                          " <cluster 1> <cluster 2> <dist> <cluster 3>"
                          "\nException: {}", e.what());
            exit(1);
        }
    }

    return linkage;
}

template <>
std::unique_ptr<MultiBRWTAnnotator>
convert_to_BRWT<MultiBRWTAnnotator>(const std::vector<std::string> &annotation_files,
                                    const std::string &linkage_matrix_file,
                                    size_t num_parallel_nodes,
                                    size_t num_threads,
                                    const std::filesystem::path &tmp_path) {
    std::vector<std::pair<uint64_t, std::string>> column_names;
    std::mutex mu;

    auto get_columns = [&](const BRWTBottomUpBuilder::CallColumn &call_column) {
        bool success = ColumnCompressed<>::merge_load(
            annotation_files,
            [&](uint64_t column_index,
                    const std::string &label,
                    std::unique_ptr<bit_vector>&& column) {
                call_column(column_index, std::move(column));
                std::lock_guard<std::mutex> lock(mu);
                column_names.emplace_back(column_index, label);
            },
            num_threads
        );
        if (!success) {
            logger->error("Can't load annotation columns");
            exit(1);
        }
    };

    auto linkage = parse_linkage_matrix(linkage_matrix_file);

    auto matrix = std::make_unique<BRWT>(
        BRWTBottomUpBuilder::build(get_columns, linkage, tmp_path,
                                   num_parallel_nodes, num_threads));

    std::sort(column_names.begin(), column_names.end(), utils::LessFirst());
    column_names.erase(std::unique(column_names.begin(), column_names.end()),
                       column_names.end());

    assert(matrix->num_columns() == column_names.size());

    LEncoder label_encoder;
    for (const auto &[i, label] : column_names) {
        label_encoder.insert_and_encode(label);
    }

    return std::make_unique<MultiBRWTAnnotator>(std::move(matrix), label_encoder);
}

template <>
void relax_BRWT<MultiBRWTAnnotator>(MultiBRWTAnnotator *annotation,
                                    size_t relax_max_arity,
                                    size_t num_threads) {
    if (relax_max_arity > 1)
        BRWTOptimizer::relax(const_cast<BRWT*>(&annotation->get_matrix()),
                             relax_max_arity,
                             num_threads);
}


template <>
std::unique_ptr<BinRelWT_sdslAnnotator>
convert<BinRelWT_sdslAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT_sdsl>(
        [&](auto callback) {
            utils::RowsFromColumnsTransformer(annotator.get_matrix().data()).call_rows(callback);
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWT_sdslAnnotator>(std::move(matrix),
                                                    annotator.get_label_encoder());
}

template <>
std::unique_ptr<BinRelWTAnnotator>
convert<BinRelWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    auto &columns = const_cast<std::vector<std::unique_ptr<bit_vector>>&>(
        annotator.get_matrix().data()
    );

    return std::make_unique<BinRelWTAnnotator>(
        std::make_unique<BinRelWT>(std::move(columns)),
        annotator.get_label_encoder()
    );
}

template <typename Label>
void merge_rows(const std::vector<const LabelEncoder<Label> *> &label_encoders,
                std::function<const BinaryMatrix::SetBitPositions(uint64_t)> get_next_row,
                uint64_t num_rows,
                const std::string &outfile) {
    LabelEncoder<Label> merged_label_enc;
    std::vector<std::vector<uint64_t> > label_mappings;

    for (const auto *label_encoder : label_encoders) {
        merged_label_enc.merge(*label_encoder);

        std::vector<uint64_t> v;
        for (size_t j = 0; j < label_encoder->size(); ++j) {
            v.push_back(merged_label_enc.encode(label_encoder->decode(j)));
        }
        label_mappings.push_back(v);
    }

    RowCompressed<Label>::serialize(
        outfile,
        merged_label_enc,
        [&](BinaryMatrix::RowCallback write_row) {
            std::set<uint64_t> indexes;
            for (uint64_t r = 0; r < num_rows; ++r) {
                for (uint64_t a = 0; a < label_encoders.size(); ++a) {
                    const auto &old_index_to_new = label_mappings.at(a);
                    for (auto old_i : get_next_row(a)) {
                        indexes.insert(old_index_to_new[old_i]);
                    }
                }
                BinaryMatrix::SetBitPositions merged_row(indexes.begin(), indexes.end());
                write_row(merged_row);

                indexes.clear();
            }
        }
    );
}

// TODO: move row iterators to BinaryMatrix
template <class T>
class Iterate {
  public:
    virtual ~Iterate() {}
    virtual T next() = 0;
};

template <class Annotator>
class IterateRows : public Iterate<BinaryMatrix::SetBitPositions> {
  public:
    IterateRows(const Annotator &annotator) : matrix_(annotator.get_matrix()) {};

    BinaryMatrix::SetBitPositions next() override {
        return matrix_.get_row(i_++);
    };

  private:
    typename BinaryMatrix::Row i_ = 0;
    const BinaryMatrix &matrix_;
};

// TODO: remove this if not crucial.
// Support only merge with streaming and only converted to row-major.
template <>
class IterateRows<ColumnCompressed<std::string>>
            : public Iterate<BinaryMatrix::SetBitPositions> {
  public:
    IterateRows(const ColumnCompressed<std::string> &annotator)
          : row_iterator_(std::make_unique<utils::RowsFromColumnsTransformer>(annotator.get_matrix().data())) {}

    BinaryMatrix::SetBitPositions next() override {
        return row_iterator_.next_row<BinaryMatrix::SetBitPositions>();
    }

  private:
    utils::RowsFromColumnsIterator row_iterator_;
};

template <class ToAnnotation, typename Label>
void merge(std::vector<std::unique_ptr<MultiLabelEncoded<Label>>>&& annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile) {
    static_assert(std::is_same_v<typename ToAnnotation::Label, Label>);

    assert((annotators.size() || filenames.size()) && "nothing to merge");

    uint64_t num_rows;
    uint64_t num_relations;

    if (annotators.size()) {
        num_rows = annotators.at(0)->num_objects();
    } else {
        RowCompressed<Label>::load_shape(filenames.at(0), &num_rows, &num_relations);
    }
    assert(num_rows);

    std::vector<const LEncoder*> label_encoders;

    std::vector<std::unique_ptr<Iterate<BinaryMatrix::SetBitPositions>>> annotator_row_iterators;
    for (const auto &annotator : annotators) {
        if (annotator->num_objects() != num_rows)
            throw std::runtime_error("Annotators have different number of rows");

        label_encoders.push_back(&annotator->get_label_encoder());
        annotator_row_iterators.push_back(
            std::make_unique<IterateRows<MultiLabelEncoded<Label>>>(*annotator)
        );
    }

    std::vector<std::unique_ptr<const LEncoder> > loaded_label_encoders;
    std::vector<std::unique_ptr<StreamRows<>>> streams;
    for (auto filename : filenames) {
        if (utils::ends_with(filename, RowCompressed<Label>::kExtension)) {

            auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);
            label_encoders.push_back(label_encoder.get());
            loaded_label_encoders.push_back(std::move(label_encoder));

            streams.emplace_back(new StreamRows<>(RowCompressed<Label>::get_row_streamer(filename)));

        } else {
            throw std::runtime_error("streaming only supported for rowcompressed annotator");
        }
    }

    merge_rows(
        label_encoders,
        [&](uint64_t annotator_idx) -> const BinaryMatrix::SetBitPositions {
            if (annotator_idx < annotators.size()) {
                return annotator_row_iterators.at(annotator_idx)->next();
            } else {
                return *streams[annotator_idx-annotators.size()]->next_row();
            }
        },
        num_rows,
        outfile
    );

    if constexpr(!std::is_same_v<RowCompressed<Label>, ToAnnotation>) {
        auto out_annotator = convert<ToAnnotation>(outfile);
        out_annotator->serialize(outfile);
    }
}

#define INSTANTIATE_MERGE(A, L) \
            template void \
            merge<A, L>(std::vector<std::unique_ptr<MultiLabelEncoded<L>>>&&, \
                        const std::vector<std::string>&, \
                        const std::string&);
INSTANTIATE_MERGE(RowFlatAnnotator, std::string);
INSTANTIATE_MERGE(RainbowfishAnnotator, std::string);
INSTANTIATE_MERGE(BinRelWTAnnotator, std::string);
INSTANTIATE_MERGE(BinRelWT_sdslAnnotator, std::string);
INSTANTIATE_MERGE(RowCompressed<>, std::string);


template<>
void merge<MultiBRWTAnnotator, std::string>(
        std::vector<std::unique_ptr<MultiLabelEncoded<std::string>>>&& annotators,
        const std::vector<std::string> &filenames,
        const std::string &outfile) {

    assert((annotators.size() || filenames.size()) && "nothing to merge");

    if (filenames.size()) {
        throw std::runtime_error("streaming only supported for rowcompressed annotator");
    }

    uint64_t num_rows = annotators.at(0)->num_objects();

    LEncoder label_encoder;

    std::vector<BRWT> brwts;

    for (auto&& annotator : annotators) {
        if (!dynamic_cast<MultiBRWTAnnotator*>(annotator.get()))
            throw std::runtime_error("merging of arbitrary annotations into BRWT is not implemented");

        if (annotator->num_objects() != num_rows)
            throw std::runtime_error("Annotators have different number of rows");

        for (const auto &label : annotator->get_label_encoder().get_labels()) {
            if (label_encoder.label_exists(label))
                throw std::runtime_error("merging of BRWT with same labels is not implemented");

            label_encoder.insert_and_encode(label);
        }

        brwts.push_back(std::move(const_cast<BRWT&>(
            dynamic_cast<MultiBRWTAnnotator&>(*annotator).get_matrix()
        )));

        annotator.reset();
    }

    MultiBRWTAnnotator annotation(
        std::make_unique<BRWT>(BRWTBottomUpBuilder::merge(
            std::move(brwts),
            BRWTBottomUpBuilder::get_basic_partitioner(-1),
            1,
            get_num_threads()
        )),
        label_encoder
    );

    annotation.serialize(outfile);
}

template <typename Label>
void convert_to_row_annotator(const ColumnCompressed<Label> &annotator,
                              const std::string &outfbase,
                              size_t num_threads) {
    uint64_t num_rows = annotator.num_objects();

    ProgressBar progress_bar(num_rows, "Serialize rows",
                             std::cerr, !common::get_verbose());

    RowCompressed<Label>::serialize(
        outfbase,
        annotator.get_label_encoder(),
        [&](BinaryMatrix::RowCallback write_row) {

            #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
            for (uint64_t i = 0; i < num_rows; i += kNumRowsInBlock) {

                uint64_t begin = i;
                uint64_t end = std::min(i + kNumRowsInBlock, num_rows);

                std::vector<BinaryMatrix::SetBitPositions> rows(end - begin);

                assert(begin <= end);
                assert(end <= num_rows);

                // TODO: use RowsFromColumnsTransformer
                for (const auto &label : annotator.get_all_labels()) {
                    size_t j = annotator.get_label_encoder().encode(label);
                    annotator.get_column(label).call_ones_in_range(begin, end,
                        [&](uint64_t idx) { rows[idx - begin].push_back(j); }
                    );
                }

                #pragma omp ordered
                {
                    for (const auto &row : rows) {
                        write_row(row);
                        ++progress_bar;
                    }
                }
            }
        }
    );
}

template
void convert_to_row_annotator(const ColumnCompressed<std::string> &annotator,
                              const std::string &outfbase,
                              size_t num_threads);

template <typename Label>
void convert_to_row_annotator(const ColumnCompressed<Label> &source,
                              RowCompressed<Label> *annotator,
                              size_t num_threads) {
    assert(annotator);
    assert(source.num_objects() == annotator->get_matrix().num_rows());

    uint64_t num_rows = source.num_objects();

    ProgressBar progress_bar(num_rows, "Rows processed",
                             std::cerr, !common::get_verbose());

    const_cast<LabelEncoder<Label>&>(annotator->get_label_encoder())
            = source.get_label_encoder();

    if (num_threads <= 1) {
        add_labels(source, 0, num_rows,
                   const_cast<BinaryMatrixRowDynamic*>(&annotator->get_matrix()),
                   &progress_bar);
        return;
    }

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t i = 0; i < num_rows; i += kNumRowsInBlock) {
        add_labels(source, i, std::min(i + kNumRowsInBlock, num_rows),
                   const_cast<BinaryMatrixRowDynamic*>(&annotator->get_matrix()),
                   &progress_bar);
    }
}

template <typename Label>
void add_labels(const ColumnCompressed<Label> &source,
                uint64_t begin, uint64_t end,
                BinaryMatrixRowDynamic *matrix,
                ProgressBar *progress_bar) {
    assert(matrix);
    assert(begin <= end);
    assert(end <= matrix->num_rows());

    // TODO: use RowsFromColumnsTransformer
    for (const auto &label : source.get_all_labels()) {
        size_t j = source.get_label_encoder().encode(label);
        source.get_column(label).call_ones_in_range(begin, end,
            [&](uint64_t idx) { matrix->set(idx, j); }
        );
    }
    if (progress_bar)
        *progress_bar += end - begin;
}

template
void convert_to_row_annotator(const ColumnCompressed<std::string> &source,
                              RowCompressed<std::string> *annotator,
                              size_t num_threads);

} // namespace annotate
