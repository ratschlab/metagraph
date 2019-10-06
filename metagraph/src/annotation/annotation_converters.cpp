#include "annotation_converters.hpp"

#include <cassert>
#include <vector>
#include <functional>
#include <filesystem>
#include <progress_bar.hpp>

#include "annotate.hpp"
#include "static_annotators_def.hpp"
#include "annotate_column_compressed.hpp"
#include "BRWT_builders.hpp"
#include "partitionings.hpp"
#include "utils.hpp"
#include "binary_matrix.hpp"
#include "annotate_row_compressed.hpp"
#include "vector_row_binmat.hpp"


namespace annotate {

typedef LabelEncoder<std::string> LEncoder;


// RowCompressed -> other

template <>
std::unique_ptr<RowFlatAnnotator>
convert<RowFlatAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    ProgressBar progress_bar(num_rows, "Processing rows", std::cerr, !utils::get_verbose());

    if (dynamic_cast<VectorRowBinMat*>(annotator.matrix_.get()))
        dynamic_cast<VectorRowBinMat&>(*annotator.matrix_).standardize_rows();

    auto matrix = std::make_unique<RowConcatenated<>>(
        [&](auto callback) {
            utils::call_rows(
                [&](const auto &row) {
                    assert(std::is_sorted(row.begin(), row.end()));
                    callback(row);
                    ++progress_bar;
                },
                dynamic_cast<const BinaryMatrix &>(*annotator.matrix_)
            );
        },
        num_columns,
        num_rows,
        num_set_bits
    );

    return std::make_unique<RowFlatAnnotator>(std::move(matrix),
                                              annotator.label_encoder_);
}

template <>
std::unique_ptr<RainbowfishAnnotator>
convert<RainbowfishAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<Rainbowfish>([&](auto callback) {
        utils::call_rows(callback, dynamic_cast<const BinaryMatrix &>(*annotator.matrix_));
    }, num_columns);

    return std::make_unique<RainbowfishAnnotator>(std::move(matrix),
                                                  annotator.label_encoder_);
}

template <>
std::unique_ptr<BinRelWT_sdslAnnotator>
convert<BinRelWT_sdslAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT_sdsl>(
        [&](auto callback) {
            utils::call_rows(callback, dynamic_cast<const BinaryMatrix &>(*annotator.matrix_));
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWT_sdslAnnotator>(std::move(matrix),
                                                    annotator.label_encoder_);
}

template <>
std::unique_ptr<BinRelWTAnnotator>
convert<BinRelWTAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT>(
        [&](auto callback) {
            utils::call_rows(callback, dynamic_cast<const BinaryMatrix &>(*annotator.matrix_));
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWTAnnotator>(std::move(matrix),
                                               annotator.label_encoder_);
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
    RowCompressed<Label>::stream_counts(filename, &num_rows, &num_relations);

    constexpr size_t num_passes = std::is_same_v<MatrixType, Rainbowfish> ? 2u : 1u;
    ProgressBar progress_bar(num_rows * num_passes, "Processing rows", std::cerr, !utils::get_verbose());

    auto call_rows = [&](BinaryMatrix::RowCallback callback) {
        typename RowCompressed<Label>::template StreamRows<> row_streamer(filename);
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
    annotator.flush();

    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowConcatenated<>>([&](auto callback) {
        utils::call_rows(callback, annotator.bitmatrix_);
    }, num_columns, num_rows, num_set_bits);

    return std::make_unique<RowFlatAnnotator>(std::move(matrix),
                                              annotator.label_encoder_);
}

template <>
std::unique_ptr<RainbowfishAnnotator>
convert<RainbowfishAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    annotator.flush();

    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<Rainbowfish>([&](auto callback) {
        utils::call_rows(callback, annotator.bitmatrix_);
    }, num_columns);

    return std::make_unique<RainbowfishAnnotator>(std::move(matrix),
                                                  annotator.label_encoder_);
}

template <class StaticAnnotation, typename Label, class Partitioning>
typename std::unique_ptr<StaticAnnotation>
convert_to_BRWT(ColumnCompressed<Label>&& annotator,
                Partitioning partitioning,
                size_t num_threads) {
    annotator.flush();

    std::vector<std::unique_ptr<bit_vector>> columns;
    for (size_t j = 0; j < annotator.bitmatrix_.size(); ++j) {
        columns.emplace_back(std::move(annotator.bitmatrix_[j]));
    }

    auto matrix = std::make_unique<BRWT>(
        BRWTBottomUpBuilder::build(std::move(columns), partitioning, num_threads)
    );

    return std::make_unique<StaticAnnotation>(std::move(matrix),
                                              annotator.label_encoder_);
}

template <>
std::unique_ptr<BRWTCompressed<>>
convert_to_greedy_BRWT<BRWTCompressed<>, std::string>(ColumnCompressed<std::string>&& annotation,
                                                      size_t num_threads) {
    return convert_to_BRWT<BRWTCompressed<>>(
        std::move(annotation),
        get_parallel_binary_grouping_greedy(num_threads),
        num_threads
    );
}

template <>
std::unique_ptr<BRWTCompressed<>>
convert_to_simple_BRWT<BRWTCompressed<>, std::string>(ColumnCompressed<std::string>&& annotation,
                                                      size_t grouping_arity,
                                                      size_t num_threads) {
    return convert_to_BRWT<BRWTCompressed<>>(
        std::move(annotation),
        BRWTBottomUpBuilder::get_basic_partitioner(grouping_arity),
        num_threads
    );
}

template <>
void relax_BRWT<BRWTCompressed<>>(BRWTCompressed<> *annotation,
                                  size_t relax_max_arity,
                                  size_t num_threads) {
    if (relax_max_arity > 1)
        BRWTOptimizer::relax(const_cast<BRWT*>(&annotation->data()),
                             relax_max_arity,
                             num_threads);
}


template <>
std::unique_ptr<BinRelWT_sdslAnnotator>
convert<BinRelWT_sdslAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    annotator.flush();

    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT_sdsl>(
        [&](auto callback) {
            utils::call_rows(callback, annotator.bitmatrix_);
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWT_sdslAnnotator>(std::move(matrix),
                                                    annotator.label_encoder_);
}

template <>
std::unique_ptr<BinRelWTAnnotator>
convert<BinRelWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    annotator.flush();

    return std::make_unique<BinRelWTAnnotator>(
        std::make_unique<BinRelWT>(std::move(annotator.bitmatrix_)),
        annotator.label_encoder_
    );
}

template <typename Label>
void merge_rows(const std::vector<const LabelEncoder<Label> *> &label_encoders,
                std::function<const std::vector<uint64_t>(uint64_t)> get_next_row,
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

    RowCompressed<Label>::write_rows(
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
                std::vector<uint64_t> merged_row(indexes.begin(), indexes.end());
                write_row(merged_row);

                indexes.clear();
            }
        }
    );
}


template <class ToAnnotation, typename Label>
void merge(const std::vector<const MultiLabelEncoded<uint64_t, Label>*> &annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile) {
    static_assert(std::is_same_v<typename ToAnnotation::Label, Label>);

    assert((annotators.size() || filenames.size()) && "nothing to merge");

    uint64_t num_rows;
    uint64_t num_relations;

    if (annotators.size()) {
        num_rows = annotators.at(0)->num_objects();
    } else {
        RowCompressed<Label>::stream_counts(filenames.at(0), &num_rows, &num_relations);
    }
    assert(num_rows);

    std::vector<const LEncoder*> label_encoders;
    std::vector<std::unique_ptr<IterateRows>> annotator_row_iterators;
    for (const auto *annotator : annotators) {
        if (annotator->num_objects() != num_rows)
            throw std::runtime_error("Annotators have different number of rows");

        label_encoders.push_back(&annotator->get_label_encoder());
        annotator_row_iterators.push_back(annotator->iterator());
    }

    std::vector<std::unique_ptr<const LEncoder> > loaded_label_encoders;
    std::vector<std::unique_ptr<typename RowCompressed<Label>::template StreamRows<>>> streams;
    for (auto filename : filenames) {
        if (utils::ends_with(filename, RowCompressed<Label>::kExtension)) {

            auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);
            label_encoders.push_back(label_encoder.get());
            loaded_label_encoders.push_back(std::move(label_encoder));

            auto annotator = std::make_unique<typename RowCompressed<Label>::template StreamRows<>>(filename);
            streams.push_back(std::move(annotator));

        } else {
            throw std::runtime_error("streaming only supported for rowcompressed annotator");
        }
    }

    merge_rows(
        label_encoders,
        [&](uint64_t annotator_idx) -> const std::vector<uint64_t> {
            if (annotator_idx < annotators.size()) {
                return annotator_row_iterators.at(annotator_idx)->next_row();
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
            merge<A, L>(const std::vector<const MultiLabelEncoded<uint64_t, L>*>&, \
                        const std::vector<std::string>&, \
                        const std::string&);
INSTANTIATE_MERGE(RowFlatAnnotator, std::string);
INSTANTIATE_MERGE(RainbowfishAnnotator, std::string);
INSTANTIATE_MERGE(BinRelWTAnnotator, std::string);
INSTANTIATE_MERGE(BinRelWT_sdslAnnotator, std::string);
INSTANTIATE_MERGE(RowCompressed<>, std::string);


} // namespace annotate
