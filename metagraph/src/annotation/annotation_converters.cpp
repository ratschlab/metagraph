#include "annotation_converters.hpp"

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


template <>
std::unique_ptr<RowFlatAnnotator>
convert<RowFlatAnnotator, std::string>(RowCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    ProgressBar progress_bar(num_rows, "Processing rows");

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
std::unique_ptr<RowFlatAnnotator>
convert<RowCompressed<>, RowFlatAnnotator, std::string, false>(const std::string &filename) {
    uint64_t num_rows;
    uint64_t num_relations;
    RowCompressed<>::stream_counts(filename, &num_rows, &num_relations, false);

    auto label_encoder = RowCompressed<std::string>::load_label_encoder(filename);

    auto annotator = new RowCompressed<>::StreamRows(filename, false);

    ProgressBar progress_bar(num_rows, "Processing rows");

    auto matrix = std::make_unique<RowConcatenated<>>(
        [&](auto callback) {
            for (uint64_t r = 0; r < num_rows; ++r) {
                auto row = annotator->next_row();
                std::sort(row->begin(), row->end());
                callback(*row);
                ++progress_bar;
            }
        },
        label_encoder->size(),
        num_rows,
        num_relations
    );

    return std::make_unique<RowFlatAnnotator>(std::move(matrix), *label_encoder);
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

uint64_t merge_rows(
    const std::vector<const LEncoder*> &label_encoders,
    const std::function<const std::vector<uint64_t>(const uint64_t, const uint64_t)> &get_row,
    const uint64_t &num_rows,
    const std::function<uint64_t(LEncoder&, const std::function<void (BinaryMatrix::RowCallback&)>)> &write_rows
) {
    LEncoder merged_label_enc;
    merged_label_enc.merge(label_encoders);

    std::vector<std::vector<uint64_t> > label_mappings;
    for (size_t i = 0; i < label_encoders.size(); ++i) {
        std::vector<uint64_t> v;
        for (size_t j = 0; j < label_encoders.at(i)->size(); ++j) {
            v.push_back(merged_label_enc.encode(label_encoders.at(i)->decode(j)));
        }
        label_mappings.push_back(v);
    }

    return write_rows(merged_label_enc, [&](BinaryMatrix::RowCallback &write_row) {
        std::set<uint64_t> label_set;
        for (uint64_t r = 0; r < num_rows; ++r) {
            label_set.clear();
            for (uint64_t a = 0; a < label_encoders.size(); ++a) {
                for (auto label : get_row(r, a)) {
                    label_set.insert(label_mappings.at(a)[label]);
                }
            }
            std::vector<uint64_t> merged_row(label_set.begin(), label_set.end());
            write_row(merged_row);
        }
    });
}

template <>
uint64_t
merge<RowCompressed<>, RowCompressed<>, std::string, false>(const std::vector<std::string> &filenames,
                                                            const std::string &outfile) {
    assert(filenames.size()>0);
    uint64_t num_rows;
    uint64_t num_relations;
    RowCompressed<>::stream_counts(filenames.at(0), &num_rows, &num_relations, false);

    std::vector<std::unique_ptr<const LEncoder> > label_encoders_;
    std::vector<const LEncoder*> label_encoders;
    for (auto filename : filenames) {
        auto label_encoder = RowCompressed<std::string>::load_label_encoder(filename);
        label_encoders.push_back(label_encoder.get());
        label_encoders_.push_back(std::move(label_encoder));
    }

    std::vector<std::unique_ptr<RowCompressed<>::StreamRows> > annotators;
    for (size_t i = 0; i < filenames.size(); ++i) {
        auto annotator = std::make_unique<RowCompressed<>::StreamRows>(filenames.at(i), false);
        annotators.push_back(std::move(annotator));
    }

    return merge_rows(
        label_encoders,
        [&](const uint64_t row_idx, const uint64_t annotator_idx) -> const std::vector<uint64_t> {
            return *annotators[annotator_idx]->next_row();
        },
        num_rows,
        [&](LEncoder &merged_label_enc, const std::function<void (BinaryMatrix::RowCallback&)> &callback) {
            return RowCompressed<>::write_rows(outfile,
                                               merged_label_enc,
                                               callback,
                                               false);
        }
    );
}

template <>
uint64_t
merge<RowFlatAnnotator, RowCompressed<>, std::string, false>(const std::vector<const RowFlatAnnotator*> &annotators,
                                                             const std::string &outfile) {
    assert(annotators.size()>0);
    const uint64_t num_rows = annotators.at(0)->num_objects();

    std::vector<const LEncoder*> label_encoders;
    for (auto annotator : annotators) {
        label_encoders.push_back(&annotator->label_encoder_);
    }

    return merge_rows(
        label_encoders,
        [&](const uint64_t row_idx, const uint64_t annotator_idx) -> const std::vector<uint64_t> {
            return annotators.at(annotator_idx)->matrix_->get_row(row_idx);
        },
        num_rows,
        [&](LEncoder &merged_label_enc, const std::function<void (BinaryMatrix::RowCallback&)> &callback) {
            return RowCompressed<>::write_rows(outfile,
                                               merged_label_enc,
                                               callback,
                                               false);
        }
    );
}

} // namespace annotate
