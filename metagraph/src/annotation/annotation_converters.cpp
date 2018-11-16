#include "annotation_converters.hpp"

#include "static_annotators_def.hpp"
#include "annotate_column_compressed.hpp"
#include "BRWT_builders.hpp"
#include "partitionings.hpp"
#include "utils.hpp"


namespace annotate {


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

} // namespace annotate
