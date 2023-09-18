#include "annotation_converters.hpp"

#include <cassert>
#include <vector>
#include <functional>
#include <filesystem>

#include <ips4o.hpp>
#include <progress_bar.hpp>
#include <tsl/hopscotch_map.h>

#include "row_diff_builder.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/hashers/hash.hpp"
#include "common/sorted_sets/sorted_multiset.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/transpose.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "binary_matrix/row_vector/vector_row_binmat.hpp"
#include "binary_matrix/multi_brwt/brwt_builders.hpp"
#include "binary_matrix/multi_brwt/clustering.hpp"
#include "int_matrix/rank_extended/tuple_csc_matrix.hpp"
#include "representation/annotation_matrix/static_annotators_def.hpp"
#include "representation/column_compressed/annotate_column_compressed.hpp"
#include "representation/row_compressed/annotate_row_compressed.hpp"


namespace mtg {
namespace annot {

using namespace mtg::annot::matrix;

using mtg::common::logger;

namespace fs = std::filesystem;

typedef LabelEncoder<std::string> LEncoder;
typedef TupleCSCMatrix<ColumnMajor> TupleCSC;
typedef TupleCSCMatrix<BRWT> TupleBRWT;

const size_t kNumRowsInBlock = 50'000;
const uint64_t ROW_DIFF_BUFFER_BYTES = 8'000'000;


// RowCompressed -> other

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation> convert(const std::string &filename) {
    using MatrixType = typename StaticAnnotation::binary_matrix_type;
    using Label = typename StaticAnnotation::Label;

    auto label_encoder = RowCompressed<Label>::read_label_encoder(filename);

    uint64_t num_rows;
    uint64_t num_relations;
    RowCompressed<Label>::read_shape(filename, &num_rows, &num_relations);

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

    if constexpr(std::is_same_v<MatrixType, RowFlat<>>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder.size(), num_rows, num_relations);

    } else if constexpr(std::is_same_v<MatrixType, Rainbowfish>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder.size());

    } else if constexpr(std::is_same_v<MatrixType, BinRelWT>) {
        matrix = std::make_unique<MatrixType>(call_rows, num_relations, label_encoder.size());

    } else if constexpr(std::is_same_v<MatrixType, RowSparse>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder.size(), num_rows, num_relations);

    } else {
        static_assert(utils::dependent_false<StaticAnnotation>::value);
    }

    return std::make_unique<StaticAnnotation>(std::move(matrix), label_encoder);
}

template std::unique_ptr<RowFlatAnnotator> convert(const std::string &filename);
template std::unique_ptr<RowSparseAnnotator> convert(const std::string &filename);
template std::unique_ptr<RainbowfishAnnotator> convert(const std::string &filename);
template std::unique_ptr<BinRelWTAnnotator> convert(const std::string &filename);


// ColumnCompressed -> other

template <>
std::unique_ptr<RowFlatAnnotator>
convert<RowFlatAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowFlat<>>([&](auto callback) {
        utils::call_rows(annotator.get_matrix().data(), callback, annotator.num_objects());
    }, num_columns, num_rows, num_set_bits);

    return std::make_unique<RowFlatAnnotator>(std::move(matrix),
                                              annotator.get_label_encoder());
}

template <>
void convert<RowFlatAnnotator, std::string>(ColumnCompressed<std::string>&& annotator,
                                            const std::string &outfbase) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto call_rows = [&](RowFlat<>::RowCallback callback) {
        utils::call_rows(annotator.get_matrix().data(), callback, annotator.num_objects());
    };

    const auto &fname = utils::make_suffix(outfbase, RowFlatAnnotator::kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + fname);

    annotator.get_label_encoder().serialize(out);
    out.close();

    RowFlat<>::serialize(call_rows, num_columns, num_rows, num_set_bits, fname, true);
}

template <>
std::unique_ptr<RowSparseAnnotator>
convert<RowSparseAnnotator, std::string>(ColumnCompressed<std::string> &&annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowSparse>(
        [&](auto callback) {
            utils::call_rows(annotator.get_matrix().data(), callback, annotator.num_objects());
        },
        num_columns, num_rows, num_set_bits
    );
    return std::make_unique<RowSparseAnnotator>(std::move(matrix),
                                                annotator.get_label_encoder());
}

template <>
std::unique_ptr<RainbowfishAnnotator>
convert<RainbowfishAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<Rainbowfish>([&](auto callback) {
        utils::call_rows(annotator.get_matrix().data(), callback);
    }, num_columns);

    return std::make_unique<RainbowfishAnnotator>(std::move(matrix),
                                                  annotator.get_label_encoder());
}

template <>
std::unique_ptr<UniqueRowAnnotator>
convert<UniqueRowAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<UniqueRowBinmat>([&](auto callback) {
        utils::call_rows(annotator.get_matrix().data(), callback);
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

std::unique_ptr<RowDiffBRWTAnnotator>
convert_row_diff_to_BRWT(RowDiffColumnAnnotator &&annotator,
                         BRWTBottomUpBuilder::Partitioner partitioning,
                         size_t num_parallel_nodes,
                         size_t num_threads) {
    const graph::DBGSuccinct* graph = annotator.get_matrix().graph();

    auto matrix = std::make_unique<BRWT>(
            BRWTBottomUpBuilder::build(std::move(annotator.release_matrix()->diffs().data()),
                                       partitioning,
                                       num_parallel_nodes,
                                       num_threads)
    );

    return std::make_unique<RowDiffBRWTAnnotator>(
            std::make_unique<RowDiff<BRWT>>(graph, std::move(*matrix)),
            annotator.get_label_encoder());
}

std::unique_ptr<MultiBRWTAnnotator>
convert_to_greedy_BRWT(ColumnCompressed<std::string> &&annotation,
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

std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_greedy_BRWT(RowDiffColumnAnnotator &&annotation,
                       size_t num_parallel_nodes,
                       size_t num_threads,
                       uint64_t num_rows_subsampled) {
    return convert_row_diff_to_BRWT(
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

std::unique_ptr<MultiBRWTAnnotator>
convert_to_simple_BRWT(ColumnCompressed<std::string>&& annotation,
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

std::unique_ptr<RowDiffBRWTAnnotator> convert_to_simple_BRWT(RowDiffColumnAnnotator&& annotation,
                                                             size_t grouping_arity,
                                                             size_t num_parallel_nodes,
                                                             size_t num_threads) {
    return convert_row_diff_to_BRWT(
            std::move(annotation),
            BRWTBottomUpBuilder::get_basic_partitioner(grouping_arity),
            num_parallel_nodes, num_threads);
}

std::pair<std::string, std::string> get_anchors_and_fork_fnames(const std::string &fbase) {
    std::string anchors_file = fbase + kRowDiffAnchorExt;
    if (!std::filesystem::exists(anchors_file)) {
        logger->error("Anchor bitmap {} does not exist. Run the row_diff"
                      " transform followed by anchor optimization.", anchors_file);
        std::exit(1);
    }
    std::string fork_succ_file = fbase + kRowDiffForkSuccExt;
    if (!std::filesystem::exists(fork_succ_file)) {
        logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
        std::exit(1);
    }
    return std::make_pair(anchors_file, fork_succ_file);
}

std::tuple<std::vector<std::unique_ptr<bit_vector>>, LabelEncoder<std::string>, uint64_t>
load_row_diff_columns(const std::vector<std::string> &files) {
    std::vector<std::unique_ptr<bit_vector>> columns;
    LabelEncoder<std::string> label_encoder;

    std::mutex mu;
    uint64_t total_num_set_bits = 0;
    bool load_successful = merge_load_row_diff(
        files,
        [&](uint64_t, const std::string &label, std::unique_ptr<bit_vector> &&column) {
            uint64_t num_set_bits = column->num_set_bits();
            logger->trace("RowDiff column: {}, Density: {}, Set bits: {}", label,
                          static_cast<double>(num_set_bits) / column->size(),
                          num_set_bits);

            std::lock_guard<std::mutex> lock(mu);

            if (columns.size() && column->size() != columns.back()->size()) {
                logger->error("Column {} has {} rows, previous column has {} rows",
                              label, column->size(), columns.back()->size());
                exit(1);
            }
            size_t col = label_encoder.insert_and_encode(label);
            if (col != columns.size()) {
                logger->error("Duplicate columns {}", label);
                exit(1);
            }
            columns.push_back(std::move(column));
            total_num_set_bits += num_set_bits;
        },
        get_num_threads()
    );

    if (!load_successful) {
        logger->error("Error while loading row-diff columns");
        exit(1);
    }

    return std::make_tuple(std::move(columns), std::move(label_encoder), total_num_set_bits);
}

template <>
void convert_to_row_diff<RowDiffRowFlatAnnotator>(
            const std::vector<std::string> &files,
            const std::string &anchors_file_fbase,
            const std::string &outfbase,
            size_t /*num_threads*/,
            size_t /*mem_bytes*/) {

    std::string anchors_file;
    std::string fork_succ_file;
    std::tie(anchors_file, fork_succ_file) = get_anchors_and_fork_fnames(anchors_file_fbase);

    logger->trace("Loading annotation from disk...");

    std::vector<std::unique_ptr<bit_vector>> columns;
    LabelEncoder<std::string> label_encoder;
    uint64_t num_set_bits;
    std::tie(columns, label_encoder, num_set_bits) = load_row_diff_columns(files);

    const auto &fname = utils::make_suffix(outfbase, RowDiffRowFlatAnnotator::kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + fname);

    label_encoder.serialize(out);

    // serialize RowDiff<RowFlat<>>
    out.write("v2.0", 4);
    out.close();
    utils::append_file(anchors_file, fname);
    utils::append_file(fork_succ_file, fname);

    uint64_t num_rows = columns.at(0)->size();
    uint64_t num_columns = columns.size();

    RowFlat<>::serialize([&](auto callback) { utils::call_rows(columns, callback); },
                         num_columns, num_rows, num_set_bits, fname, true);
    logger->trace("Annotation converted");
}

template <>
void convert_to_row_diff<RowDiffRowSparseAnnotator>(
            const std::vector<std::string> &files,
            const std::string &anchors_file_fbase,
            const std::string &outfbase,
            size_t /*num_threads*/,
            size_t /*mem_bytes*/) {

    std::string anchors_file;
    std::string fork_succ_file;
    std::tie(anchors_file, fork_succ_file) = get_anchors_and_fork_fnames(anchors_file_fbase);

    logger->trace("Loading annotation from disk...");

    std::vector<std::unique_ptr<bit_vector>> columns;
    LabelEncoder<std::string> label_encoder;
    uint64_t num_set_bits;
    std::tie(columns, label_encoder, num_set_bits) = load_row_diff_columns(files);

    uint64_t num_rows = columns.at(0)->size();
    uint64_t num_columns = columns.size();

    RowDiffRowSparseAnnotator row_sparse(label_encoder, nullptr,
        [&](auto callback) { utils::call_rows(columns, callback); },
        num_columns, num_rows, num_set_bits
    );

    const_cast<RowDiff<RowSparse> &>(row_sparse.get_matrix())
            .load_anchor(anchors_file);
    const_cast<RowDiff<RowSparse> &>(row_sparse.get_matrix())
            .load_fork_succ(fork_succ_file);

    logger->trace("Annotation converted");
    row_sparse.serialize(outfbase);
}


std::unique_ptr<MultiBRWTAnnotator>
convert_to_BRWT(
        const std::vector<std::vector<uint64_t>> &linkage,
        size_t num_parallel_nodes,
        size_t num_threads,
        const fs::path &tmp_path,
        const std::function<void(const BRWTBottomUpBuilder::CallColumn &)> &get_columns,
        std::vector<std::pair<uint64_t, std::string>> &&column_names) {
    if (!linkage.size())
        logger->warn("Parsed empty linkage tree");

    auto matrix = std::make_unique<BRWT>(
            BRWTBottomUpBuilder::build(get_columns, linkage, tmp_path,
                                       num_parallel_nodes, num_threads));

    std::sort(column_names.begin(), column_names.end(), utils::LessFirst());
    column_names.erase(std::unique(column_names.begin(), column_names.end()),
                       column_names.end());

    assert(matrix->num_columns() == column_names.size());

    LEncoder label_encoder;
    for (const auto &[i, label] : column_names) {
        if (label_encoder.label_exists(label)) {
            logger->error("Label '{}' is found in multiple annotations\nMerging labels in not supported"
                          " in the Multi-BRWT convertor. Consider merging columns before the conversion.",
                          label);
            exit(1);
        }
        label_encoder.insert_and_encode(label);
    }

    return std::make_unique<MultiBRWTAnnotator>(std::move(matrix), label_encoder);
}

template <>
std::unique_ptr<MultiBRWTAnnotator> convert_to_BRWT<MultiBRWTAnnotator>(
        const std::vector<std::string> &annotation_files,
        const std::vector<std::vector<uint64_t>> &linkage,
        size_t num_parallel_nodes,
        size_t num_threads,
        const fs::path &tmp_path) {
    std::vector<std::pair<uint64_t, std::string>> column_names;
    std::mutex mu;

    auto get_columns = [&](const BRWTBottomUpBuilder::CallColumn &call_column) {
        bool success = ColumnCompressed<>::merge_load(
            annotation_files,
            [&](uint64_t j,
                    const std::string &label,
                    std::unique_ptr<bit_vector>&& column) {
                call_column(j, std::move(column));
                std::lock_guard<std::mutex> lock(mu);
                column_names.emplace_back(j, label);
            },
            num_threads
        );
        if (!success) {
            logger->error("Can't load annotation columns");
            exit(1);
        }
    };
    return convert_to_BRWT(linkage, num_parallel_nodes, num_threads,
                           tmp_path, get_columns, std::move(column_names));
}

template<>
std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_BRWT<RowDiffBRWTAnnotator>(const std::vector<std::string> &annotation_files,
                         const std::vector<std::vector<uint64_t>> &linkage,
                         size_t num_parallel_nodes,
                         size_t num_threads,
                         const fs::path &tmp_path) {
    std::vector<std::pair<uint64_t, std::string>> column_names;
    std::mutex mu;

    auto get_columns = [&](const BRWTBottomUpBuilder::CallColumn &call_column) {
        bool success = merge_load_row_diff(
            annotation_files,
            [&](uint64_t j,
                    const std::string &label,
                    std::unique_ptr<bit_vector>&& column) {
                call_column(j, std::move(column));
                std::lock_guard<std::mutex> lock(mu);
                column_names.emplace_back(j, label);
            },
            num_threads
        );
        if (!success) {
            logger->error("Can't load annotation columns");
            exit(1);
        }
    };

    std::unique_ptr<MultiBRWTAnnotator> annotator
            = convert_to_BRWT(linkage, num_parallel_nodes, num_threads,
                              tmp_path, get_columns, std::move(column_names));

    return std::make_unique<RowDiffBRWTAnnotator>(
            std::make_unique<RowDiff<BRWT>>(nullptr, std::move(*annotator->release_matrix())),
            annotator->get_label_encoder());
}

void relax_BRWT(BRWT *annotation, size_t relax_max_arity, size_t num_threads) {
    if (relax_max_arity > 1)
        BRWTOptimizer::relax(annotation, relax_max_arity, num_threads);
}

using CallColumn = std::function<void(std::unique_ptr<bit_vector>&&)>;

std::vector<uint64_t>
get_row_classes(const std::function<void(const CallColumn &)> &call_columns,
                size_t num_columns) {
    std::vector<uint64_t> row_classes;
    uint64_t max_class = 0;

    ProgressBar progress_bar(num_columns, "Iterate columns",
                             std::cerr, !common::get_verbose());

    tsl::hopscotch_map<uint64_t, uint64_t> new_class;

    call_columns([&](const auto &col_ptr) {
        new_class.clear();

        if (row_classes.empty())
            row_classes.assign(col_ptr->size(), 0);

        const size_t batch_size = 5'000'000;
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (uint64_t begin = 0; begin < col_ptr->size(); begin += batch_size) {
            uint64_t end = std::min(begin + batch_size, col_ptr->size());

            std::vector<uint64_t> pos;
            pos.reserve(batch_size);

            col_ptr->call_ones_in_range(begin, end, [&](uint64_t i) {
                pos.push_back(i);
            });

            #pragma omp critical
            {
                for (uint64_t i : pos) {
                    uint64_t &c = row_classes[i];

                    auto [it, inserted] = new_class.try_emplace(c, max_class + 1);
                    if (inserted) {
                        ++max_class;

                        if (max_class == std::numeric_limits<uint64_t>::max()) {
                            logger->error("Too many distinct rows");
                            exit(1);
                        }
                    }

                    c = it->second;
                }
            }
        }

        ++progress_bar;
    });

    return row_classes;
}

std::unique_ptr<Rainbow<BRWT>>
convert_to_RainbowBRWT(const std::function<void(const CallColumn &)> &call_columns,
                       size_t max_brwt_arity = 1) {
    uint64_t num_columns = 0;
    uint64_t num_ones = 0;
    call_columns([&](const auto &col_ptr) {
        num_columns++;
        num_ones += col_ptr->num_set_bits();
    });

    if (!num_columns)
        return std::make_unique<Rainbow<BRWT>>();

    logger->trace("Identifying unique rows");
    std::vector<uint64_t> row_classes = get_row_classes(call_columns, num_columns);

    size_t batch_size = 100'000'000;
    common::SortedMultiset<uint64_t, uint64_t> row_counter(get_num_threads(), batch_size);
    {
        ProgressBar progress_bar(row_classes.size(),
                                 "Counting row classes",
                                 std::cerr, !common::get_verbose());
        auto it = row_classes.begin();
        while (it + batch_size < row_classes.end()) {
            row_counter.insert(it, it + batch_size);
            it += batch_size;
            progress_bar += batch_size;
        }
        row_counter.insert(it, row_classes.end());
        progress_bar += row_classes.end() - it;
    }
    auto &pairs = row_counter.data();
    logger->trace("Number of unique rows: {}", pairs.size());

    ips4o::parallel::sort(pairs.begin(), pairs.end(), utils::GreaterSecond(),
                          get_num_threads());

    if (common::get_verbose()) {
        // print historgram of row multiplicities
        std::stringstream mult_log_message;
        uint64_t multiplicity = pairs[0].second;
        uint64_t num_unique_rows = 0;
        for (const auto &[row_class, count] : pairs) {
            assert(count <= multiplicity);
            if (count == multiplicity) {
                num_unique_rows++;
            } else {
                mult_log_message << multiplicity << ": " << num_unique_rows << ", ";
                multiplicity = count;
                num_unique_rows = 1;
            }
        }
        mult_log_message << multiplicity << ": " << num_unique_rows;
        logger->trace("<Row multiplicity: num unique rows>: {}", mult_log_message.str());
    }

    tsl::hopscotch_map<uint64_t, uint64_t> class_to_code;
    for (size_t i = 0; i < pairs.size(); ++i) {
        class_to_code.emplace(pairs[i].first, i);
    }

    // maps unique row ids to the respective rows of the original matrix
    std::vector<uint64_t> row_pointers(class_to_code.size());

    uint64_t total_code_length = 0;
    for (size_t i = 0; i < row_classes.size(); ++i) {
        uint64_t &c = row_classes[i];
        c = class_to_code[c];
        row_pointers[c] = i;
        total_code_length += sdsl::bits::hi(c) + 1;
    }
    class_to_code.clear();

    logger->trace("Started matrix reduction");

    ProgressBar progress_bar(num_columns, "Reduce columns",
                             std::cerr, !common::get_verbose());

    std::vector<std::unique_ptr<bit_vector>> columns(num_columns);
    size_t j = 0;
    ThreadPool thread_pool(get_num_threads());
    call_columns([&](auto &&col_ptr) {
        thread_pool.enqueue([&](size_t j, const auto &col_ptr) {
            sdsl::bit_vector reduced_column(row_pointers.size(), false);
            for (size_t r = 0; r < row_pointers.size(); ++r) {
                reduced_column[r] = (*col_ptr)[row_pointers[r]];
            }
            columns[j] = std::make_unique<bit_vector_smart>(std::move(reduced_column));
            ++progress_bar;
        }, j++, std::move(col_ptr));
    });
    thread_pool.join();

    logger->trace("Start compressing the assignment vector");

    sdsl::bit_vector code_bv(total_code_length);
    sdsl::bit_vector boundary_bv(total_code_length, false);

    uint64_t pos = 0;
    for (uint64_t c : row_classes) {
        uint8_t code_len = sdsl::bits::hi(c) + 1;
        code_bv.set_int(pos, c, code_len);
        boundary_bv[pos + code_len - 1] = true;
        pos += code_len;
    }
    assert(pos == code_bv.size());

    bit_vector_rrr<> row_code_delimiters(std::move(boundary_bv));
    logger->trace("Assignment vector constructed"
                  "\nRow codes: {} bits\nBoundary bitmap: {} bits, {} set",
                  total_code_length, total_code_length, row_code_delimiters.num_set_bits());

    logger->trace("Building Multi-BRWT");

    size_t num_threads = get_num_threads();
    size_t num_rows_subsampled = 1'000'000;
    size_t num_parallel_nodes = num_threads;

    auto partitioning = [num_threads,num_rows_subsampled](const auto &columns) {
        std::vector<sdsl::bit_vector> subvectors
            = random_submatrix(columns, num_rows_subsampled, num_threads);
        return greedy_matching(subvectors, num_threads);
    };

    BRWT reduced_matrix = BRWTBottomUpBuilder::build(std::move(columns),
                                                     partitioning,
                                                     num_parallel_nodes,
                                                     num_threads);

    if (max_brwt_arity > 1)
        BRWTOptimizer::relax(&reduced_matrix, max_brwt_arity, num_threads);

    return std::make_unique<Rainbow<BRWT>>(std::move(reduced_matrix),
                                           std::move(code_bv),
                                           std::move(row_code_delimiters),
                                           num_ones);
}

template <>
std::unique_ptr<RbBRWTAnnotator>
convert_to_RbBRWT<RbBRWTAnnotator>(const std::vector<std::string> &annotation_files,
                                   size_t max_brwt_arity) {
    LEncoder label_encoder;

    auto call_columns = [&](const CallColumn &call_column) {
        label_encoder.clear();

        bool success = ColumnCompressed<>::merge_load(
            annotation_files,
            [&](uint64_t /*j*/,
                    const std::string &label,
                    std::unique_ptr<bit_vector>&& column) {
                call_column(std::move(column));
                label_encoder.insert_and_encode(label);
            },
            0
        );
        if (!success) {
            logger->error("Can't load annotation columns");
            exit(1);
        }
    };

    auto matrix = convert_to_RainbowBRWT(call_columns, max_brwt_arity);

    assert(matrix->num_columns() == label_encoder.size());

    return std::make_unique<RbBRWTAnnotator>(std::move(matrix), label_encoder);
}

void convert_batch_to_row_disk(
        const std::function<bool(const std::vector<std::string> &, const ColumnCallback &, size_t)>
                &get_cols,
        const std::vector<std::string> &source_files,
        const std::string &outfname,
        const std::function<void(std::ofstream &)> &decorate,
        uint64_t num_rows,
        size_t num_threads) {
    std::vector<std::string> col_names;
    std::vector<std::unique_ptr<bit_vector>> columns;
    uint64_t num_set_bits = 0;
    std::mutex mu;
    bool success = get_cols(
            source_files,
            [&](uint64_t j,
                    const std::string &label,
                    std::unique_ptr<bit_vector>&& column) {
                if (num_rows != column->size()) {
                    logger->error("Annotations have different number of rows and can't be merged");
                    std::exit(1);
                }

                std::lock_guard<std::mutex> lock(mu);

                num_set_bits += column->num_set_bits();
                while (columns.size() <= j) {
                    columns.emplace_back();
                    col_names.emplace_back();
                }
                columns[j] = std::move(column);
                col_names[j] = label;
            },
            num_threads
    );

    if (!success) {
        logger->error("Can't load annotation columns");
        exit(1);
    }

    // this must be done after loading all columns
    // to keep their order correct
    // (label_encoder.insert_and_encode cannot be inside get_cols callback)
    LEncoder label_encoder;
    for (const auto &label : col_names) {
        size_t col = label_encoder.insert_and_encode(label);
        if (col + 1 != label_encoder.size()) {
            logger->error("Duplicate columns {}", label);
            exit(1);
        }
    }

    std::ofstream out = utils::open_new_ofstream(outfname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + outfname);

    label_encoder.serialize(out);

    decorate(out);

    out.close();

    RowDisk::serialize(outfname,
        [&](BinaryMatrix::RowCallback write_row) {
            utils::call_rows(columns, write_row);
        },
        label_encoder.size(), num_rows, num_set_bits
    );
}

uint64_t get_num_rows_from_column_anno(const std::string &fname) {
    std::ifstream in(fname, std::ios::binary);
    if (!in.good())
        throw std::ifstream::failure("can't open file");

    return load_number(in);
}


void merge_row_disk_annotations(const std::vector<std::string> &files,
                                const std::string &outfname,
                                const std::function<void(std::ofstream &)> &decorate,
                                uint64_t num_rows,
                                size_t num_threads,
                                size_t mem_bytes) {
    if (files.size() == 1) {
        fs::rename(files[0], outfname);
        return;
    }

    logger->trace("Merge {} row_disk annotations", files.size());

    LEncoder label_encoder;
    std::vector<std::unique_ptr<RowDisk>> matrices(files.size());
    std::vector<uint64_t> offsets(files.size() + 1, 0);

    uint64_t num_set_bits = 0;
    size_t tot_boundary_bytes = 0;
    for (size_t j = 0; j < files.size(); ++j) {
        LEncoder le;
        matrices[j] = std::make_unique<RowDisk>();

        sdsl::mmap_ifstream in(files[j], std::ios::binary);
        if (!le.load(in)) {
            logger->error("Can't load label encoder from {}", files[j]);
            exit(1);
        }

        for (const auto &label : le.get_labels()) {
            size_t col = label_encoder.insert_and_encode(label);
            if (col + 1 != label_encoder.size()) {
                logger->error("Duplicate columns {}", label);
                exit(1);
            }
        }

        offsets[j + 1] = le.size();
        if (!matrices[j]->load(in)) {
            logger->error("Can't load matrix from {}", files[j]);
            exit(1);
        }

        num_set_bits += matrices[j]->num_relations();
        tot_boundary_bytes += matrices[j]->get_boundary().size() / 8;
    }

    size_t predicted_serialize_memory
            = bit_vector_small::predict_size(num_set_bits + num_rows, num_rows) / 8 // boundary
                + num_threads
                    * (2 * kNumRowsInBlock * sizeof(uint64_t) * num_set_bits / num_rows
                       + kNumRowsInBlock * sizeof(BinaryMatrix::Row));

    size_t tot_predicted_mem = predicted_serialize_memory + tot_boundary_bytes;
    if (tot_predicted_mem > mem_bytes) {
        logger->error(
                "{} MB is not enough to complete the conversion with {} threads."
                " > {} MB is requred", mem_bytes / 1e6, num_threads, tot_predicted_mem / 1e6);
        exit(1);
    }

    std::ofstream out = utils::open_new_ofstream(outfname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + outfname);

    label_encoder.serialize(out);
    decorate(out);
    out.close();

    // compute global offsets (partial sums)
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());

    ProgressBar progress_bar(num_rows, "Merge rows", std::cerr, !common::get_verbose());

    RowDisk::serialize(
            outfname,
            [&](BinaryMatrix::RowCallback write_row) {
                // maybe better if single thread reads single annotation?
                #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
                for (uint64_t begin = 0; begin < num_rows; begin += kNumRowsInBlock) {
                    uint64_t end = std::min(begin + kNumRowsInBlock, num_rows);

                    assert(begin <= end);
                    assert(end <= num_rows);

                    std::vector<BinaryMatrix::SetBitPositions> rows(end - begin);

                    std::vector<BinaryMatrix::Row> query_row_id(end - begin);
                    std::iota(query_row_id.begin(), query_row_id.end(), begin);

                    for (size_t b = 0; b < files.size(); ++b) {
                        auto input_rows = matrices[b]->get_rows(query_row_id);
                        for (size_t r = 0; r < input_rows.size(); ++r) {
                            for (auto j : input_rows[r]) {
                                rows[r].push_back(j + offsets[b]);
                            }
                        }
                    }

                    #pragma omp ordered
                    {
                        for (const auto &row : rows) {
                            write_row(row);
                            ++progress_bar;
                        }
                    }
                }
            },
            label_encoder.size(), num_rows, num_set_bits);

    for (const auto &fname : files) {
        fs::remove(fname);
    }
}

// works for cols -> row_disk and row_diff -> row_diff_disk
void convert_to_row_disk(
        const std::vector<std::string> &files,
        const std::string &outfbase,
        size_t num_threads,
        size_t mem_bytes,
        const std::string &result_fname,
        const std::string &file_extension,
        uint64_t num_rows,
        const std::function<void(std::ofstream &)> &decorate,
        const std::function<bool(const std::vector<std::string> &, const ColumnCallback &, size_t)>
                &get_cols) {
    if (!files.size())
        return;

    assert(utils::ends_with(result_fname, file_extension));

    double density_prediction = 0.01; // TODO: predict it better?

    size_t conversion_overhead_per_column = num_threads * kNumRowsInBlock * 8 * 2
            * density_prediction; // 8 bytes per stored position, times 2 because of capacity

    std::vector<std::string> parts;

    // load as many columns as we can fit in memory, and convert them
    for (uint32_t i = 0; i < files.size(); ) {
        logger->trace("Loading columns for batch-conversion...");
        size_t columns_size = 0;
        std::vector<std::string> file_batch;

        for (; i < files.size(); ++i) {
            uint64_t file_size = utils::with_mmap() ? 0 : fs::file_size(files[i]);

            if (file_size > mem_bytes) {
                logger->error("Not enough memory to load {}, requires at least {} MB",
                              files[i], file_size / 1e6);
                exit(1);
            }

            // TODO: construct with disk swap (sd_vector_swap), so this won't need to reserve anything
            size_t boundary_size
                    = bit_vector_small::predict_size(num_rows * (file_batch.size() + 1)
                                                                     * density_prediction
                                                             + num_rows,
                                                     num_rows)
                    / 8;

            size_t conversion_overhead
                    = conversion_overhead_per_column * (file_batch.size() + 1);

            if (columns_size + file_size + boundary_size + conversion_overhead > mem_bytes) {
                if (file_batch.empty()) {
                    logger->error("{} MB is not enough to perform the conversion. > {} MB requred",
                                  mem_bytes / 1e6,
                                  (columns_size + file_size + boundary_size + conversion_overhead) / 1e6);
                    exit(1);
                }
                break;
            }

            columns_size += file_size;

            file_batch.push_back(files[i]);
        }

        parts.push_back(utils::make_suffix(fmt::format("{}-{}", outfbase, parts.size()),
                                           file_extension));
        Timer timer;
        logger->trace("Annotations in batch: {}", file_batch.size());

        convert_batch_to_row_disk(get_cols, file_batch, parts.back(),
                                  file_batch.size() == files.size() ? decorate
                                                                    : [](std::ofstream &) {},
                                  num_rows, num_threads);

        logger->trace("Batch processed in {} sec", timer.elapsed());
    }

    // merge the temp row sparse disk annotations
    merge_row_disk_annotations(parts, result_fname, decorate, num_rows,
                               num_threads, mem_bytes);
}


uint64_t get_num_rows_from_row_diff_anno(const std::string &fname) {
    std::unique_ptr<std::ifstream> in = utils::open_ifstream(fname);
    if (!in->good())
        throw std::ifstream::failure("can't open file");

    LabelEncoder<std::string> label_encoder;
    RowDiff<ColumnMajor> matrix;

    if (!label_encoder.load(*in) || !matrix.load(*in)) {
        logger->error("Can't load {}", fname);
        exit(1);
    }

    return matrix.num_rows();
}

template <>
void convert_to_row_diff<RowDiffDiskAnnotator>(
            const std::vector<std::string> &files,
            const std::string &anchors_file_fbase,
            const std::string &outfbase,
            size_t num_threads,
            size_t mem_bytes) {

    std::string anchors_file;
    std::string fork_succ_file;
    std::tie(anchors_file, fork_succ_file) = get_anchors_and_fork_fnames(anchors_file_fbase);
    auto write_anchors = [&](std::ofstream &out) {
        IRowDiff row_diff;
        row_diff.load_anchor(anchors_file);
        row_diff.load_fork_succ(fork_succ_file);
        out.write("v2.0", 4);
        row_diff.anchor().serialize(out);
        row_diff.fork_succ().serialize(out);
    };

    convert_to_row_disk(
            files, outfbase, num_threads, mem_bytes,
            utils::make_suffix(outfbase, RowDiffDiskAnnotator::kExtension),
            RowDiffDiskAnnotator::kExtension,
            get_num_rows_from_row_diff_anno(files[0]),
            write_anchors, merge_load_row_diff
    );
}

template <>
void convert_to_row_diff<IntRowDiffDiskAnnotator>(
            const std::vector<std::string> &files,
            const std::string &anchors_file_fbase,
            const std::string &outfbase,
            size_t num_threads,
            size_t /*mem_bytes*/) {

    std::string anchors_file;
    std::string fork_succ_file;
    std::tie(anchors_file, fork_succ_file) = get_anchors_and_fork_fnames(anchors_file_fbase);

    auto outfname = utils::make_suffix(outfbase, IntRowDiffDiskAnnotator::kExtension);

    uint64_t num_rows = get_num_rows_from_column_anno(files[0]);
    size_t num_set_bits = 0;

    std::vector<std::string> col_names;
    std::vector<std::unique_ptr<bit_vector>> columns;
    std::vector<sdsl::int_vector<>> col_values;

    std::mutex mu;

    size_t max_val = 0;

    // TODO: convert with streaming
    ColumnCompressed<>::load_columns_and_values(
            files,
            [&](size_t j,
                    const std::string &label,
                    std::unique_ptr<bit_vector> &&column,
                    sdsl::int_vector<> &&values) {
                size_t local_max_val = 0;
                if (values.size())
                    local_max_val = *std::max_element(values.begin(), values.end());

                std::lock_guard<std::mutex> lock(mu);

                if (local_max_val > max_val)
                    max_val = local_max_val;

                num_set_bits += column->num_set_bits();

                while (columns.size() <= j) {
                    columns.emplace_back();
                    col_values.emplace_back();
                    col_names.emplace_back();
                }
                columns[j] = std::move(column);
                col_values[j] = std::move(values);
                col_names[j] = label;
            },
            num_threads);

    // this must be done after loading all columns
    // to keep their order correct
    // (label_encoder.insert_and_encode cannot be inside get_cols callback)
    LEncoder label_encoder;
    for (const auto &label : col_names) {
        size_t col = label_encoder.insert_and_encode(label);
        if (col + 1 != label_encoder.size()) {
            logger->error("Duplicate columns {}", label);
            exit(1);
        }
    }

    std::ofstream out = utils::open_new_ofstream(outfname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + outfname);

    label_encoder.serialize(out);

    IRowDiff row_diff;
    row_diff.load_anchor(anchors_file);
    row_diff.load_fork_succ(fork_succ_file);
    out.write("v2.0", 4);
    row_diff.anchor().serialize(out);
    row_diff.fork_succ().serialize(out);

    out.close();

    ProgressBar progress_bar(num_rows, "Serialize rows", std::cerr, !common::get_verbose());

    IntRowDisk::serialize(
            outfname,
            [&](std::function<void(const IntMatrix::RowValues &)> write_row_with_values) {
                #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
                for (uint64_t begin = 0 ; begin < num_rows; begin += kNumRowsInBlock) {
                    uint64_t end = std::min(begin + kNumRowsInBlock, num_rows);

                    assert(begin <= end);
                    assert(end <= num_rows);

                    std::vector<IntMatrix::RowValues> rows(end - begin);

                    for (size_t j = 0 ; j < columns.size() ; ++j) {
                        size_t val_idx = begin ? columns[j]->rank1(begin - 1) : 0;
                        columns[j]->call_ones_in_range(begin, end,
                            [&](uint64_t i) {
                                rows[i - begin].emplace_back(j, col_values[j][val_idx++]);
                            }
                        );
                    }

                    #pragma omp ordered
                    {
                        for (const auto &row_with_values : rows) {
                            write_row_with_values(row_with_values);
                            ++progress_bar;
                        }
                    }
                }
            },
            columns.size(), num_rows, num_set_bits, max_val);
}

template <>
void convert_to_row_diff<RowDiffDiskCoordAnnotator>(
            const std::vector<std::string> &files,
            const std::string &anchors_file_fbase,
            const std::string &outfbase,
            size_t num_threads,
            size_t /*mem_bytes*/) {

    std::string anchors_file;
    std::string fork_succ_file;
    std::tie(anchors_file, fork_succ_file) = get_anchors_and_fork_fnames(anchors_file_fbase);

    auto outfname = utils::make_suffix(outfbase, RowDiffDiskCoordAnnotator::kExtension);

    uint64_t num_rows = get_num_rows_from_column_anno(files[0]);
    size_t num_set_bits = 0;

    size_t num_values = 0;

    std::vector<std::string> col_names;
    std::vector<std::unique_ptr<bit_vector>> columns;
    std::vector<sdsl::int_vector<>> col_values;
    std::vector<bit_vector_small> col_delims;

    std::mutex mu;

    size_t max_tuple_size = 0;
    size_t max_val = 0;

    // TODO: convert with streaming
    ColumnCompressed<>::load_columns_delims_and_values(
            files,
            [&](size_t j, const std::string &label,
                    std::unique_ptr<bit_vector> &&column,
                    bit_vector_small&& delims,
                    sdsl::int_vector<> &&values) {
                size_t local_max_val = 0;
                if (values.size())
                    local_max_val = *std::max_element(values.begin(), values.end());

                std::lock_guard<std::mutex> lock(mu);

                if (local_max_val > max_val)
                    max_val = local_max_val;

                if (values.size() > max_tuple_size)
                    max_tuple_size = values.size();

                num_set_bits += column->num_set_bits();

                num_values += values.size();

                while (columns.size() <= j) {
                    columns.emplace_back();
                    col_values.emplace_back();
                    col_delims.emplace_back();
                    col_names.emplace_back();
                }
                columns[j] = std::move(column);
                col_values[j] = std::move(values);
                col_delims[j] = std::move(delims);
                col_names[j] = label;
            },
            num_threads);

    // this must be done after loading all columns
    // to keep their order correct
    // (label_encoder.insert_and_encode cannot be inside get_cols callback)
    LEncoder label_encoder;
    for (const auto &label : col_names) {
        size_t col = label_encoder.insert_and_encode(label);
        if (col + 1 != label_encoder.size()) {
            logger->error("Duplicate columns {}", label);
            exit(1);
        }
    }

    std::ofstream out = utils::open_new_ofstream(outfname);
    if (!out.good())
        throw std::ofstream::failure("Can't write to " + outfname);

    label_encoder.serialize(out);

    IRowDiff row_diff;
    row_diff.load_anchor(anchors_file);
    row_diff.load_fork_succ(fork_succ_file);
    out.write("v2.0", 4);
    row_diff.anchor().serialize(out);
    row_diff.fork_succ().serialize(out);

    out.close();

    ProgressBar progress_bar(num_rows, "Serialize rows", std::cerr, !common::get_verbose());

    CoordRowDisk::serialize(
            outfname,
            [&](std::function<void(const MultiIntMatrix::RowTuples &)> write_row_with_tuples) {
                #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
                for (uint64_t begin = 0 ; begin < num_rows; begin += kNumRowsInBlock) {
                    uint64_t end = std::min(begin + kNumRowsInBlock, num_rows);

                    assert(begin <= end);
                    assert(end <= num_rows);

                    std::vector<MultiIntMatrix::RowTuples> rows(end - begin);

                    for (size_t j = 0 ; j < columns.size() ; ++j) {
                        size_t r = begin ? columns[j]->rank1(begin - 1) + 1 : 1;
                        if (r > columns[j]->num_set_bits())
                            continue;
                        size_t tb = col_delims[j].select1(r) + 1 - r;

                        columns[j]->call_ones_in_range(begin, end, [&](uint64_t i) {
                            size_t te = col_delims[j].select1(r + 1) - r;
                            rows[i - begin].emplace_back(j,
                                    MultiIntMatrix::Tuple(col_values[j].begin() + tb,
                                                          col_values[j].begin() + te));
                            tb = te;
                            ++r;
                        });
                    }

                    #pragma omp ordered
                    {
                        for (const auto& row_with_tuples : rows) {
                            write_row_with_tuples(row_with_tuples);
                            ++progress_bar;
                        }
                    }
                }
            },
            columns.size(), num_rows, num_set_bits, num_values, max_val, max_tuple_size);
}

template <>
std::unique_ptr<RbBRWTAnnotator>
convert<RbBRWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    auto &columns = const_cast<std::vector<std::unique_ptr<bit_vector>>&>(
        annotator.get_matrix().data()
    );
    auto matrix = convert_to_RainbowBRWT(
        [&](const auto &callback) {
            for (auto &column : columns) {
                callback(std::move(column));
            }
        }
    );
    return std::make_unique<RbBRWTAnnotator>(std::move(matrix),
                                             annotator.get_label_encoder());
}

template <>
std::unique_ptr<BinRelWTAnnotator>
convert<BinRelWTAnnotator, std::string>(ColumnCompressed<std::string>&& annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<BinRelWT>(
        [&](auto callback) {
            utils::call_rows(annotator.get_matrix().data(), callback);
        },
        num_set_bits,
        num_columns
    );

    return std::make_unique<BinRelWTAnnotator>(std::move(matrix),
                                                    annotator.get_label_encoder());
}

template <typename Label>
void merge_rows(const std::vector<LabelEncoder<Label>> &label_encoders,
                std::function<const BinaryMatrix::SetBitPositions(uint64_t)> get_next_row,
                uint64_t num_rows,
                const std::string &outfile) {
    LabelEncoder<Label> merged_label_enc;
    std::vector<std::vector<uint64_t> > label_mappings;

    for (const auto &label_encoder : label_encoders) {
        merged_label_enc.merge(label_encoder);

        std::vector<uint64_t> v;
        for (size_t j = 0; j < label_encoder.size(); ++j) {
            v.push_back(merged_label_enc.encode(label_encoder.decode(j)));
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

void merge_row_compressed(const std::vector<std::string> &filenames,
                          const std::string &outfile) {
    assert(filenames.size() && "nothing to merge");

    uint64_t num_rows;
    uint64_t num_relations;

    RowCompressed<>::read_shape(filenames.at(0), &num_rows, &num_relations);
    assert(num_rows);

    std::vector<LEncoder> label_encoders;
    std::vector<std::unique_ptr<StreamRows<>>> streams;

    for (auto filename : filenames) {
        if (!utils::ends_with(filename, RowCompressed<>::kExtension))
            throw std::runtime_error("Can't merge annotations of mixed types");

        label_encoders.emplace_back(RowCompressed<>::read_label_encoder(filename));

        streams.emplace_back(new StreamRows<>(RowCompressed<>::get_row_streamer(filename)));
    }

    merge_rows(
        label_encoders,
        [&](uint64_t annotator_idx) -> const BinaryMatrix::SetBitPositions {
            return *streams[annotator_idx]->next_row();
        },
        num_rows,
        outfile
    );
}

void merge_brwt(const std::vector<std::string> &filenames,
                const std::string &outfile) {
    assert(filenames.size() && "nothing to merge");

    uint64_t num_rows = 0;

    LEncoder label_encoder;

    std::vector<BRWT> brwts;

    for (auto filename : filenames) {
        MultiBRWTAnnotator annotator;
        if (!annotator.load(filename)) {
            logger->error("Cannot load annotations from file '{}'", filename);
            exit(1);
        }

        if (filename == filenames[0])
            num_rows = annotator.num_objects();

        if (annotator.num_objects() != num_rows)
            throw std::runtime_error("Annotators have different number of rows");

        for (const auto &label : annotator.get_label_encoder().get_labels()) {
            if (label_encoder.label_exists(label))
                throw std::runtime_error("merging of BRWT with same labels is not implemented");

            label_encoder.insert_and_encode(label);
        }

        brwts.push_back(std::move(const_cast<BRWT&>(annotator.get_matrix())));
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
                              const std::string &outfbase) {
    RowCompressed<Label>::serialize(outfbase, annotator.get_label_encoder(),
        [&](BinaryMatrix::RowCallback write_row) {
            utils::call_rows(annotator.get_matrix().data(), write_row, annotator.num_objects());
        }
    );
}

template
void convert_to_row_annotator(const ColumnCompressed<std::string> &annotator,
                              const std::string &outfbase);

fs::path get_count_filename(fs::path vector_fname,
                            const std::string &extension,
                            fs::path out_dir,
                            fs::path fallback_filename) {
    if (vector_fname.extension() != extension) {
        if (!vector_fname.empty()) {
            logger->warn(
                "The output filename {} does not have the expected extension"
                " '{}'. The name of the row count vector will be derived"
                " from the first column in batch.", vector_fname, extension);
        }

        // fix the name, set to the fallback value
        vector_fname = out_dir/fs::path(fallback_filename)
                                .filename()
                                .replace_extension()
                                .replace_extension(extension);

        logger->info("The name of the row count vector was automatically"
                     " derived from the first column in batch and set to {}",
                     vector_fname);

        if (fs::exists(vector_fname)) {
            logger->warn("Found row count vector {}, will be overwritten",
                         vector_fname);
            fs::remove(vector_fname);
        }
    }

    return vector_fname;
}

void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string &graph_fname,
                         size_t mem_bytes,
                         uint32_t max_path_length,
                         fs::path out_dir,
                         fs::path swap_dir,
                         RowDiffStage construction_stage,
                         fs::path count_vector_fname,
                         bool with_values,
                         bool with_coordinates,
                         size_t num_coords_per_seq) {
    assert(!with_values || !with_coordinates);

    if (out_dir.empty())
        out_dir = "./";

    if (construction_stage != RowDiffStage::COUNT_LABELS)
        build_pred_succ(graph_fname, graph_fname, out_dir,
                        ".row_count", get_num_threads());

    if (construction_stage == RowDiffStage::CONVERT) {
        assign_anchors(graph_fname, graph_fname, out_dir, max_path_length,
                       ".row_reduction", get_num_threads());

        const std::string anchors_fname = graph_fname + kRowDiffAnchorExt;
        if (!fs::exists(anchors_fname)) {
            logger->error("Can't find anchors bitmap at {}", anchors_fname);
            exit(1);
        }
        if (!utils::with_mmap()) { // only reserve space for anchors with no mmap
            uint64_t anchor_size = fs::file_size(anchors_fname);
            if (anchor_size > mem_bytes) {
                logger->warn("Anchor bitmap ({} MB) is larger than the memory"
                             " allocated ({} MB). Reserve more RAM.",
                             anchor_size / 1e6, mem_bytes / 1e6);
                return;
            }
            mem_bytes -= anchor_size;
        }
    }

    if (!files.size())
        return;

    // load as many columns as we can fit in memory, and convert them
    for (uint32_t i = 0; i < files.size(); ) {
        logger->trace("Loading columns for batch-conversion...");
        size_t mem_bytes_left = mem_bytes;
        std::vector<std::string> file_batch;
        for ( ; i < files.size(); ++i) {
            // also add some space for buffers for each column
            uint64_t file_size = (utils::with_mmap() ? 0 : fs::file_size(files[i])) + ROW_DIFF_BUFFER_BYTES;
            if (with_values && !utils::with_mmap()) {
                // also add k-mer counts
                try {
                    const auto &values_fname
                        = utils::remove_suffix(files[i], ColumnCompressed<>::kExtension)
                                                    + ColumnCompressed<>::kCountExtension;
                    file_size += fs::file_size(values_fname);
                } catch (...) {
                    // Count vectors may be missing for empty annotations. If a count file
                    // is missing for a non-empty annotation, the error will be thrown later
                    // in convert_batch_to_row_diff, so we skip it here in any case.
                }
            } else if (with_coordinates && !utils::with_mmap()) {
                // also add k-mer coordinates
                try {
                    const auto &coord_fname
                        = utils::remove_suffix(files[i], ColumnCompressed<>::kExtension)
                                                    + ColumnCompressed<>::kCoordExtension;
                    file_size += fs::file_size(coord_fname);
                } catch (...) {
                    // Attribute vectors may be missing for empty annotations
                }
            }
            if (file_size > mem_bytes) {
                logger->error("Not enough memory to load {}, requires {} MB",
                              files[i], file_size / 1e6);
                exit(1);
            }
            if (file_size > mem_bytes_left)
                break;

            mem_bytes_left -= file_size;
            file_batch.push_back(files[i]);

            // get a file name of the count vector or derive from the first file in batch
            if (construction_stage == RowDiffStage::COUNT_LABELS) {
                count_vector_fname = get_count_filename(count_vector_fname, ".row_count",
                                                        out_dir, file_batch.front());
            } else if (construction_stage == RowDiffStage::COMPUTE_REDUCTION) {
                count_vector_fname = get_count_filename(count_vector_fname, ".row_reduction",
                                                        out_dir, file_batch.front());
            }
        }

        Timer timer;
        logger->trace("Annotations in batch: {}", file_batch.size());

        if (construction_stage == RowDiffStage::COUNT_LABELS) {
            count_labels_per_row(file_batch, count_vector_fname, with_coordinates);
        } else {
            convert_batch_to_row_diff(graph_fname,
                    file_batch, out_dir, swap_dir, count_vector_fname, ROW_DIFF_BUFFER_BYTES,
                    construction_stage == RowDiffStage::COMPUTE_REDUCTION,
                    with_values, with_coordinates, num_coords_per_seq);
        }

        logger->trace("Batch processed in {} sec", timer.elapsed());
    }
}

void convert_row_diff_to_col_compressed(const std::vector<std::string> &files,
                                        const std::string &outfbase) {
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (uint32_t i = 0; i < files.size(); ++i) {
        std::string file = files[i];
        RowDiffColumnAnnotator input_anno;
        input_anno.load(file);

        std::string outname = utils::remove_suffix(fs::path(file).filename(),
                                                   RowDiffColumnAnnotator::kExtension);
        std::string out_path
                = (fs::path(outfbase).remove_filename() / outname).string()
                + "_row_diff" + ColumnCompressed<std::string>::kExtension;
        std::ofstream out(out_path, std::ios::binary);
        logger->trace("Transforming {} to {}", file, out_path);

        serialize_number(out, input_anno.num_objects());
        input_anno.get_label_encoder().serialize(out);
        input_anno.get_matrix().diffs().serialize(out);
        out.close();
    }
}

template <class Annotator>
StaticBinRelAnnotator<TupleCSCMatrix<typename Annotator::binary_matrix_type>, std::string>
load_coords(Annotator&& anno, const std::vector<std::string> &files) {
    std::vector<bit_vector_smart> delimiters(anno.num_labels());
    std::vector<sdsl::int_vector<>> column_values(anno.num_labels());

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < files.size(); ++i) {
        auto label_encoder = ColumnCompressed<>::read_label_encoder(files[i]);

        auto coords_fname = utils::remove_suffix(files[i], ColumnCompressed<>::kExtension)
                                                        + ColumnCompressed<>::kCoordExtension;
        std::unique_ptr<std::ifstream> in = utils::open_ifstream(coords_fname);
        size_t j = 0;
        try {
            TupleCSC::load_tuples(*in, label_encoder.size(), [&](auto&& delims, auto&& values) {
                size_t idx;
                try {
                    idx = anno.get_label_encoder().encode(label_encoder.decode(j));
                } catch (...) {
                    logger->error("Label '{}' from {} is missing in the target annotator",
                                  label_encoder.decode(j), files[i]);
                    exit(1);
                }
                if (delimiters[idx].size()) {
                    logger->error("Merging coordinate annotations with overlapping"
                                  " labels is not implemented");
                    exit(1);
                }
                delimiters[idx] = std::move(delims);
                column_values[idx] = std::move(values);
                j++;
            });
        } catch (const std::exception &e) {
            logger->error("Couldn't load coordinates from {}\nException: {}", coords_fname, e.what());
            exit(1);
        } catch (...) {
            logger->error("Couldn't load coordinates from {}", coords_fname);
            exit(1);
        }
        assert(j == label_encoder.size());
    }

    auto label_encoder = anno.get_label_encoder();

    return StaticBinRelAnnotator<TupleCSCMatrix<typename Annotator::binary_matrix_type>, std::string>(
            std::make_unique<TupleCSCMatrix<typename Annotator::binary_matrix_type>>(
                    std::move(*anno.release_matrix()),
                    std::move(delimiters),
                    std::move(column_values)),
            std::move(label_encoder));
}

template
StaticBinRelAnnotator<TupleCSC, std::string>
load_coords<ColumnCompressed<>>(ColumnCompressed<>&&, const std::vector<std::string> &);

template
StaticBinRelAnnotator<TupleBRWT, std::string>
load_coords<MultiBRWTAnnotator>(MultiBRWTAnnotator&&, const std::vector<std::string> &);

} // namespace annot
} // namespace mtg
