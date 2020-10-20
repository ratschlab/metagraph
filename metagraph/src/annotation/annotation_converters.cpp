#include "annotation_converters.hpp"

#include <cassert>
#include <vector>
#include <functional>
#include <filesystem>

#include <ips4o.hpp>
#include <progress_bar.hpp>
#include <tsl/hopscotch_map.h>

#include "annotation/row_diff_builder.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/hashers/hash.hpp"
#include "common/sorted_sets/sorted_multiset.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "binary_matrix/row_vector/vector_row_binmat.hpp"
#include "binary_matrix/multi_brwt/brwt_builders.hpp"
#include "binary_matrix/multi_brwt/clustering.hpp"
#include "representation/annotation_matrix/static_annotators_def.hpp"
#include "representation/column_compressed/annotate_column_compressed.hpp"
#include "representation/row_compressed/annotate_row_compressed.hpp"


namespace mtg {
namespace annot {

using namespace mtg::annot::binmat;

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

    } else if constexpr(std::is_same_v<MatrixType, RowSparse>) {
        matrix = std::make_unique<MatrixType>(call_rows, label_encoder->size(), num_rows, num_relations);

    } else {
        static_assert(utils::dependent_false<StaticAnnotation>::value);
    }

    return std::make_unique<StaticAnnotation>(std::move(matrix), *label_encoder);
}

template std::unique_ptr<RowFlatAnnotator> convert(const std::string &filename);
template std::unique_ptr<RowSparseAnnotator> convert(const std::string &filename);
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
std::unique_ptr<RowSparseAnnotator>
convert<RowSparseAnnotator, std::string>(ColumnCompressed<std::string> &&annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowSparse>(
            [&](auto callback) {
                utils::RowsFromColumnsTransformer(annotator.get_matrix().data()).call_rows(callback);
            },
            num_columns, num_rows, num_set_bits);

    return std::make_unique<RowSparseAnnotator>(std::move(matrix),
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

std::unique_ptr<RowDiffBRWTAnnotator>
convert_row_diff_to_BRWT(RowDiffAnnotator &&annotator,
                         BRWTBottomUpBuilder::Partitioner partitioning,
                         size_t num_parallel_nodes,
                         size_t num_threads) {
    // we are going to take the columns from the annotator and thus
    // have to replace them with empty columns to keep the structure valid
    const graph::DBGSuccinct* graph = annotator.get_matrix().graph();
    std::vector<std::unique_ptr<bit_vector>> columns
            = annotator.release_matrix()->diffs().release_columns();

    auto matrix = std::make_unique<BRWT>(
            BRWTBottomUpBuilder::build(std::move(columns),
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
convert_to_greedy_BRWT(RowDiffAnnotator &&annotation,
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

std::unique_ptr<RowDiffBRWTAnnotator> convert_to_simple_BRWT(RowDiffAnnotator&& annotation,
                                                             size_t grouping_arity,
                                                             size_t num_parallel_nodes,
                                                             size_t num_threads) {
    return convert_row_diff_to_BRWT(
            std::move(annotation),
            BRWTBottomUpBuilder::get_basic_partitioner(grouping_arity),
            num_parallel_nodes, num_threads);
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
                          " Each line must contain exactly 4 values:"
                          " <cluster 1> <cluster 2> <dist> <cluster 3>"
                          "\nException: {}", e.what());
            exit(1);
        }
    }

    return linkage;
}

std::unique_ptr<RowSparseAnnotator> convert(const RowDiffAnnotator &annotator) {
    uint64_t num_set_bits = annotator.num_relations();
    uint64_t num_rows = annotator.num_objects();
    uint64_t num_columns = annotator.num_labels();

    auto matrix = std::make_unique<RowSparse>(
            [&](auto callback) {
                utils::RowsFromColumnsTransformer(annotator.get_matrix().diffs().data())
                        .call_rows(callback);
            },
            num_columns, num_rows, num_set_bits);

    return std::make_unique<RowSparseAnnotator>(std::move(matrix),
                                                annotator.get_label_encoder());
}


std::unique_ptr<MultiBRWTAnnotator>
convert_to_BRWT(
        const std::string &linkage_matrix_file,
        size_t num_parallel_nodes,
        size_t num_threads,
        const std::filesystem::path &tmp_path,
        const std::function<void(const BRWTBottomUpBuilder::CallColumn &)> &get_columns,
        std::vector<std::pair<uint64_t, std::string>> &&column_names) {
    auto linkage = parse_linkage_matrix(linkage_matrix_file);

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
        label_encoder.insert_and_encode(label);
    }

    return std::make_unique<MultiBRWTAnnotator>(std::move(matrix), label_encoder);
}

template <>
std::unique_ptr<MultiBRWTAnnotator> convert_to_BRWT<MultiBRWTAnnotator>(
        const std::vector<std::string> &annotation_files,
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
    return convert_to_BRWT(linkage_matrix_file, num_parallel_nodes, num_threads, tmp_path,
                    get_columns, std::move(column_names));
}

template<>
std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_BRWT<RowDiffBRWTAnnotator>(const std::vector<std::string> &annotation_files,
                         const std::string &linkage_matrix_file,
                         size_t num_parallel_nodes,
                         size_t num_threads,
                         const std::filesystem::path &tmp_path) {
    std::vector<std::pair<uint64_t, std::string>> column_names;
    std::mutex mu;

    auto get_columns = [&](const BRWTBottomUpBuilder::CallColumn &call_column) {
        for (const auto &fname : annotation_files) {
            RowDiffAnnotator annotator;
            if (!annotator.merge_load({fname})) {
                logger->error("Could not load {}", fname);
                std::exit(1);
            }
            std::vector<std::unique_ptr<bit_vector>> cols
                    = std::move(annotator.release_matrix()->diffs().release_columns());
            for (uint32_t idx = 0; idx < cols.size(); ++idx) {
                std::string label = annotator.get_label_encoder().get_labels()[idx];
                call_column(idx, std::move(cols[idx]));

                std::lock_guard<std::mutex> lock(mu);
                column_names.emplace_back(idx, label);
            };
        }
    };

    std::unique_ptr<MultiBRWTAnnotator> annotator
            = convert_to_BRWT(linkage_matrix_file, num_parallel_nodes, num_threads,
                              tmp_path, get_columns, std::move(column_names));

    return std::make_unique<RowDiffBRWTAnnotator>(
            std::make_unique<RowDiff<BRWT>>(nullptr, std::move(*annotator->release_matrix())),
            annotator->get_label_encoder());
}

void relax_BRWT(binmat::BRWT *annotation, size_t relax_max_arity, size_t num_threads) {
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
            [&](uint64_t /*column_index*/,
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

void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string& graph_fname,
                         size_t mem_bytes,
                         uint32_t max_path_length,
                         const std::filesystem::path &dest_dir) {
    logger->trace("Loading graph...");
    graph::DBGSuccinct graph(2);
    bool result = graph.load(graph_fname);
    if (!result) {
        logger->error("Cannot load graph from {}", graph_fname);
        std::exit(1);
    }

    mem_bytes -= std::filesystem::file_size(graph_fname);
    // load as many columns as we can fit in memory, and convert them
    for (uint32_t i = 0; i < files.size();) {
        logger->trace("Loading columns for batch-conversion...");
        size_t cur_mem_bytes = mem_bytes;
        std::vector<std::string> file_batch;
        for (; i < files.size(); ++i) {
            // *2 in order to account for constructing the sparsified column
            size_t file_size = 2 * std::filesystem::file_size(files[i]);
            if (file_size > mem_bytes) {
                logger->warn(
                        "Not enough memory to process {}, requires {} MB",
                        files[i], file_size/1e6);
                continue;
            }
            if (file_size > cur_mem_bytes || file_batch.size() >= 15'000)
                break;

            cur_mem_bytes -= file_size;
            file_batch.push_back(files[i]);
        }

        Timer timer;
        logger->trace("Starting converting column-batch with {} columns ...",
                      file_batch.size());
        convert_batch_to_row_diff(graph, graph_fname, file_batch, dest_dir, max_path_length);
        logger->trace("Column-batch converted in {} sec", timer.elapsed());
    }

}

void convert_row_diff_to_col_compressed(const std::vector<std::string> &files,
                                        const std::string &outfbase) {
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (uint32_t i = 0; i < files.size(); ++i) {
        std::string file = files[i];
        RowDiffAnnotator input_anno;
        input_anno.load(file);

        std::string outname = utils::remove_suffix(std::filesystem::path(file).filename(), RowDiffAnnotator::kExtension);
        std::string out_path
                = (std::filesystem::path(outfbase).remove_filename() / outname).string()
                + "_row_diff" + ColumnCompressed<std::string>::kExtension;
        std::ofstream outstream(out_path, std::ios::binary);
        logger->trace("Transforming {} to {}", file, out_path);

        serialize_number(outstream, input_anno.num_objects());
        input_anno.get_label_encoder().serialize(outstream);
        input_anno.get_matrix().diffs().serialize(outstream);
        outstream.close();
    }
}

template <class BinaryMatrix>
void wrap_in_row_diff(const std::string &anno_file,
                      const std::string &graph_file,
                      const std::string &out_dir) {
    if (!std::filesystem::is_directory(out_dir)) {
        logger->error("Output path must be a directory, not a file: {}", out_dir);
        std::exit(1);
    }
    using std::filesystem::path;
    path file_name = path(anno_file).filename().replace_extension("");
    std::string old_extension = file_name.extension();
    path out_file = path(out_dir)/ (file_name.replace_extension("").string() + "row_diff_" + old_extension
               + ".annodbg");
    std::ofstream out(out_file, ios::binary);
    uint64_t v = 0;
    out.write(reinterpret_cast<char *>(&v), sizeof(uint64_t));
    out.close();
    std::system(("cat " + anno_file + " >> out_file").c_str());
}


} // namespace annot
} // namespace mtg
