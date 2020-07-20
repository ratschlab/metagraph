#ifndef __METHOD_CONSTRUCTORS_HPP__
#define __METHOD_CONSTRUCTORS_HPP__

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "annotation/binary_matrix/multi_brwt/clustering.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_vector/vector_row_binmat.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "common/data_generation.hpp"


namespace mtg {
namespace experiments {

template <typename T>
using UniquePtrs = typename std::vector<std::unique_ptr<T>>;


enum class MatrixType {
    ROW = 0,
    COLUMN,
    BRWT,
    BRWT_EXTRA,
    BIN_REL_WT_SDSL,
    BIN_REL_WT,
    ROW_FLAT,
    RAINBOWFISH
};

MatrixType string_to_matrix_type(const std::string &string) {
    if (string == "column") {
        return MatrixType::COLUMN;

    } else if (string == "row") {
        return MatrixType::ROW;

    } else if (string == "brwt") {
        return MatrixType::BRWT;

    } else if (string == "brwt_extra") {
        return MatrixType::BRWT_EXTRA;

    } else if (string == "bin_rel_wt_sdsl") {
        return MatrixType::BIN_REL_WT_SDSL;

    } else if (string == "bin_rel_wt") {
        return MatrixType::BIN_REL_WT;

    } else if (string == "flat") {
        return MatrixType::ROW_FLAT;

    } else if (string == "rbfish") {
        return MatrixType::RAINBOWFISH;

    } else {
        std::cerr << "Error: unknown annotation representation" << std::endl;
        exit(1);
    }
}

std::unique_ptr<anno::binmat::BinaryMatrix>
matrix_type_to_data(const std::string &file, MatrixType type) {
    std::unique_ptr<anno::binmat::BinaryMatrix> matrix_ptr;
    if (type == MatrixType::COLUMN) {
        matrix_ptr.reset(new anno::binmat::ColumnMajor());
    } else if (type == MatrixType::ROW) {
        matrix_ptr.reset(new anno::binmat::VectorRowBinMat());
    } else if (type == MatrixType::BRWT) {
        matrix_ptr.reset(new anno::binmat::BRWT());
    } else if (type == MatrixType::BRWT_EXTRA) {
        matrix_ptr.reset(new anno::binmat::BRWT());
    } else if (type == MatrixType::BIN_REL_WT_SDSL) {
        matrix_ptr.reset(new anno::binmat::BinRelWT_sdsl());
    } else if (type == MatrixType::BIN_REL_WT) {
        matrix_ptr.reset(new anno::binmat::BinRelWT());
    } else if (type == MatrixType::ROW_FLAT) {
        matrix_ptr.reset(new anno::binmat::RowConcatenated<>());
    } else if (type == MatrixType::RAINBOWFISH) {
        matrix_ptr.reset(new anno::binmat::Rainbowfish());
    } else {
        std::cerr << "Error: invalid matrix type" << std::endl;
        exit(1);
    }
    std::ifstream in(file, std::ios::binary);
    try {
        matrix_ptr->load(in);
    } catch (...) {
        std::cerr << "Error: invalid file or mismatched matrix type" << std::endl;
        exit(1);
    }

    return matrix_ptr;
}


template <typename BitVector>
UniquePtrs<bit_vector> convert_to(UniquePtrs<bit_vector>&& input) {
    UniquePtrs<bit_vector> result;
    result.reserve(input.size());

    for (auto &vector : input) {
        result.emplace_back(
            new BitVector(vector->convert_to<BitVector>())
        );
        vector.reset();
    }
    return result;
}

std::unique_ptr<anno::binmat::BinaryMatrix>
generate_brwt_from_rows(std::vector<std::unique_ptr<bit_vector>>&& columns,
                        size_t arity, bool greedy, size_t relax_arity_limit) {
    std::unique_ptr<anno::binmat::BRWT> binary_matrix;

    if (greedy) {
        binary_matrix = std::make_unique<anno::binmat::BRWT>(
            anno::binmat::BRWTBottomUpBuilder::build(std::move(columns),
                [](const auto &columns) {
                    std::vector<sdsl::bit_vector> subvectors
                        = anno::binmat::random_submatrix(columns, 1'000'000);
                    return anno::binmat::greedy_matching(subvectors);
                }
            )
        );
    } else {
        binary_matrix = std::make_unique<anno::binmat::BRWT>(
            anno::binmat::BRWTBottomUpBuilder::build(
                std::move(columns),
                anno::binmat::BRWTBottomUpBuilder::get_basic_partitioner(arity)
            )
        );
    }

    if (relax_arity_limit > 1)
        anno::binmat::BRWTOptimizer::relax(binary_matrix.get(), relax_arity_limit);

    return binary_matrix;
}

template <class... Args>
std::unique_ptr<anno::binmat::BinaryMatrix>
generate_from_rows(std::vector<std::unique_ptr<bit_vector>>&& columns,
                   MatrixType type, Args&&... args) {
    std::unique_ptr<anno::binmat::BinaryMatrix> binary_matrix;

    switch (type) {
        case MatrixType::ROW: {
            anno::binmat::VectorRowBinMat mat(columns[0]->size());
            for (size_t i = 0; i < columns.size(); ++i) {
                columns[i]->call_ones([&mat,i](auto p) { mat.set(p, i); });
            }
            binary_matrix.reset(new anno::binmat::VectorRowBinMat(std::move(mat)));
            break;
        }
        case MatrixType::COLUMN: {
            binary_matrix.reset(new anno::binmat::ColumnMajor(
                convert_to<bit_vector_sd>(std::move(columns))
            ));
            break;
        }
        case MatrixType::BRWT: {
            binary_matrix.reset(new anno::binmat::BRWT(
                anno::binmat::BRWTBottomUpBuilder::build(std::move(columns))
            ));
            break;
        }
        case MatrixType::BRWT_EXTRA: {
            binary_matrix.reset(new anno::binmat::BRWT(
                anno::binmat::BRWTBottomUpBuilder::build(std::move(columns),
                    [](const auto &columns) {
                        std::vector<sdsl::bit_vector> subvectors
                            = anno::binmat::random_submatrix(columns, 1'000'000);
                        return anno::binmat::greedy_matching(subvectors);
                    }
                )
            ));
            break;
        }
        case MatrixType::BIN_REL_WT_SDSL: {
            uint64_t num_set_bits = 0;
            uint64_t num_columns = columns.size();

            for (const auto &vector_ptr : columns) {
                num_set_bits += vector_ptr->num_set_bits();
            }

            binary_matrix.reset(new anno::binmat::BinRelWT_sdsl(
                [&](const auto &callback) {
                    utils::RowsFromColumnsTransformer(columns).call_rows(callback);
                },
                num_set_bits, num_columns
            ));

            break;
        }
        case MatrixType::BIN_REL_WT: {
            binary_matrix.reset(new anno::binmat::BinRelWT(std::move(columns)));
            break;
        }
        case MatrixType::ROW_FLAT: {
            const auto num_columns = columns.size();
            const auto num_rows = columns[0]->size();

            uint64_t num_set_bits = 0;
            for (const auto &vector_ptr : columns) {
                num_set_bits += vector_ptr->num_set_bits();
            }

            binary_matrix.reset(new anno::binmat::RowConcatenated<>(
                [&](const auto &callback) {
                    utils::RowsFromColumnsTransformer(columns).call_rows(callback);
                },
                num_columns,
                num_rows,
                num_set_bits,
                std::forward<Args>(args)...
            ));

            break;
        }
        case MatrixType::RAINBOWFISH: {
            const auto num_columns = columns.size();

            binary_matrix.reset(new anno::binmat::Rainbowfish(
                [&](const auto &callback) {
                    utils::RowsFromColumnsTransformer(columns).call_rows(callback);
                },
                num_columns,
                std::forward<Args>(args)...
            ));
            break;
        }
        default: {
            throw std::runtime_error("ERROR: Unsupported compression method.");
        }
    }

    return binary_matrix;
}

void improve_linear_order(std::vector<std::unique_ptr<bit_vector>> *columns_ptr) {
    assert(columns_ptr);
    auto &columns = *columns_ptr;
    assert(columns[0]->size());
    uint64_t num_rows = columns[0]->size();

    // measure Hamming distance matrix
    std::vector<std::vector<uint64_t>> matches;
    std::vector<std::pair<uint64_t, uint64_t>> best_match;
    best_match.insert(best_match.end(), columns.size(), std::make_pair(0, -1));
    matches.reserve(columns.size());
    for (uint64_t i = 0; i < columns.size(); ++i) {
        matches.emplace_back();
        matches.back().reserve(columns.size() - i);
        matches.back().emplace_back(num_rows);

        auto bv_i = columns[i]->copy_to<sdsl::bit_vector>();

        for (uint64_t j = i + 1; j < columns.size(); ++j) {
            auto bv = bv_i;
            bv ^= columns[j]->copy_to<sdsl::bit_vector>();
            matches.back().emplace_back(bv.size() - sdsl::util::cnt_one_bits(bv));
            if (matches.back().back() < best_match[i].second) {
                best_match[i].first = j;
                best_match[i].second = matches.back().back();

                best_match[j].first = i;
                best_match[j].second = matches.back().back();
            }
        }
    }

    std::vector<std::unique_ptr<bit_vector>> columns_new;
    sdsl::bit_vector match_chosen(columns.size(), true);
    for (uint64_t i = 0; i < best_match.size(); ++i) {
        if (!match_chosen[i])
            continue;

        columns_new.emplace_back(std::move(columns[i]));
        match_chosen[i] = false;
        if (match_chosen[best_match[i].first]) {
            columns_new.emplace_back(std::move(columns[best_match[i].first]));
            match_chosen[best_match[i].first] = false;
        }
    }

    std::swap(*columns_ptr, columns_new);
}

template <class BitVector>
std::vector<std::unique_ptr<bit_vector>>
subsample_rows(const std::vector<std::unique_ptr<BitVector>> &source_columns,
               uint64_t n,
               uint64_t num_selected_rows,
               uint64_t num_threads = 1) {
    DataGenerator generator;
    generator.set_seed(42);
    auto selector = generator.generate_random_column_fixed_size(n, num_selected_rows);
    assert(selector->num_set_bits() == num_selected_rows);

    std::vector<std::unique_ptr<bit_vector>> columns(source_columns.size());

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t i = 0; i < columns.size(); ++i) {
        sdsl::bit_vector column(num_selected_rows, false);
        selector->call_ones(
            [&](uint64_t pos) {
                if ((*source_columns[i])[pos])
                    column[selector->rank1(pos)] = true;
            }
        );
        columns[i].reset(new bit_vector_small(std::move(column)));
    }

    return columns;
}

} // namespace experiments
} // namespace mtg

#endif // __METHOD_CONSTRUCTORS_HPP__
