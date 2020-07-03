#include <random>

#include "gtest/gtest.h"

#include "test_matrix_helpers.hpp"

#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

template <typename BinMat>
class BinaryMatrixTest : public ::testing::Test { };
typedef ::testing::Types<BRWT,
                         BRWTOptimized,
                         ColumnMajor,
                         BinRelWT,
                         BinRelWT_sdsl,
                         RowConcatenated<>,
                         UniqueRowBinmat,
                         Rainbow<BRWT>,
                         Rainbowfish> BinMatTypes;
TYPED_TEST_SUITE(BinaryMatrixTest, BinMatTypes);


TYPED_TEST(BinaryMatrixTest, DefaultConstructor) {
    TypeParam matrix;
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TYPED_TEST(BinaryMatrixTest, DefaultConstructorFromEmptyList) {
    TypeParam matrix {};
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TYPED_TEST(BinaryMatrixTest, BuildEmpty) {
    {
        auto matrix = build_matrix_from_columns<TypeParam>();
        EXPECT_EQ(0u, matrix.num_columns());
        EXPECT_EQ(0u, matrix.num_rows());
    }
}

TYPED_TEST(BinaryMatrixTest, BuildOneColumn) {
    BitVectorPtrArray columns, copy;
    columns.emplace_back(new bit_vector_stat(10, true));
    copy.emplace_back(new bit_vector_stat(10, true));

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), 10), columns);
}

TYPED_TEST(BinaryMatrixTest, BuildOneBigColumn) {
    BitVectorPtrArray columns, copy;
    columns.emplace_back(new bit_vector_stat(100'000, true));
    copy.emplace_back(new bit_vector_stat(100'000, true));

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), 100'000), columns);
}

TYPED_TEST(BinaryMatrixTest, AllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllOne100Rows) {
    for (size_t num_rows = 1; num_rows < 100; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 5; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed1) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {

                sdsl::bit_vector bv(num_rows, 0);

                for (size_t i = 0; i < num_rows; ++i) {
                    bv[i] = (i + j) % 2;
                }

                columns.emplace_back(new bit_vector_stat(std::move(bv)));

                copy.push_back(columns.back()->copy());
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed1ManyRowsOdd) {
    for (size_t num_rows = 1; num_rows < 100; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 8; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {

                sdsl::bit_vector bv(num_rows, 0);

                for (size_t i = 0; i < num_rows; ++i) {
                    bv[i] = (i + j * i) % 2;
                }

                columns.emplace_back(new bit_vector_stat(std::move(bv)));

                copy.push_back(columns.back()->copy());
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed1ManyRowsEven) {
    for (size_t num_rows = 1; num_rows < 100; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 8; ++num_columns) {
            BitVectorPtrArray columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {

                sdsl::bit_vector bv(num_rows, 0);

                for (size_t i = 0; i < num_rows; ++i) {
                    bv[i] = !((i + j * i) % 2);
                }

                columns.emplace_back(new bit_vector_stat(std::move(bv)));

                copy.push_back(columns.back()->copy());
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed2) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy;

            columns.emplace_back(new bit_vector_stat(num_rows));

            for (size_t i = 1; i < num_columns; ++i) {
                sdsl::bit_vector bv(num_rows, 0);
                for (size_t j = 0; j < num_rows; ++j) {
                    bv[j] = (i + j) % 2;
                }
                columns.emplace_back(new bit_vector_stat(std::move(bv)));
            }

            for (const auto &column : columns) {
                copy.push_back(column->copy());
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed3) {
    size_t num_rows = 7;
    std::vector<sdsl::bit_vector> bvs(2, sdsl::bit_vector(num_rows, 0));
    BitVectorPtrArray columns, copy;

    bvs[0][1] = true;
    bvs[0][2] = true;
    bvs[0][5] = true;
    bvs[1][1] = true;
    bvs[1][2] = true;
    bvs[1][3] = true;
    bvs[1][4] = true;
    bvs[1][5] = true;

    columns.emplace_back(new bit_vector_stat(std::move(bvs[0])));
    columns.emplace_back(new bit_vector_stat(std::move(bvs[1])));

    for (const auto &column : columns) {
        copy.push_back(column->copy());
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed4) {
    size_t num_rows = 4;
    std::vector<sdsl::bit_vector> bvs(2, sdsl::bit_vector(num_rows, 0));
    BitVectorPtrArray columns, copy;

    bvs[0][0] = true;
    bvs[0][1] = true;
    bvs[0][2] = true;
    bvs[0][3] = true;
    bvs[1][0] = true;
    bvs[1][1] = true;
    bvs[1][2] = true;
    bvs[1][3] = true;

    columns.emplace_back(new bit_vector_stat(std::move(bvs[0])));
    columns.emplace_back(new bit_vector_stat(std::move(bvs[1])));

    for (const auto &column : columns) {
        copy.push_back(column->copy());
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows),columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed5) {
    size_t num_rows = 4;
    std::vector<sdsl::bit_vector> bvs(2, sdsl::bit_vector(num_rows, 0));
    BitVectorPtrArray columns, copy;

    bvs[0][0] = true;
    bvs[0][1] = true;
    bvs[0][2] = true;
    bvs[0][3] = true;

    columns.emplace_back(new bit_vector_stat(std::move(bvs[0])));
    columns.emplace_back(new bit_vector_stat(std::move(bvs[1])));

    for (const auto &column : columns) {
        copy.push_back(column->copy());
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed6) {
    size_t num_rows = 4;
    std::vector<sdsl::bit_vector> bvs(8, sdsl::bit_vector(num_rows, 0));
    BitVectorPtrArray columns, copy;

    bvs[0][0] = true;
    bvs[0][1] = true;
    bvs[0][2] = true;
    bvs[0][3] = true;
    bvs[7][3] = true;

    for (size_t j = 0; j < bvs.size(); ++j) {
        columns.emplace_back(new bit_vector_stat(std::move(bvs[j])));
    }

    for (const auto &column : columns) {
        copy.push_back(column->copy());
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy), num_rows), columns);
}

} // namespace
