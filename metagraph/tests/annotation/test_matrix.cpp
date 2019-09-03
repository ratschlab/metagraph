#include <random>

#include "gtest/gtest.h"

#include "test_matrix_helpers.hpp"

#include "bin_rel_wt.hpp"
#include "bin_rel_wt_sdsl.hpp"
#include "column_major.hpp"


template <typename BinMat>
class BinaryMatrixTest : public ::testing::Test { };
typedef ::testing::Types<BRWT,
                         BRWTOptimized,
                         ColMajorCompressed,
                         BinRelWT,
                         BinRelWT_sdsl,
                         RowConcatenated<>,
                         Rainbowfish> BinMatTypes;
TYPED_TEST_CASE(BinaryMatrixTest, BinMatTypes);


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
    {
        auto matrix = build_matrix_from_rows<TypeParam>();
        EXPECT_EQ(0u, matrix.num_columns());
        EXPECT_EQ(0u, matrix.num_rows());
    }
}

TYPED_TEST(BinaryMatrixTest, BuildOneColumn) {
    BitVectorPtrArray columns, copy1, copy2;
    columns.emplace_back(new bit_vector_stat(10, true));
    copy1.emplace_back(new bit_vector_stat(10, true));
    copy2.emplace_back(new bit_vector_stat(10, true));

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), 10), columns);
    test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), 10), columns);
}

TYPED_TEST(BinaryMatrixTest, AllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy1, copy2;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
                copy1.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
                copy2.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
            test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy1, copy2;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
                copy1.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
                copy2.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
            test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed1) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy1, copy2;

            for (size_t j = 0; j < num_columns; ++j) {

                columns.emplace_back(new bit_vector_stat(num_rows));

                for (size_t i = 0; i < num_rows; ++i) {
                    columns.back()->set(i, (i + 2 * j) % 2);
                }
                copy1.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
                copy2.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
            test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed2) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BitVectorPtrArray columns, copy1, copy2;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows));
            }

            for (size_t j = 0; j < num_rows; ++j) {
                for (size_t i = 1; i < num_columns - 1; ++i) {
                    columns.at(i)->set(j, (i + j) % 2);
                }
            }

            for (const auto &column : columns) {
                copy1.emplace_back(new bit_vector_stat(column->to_vector()));
                copy2.emplace_back(new bit_vector_stat(column->to_vector()));
            }

            test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
            test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
        }
    }
}

TYPED_TEST(BinaryMatrixTest, AllMixed3) {
    size_t num_rows = 7;
    size_t num_columns = 2;
    BitVectorPtrArray columns, copy1, copy2;

    for (size_t j = 0; j < num_columns; ++j) {
        columns.emplace_back(new bit_vector_stat(num_rows));
    }

    columns.at(0)->set(1, true);
    columns.at(0)->set(2, true);
    columns.at(0)->set(5, true);
    columns.at(1)->set(1, true);
    columns.at(1)->set(2, true);
    columns.at(1)->set(3, true);
    columns.at(1)->set(4, true);
    columns.at(1)->set(5, true);

    for (const auto &column : columns) {
        copy1.emplace_back(new bit_vector_stat(column->to_vector()));
        copy2.emplace_back(new bit_vector_stat(column->to_vector()));
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
    test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed4) {
    size_t num_rows = 4;
    size_t num_columns = 2;
    BitVectorPtrArray columns, copy1, copy2;

    for (size_t j = 0; j < num_columns; ++j) {
        columns.emplace_back(new bit_vector_stat(num_rows));
    }

    columns.at(0)->set(0, true);
    columns.at(0)->set(1, true);
    columns.at(0)->set(2, true);
    columns.at(0)->set(3, true);
    columns.at(1)->set(0, true);
    columns.at(1)->set(1, true);
    columns.at(1)->set(2, true);
    columns.at(1)->set(3, true);

    for (const auto &column : columns) {
        copy1.emplace_back(new bit_vector_stat(column->to_vector()));
        copy2.emplace_back(new bit_vector_stat(column->to_vector()));
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows),columns);
    test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed5) {
    size_t num_rows = 4;
    size_t num_columns = 2;
    BitVectorPtrArray columns, copy1, copy2;

    for (size_t j = 0; j < num_columns; ++j) {
        columns.emplace_back(new bit_vector_stat(num_rows));
    }

    columns.at(0)->set(0, true);
    columns.at(0)->set(1, true);
    columns.at(0)->set(2, true);
    columns.at(0)->set(3, true);

    for (const auto &column : columns) {
        copy1.emplace_back(new bit_vector_stat(column->to_vector()));
        copy2.emplace_back(new bit_vector_stat(column->to_vector()));
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
    test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
}

TYPED_TEST(BinaryMatrixTest, AllMixed6) {
    size_t num_rows = 4;
    size_t num_columns = 8;
    BitVectorPtrArray columns, copy1, copy2;

    for (size_t j = 0; j < num_columns; ++j) {
        columns.emplace_back(new bit_vector_stat(num_rows));
    }

    columns.at(0)->set(0, true);
    columns.at(0)->set(1, true);
    columns.at(0)->set(2, true);
    columns.at(0)->set(3, true);
    columns.at(7)->set(3, true);

    for (const auto &column : columns) {
        copy1.emplace_back(new bit_vector_stat(column->to_vector()));
        copy2.emplace_back(new bit_vector_stat(column->to_vector()));
    }

    test_matrix(build_matrix_from_columns<TypeParam>(std::move(copy1), num_rows), columns);
    test_matrix(build_matrix_from_rows<TypeParam>(std::move(copy2), num_rows), columns);
}
