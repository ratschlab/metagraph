#include <random>

#include "gtest/gtest.h"

#include "BRWT.hpp"
#include "BRWT_builders.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_BRWT";


TEST(BRWTOptimizer, EmptyConstructor) {
    BRWT matrix;
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(BRWTOptimizer, BuildBottomUPEmpty) {
    BRWT matrix = BRWTBottomUpBuilder::build({});
    BRWTOptimizer::relax(&matrix);

    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(BRWTOptimizer, BuildBottomUPOneColumn) {
    BRWTBottomUpBuilder::VectorsPtr columns;
    columns.emplace_back(new bit_vector_stat(10, true));

    BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));
    BRWTOptimizer::relax(&matrix);

    EXPECT_EQ(1u, matrix.num_columns());
    EXPECT_EQ(10u, matrix.num_rows());
}

TEST(BRWTOptimizer, ArityEmpty) {
    BRWT matrix;
    // empty root
    EXPECT_EQ(1u, matrix.num_nodes());
    ASSERT_EQ(0u, matrix.avg_arity());
}

TEST(BRWTOptimizer, ArityOneCol) {
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, true));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));
        BRWTOptimizer::relax(&matrix);

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, false));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));
        BRWTOptimizer::relax(&matrix);

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
}

TEST(BRWTOptimizer, ArityTwoCol) {
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, true));
        columns.emplace_back(new bit_vector_stat(10, false));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));
        BRWTOptimizer::relax(&matrix);

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, false));
        columns.emplace_back(new bit_vector_stat(10, true));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));
        BRWTOptimizer::relax(&matrix);

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
}

void test_brwt(const BRWT &matrix,
               const BRWTBottomUpBuilder::VectorsPtr &columns);

TEST(BRWTOptimizer, BuildBottomUPAllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            BRWT matrix = BRWTBottomUpBuilder::build(std::move(copy));
            BRWTOptimizer::relax(&matrix);

            EXPECT_EQ(0u, matrix.num_relations());

            test_brwt(matrix, columns);
        }
    }
}

TEST(BRWTOptimizer, BuildBottomUPAllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            BRWT matrix = BRWTBottomUpBuilder::build(std::move(copy));
            BRWTOptimizer::relax(&matrix);

            EXPECT_EQ(num_rows * num_columns, matrix.num_relations());

            test_brwt(matrix, columns);
        }
    }
}

TEST(BRWTOptimizer, BuildBottomUPAllMixed) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {

                columns.emplace_back(new bit_vector_stat(num_rows));

                for (size_t i = 0; i < num_rows; ++i) {
                    columns.back()->set(i, (i + 2 * j) % 2);
                }
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            BRWT matrix = BRWTBottomUpBuilder::build(std::move(copy));
            BRWTOptimizer::relax(&matrix);

            test_brwt(matrix, columns);
        }
    }
}

void test_serialization(const BRWT &matrix);

TEST(BRWTOptimizer, SerializationEmpty) {
    BRWT matrix;
    test_serialization(matrix);
}

TEST(BRWTOptimizer, SerializationAllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
            }

            auto matrix = BRWTBottomUpBuilder::build(std::move(columns));
            BRWTOptimizer::relax(&matrix);
            test_serialization(matrix);
        }
    }
}

TEST(BRWTOptimizer, SerializationAllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
            }

            auto matrix = BRWTBottomUpBuilder::build(std::move(columns));
            BRWTOptimizer::relax(&matrix);
            test_serialization(matrix);
        }
    }
}

TEST(BRWTOptimizer, SerializationMixed) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {

                columns.emplace_back(new bit_vector_stat(num_rows));

                for (size_t i = 0; i < num_rows; ++i) {
                    columns.back()->set(i, (i + 2 * j) % 2);
                }
            }

            auto matrix = BRWTBottomUpBuilder::build(std::move(columns));
            BRWTOptimizer::relax(&matrix);
            test_serialization(matrix);
        }
    }
}
