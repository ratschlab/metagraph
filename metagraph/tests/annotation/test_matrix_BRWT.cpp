#include <random>

#include "gtest/gtest.h"

#include "BRWT.hpp"
#include "BRWT_builders.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_BRWT";


TEST(BRWT, EmptyConstructor) {
    BRWT matrix;
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(BRWT, BuildBottomUPEmpty) {
    BRWT matrix = BRWTBottomUpBuilder::build({});
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(BRWT, BuildBottomUPOneColumn) {
    BRWTBottomUpBuilder::VectorsPtr columns;
    columns.emplace_back(new bit_vector_stat(10, true));

    BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));

    EXPECT_EQ(1u, matrix.num_columns());
    EXPECT_EQ(10u, matrix.num_rows());
}

TEST(BRWT, ArityEmpty) {
    BRWT matrix;
    // empty root
    EXPECT_EQ(1u, matrix.num_nodes());
    ASSERT_EQ(0u, matrix.avg_arity());
}

TEST(BRWT, ArityOneCol) {
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, true));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, false));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));

        // only root
        EXPECT_EQ(1u, matrix.num_nodes());
        EXPECT_EQ(0u, matrix.avg_arity());
    }
}

TEST(BRWT, ArityTwoCol) {
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, true));
        columns.emplace_back(new bit_vector_stat(10, false));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
    {
        BRWTBottomUpBuilder::VectorsPtr columns;
        columns.emplace_back(new bit_vector_stat(10, false));
        columns.emplace_back(new bit_vector_stat(10, true));

        BRWT matrix = BRWTBottomUpBuilder::build(std::move(columns));

        // root + 2 leaves
        EXPECT_EQ(3u, matrix.num_nodes());
        EXPECT_EQ(2u, matrix.avg_arity());
    }
}

void test_brwt(const BRWT &matrix,
               const BRWTBottomUpBuilder::VectorsPtr &columns) {
    // check if the number of columns is the same
    ASSERT_EQ(columns.size(), matrix.num_columns());
    // if no columns, check if the number of rows is zero as well
    if (!columns.size()) {
        ASSERT_EQ(0u, matrix.num_rows());
        return;
    }

    // check if the number of rows is the same
    ASSERT_EQ(columns.at(0)->size(), matrix.num_rows());

    // check get_column
    for (size_t j = 0; j < matrix.num_columns(); ++j) {
        const auto &col = *columns[j];

        assert(col.size() == matrix.num_rows());

        auto col_set_bits = matrix.get_column(j);

        // make sure all returned indexes are unique
        ASSERT_EQ(col_set_bits.size(), convert_to_set(col_set_bits).size());

        ASSERT_EQ(col.num_set_bits(), col_set_bits.size())
            << "Column: " << j;

        for (auto i : col_set_bits) {
            ASSERT_TRUE(i < matrix.num_rows());
            ASSERT_TRUE(col[i]);
        }
    }

    // check get_row
    for (size_t i = 0; i < matrix.num_rows(); ++i) {
        auto row_set_bits = matrix.get_row(i);

        // make sure all returned indexes are unique
        ASSERT_EQ(row_set_bits.size(), convert_to_set(row_set_bits).size());

        for (auto j : row_set_bits) {
            ASSERT_TRUE(j < matrix.num_columns());
            EXPECT_TRUE((*columns[j])[i]);
        }

        auto set_bits = convert_to_set(row_set_bits);
        for (size_t j = 0; j < columns.size(); ++j) {
            EXPECT_EQ((*columns[j])[i], set_bits.count(j));
        }
    }

    // check get
    for (size_t i = 0; i < matrix.num_rows(); ++i) {
        for (size_t j = 0; j < matrix.num_columns(); ++j) {
            EXPECT_EQ(columns[j]->operator[](i), matrix.get(i, j))
                << i << " " << j;
        }
    }
}

TEST(BRWT, BuildBottomUPAllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            BRWT matrix = BRWTBottomUpBuilder::build(std::move(copy));

            EXPECT_EQ(0u, matrix.num_relations());

            test_brwt(matrix, columns);
        }
    }
}

TEST(BRWT, BuildBottomUPAllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
                copy.emplace_back(new bit_vector_stat(columns.back()->to_vector()));
            }

            BRWT matrix = BRWTBottomUpBuilder::build(std::move(copy));

            EXPECT_EQ(num_rows * num_columns, matrix.num_relations());

            test_brwt(matrix, columns);
        }
    }
}

TEST(BRWT, BuildBottomUPAllMixed) {
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

            EXPECT_TRUE(matrix.avg_arity() <= 2) << matrix.avg_arity();

            test_brwt(matrix, columns);
        }
    }
}

void test_serialization(const BRWT &matrix) {
    {
        std::ofstream out(test_dump_basename_vec_good, std::ios::binary);
        matrix.serialize(out);
        out.close();
    }

    {
        BRWT loaded;
        std::ifstream in(test_dump_basename_vec_bad, std::ios::binary);
        ASSERT_FALSE(loaded.load(in));
    }

    BRWT loaded;
    {
        std::ifstream in(test_dump_basename_vec_good, std::ios::binary);
        ASSERT_TRUE(loaded.load(in));
        ASSERT_TRUE(in.good());
    }

    ASSERT_EQ(matrix.num_columns(), loaded.num_columns());
    ASSERT_EQ(matrix.num_rows(), loaded.num_rows());
    for (size_t j = 0; j < loaded.num_columns(); ++j) {
        EXPECT_EQ(matrix.get_column(j), loaded.get_column(j));
    }
}

TEST(BRWT, SerializationEmpty) {
    BRWT matrix;
    test_serialization(matrix);
}

TEST(BRWT, SerializationAllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 0));
            }

            test_serialization(BRWTBottomUpBuilder::build(std::move(columns)));
        }
    }
}

TEST(BRWT, SerializationAllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_stat(num_rows, 1));
            }

            test_serialization(BRWTBottomUpBuilder::build(std::move(columns)));
        }
    }
}

TEST(BRWT, SerializationMixed) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            BRWTBottomUpBuilder::VectorsPtr columns;

            for (size_t j = 0; j < num_columns; ++j) {

                columns.emplace_back(new bit_vector_stat(num_rows));

                for (size_t i = 0; i < num_rows; ++i) {
                    columns.back()->set(i, (i + 2 * j) % 2);
                }
            }

            test_serialization(BRWTBottomUpBuilder::build(std::move(columns)));
        }
    }
}
