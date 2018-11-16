#include <random>

#include "gtest/gtest.h"

#include "column_major.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column";


TEST(ColMajorCompressed, EmptyConstructor) {
    ColMajorCompressed matrix;
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(ColMajorCompressed, Empty) {
    ColMajorCompressed matrix {};
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(ColMajorCompressed, OneColumn) {
    std::vector<std::unique_ptr<bit_vector_sd>> columns;
    columns.emplace_back(new bit_vector_sd(10, true));

    ColMajorCompressed matrix(std::move(columns));

    EXPECT_EQ(1u, matrix.num_columns());
    EXPECT_EQ(10u, matrix.num_rows());
}

void test_column_major(const ColMajorCompressed &matrix,
                       const std::vector<std::unique_ptr<bit_vector_sd>> &columns) {
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
            EXPECT_TRUE(col[i]);
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

TEST(ColMajorCompressed, AllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_sd(num_rows, 0));
                copy.emplace_back(new bit_vector_sd(columns.back()->to_vector()));
            }

            ColMajorCompressed matrix_copy(columns);

            test_column_major(matrix_copy, columns);


            ColMajorCompressed matrix(std::move(copy));

            EXPECT_EQ(0u, matrix.num_relations());

            test_column_major(matrix, columns);
        }
    }
}

TEST(ColMajorCompressed, AllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_sd(num_rows, 1));
                copy.emplace_back(new bit_vector_sd(columns.back()->to_vector()));
            }

            ColMajorCompressed matrix_copy(columns);

            test_column_major(matrix_copy, columns);


            ColMajorCompressed matrix(std::move(copy));

            EXPECT_EQ(num_rows * num_columns, matrix.num_relations());

            test_column_major(matrix, columns);
        }
    }
}

TEST(ColMajorCompressed, AllMixed) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns, copy;

            for (size_t j = 0; j < num_columns; ++j) {
                sdsl::bit_vector bv(num_rows, false);
                for (size_t i = 0; i < num_rows; ++i) {
                    bv[i] = (i + 2 * j) % 2;
                }
                columns.emplace_back(new bit_vector_sd(std::move(bv)));
                copy.emplace_back(new bit_vector_sd(columns.back()->to_vector()));
            }

            ColMajorCompressed matrix_copy(columns);

            test_column_major(matrix_copy, columns);


            ColMajorCompressed matrix(std::move(copy));

            test_column_major(matrix, columns);
        }
    }
}

void test_serialization(const ColMajorCompressed &matrix) {
    {
        std::ofstream out(test_dump_basename_vec_good, std::ios::binary);
        matrix.serialize(out);
        out.close();
    }

    {
        ColMajorCompressed loaded;
        std::ifstream in(test_dump_basename_vec_bad, std::ios::binary);
        ASSERT_FALSE(loaded.load(in));
    }

    ColMajorCompressed loaded;
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

TEST(ColMajorCompressed, SerializationEmpty) {
    ColMajorCompressed matrix;
    test_serialization(matrix);
}

TEST(ColMajorCompressed, SerializationAllZero) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_sd(num_rows, 0));
            }

            test_serialization(ColMajorCompressed(std::move(columns)));
        }
    }
}

TEST(ColMajorCompressed, SerializationAllOne) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns;

            for (size_t j = 0; j < num_columns; ++j) {
                columns.emplace_back(new bit_vector_sd(num_rows, 1));
            }

            test_serialization(ColMajorCompressed(std::move(columns)));
        }
    }
}

TEST(ColMajorCompressed, SerializationMixed) {
    for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
        for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
            std::vector<std::unique_ptr<bit_vector_sd>> columns;

            for (size_t j = 0; j < num_columns; ++j) {
                sdsl::bit_vector bv(num_rows, false);
                for (size_t i = 0; i < num_rows; ++i) {
                    bv[i] = (i + 2 * j) % 2;
                }
                columns.emplace_back(new bit_vector_sd(std::move(bv)));
            }

            test_serialization(ColMajorCompressed(std::move(columns)));
        }
    }
}
