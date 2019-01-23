#include <random>

#include "gtest/gtest.h"

#include "bin_rel_wt.hpp"
#include "column_major.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column";

TEST(BinRelWT, EmptyConstructorBinRelWT) {
    BinRelWT matrix;
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

TEST(BinRelWT, EmptyBinRelWT) {
    BinRelWT matrix {};
    EXPECT_EQ(0u, matrix.num_columns());
    EXPECT_EQ(0u, matrix.num_rows());
}

void test_bin_rel_wt(const BinRelWT &matrix,
                       const std::vector<std::unique_ptr<bit_vector_sd>> &rows) {
    // check if the number of rows is the same
    ASSERT_EQ(rows.size(), matrix.num_rows());
    // if no rows, check if the number of columnss is zero as well
    if (!rows.size()) {
        ASSERT_EQ(0u, matrix.num_columns());
        return;
    }

    // check if the number of columns is the same
    ASSERT_EQ(rows.at(0)->size(), matrix.num_columns());
    // check get_row
    for (size_t j = 0; j < matrix.num_rows(); ++j) {
        const auto &row = *(rows[j]);
        ASSERT_TRUE(row.size() >= matrix.num_columns());
        auto row_set_bits = matrix.get_row(j);

        // make sure all returned indexes are unique
        ASSERT_EQ(row_set_bits.size(), convert_to_set(row_set_bits).size());

        ASSERT_EQ(row.num_set_bits(), row_set_bits.size())
            << "Row: " << j;

        for (auto i : row_set_bits) {
            ASSERT_TRUE(i < matrix.num_columns()) << "i: " << i << " num_cols: " << matrix.num_columns();
            EXPECT_TRUE(row[i]);
        }
    }

    // check get_column
    for (size_t i = 0; i < matrix.num_columns(); ++i) {
        auto col_set_bits = matrix.get_column(i);

        // make sure all returned indexes are unique
        ASSERT_EQ(col_set_bits.size(), convert_to_set(col_set_bits).size());

        for (auto j : col_set_bits) {
            ASSERT_TRUE(j < matrix.num_rows());
            EXPECT_TRUE((*rows.at(j))[i]);
        }

        auto set_bits = convert_to_set(col_set_bits);
        for (size_t j = 0; j < rows.size(); ++j) {
            EXPECT_EQ((*rows.at(j))[i], set_bits.count(j));
        }
    }

    // check get
    for (size_t i = 0; i < matrix.num_rows(); ++i) {
        for (size_t j = 0; j < matrix.num_columns(); ++j) {
            EXPECT_EQ((*rows.at(i))[j], matrix.get(i, j))
                << i << " " << j;
        }
    }
}

TEST(BinRelWT, AllZero) {
    for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;

            for (size_t j = 0; j < num_rows; ++j) {
               rows.emplace_back(new bit_vector_sd(num_columns, 0));
            }
            size_t num_set_bits = 0;
            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows[r]->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);
            
            test_bin_rel_wt(matrix, rows);
        }
    }
}

TEST(BinRelWT, AllOne) {
    for (size_t num_columns = 2; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 2; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;

            for (size_t j = 0; j < num_rows; ++j) {
                rows.emplace_back(new bit_vector_sd(num_columns, 1));
            }
            size_t num_set_bits = num_rows * num_columns;

            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows.at(r)->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);

            test_bin_rel_wt(matrix, rows);
        }
    }
}

TEST(BinRelWT, AllMixed) {
    for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;
            size_t num_set_bits = 0;

            for (size_t j = 0; j < num_rows; ++j) {
                sdsl::bit_vector bv(num_columns, 0);
                // Let the first and last column be empty.
                for (size_t i = 1; i < num_columns - 1; ++i) {
                    bv[i] = (i + j) % 2;
                }
                rows.emplace_back(new bit_vector_sd(std::move(bv)));
                num_set_bits += rows[rows.size() - 1]->num_set_bits();
            }

            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows[r]->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);

            test_bin_rel_wt(matrix, rows);
       }
    }
}

void test_serialization(const BinRelWT &matrix) {
    {
        std::ofstream out(test_dump_basename_vec_good, std::ios::binary);
        matrix.serialize(out);
        out.close();
    }

    {
        BinRelWT loaded;
        std::ifstream in(test_dump_basename_vec_bad, std::ios::binary);
        ASSERT_FALSE(loaded.load(in));
    }

    BinRelWT loaded;
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

TEST(BinRelWT, SerializationEmpty) {
    BinRelWT matrix;
    test_serialization(matrix);
}

TEST(BinRelWT, SerializationAllZero) {
    for (size_t num_columns = 1; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 1; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;

            for (size_t j = 0; j < num_rows; ++j) {
                rows.emplace_back(new bit_vector_sd(num_columns, 0));
            }

            size_t num_set_bits = 0;
            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows[r]->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);

            test_serialization(matrix);
        }
    }
}

TEST(BinRelWT, SerializationAllOne) {
    for (size_t num_columns = 2; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 2; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;

            for (size_t j = 0; j < num_rows; ++j) {
                rows.emplace_back(new bit_vector_sd(num_columns, 1));
            }
            size_t num_set_bits = num_rows * num_columns;

            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows[r]->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);
 
            test_serialization(matrix);
        }
    }
}

TEST(BinRelWT, SerializationMixed) {
    for (size_t num_columns = 2; num_columns < 20; ++num_columns) {
        for (size_t num_rows = 2; num_rows < 20; ++num_rows) {
            std::vector<std::unique_ptr<bit_vector_sd>> rows;
            size_t num_set_bits = 0;

            for (size_t j = 0; j < num_rows; ++j) {
               sdsl::bit_vector bv(num_columns, 0);
                for (size_t i = 0; i < num_columns; ++i) {
                    bv[i] = (i + j) % 2;
                }
                rows.emplace_back(new bit_vector_sd(std::move(bv)));
                num_set_bits += rows[rows.size() - 1]->num_set_bits();
            }

            auto matrix = BinRelWT([&](auto callback) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    std::vector<BinRelWT::Column> row_set_indices;
                    for (size_t c = 0; c < rows[r]->size(); ++c) {
                        if ((*rows[r])[c])
                            row_set_indices.push_back(c);
                    }
                    callback(row_set_indices);
                }
            }, num_set_bits, num_columns);
            
            test_serialization(matrix);
        }
    }
}
