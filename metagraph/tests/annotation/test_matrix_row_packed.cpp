#include "gtest/gtest.h"

#include "rainbowfish.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "utils.hpp"
#include "bit_vector.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

const std::vector<sdsl::bit_vector> vectors {
    { 0, 1, 1, 0, 0, 1, 0 },
    { 0, 1, 1, 1, 1, 1, 0 }
};

const std::vector<sdsl::bit_vector> matrix {
    { 0, 0 },
    { 1, 1 },
    { 1, 1 },
    { 0, 1 },
    { 0, 1 },
    { 1, 1 },
    { 0, 0 }
};

const std::vector<std::vector<uint64_t>> matrix_rows {
    { },
    { 0, 1 },
    { 0, 1 },
    { 1 },
    { 1 },
    { 0, 1 },
    { }
};

const std::vector<std::vector<uint64_t>> matrix_columns {
    { 1, 2, 5 },
    { 1, 2, 3, 4, 5 }
};

const std::vector<std::vector<uint64_t>> indices {
    { 1, 0 }, { 1, 1 },
    { 2, 0 }, { 2, 1 },
    { 3, 1 },
    { 4, 1 },
    { 5, 0 }, { 5, 1 }
};

using Row = BinaryMatrix::Row;
using Column = BinaryMatrix::Column;


template <class BinMat>
void test_serialize_load_test(const BinMat &binary_matrix,
                              uint64_t num_rows = matrix.size(),
                              uint64_t num_columns = matrix[0].size(),
                              uint64_t num_relations = indices.size()) {
    test_binary_matrix(binary_matrix, num_rows, num_columns, num_relations);

    std::ofstream out(test_dump_basename, std::ios::binary);
    binary_matrix.serialize(out);
    out.close();

    BinMat binary_matrix_loaded;
    std::ifstream in(test_dump_basename, std::ios::binary);
    binary_matrix_loaded.load(in);

    test_binary_matrix(binary_matrix_loaded, num_rows, num_columns, num_relations);
}

void test_binary_matrix(const BinaryMatrix &binary_matrix,
                        uint64_t num_rows,
                        uint64_t num_columns,
                        uint64_t num_relations) {
    ASSERT_EQ(num_rows, binary_matrix.num_rows());
    ASSERT_EQ(num_columns, binary_matrix.num_columns());
    ASSERT_EQ(num_relations, binary_matrix.num_relations());

    for (Row i = 0; i < binary_matrix.num_rows(); ++i) {
        for (Column j = 0; j < binary_matrix.num_columns(); ++j) {
            EXPECT_EQ(matrix[i][j], binary_matrix.get(i, j)) << i << " " << j;
        }
        EXPECT_EQ(matrix_rows[i], binary_matrix.get_row(i)) << i;
    }

    for (Column j = 0; j < binary_matrix.num_columns(); ++j) {
        std::vector<Column> col_indices;
        for (const auto &a : matrix_columns[j]) {
            if (a >= num_rows)
                break;

            col_indices.push_back(a);
        }
        EXPECT_EQ(col_indices, binary_matrix.get_column(j)) << j;
    }
}

TEST(RowPacked, RowConcatenated) {
    std::vector<std::unique_ptr<bit_vector_small>> small_vectors;
    small_vectors.emplace_back(new bit_vector_small(vectors[0]));
    small_vectors.emplace_back(new bit_vector_small(vectors[1]));
    test_serialize_load_test(RowConcatenated<>(
        [&](auto vector_callback) {
            utils::call_rows(vector_callback, small_vectors);
        }, vectors.size(), matrix.size(), indices.size()
    ));
}

TEST(RowPacked, RainbowfishNoBuffer) {
    std::vector<std::unique_ptr<bit_vector_small>> small_vectors;
    small_vectors.emplace_back(new bit_vector_small(vectors[0]));
    small_vectors.emplace_back(new bit_vector_small(vectors[1]));
    test_serialize_load_test(Rainbowfish(
        [&](auto vector_callback) {
            utils::call_rows(vector_callback, small_vectors);
        }, vectors.size()
    ));
}

TEST(RowPacked, TestConstructor) {
    std::unique_ptr<Rainbowfish> rbf;

    rbf.reset(new Rainbowfish(
        [&](auto row_callback) {
            row_callback(std::vector<Rainbowfish::Column> { 0, 1, 2, 3 });
            row_callback(std::vector<Rainbowfish::Column> { 0, 1, 2, 3 });
        },
        4
    ));

    rbf.reset(new Rainbowfish(
        [&](auto row_callback) {
            row_callback(std::vector<Rainbowfish::Column> { 0, 1, 2, 3 });
            row_callback(std::vector<Rainbowfish::Column> {  });
        },
        4
    ));

    rbf.reset(new Rainbowfish(
        [&](auto row_callback) {
            row_callback(std::vector<Rainbowfish::Column> { 0, 1, 2, 3 });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> {  });
            row_callback(std::vector<Rainbowfish::Column> { 3 });
        },
        4
    ));

    rbf.reset(new Rainbowfish(
        [&](auto row_callback) {
            for (size_t i = 0; i < 1000; ++i) {
                row_callback(std::vector<Rainbowfish::Column> { 0, 1, 2, 3 });
            }
        },
        4
    ));

    rbf.reset(new Rainbowfish(
        [&](auto row_callback) {
            for (size_t i = 0; i < 1000; ++i) {
                for (size_t j = 0; j < 100; ++j) {
                    row_callback(std::vector<Rainbowfish::Column> { j });
                }
            }
        },
        1000
    ));

    ASSERT_TRUE(true);
}

TEST(RowPacked, RainbowfishBufferAllSizes) {
    std::vector<std::vector<bool>> columns(vectors.size());
    uint64_t num_relations = 0;
    for (uint64_t i = 0; i < 8; ++i) {
        std::vector<bit_vector_small> columns_small;
        for (auto &a : columns) {
            columns_small.emplace_back(to_sdsl(a));
        }

        std::vector<const bit_vector_small*> columns_small_ptr;
        for (auto &a : columns_small) {
            columns_small_ptr.emplace_back(&a);
        }

        for (uint64_t j = 1; j < i; ++j) {
            test_serialize_load_test(Rainbowfish(
                [&](auto vector_callback) {
                    utils::call_rows(vector_callback, columns_small_ptr);
                }, 2, j),
            i, columns.size(), num_relations);
        }
        if (i < matrix.size()) {
            for (uint64_t j = 0; j < columns.size(); ++j) {
                columns[j].push_back(matrix[i][j]);
                if (matrix[i][j])
                    num_relations++;
            }
        }
    }
}

const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_row_flat";

//TEST(RowFlat, MergeLoadDisjoint) {
//    {
//        annotate::RowFlatAnnotator annotation;
//        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
//        annotation.set_labels(2, { "Label1", "Label2" });
//        annotation.set_labels(4, { "Label8" });
//
//        annotation.serialize(test_dump_basename_vec_good + "_1");
//    }
//    {
//        annotate::RowFlatAnnotator annotation;
//        annotation.set_labels(1, { "2_Label0", "2_Label2", "2_Label8" });
//        annotation.set_labels(2, { "2_Label1", "2_Label9", "2_Label0" });
//        annotation.set_labels(3, { "2_Label8" });
//
//        annotation.serialize(test_dump_basename_vec_good + "_2");
//    }
//    {
//        annotate::RowFlatAnnotator annotation;
//        ASSERT_TRUE(annotation.merge_load({ test_dump_basename_vec_good + "_1",
//                                            test_dump_basename_vec_good + "_2" }));
//
//        EXPECT_EQ(5u, annotation.num_objects());
//        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
//        EXPECT_EQ(convert_to_set({ "2_Label0", "2_Label2", "2_Label8" }), convert_to_set(annotation.get(1)));
//        EXPECT_EQ(convert_to_set({ "Label1", "2_Label1", "Label2", "2_Label9", "2_Label0" }),
//                    convert_to_set(annotation.get(2)));
//        EXPECT_EQ(convert_to_set({ "2_Label8" }), convert_to_set(annotation.get(3)));
//        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
//    }
//}

TEST(RowFlat, MergeLoad) {
    {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(4, {"Label8"});

        auto annotation = annotate::convert<annotate::RowFlatAnnotator>(std::move(column_annotator));
        annotation->serialize(test_dump_basename_vec_good + "_1");
    }
    {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.set_labels(1, { "Label0", "Label2", "Label8" });
        column_annotator.set_labels(2, { "Label1", "Label9", "Label0" });
        column_annotator.set_labels(3, { "Label8" });

        auto annotation = annotate::convert<annotate::RowFlatAnnotator>(std::move(column_annotator));
        annotation->serialize(test_dump_basename_vec_good + "_2");
    }
    {
        annotate::RowFlatAnnotator annotation;
        ASSERT_TRUE(annotation.merge_load({ test_dump_basename_vec_good + "_1",
                                            test_dump_basename_vec_good + "_2" }));

        EXPECT_EQ(5u, annotation.num_objects());
        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label9", "Label0" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}
