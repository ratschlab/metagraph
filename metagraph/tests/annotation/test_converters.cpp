#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/annotation_converters.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_row_compressed_merge = test_dump_basename + "_row_compressed_merge";
const std::string test_dump_basename_rowflat_merge = test_dump_basename + "_rowflat_merge";
const std::string test_dump_basename_row_compressed_to_rowflat = test_dump_basename + "_row_compressed_to_rowflat";


class ConvertFromRowCompressed : public ::testing::Test {
  protected:
    static const uint64_t num_rows = 5;
    static anno::RowCompressed<> *initial_annotation;
    static anno::MultiLabelEncoded<std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new anno::RowCompressed<>(num_rows);
        initial_annotation->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        initial_annotation->add_labels({ 2 }, {"Label1", "Label2"});
        initial_annotation->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        initial_annotation->add_labels({ 4 }, {"Label2"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
                  convert_to_set(annotation->get(0)));
        EXPECT_EQ(std::vector<std::string>({}),
                  annotation->get(1));
        EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
                  convert_to_set(annotation->get(2)));
        EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
                  convert_to_set(annotation->get(3)));
        EXPECT_EQ(convert_to_set({"Label2"}),
                  convert_to_set(annotation->get(4)));

        delete initial_annotation;
        delete annotation;
    }
};

const uint64_t ConvertFromRowCompressed::num_rows;
anno::RowCompressed<> *ConvertFromRowCompressed::initial_annotation = nullptr;
anno::MultiLabelEncoded<std::string> *ConvertFromRowCompressed::annotation = nullptr;

class MergeAnnotators : public ::testing::Test {
  protected:
    static const uint64_t num_rows = 5;
    static anno::RowCompressed<> *input_annotation_1;
    static anno::RowCompressed<> *input_annotation_2;
    static anno::RowCompressed<> *merged_annotation_expected;
    static anno::MultiLabelEncoded<std::string> *merged_annotation;

    virtual void SetUp() {

        input_annotation_1 = new anno::RowCompressed<>(num_rows);
        input_annotation_1->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        input_annotation_1->add_labels({ 2 }, {"Label1", "Label2"});
        input_annotation_1->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        input_annotation_1->add_labels({ 4 }, {"Label2"});

        input_annotation_2 = new anno::RowCompressed<>(num_rows);
        input_annotation_2->add_labels({ 1 }, {"Label0", "Label3"});
        input_annotation_2->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        input_annotation_2->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});

        merged_annotation_expected = new anno::RowCompressed<>(num_rows);
        merged_annotation_expected->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        merged_annotation_expected->add_labels({ 1 }, {"Label0", "Label3"});
        merged_annotation_expected->add_labels({ 2 }, {"Label0", "Label2", "Label1", "Label9", "Label7"});
        merged_annotation_expected->add_labels({ 3 }, {"Label2", "Label8", "Label1"});
        merged_annotation_expected->add_labels({ 4 }, {"Label2", "Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(merged_annotation);
        for(uint64_t i = 0; i < num_rows; ++i) {
            auto row_expected = merged_annotation_expected->get(i);
            auto row = merged_annotation->get(i);
            EXPECT_EQ(row_expected, row);
        }

        delete input_annotation_1;
        delete input_annotation_2;
        delete merged_annotation_expected;
        delete merged_annotation;
    }
};

const uint64_t MergeAnnotators::num_rows;
anno::RowCompressed<> *MergeAnnotators::input_annotation_1 = nullptr;
anno::RowCompressed<> *MergeAnnotators::input_annotation_2 = nullptr;
anno::RowCompressed<> *MergeAnnotators::merged_annotation_expected = nullptr;
anno::MultiLabelEncoded<std::string> *MergeAnnotators::merged_annotation = nullptr;


class ConvertFromColumnCompressed : public ::testing::Test {
  protected:
    static anno::ColumnCompressed<> *initial_annotation;
    static anno::MultiLabelEncoded<std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new anno::ColumnCompressed<>(5);
        initial_annotation->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        initial_annotation->add_labels({ 2 }, {"Label1", "Label2"});
        initial_annotation->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        initial_annotation->add_labels({ 4 }, {"Label2"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
                  convert_to_set(annotation->get(0)));
        EXPECT_EQ(std::vector<std::string>({}),
                  annotation->get(1));
        EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
                  convert_to_set(annotation->get(2)));
        EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
                  convert_to_set(annotation->get(3)));
        EXPECT_EQ(convert_to_set({"Label2"}),
                  convert_to_set(annotation->get(4)));

        delete initial_annotation;
        delete annotation;
    }

};

anno::ColumnCompressed<> *ConvertFromColumnCompressed::initial_annotation = nullptr;
anno::MultiLabelEncoded<std::string> *ConvertFromColumnCompressed::annotation = nullptr;


// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT) {
//     anno::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::BinRelWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT) {
    annotation = anno::convert<anno::BinRelWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT_sdsl) {
//     anno::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::BinRelWT_sdslAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT_sdsl) {
    annotation = anno::convert<anno::BinRelWT_sdslAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_RowFlat) {
//     anno::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RowFlat) {
    annotation = anno::convert<anno::RowFlatAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_Rainbowfish) {
//     anno::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RainbowfishAnnotator) {
    annotation = anno::convert<anno::RainbowfishAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_GreedyBRWT) {
//     anno::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert_to_greedy_BRWT<anno::MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_GreedyBRWT) {
    annotation = anno::convert_to_greedy_BRWT<anno::MultiBRWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}


TEST(ConvertFromRowCompressedEmpty, to_BinRelWT) {
    anno::RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = anno::convert<anno::BinRelWTAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT) {
    annotation = anno::convert<anno::BinRelWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

TEST(ConvertFromRowCompressedEmpty, to_BinRelWT_sdsl) {
    anno::RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = anno::convert<anno::BinRelWT_sdslAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT_sdsl) {
    annotation = anno::convert<anno::BinRelWT_sdslAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_RowFlat) {
//     anno::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RowFlat) {
    annotation = anno::convert<anno::RowFlatAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat) {
    initial_annotation->serialize(test_dump_basename_row_compressed_to_rowflat);

    annotation = anno::convert<anno::RowFlatAnnotator>(
        test_dump_basename_row_compressed_to_rowflat
    ).release();
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat2) {
    const auto rowflat_filename = test_dump_basename_row_compressed_to_rowflat
                                                + anno::RowCompressed<>::kExtension;

    initial_annotation->serialize(rowflat_filename);

    annotation = anno::convert<anno::RowFlatAnnotator>(
        rowflat_filename
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_Rainbowfish) {
//     anno::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert<anno::RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RainbowfishAnnotator) {
    annotation = anno::convert<anno::RainbowfishAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_GreedyBRWT) {
//     anno::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = anno::convert_to_greedy_BRWT<anno::MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

// TEST_F(ConvertFromRowCompressed, to_GreedyBRWT) {
//     annotation = anno::convert_to_greedy_BRWT<anno::MultiBRWTAnnotator>(
//         std::move(*initial_annotation)
//     ).release();
// }

TEST_F(MergeAnnotators, RowCompressed) {
    std::vector<std::string> filenames;
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_1";
        input_annotation_1->serialize(filename);
        filenames.push_back(filename + anno::RowCompressed<>::kExtension);
    }
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_2";
        input_annotation_2->serialize(filename);
        filenames.push_back(filename + anno::RowCompressed<>::kExtension);
    }

    anno::merge<anno::RowCompressed<>, std::string>(
        {}, filenames, test_dump_basename_row_compressed_merge + "_merged"
    );

    merged_annotation = new anno::RowCompressed<>(num_rows);
    merged_annotation->merge_load({ test_dump_basename_row_compressed_merge + "_merged" });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowCompressed) {
    std::vector<std::unique_ptr<anno::MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(anno::convert<anno::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(anno::convert<anno::RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowcompressed";
    anno::merge<anno::RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new anno::RowCompressed<>(num_rows);
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowFlat) {
    std::vector<std::unique_ptr<anno::MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(anno::convert<anno::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(anno::convert<anno::RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowflat";
    anno::merge<anno::RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new anno::RowFlatAnnotator();
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, Mixed_to_RowFlat) {
    std::vector<std::unique_ptr<anno::MultiLabelEncoded<std::string>>> annotators;
    std::vector<std::string> filenames;
    {
        auto annotator = anno::convert<anno::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        );
        annotators.push_back(std::move(annotator));
    }
    {
        auto annotator = std::make_unique<anno::ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    //TODO
    {
        auto annotator = std::make_unique<anno::ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    {
        //TODO: move into fixture as input_annotation_3 and make non-overlapping
        const std::string filename = test_dump_basename_row_compressed_merge + "_mixed_2";
        auto annotation = std::make_unique<anno::RowCompressed<> >(num_rows);
        annotation->add_labels({ 0 }, {"Label0"});
        annotation->add_labels({ 2 }, {"Label1"});
        annotation->add_labels({ 3 }, {"Label1"});
        annotation->add_labels({ 4 }, {"Label2"});
        annotation->serialize(filename);
        filenames.push_back(filename + anno::RowCompressed<>::kExtension);
    }

    const auto outfile = test_dump_basename_rowflat_merge + "_mixed_to_rowflat";
    anno::merge<anno::RowFlatAnnotator, std::string>(std::move(annotators), filenames, outfile);

    merged_annotation = new anno::RowFlatAnnotator();
    merged_annotation->merge_load({ outfile });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST(ColumnCompressed, ToRowAnnotator) {
    {
        anno::ColumnCompressed<> annotation(0);

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);
    }
    {
        anno::ColumnCompressed<> annotation(1);

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        anno::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorParallel) {
    {
        anno::ColumnCompressed<> annotation(0);

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);
    }
    {
        anno::ColumnCompressed<> annotation(1);

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        anno::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        anno::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreaming) {
    {
        anno::ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        anno::ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        anno::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreamingParallel) {
    {
        anno::ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        anno::ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        anno::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        anno::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename, 10);

        anno::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

} // namespace
