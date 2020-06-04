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
    static annotate::RowCompressed<> *initial_annotation;
    static annotate::MultiLabelEncoded<std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new annotate::RowCompressed<>(num_rows);
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
annotate::RowCompressed<> *ConvertFromRowCompressed::initial_annotation = nullptr;
annotate::MultiLabelEncoded<std::string> *ConvertFromRowCompressed::annotation = nullptr;

class MergeAnnotators : public ::testing::Test {
  protected:
    static const uint64_t num_rows = 5;
    static annotate::RowCompressed<> *input_annotation_1;
    static annotate::RowCompressed<> *input_annotation_2;
    static annotate::RowCompressed<> *merged_annotation_expected;
    static annotate::MultiLabelEncoded<std::string> *merged_annotation;

    virtual void SetUp() {

        input_annotation_1 = new annotate::RowCompressed<>(num_rows);
        input_annotation_1->add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        input_annotation_1->add_labels({ 2 }, {"Label1", "Label2"});
        input_annotation_1->add_labels({ 3 }, {"Label1", "Label2", "Label8"});
        input_annotation_1->add_labels({ 4 }, {"Label2"});

        input_annotation_2 = new annotate::RowCompressed<>(num_rows);
        input_annotation_2->add_labels({ 1 }, {"Label0", "Label3"});
        input_annotation_2->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        input_annotation_2->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});

        merged_annotation_expected = new annotate::RowCompressed<>(num_rows);
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
annotate::RowCompressed<> *MergeAnnotators::input_annotation_1 = nullptr;
annotate::RowCompressed<> *MergeAnnotators::input_annotation_2 = nullptr;
annotate::RowCompressed<> *MergeAnnotators::merged_annotation_expected = nullptr;
annotate::MultiLabelEncoded<std::string> *MergeAnnotators::merged_annotation = nullptr;


class ConvertFromColumnCompressed : public ::testing::Test {
  protected:
    static annotate::ColumnCompressed<> *initial_annotation;
    static annotate::MultiLabelEncoded<std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new annotate::ColumnCompressed<>(5);
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

annotate::ColumnCompressed<> *ConvertFromColumnCompressed::initial_annotation = nullptr;
annotate::MultiLabelEncoded<std::string> *ConvertFromColumnCompressed::annotation = nullptr;


// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::BinRelWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT) {
    annotation = annotate::convert<annotate::BinRelWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_BinRelWT_sdsl) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_BinRelWT_sdsl) {
    annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_RowFlat) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RowFlat) {
    annotation = annotate::convert<annotate::RowFlatAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_Rainbowfish) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_RainbowfishAnnotator) {
    annotation = annotate::convert<annotate::RainbowfishAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromColumnCompressedEmpty, to_GreedyBRWT) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_GreedyBRWT) {
    annotation = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}


TEST(ConvertFromRowCompressedEmpty, to_BinRelWT) {
    annotate::RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = annotate::convert<annotate::BinRelWTAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT) {
    annotation = annotate::convert<annotate::BinRelWTAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

TEST(ConvertFromRowCompressedEmpty, to_BinRelWT_sdsl) {
    annotate::RowCompressed<> empty_column_annotator(5);
    auto empty_annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(5u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

TEST_F(ConvertFromRowCompressed, to_BinRelWT_sdsl) {
    annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_RowFlat) {
//     annotate::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::RowFlatAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RowFlat) {
    annotation = annotate::convert<annotate::RowFlatAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat) {
    initial_annotation->serialize(test_dump_basename_row_compressed_to_rowflat);

    annotation = annotate::convert<annotate::RowFlatAnnotator>(
        test_dump_basename_row_compressed_to_rowflat
    ).release();
}

TEST_F(ConvertFromRowCompressed, stream_to_RowFlat2) {
    const auto rowflat_filename = test_dump_basename_row_compressed_to_rowflat
                                                + annotate::RowCompressed<>::kExtension;

    initial_annotation->serialize(rowflat_filename);

    annotation = annotate::convert<annotate::RowFlatAnnotator>(
        rowflat_filename
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_Rainbowfish) {
//     annotate::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert<annotate::RainbowfishAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromRowCompressed, to_RainbowfishAnnotator) {
    annotation = annotate::convert<annotate::RainbowfishAnnotator>(
        std::move(*initial_annotation)
    ).release();
}

// TEST(ConvertFromRowCompressedEmpty, to_GreedyBRWT) {
//     annotate::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

// TEST_F(ConvertFromRowCompressed, to_GreedyBRWT) {
//     annotation = annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
//         std::move(*initial_annotation)
//     ).release();
// }

TEST_F(MergeAnnotators, RowCompressed) {
    std::vector<std::string> filenames;
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_1";
        input_annotation_1->serialize(filename);
        filenames.push_back(filename + annotate::RowCompressed<>::kExtension);
    }
    {
        const std::string filename = test_dump_basename_row_compressed_merge + "_2";
        input_annotation_2->serialize(filename);
        filenames.push_back(filename + annotate::RowCompressed<>::kExtension);
    }

    annotate::merge<annotate::RowCompressed<>, std::string>(
        {}, filenames, test_dump_basename_row_compressed_merge + "_merged"
    );

    merged_annotation = new annotate::RowCompressed<>(num_rows);
    merged_annotation->merge_load({ test_dump_basename_row_compressed_merge + "_merged" });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowCompressed) {
    std::vector<std::unique_ptr<annotate::MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(annotate::convert<annotate::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(annotate::convert<annotate::RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowcompressed";
    annotate::merge<annotate::RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new annotate::RowCompressed<>(num_rows);
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, RowFlat_to_RowFlat) {
    std::vector<std::unique_ptr<annotate::MultiLabelEncoded<std::string>>> row_flat_annotators;
    {
        row_flat_annotators.push_back(annotate::convert<annotate::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        ));
    }
    {
        row_flat_annotators.push_back(annotate::convert<annotate::RowFlatAnnotator>(
            std::move(*input_annotation_2)
        ));
    }

    const auto filename = test_dump_basename_rowflat_merge + "_to_rowflat";
    annotate::merge<annotate::RowFlatAnnotator, std::string>(std::move(row_flat_annotators), {}, filename);

    merged_annotation = new annotate::RowFlatAnnotator();
    merged_annotation->merge_load({ filename });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST_F(MergeAnnotators, Mixed_to_RowFlat) {
    std::vector<std::unique_ptr<annotate::MultiLabelEncoded<std::string>>> annotators;
    std::vector<std::string> filenames;
    {
        auto annotator = annotate::convert<annotate::RowFlatAnnotator>(
            std::move(*input_annotation_1)
        );
        annotators.push_back(std::move(annotator));
    }
    {
        auto annotator = std::make_unique<annotate::ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    //TODO
    {
        auto annotator = std::make_unique<annotate::ColumnCompressed<> >(5);
        annotator->add_labels({ 1 }, {"Label0", "Label3"});
        annotator->add_labels({ 2 }, {"Label0", "Label9", "Label7"});
        annotator->add_labels({ 4 }, {"Label1", "Label3", "Label9", "Label10", "Label5", "Label6", "Label11", "Label12", "Label13", "Label14", "Label15", "Label16"});
        annotators.push_back(std::move(annotator));
    }
    {
        //TODO: move into fixture as input_annotation_3 and make non-overlapping
        const std::string filename = test_dump_basename_row_compressed_merge + "_mixed_2";
        auto annotation = std::make_unique<annotate::RowCompressed<> >(num_rows);
        annotation->add_labels({ 0 }, {"Label0"});
        annotation->add_labels({ 2 }, {"Label1"});
        annotation->add_labels({ 3 }, {"Label1"});
        annotation->add_labels({ 4 }, {"Label2"});
        annotation->serialize(filename);
        filenames.push_back(filename + annotate::RowCompressed<>::kExtension);
    }

    const auto outfile = test_dump_basename_rowflat_merge + "_mixed_to_rowflat";
    annotate::merge<annotate::RowFlatAnnotator, std::string>(std::move(annotators), filenames, outfile);

    merged_annotation = new annotate::RowFlatAnnotator();
    merged_annotation->merge_load({ outfile });
    EXPECT_EQ(num_rows, merged_annotation->num_objects());
}

TEST(ColumnCompressed, ToRowAnnotator) {
    {
        annotate::ColumnCompressed<> annotation(0);

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        annotate::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorParallel) {
    {
        annotate::ColumnCompressed<> annotation(0);

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        annotate::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }

        annotate::RowCompressed<> row_annotator(annotation.num_objects());
        convert_to_row_annotator(annotation, &row_annotator, 10);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreaming) {
    {
        annotate::ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        annotate::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreamingParallel) {
    {
        annotate::ColumnCompressed<> annotation(0);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});
        annotation.add_labels({ 2 }, {"Label1", "Label2"});
        annotation.add_labels({ 3 }, {});
        annotation.add_labels({ 4 }, {"Label1", "Label2", "Label8"});
        annotation.add_labels({ 5 }, {"Label2"});

        convert_to_row_annotator(annotation, test_dump_basename, 10);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        const size_t num_rows = 20'000'000;
        annotate::ColumnCompressed<> annotation(num_rows, 5);
        for (size_t i = 0; i < num_rows; ++i) {
            annotation.add_labels({ 0 }, {"Label0", "Label2"});
        }
        convert_to_row_annotator(annotation, test_dump_basename, 10);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

} // namespace
