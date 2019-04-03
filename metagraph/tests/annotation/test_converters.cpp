#include <random>

#include "gtest/gtest.h"

#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "utils.hpp"


class ConvertFromRowCompressed : public ::testing::Test {
  protected:
    static annotate::RowCompressed<> *initial_annotation;
    static annotate::MultiLabelEncoded<uint64_t, std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new annotate::RowCompressed<>(5);
        initial_annotation->add_labels(0, {"Label0", "Label2", "Label8"});
        initial_annotation->add_labels(2, {"Label1", "Label2"});
        initial_annotation->add_labels(3, {"Label1", "Label2", "Label8"});
        initial_annotation->add_labels(4, {"Label2"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
                  convert_to_set(annotation->get_labels(0)));
        EXPECT_EQ(std::vector<std::string>({}),
                  annotation->get_labels(1));
        EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
                  convert_to_set(annotation->get_labels(2)));
        EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
                  convert_to_set(annotation->get_labels(3)));
        EXPECT_EQ(convert_to_set({"Label2"}),
                  convert_to_set(annotation->get_labels(4)));

        delete initial_annotation;
        delete annotation;
    }

};

annotate::RowCompressed<> *ConvertFromRowCompressed::initial_annotation = nullptr;
annotate::MultiLabelEncoded<uint64_t, std::string> *ConvertFromRowCompressed::annotation = nullptr;


class ConvertFromColumnCompressed : public ::testing::Test {
  protected:
    static annotate::ColumnCompressed<> *initial_annotation;
    static annotate::MultiLabelEncoded<uint64_t, std::string> *annotation;

    virtual void SetUp() {
        initial_annotation = new annotate::ColumnCompressed<>(5);
        initial_annotation->add_labels(0, {"Label0", "Label2", "Label8"});
        initial_annotation->add_labels(2, {"Label1", "Label2"});
        initial_annotation->add_labels(3, {"Label1", "Label2", "Label8"});
        initial_annotation->add_labels(4, {"Label2"});
    }

    virtual void TearDown() {
        ASSERT_TRUE(annotation);
        ASSERT_EQ(4u, annotation->num_labels());
        ASSERT_EQ(5u, annotation->num_objects());
        ASSERT_EQ(9u, annotation->num_relations());

        EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
                  convert_to_set(annotation->get_labels(0)));
        EXPECT_EQ(std::vector<std::string>({}),
                  annotation->get_labels(1));
        EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
                  convert_to_set(annotation->get_labels(2)));
        EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
                  convert_to_set(annotation->get_labels(3)));
        EXPECT_EQ(convert_to_set({"Label2"}),
                  convert_to_set(annotation->get_labels(4)));

        delete initial_annotation;
        delete annotation;
    }

};

annotate::ColumnCompressed<> *ConvertFromColumnCompressed::initial_annotation = nullptr;
annotate::MultiLabelEncoded<uint64_t, std::string> *ConvertFromColumnCompressed::annotation = nullptr;


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
    ).release();;
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
    ).release();;
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
    ).release();;
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
    ).release();;
}

// TEST(ConvertFromColumnCompressedEmpty, to_GreedyBRWT) {
//     annotate::ColumnCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

TEST_F(ConvertFromColumnCompressed, to_GreedyBRWT) {
    annotation = annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
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
    ).release();;
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
    ).release();;
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
    ).release();;
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
    ).release();;
}

// TEST(ConvertFromRowCompressedEmpty, to_GreedyBRWT) {
//     annotate::RowCompressed<> empty_column_annotator(5);
//     auto empty_annotation = annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
//         std::move(empty_column_annotator)
//     );
//     EXPECT_EQ(0u, empty_annotation->num_labels());
//     EXPECT_EQ(5u, empty_annotation->num_objects());
//     EXPECT_EQ(0u, empty_annotation->num_relations());
// }

// TEST_F(ConvertFromRowCompressed, to_GreedyBRWT) {
//     annotation = annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
//         std::move(*initial_annotation)
//     ).release();
// }
