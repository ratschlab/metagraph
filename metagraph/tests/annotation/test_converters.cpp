#include <random>

#include "gtest/gtest.h"

#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "utils.hpp"
#include "binary_matrix.hpp"


const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_row_compressed_merge = test_dump_basename + "_row_compressed_merge";


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

TEST(RowCompressed, Merge) {
    std::vector<std::string> filenames;
    {
        auto annotation = new annotate::RowCompressed<>(5);
        annotation->add_labels(0, {"Label0", "Label2", "Label8"});
        annotation->add_labels(2, {"Label1", "Label2"});
        annotation->add_labels(3, {"Label1", "Label2", "Label8"});
        annotation->add_labels(4, {"Label2"});
        const std::string filename = test_dump_basename_row_compressed_merge + "_1";
        annotation->serialize(filename);
        filenames.push_back(filename + annotate::kRowAnnotatorExtension);
    }
    {
        auto annotation = new annotate::RowCompressed<>(5);
        annotation->add_labels(1, {"Label0", "Label3"});
        annotation->add_labels(2, {"Label0", "Label9", "Label7"});
        annotation->add_labels(4, {"Label1", "Label3", "Label9", "Label10", "9", "a", "b", "c", "d", "e", "f", "10"});
        const std::string filename = test_dump_basename_row_compressed_merge + "_2";
        annotation->serialize(filename);
        filenames.push_back(filename + annotate::kRowAnnotatorExtension);
    }
    {
        auto annotation = new annotate::RowCompressed<>(5);
        annotation->add_labels(0, {"Label0", "Label2", "Label8"});
        annotation->add_labels(1, {"Label0", "Label3"});
        annotation->add_labels(2, {"Label0", "Label2", "Label1", "Label9", "Label7"});
        annotation->add_labels(3, {"Label2", "Label8", "Label1"});
        annotation->add_labels(4, {"Label2", "Label1", "Label3", "Label9", "Label10", "9", "a", "b", "c", "d", "e", "f", "10"});
        const std::string filename = test_dump_basename_row_compressed_merge + "_TEST";
        annotation->serialize(filename);
    }
    //{
    //    auto annotation = new annotate::RowCompressed<>(5);
    //    annotation->add_labels(0, {"Label0", "Label1", "Label2", "Label3","Label7", "Label8", "Label9", "Label10"});
    //    annotation->add_labels(1, {"Label0", "Label1", "Label2", "Label3","Label7", "Label8", "Label9", "Label10"});
    //    annotation->add_labels(2, {"Label0", "Label1", "Label2", "Label3","Label7", "Label8", "Label9", "Label10"});
    //    annotation->add_labels(3, {"Label0", "Label1", "Label2", "Label3","Label7", "Label8", "Label9", "Label10"});
    //    annotation->add_labels(4, {"Label0", "Label1", "Label2", "Label3","Label7", "Label8", "Label9", "Label10"});
    //    const std::string filename = test_dump_basename_row_compressed_merge + "_TEST";
    //    annotation->serialize(filename);
    //    filenames.push_back(filename + annotate::kRowAnnotatorExtension);
    //}
    typedef annotate::LabelEncoder<std::string> LEncoder;
    std::vector<LEncoder*> label_encoders;

    for(auto filename : filenames) {
        LEncoder* label_encoder = annotate::RowCompressed<std::string>::load_label_encoder(filename);
        label_encoders.push_back(label_encoder);
        std::cout << "-------" << std::endl;
        for(size_t i = 0; i < label_encoder->size(); ++i) {
            std::cout << label_encoder->decode(i) << std::endl;
        }
    }
    LEncoder* merged_label_enc { new LEncoder() };
    merged_label_enc->merge(label_encoders);
    std::cout << "-------" << std::endl;
    for(size_t i = 0; i < merged_label_enc->size(); ++i) {
        std::cout << merged_label_enc->decode(i) << std::endl;
    }

    std::cout << "=======" << std::endl;
    //uint64_t num_rows = 0;
    //for (int file_idx = 0; file_idx <= 1; ++file_idx) {
    //    num_rows = 0;
    //    auto callback = [&](const std::vector<uint64_t> &row) {
    //        std::for_each(row.begin(), row.end(), [&](uint64_t i) { std::cout << label_encoders.at(file_idx)->decode(i) << std::endl; });
    //        std::cout << "," << std::endl;
    //        ++num_rows;
    //    };
    //    ASSERT_TRUE(annotate::RowCompressed<>::stream_rows(filenames.at(file_idx), callback, false));
    //    std::cout << "#######" << std::endl;
    //}

    std::cout << "#######" << std::endl;
    // auto callback2 = []() -> const std::vector<uint64_t> {
    //     std::vector<uint64_t> v;
    //     v.push_back(0);
    //     v.push_back(1);
    //     v.push_back(2);
    //     v.push_back(3);
    //     v.push_back(4);
    //     v.push_back(5);
    //     v.push_back(6);
    //     v.push_back(7);
    //     return v;
    // };
    // annotate::RowCompressed<>::write_rows(test_dump_basename_row_compressed_merge + "_merged",
    //                                       *merged_label_enc,
    //                                       num_rows,
    //                                       callback2,
    //                                       false);
    std::cout << "%%%%%%%" << std::endl;
    // auto callback3 = [&](const std::function<void(void)> &init, const std::function<void (const std::vector<uint64_t> &)> &row_callback, const std::function<void(void)> &finalize) {
    //     init();
    //     for(size_t i = 0; i < filenames.size(); ++i) {
    //         std::set<uint64_t> label_set;

    //         annotate::RowCompressed<>
    //         ::stream_rows(filenames.at(i), [&](const std::vector<uint64_t> &row) {
    //             for(auto label : row)
    //                 label_set.insert(merged_label_enc->encode(label_encoders.at(i)->decode(label)));
    //         }, false);
    //         std::vector<uint64_t> merged_row(label_set.begin(), label_set.end());
    //         row_callback(merged_row);
    //     }
    //     finalize();
    // };
    // annotate::RowCompressed<>::write_rows(test_dump_basename_row_compressed_merge + "_merged",
    //                                       *merged_label_enc,
    //                                       num_rows,
    //                                       callback3,
    //                                       false);
    // annotate::RowCompressed<> written(num_rows);
    // written.merge_load({ test_dump_basename_row_compressed_merge + "_merged" });
    // //std::cout << "a" << std::endl;
    // for(uint64_t i = 0; i < num_rows; ++i) {
    //     auto row = written.get_labels(i);
    //     std::for_each(row.begin(), row.end(), [](const std::string s) { std::cout << s << std::endl; });
    //     std::cout << "," << std::endl;
    // }
    std::vector<annotate::RowCompressed<>::StreamRows*> annotators;
    std::vector<std::vector<uint64_t> > reencode;

    for(size_t i = 0; i < filenames.size(); ++i) {

        auto annotator = new annotate::RowCompressed<>::StreamRows(filenames.at(i), false);
        annotators.push_back(annotator);

        std::vector<uint64_t> v;
        for(size_t j = 0; j < label_encoders.at(i)->size(); ++j)
            v.push_back(merged_label_enc->encode(label_encoders.at(i)->decode(j)));
        reencode.push_back(v);
    }

    auto callback3 = [&](const std::function<void (const std::vector<uint64_t> &)> &row_callback) {
        std::set<uint64_t> label_set;
        bool done = false;
        while (!done) {
            label_set.clear();
            for(size_t i = 0; i < annotators.size(); ++i) {
                auto row = annotators[i]->next_row();
                if(!row) {
                    done = true;
                    break;
                }
                for(auto label : *row)
                    label_set.insert(merged_label_enc->encode(label_encoders.at(i)->decode(label)));
            }
            std::vector<uint64_t> merged_row(label_set.begin(), label_set.end());
            row_callback(merged_row);
            // std::for_each(merged_row.begin(), merged_row.end(), [&](uint64_t s) { std::cout << merged_label_enc->decode(s) << std::endl; });
            // std::cout << "," << std::endl;
        }
    };

    // TODO: stream_counts to get num_rows
    const uint64_t num_rows = 5;
    annotate::RowCompressed<>::write_rows(test_dump_basename_row_compressed_merge + "_merged",
                                          *merged_label_enc,
                                          num_rows,
                                          callback3,
                                          false);
    // int i = 0;
    // std::unique_ptr<std::vector<VectorRowBinMat::Row> > row;
    // while(row = sr.next_row()) {
    //     if(++i > 20)
    //        break;
    //     std::for_each(row->begin(), row->end(), [](uint64_t s) { std::cout << s << std::endl; });
    //     std::cout << "," << std::endl;
    // }

    annotate::RowCompressed<> written(num_rows);
    written.merge_load({ test_dump_basename_row_compressed_merge + "_merged" });
    //std::cout << "a" << std::endl;
    for(uint64_t i = 0; i < num_rows; ++i) {
        auto row = written.get_labels(i);
        std::for_each(row.begin(), row.end(), [](const std::string s) { std::cout << s << std::endl; });
        std::cout << "," << std::endl;
    }
}
