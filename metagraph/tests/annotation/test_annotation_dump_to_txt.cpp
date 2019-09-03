#include "gtest/gtest.h"

#include "test_annotation.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_matrix";


template <typename Annotator>
class AnnotatorPresetDumpTest : public AnnotatorPreset2Test<Annotator> { };

typedef ::testing::Types<annotate::BinRelWTAnnotator,
                         annotate::BinRelWT_sdslAnnotator,
                         annotate::BRWTCompressed<>,
                         annotate::RainbowfishAnnotator,
                         annotate::RowFlatAnnotator,
                         annotate::ColumnCompressed<>> AnnotatorDumpTestTypes;

TYPED_TEST_CASE(AnnotatorPresetDumpTest, AnnotatorDumpTestTypes);


TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadText) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, false));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    uint64_t size, num_set_bits, pos;
    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".text.annodbg");
        ASSERT_TRUE(fin.good());

        fin >> size >> num_set_bits;

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            fin >> pos;
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadTextParallel) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, false, 3));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    uint64_t size, num_set_bits, pos;
    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".text.annodbg");
        ASSERT_TRUE(fin.good());

        fin >> size >> num_set_bits;

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            fin >> pos;
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadBinary) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, true));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".raw.annodbg",
                          std::ios::binary);
        ASSERT_TRUE(fin.good());

        uint64_t size = load_number(fin);
        uint64_t num_set_bits = load_number(fin);

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            size_t pos = load_number(fin);
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}

TYPED_TEST(AnnotatorPresetDumpTest, SerializationAndLoadBinaryParallel) {
    ASSERT_TRUE(this->annotation->dump_columns(test_dump_basename_vec_good, true, 3));

    annotate::ColumnCompressed<> loaded(this->annotation->num_objects());

    std::vector<std::string> labels;
    const auto& label_encoder = this->annotation->get_label_encoder();
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        labels.emplace_back(label_encoder.decode(i));
    }

    for (size_t i = 0; i < this->annotation->num_labels(); ++i) {
        std::ifstream fin(test_dump_basename_vec_good
                            + "." + std::to_string(i) + ".raw.annodbg",
                          std::ios::binary);
        ASSERT_TRUE(fin.good());

        uint64_t size = load_number(fin);
        uint64_t num_set_bits = load_number(fin);

        ASSERT_EQ(loaded.num_objects(), size);

        while (num_set_bits--) {
            size_t pos = load_number(fin);
            ASSERT_GT(this->annotation->num_objects(), pos);
            loaded.add_labels(pos, { labels[i] });
        }
    }

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(loaded.get(0)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }),           convert_to_set(loaded.get(2)));
    EXPECT_EQ(convert_to_set({}),                               convert_to_set(loaded.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }),                     convert_to_set(loaded.get(4)));
}
