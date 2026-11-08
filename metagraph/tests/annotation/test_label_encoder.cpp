#include <sstream>
#include <vector>
#include <string>
#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "common/serialization.hpp"
#include "common/algorithms.hpp"
#include "annotation/representation/base/annotation.hpp"


namespace {

using namespace mtg;
using namespace mtg::annot;

// Helper function to serialize LabelEncoder in the old format
// Old format: serialize_string_vector(labels), serialize_number_vector(values), serialize_string_vector(labels)
void serialize_label_encoder_old_format(std::ostream &outstream,
                                        const std::vector<std::string> &labels) {
    // First serialize labels
    serialize_string_vector(outstream, labels);
    // Then serialize values (indices 0, 1, 2, ...)
    serialize_number_vector(outstream, utils::arange<uint64_t>(0, labels.size()));
    // Finally serialize labels again
    serialize_string_vector(outstream, labels);
}

TEST(LabelEncoderBackwardCompatibility, LoadOldFormatEmpty) {
    // Create old format data for empty LabelEncoder
    std::stringstream stream;
    std::vector<std::string> empty_labels;
    serialize_label_encoder_old_format(stream, empty_labels);
    serialize_label_encoder_old_format(stream, empty_labels);
    stream.seekg(0);

    // Load it
    LabelEncoder<std::string> encoder;
    ASSERT_TRUE(encoder.load(stream));
    // Verify it's empty
    EXPECT_EQ(0u, encoder.size());

    ASSERT_EQ(std::vector<std::string>(), encoder.get_labels());

    ASSERT_TRUE(encoder.load(stream));
    ASSERT_EQ(std::vector<std::string>(), encoder.get_labels());
}

TEST(LabelEncoderBackwardCompatibility, LoadOldFormat) {
    // Create old format data for non-empty LabelEncoder
    std::stringstream stream;
    std::vector<std::string> labels = {"label1", "label2", "label3"};
    serialize_label_encoder_old_format(stream, labels);
    serialize_label_encoder_old_format(stream, labels);
    stream.seekg(0);

    // Load it
    LabelEncoder<std::string> encoder;
    ASSERT_TRUE(encoder.load(stream));

    // Verify the labels
    EXPECT_EQ(3u, encoder.size());
    EXPECT_TRUE(encoder.label_exists("label1"));
    EXPECT_TRUE(encoder.label_exists("label2"));
    EXPECT_TRUE(encoder.label_exists("label3"));

    // Verify encoding/decoding
    EXPECT_EQ(0u, encoder.encode("label1"));
    EXPECT_EQ(1u, encoder.encode("label2"));
    EXPECT_EQ(2u, encoder.encode("label3"));

    EXPECT_EQ("label1", encoder.decode(0));
    EXPECT_EQ("label2", encoder.decode(1));
    EXPECT_EQ("label3", encoder.decode(2));

    ASSERT_EQ(std::vector<std::string>({"label1", "label2", "label3"}), encoder.get_labels());

    ASSERT_TRUE(encoder.load(stream));
    ASSERT_EQ(std::vector<std::string>({"label1", "label2", "label3"}), encoder.get_labels());
}

TEST(LabelEncoderBackwardCompatibility, LoadNewFormatEmpty) {
    // Create new format data for empty LabelEncoder
    std::stringstream stream;
    LabelEncoder<std::string> empty_encoder;
    empty_encoder.serialize(stream);
    empty_encoder.serialize(stream);
    stream.seekg(0);

    // Load it
    LabelEncoder<std::string> encoder;
    ASSERT_TRUE(encoder.load(stream));
    // Verify it's empty
    EXPECT_EQ(0u, encoder.size());

    ASSERT_EQ(std::vector<std::string>(), encoder.get_labels());

    ASSERT_TRUE(encoder.load(stream));
    ASSERT_EQ(std::vector<std::string>(), encoder.get_labels());
}

TEST(LabelEncoderBackwardCompatibility, LoadNewFormat) {
    // Create new format data for non-empty LabelEncoder
    std::stringstream stream;
    LabelEncoder<std::string> original_encoder;
    original_encoder.insert_and_encode("label1");
    original_encoder.insert_and_encode("label2");
    original_encoder.insert_and_encode("label3");
    original_encoder.serialize(stream);
    original_encoder.serialize(stream);
    stream.seekg(0);

    // Load it
    LabelEncoder<std::string> encoder;
    ASSERT_TRUE(encoder.load(stream));

    // Verify the labels
    EXPECT_EQ(3u, encoder.size());
    EXPECT_TRUE(encoder.label_exists("label1"));
    EXPECT_TRUE(encoder.label_exists("label2"));
    EXPECT_TRUE(encoder.label_exists("label3"));

    // Verify encoding/decoding
    EXPECT_EQ(0u, encoder.encode("label1"));
    EXPECT_EQ(1u, encoder.encode("label2"));
    EXPECT_EQ(2u, encoder.encode("label3"));

    EXPECT_EQ("label1", encoder.decode(0));
    EXPECT_EQ("label2", encoder.decode(1));
    EXPECT_EQ("label3", encoder.decode(2));

    ASSERT_EQ(std::vector<std::string>({"label1", "label2", "label3"}), encoder.get_labels());

    ASSERT_TRUE(encoder.load(stream));
    ASSERT_EQ(std::vector<std::string>({"label1", "label2", "label3"}), encoder.get_labels());
}

// A string hash that differs from std::hash, standing in for another
// toolchain's std::hash (e.g. libc++ vs libstdc++).
struct OtherToolchainHash {
    size_t operator()(const std::string &s) const {
        return std::hash<std::string>()(s) * 0x9E3779B97F4A7C15ull + 0x632BE59BD9B4E019ull;
    }
};

TEST(LabelEncoder, LoadBucketsWrittenWithAnotherHash) {
    std::vector<std::string> labels;
    for (size_t i = 0; i < 1000; ++i) {
        labels.push_back("label" + std::to_string(i));
    }
    tsl::ordered_set<std::string, OtherToolchainHash, std::equal_to<std::string>,
                     std::allocator<std::string>, std::vector<std::string>, uint64_t>
            written(labels.begin(), labels.end());

    std::stringstream stream;
    stream.write("LE-v2.0", 7);
    Serializer serializer(stream);
    written.serialize(serializer);
    // a second record after the first one must still be readable
    LabelEncoder<std::string> second;
    second.insert_and_encode("next");
    second.serialize(stream);
    stream.seekg(0);

    LabelEncoder<std::string> encoder;
    ASSERT_TRUE(encoder.load(stream));
    ASSERT_EQ(labels, encoder.get_labels());
    for (size_t i = 0; i < labels.size(); ++i) {
        ASSERT_TRUE(encoder.label_exists(labels[i])) << labels[i];
        ASSERT_EQ(i, encoder.encode(labels[i]));
    }

    ASSERT_TRUE(encoder.load(stream));
    EXPECT_EQ(std::vector<std::string>({ "next" }), encoder.get_labels());
}

} // namespace
