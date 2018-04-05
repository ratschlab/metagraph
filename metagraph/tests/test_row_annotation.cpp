#include <random>

#include "gtest/gtest.h"

#include "annotate.hpp"
#include "annotate_row_compressed.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";


TEST(RowCompressed, EmptyConstructor) {
    annotate::RowCompressed annotation(5);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

TEST(RowCompressed, add_label) {
    annotate::RowCompressed annotation(5);
    annotation.add_label(0, "0");
    annotation.add_label(1, "0");
    annotation.add_label(2, "1");
    annotation.add_label(1, "2");

    EXPECT_EQ(std::set<std::string>({"0"}), annotation.get(0));
    EXPECT_EQ(std::set<std::string>({"0", "2"}), annotation.get(1));
    EXPECT_EQ(std::set<std::string>({"1"}), annotation.get(2));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(4));
}

TEST(RowCompressed, set_label) {
    annotate::RowCompressed annotation(5);
    annotation.set_label(0, { "Label0", "Label2", "Label8" });
    annotation.set_label(2, { "Label1", "Label2" });
    annotation.set_label(4, { "Label8" });

    EXPECT_EQ(std::set<std::string>({ "Label0", "Label2", "Label8" }), annotation.get(0));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(1));
    EXPECT_EQ(std::set<std::string>({ "Label1", "Label2" }), annotation.get(2));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
    EXPECT_EQ(std::set<std::string>({ "Label8" }), annotation.get(4));
}

TEST(RowCompressed, Serialization) {
    {
        annotate::RowCompressed annotation(5);
        annotation.set_label(0, { "Label0", "Label2", "Label8" });
        annotation.set_label(2, { "Label1", "Label2" });
        annotation.set_label(4, { "Label8" });

        annotation.serialize(test_dump_basename + "_color_compressed");
    }
    {
        annotate::RowCompressed annotation(5);
        ASSERT_FALSE(annotation.load(test_dump_basename + "_bad_file"));
        ASSERT_TRUE(annotation.load(test_dump_basename + "_color_compressed"));

        EXPECT_EQ(std::set<std::string>({ "Label0", "Label2", "Label8" }), annotation.get(0));
        EXPECT_EQ(std::set<std::string>({}), annotation.get(1));
        EXPECT_EQ(std::set<std::string>({ "Label1", "Label2" }), annotation.get(2));
        EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
        EXPECT_EQ(std::set<std::string>({ "Label8" }), annotation.get(4));
    }
}

TEST(RowCompressed, has_label) {
    std::unique_ptr<annotate::AnnotationCategory<std::set<std::string>>> annotation(
        new annotate::RowCompressed(5)
    );
    annotation->set_label(0, {"Label0", "Label2", "Label8"});
    annotation->set_label(2, {"Label1", "Label2"});
    annotation->set_label(4, {"Label8"});

    EXPECT_FALSE(annotation->has_label(0, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(0, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_TRUE(annotation->has_label(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(annotation->has_label(0, { "Label0", "Label8" }));
    EXPECT_TRUE(annotation->has_label(0, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(0, {}));

    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(1, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label8" }));
    EXPECT_FALSE(annotation->has_label(1, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(1, {}));

    EXPECT_FALSE(annotation->has_label(2, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(2, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_FALSE(annotation->has_label(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_label(2, { "Label1", "Label8" }));
    EXPECT_TRUE(annotation->has_label(2, { "Label1", "Label2" }));
    EXPECT_TRUE(annotation->has_label(2, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(2, {}));
}

TEST(RowCompressed, set_label_cache) {
    size_t graph_half_size = 500;
    annotate::RowCompressed annotation(graph_half_size * 2);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_label(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TEST(RowCompressed, set_label_random) {
    size_t graph_half_size = 500;
    annotate::RowCompressed annotation(graph_half_size * 2);

    std::vector<std::string> colors { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, colors[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}
