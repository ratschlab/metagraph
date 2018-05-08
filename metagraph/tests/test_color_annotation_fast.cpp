#include <random>

#include "gtest/gtest.h"

#include "annotate_color_compressed_fast.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_color_compressed";


TEST(FastColorCompressed, EmptyConstructor) {
    annotate::FastColorCompressed<> annotation(5);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

TEST(FastColorCompressed, add_color) {
    annotate::FastColorCompressed<> annotation(5);
    annotation.add_color(0, "0");
    annotation.add_color(1, "0");
    annotation.add_color(2, "1");
    annotation.add_color(1, "2");

    EXPECT_EQ(convert_to_set({"0"}), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({"0", "2"}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({"1"}), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(4)));
}

TEST(FastColorCompressed, set_coloring) {
    std::unique_ptr<annotate::MultiColorAnnotation<uint64_t, std::string>> annotation(
        new annotate::FastColorCompressed<>(5)
    );
    annotation->set_coloring(0, { "Label0", "Label2", "Label8" });
    annotation->set_coloring(2, { "Label1", "Label2" });
    annotation->set_coloring(4, { "Label8" });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_coloring(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_coloring(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_coloring(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_coloring(4)));
}

TEST(FastColorCompressed, Serialization) {
    {
        annotate::FastColorCompressed<> annotation(5);
        annotation.set_coloring(0, { "Label0", "Label2", "Label8" });
        annotation.set_coloring(2, { "Label1", "Label2" });
        annotation.set_coloring(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        annotate::FastColorCompressed<> annotation(5);
        ASSERT_FALSE(annotation.load(test_dump_basename_vec_bad));
        ASSERT_TRUE(annotation.load(test_dump_basename_vec_good));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(FastColorCompressed, MergeLoad) {
    {
        annotate::FastColorCompressed<> annotation(5);
        annotation.set_coloring(0, { "Label0", "Label2", "Label8" });
        annotation.set_coloring(2, { "Label1", "Label2" });
        annotation.set_coloring(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_1");
    }
    {
        annotate::FastColorCompressed<> annotation(5);
        annotation.set_coloring(1, { "Label0", "Label2", "Label8" });
        annotation.set_coloring(2, { "Label1", "Label9", "Label0" });
        annotation.set_coloring(3, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_2");
    }
    {
        annotate::FastColorCompressed<> annotation(0);
        ASSERT_TRUE(annotation.merge_load({ test_dump_basename_vec_good + "_1",
                                            test_dump_basename_vec_good + "_2" }));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label9", "Label0" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(FastColorCompressed, has_colors) {
    std::unique_ptr<annotate::MultiColorAnnotation<uint64_t, std::string>> annotation(
        new annotate::FastColorCompressed<>(5)
    );
    annotation->set(0, {"Label0", "Label2", "Label8"});
    annotation->set(2, {"Label1", "Label2"});
    annotation->set(4, {"Label8"});

    EXPECT_FALSE(annotation->has_colors(0, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_colors(0, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_TRUE(annotation->has_colors(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(annotation->has_colors(0, { "Label0", "Label8" }));
    EXPECT_TRUE(annotation->has_colors(0, { "Label2" }));
    EXPECT_TRUE(annotation->has_colors(0, {}));

    EXPECT_FALSE(annotation->has_colors(1, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_colors(1, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_colors(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_colors(1, { "Label0", "Label8" }));
    EXPECT_FALSE(annotation->has_colors(1, { "Label2" }));
    EXPECT_TRUE(annotation->has_colors(1, {}));

    EXPECT_FALSE(annotation->has_colors(2, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_colors(2, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_colors(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_colors(2, { "Label1", "Label8" }));
    EXPECT_TRUE(annotation->has_colors(2, { "Label1", "Label2" }));
    EXPECT_TRUE(annotation->has_colors(2, { "Label2" }));
    EXPECT_TRUE(annotation->has_colors(2, {}));
}

TEST(FastColorCompressed, add_color_sequential) {
    size_t graph_half_size = 1000;
    annotate::FastColorCompressed<> annotation(graph_half_size * 2);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_color(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_color(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_coloring(i).size());
    }
}

TEST(FastColorCompressed, add_color_random) {
    size_t graph_half_size = 1000;
    annotate::FastColorCompressed<> annotation(graph_half_size * 2);

    std::vector<std::string> colors { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_color(i, colors[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_coloring(i).size());
    }
}

TEST(FastColorCompressed, add_color_random_with_caching) {
    size_t graph_half_size = 1000;
    annotate::FastColorCompressed<> annotation(graph_half_size * 2, 2);

    std::vector<std::string> colors { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_color(i, colors[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_coloring(i).size());
    }
}

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector);

TEST(FastColorCompressed, get_most_frequent_colors) {
    annotate::FastColorCompressed<> annotation(5);

    annotation.add_colors(0, {"Label0", "Label2", "Label8"});
    annotation.add_colors(2, {"Label1", "Label2"});
    annotation.add_colors(3, {"Label1", "Label2", "Label8"});
    annotation.add_colors(4, {"Label2", "Label8"});

    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }));

    EXPECT_EQ(to_set(VectorCounts({ std::make_pair("Label1", 1),
                                    std::make_pair("Label2", 1) })),
              to_set(annotation.get_most_frequent_colors({ 2 })));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 1000));
}

TEST(FastColorCompressed, aggregate_colors) {
    annotate::FastColorCompressed<> annotation(5);

    annotation.add_colors(0, {"Label0", "Label2", "Label8"});
    annotation.add_colors(2, {"Label1", "Label2"});
    annotation.add_colors(3, {"Label1", "Label2", "Label8"});
    annotation.add_colors(4, {"Label2"});

    EXPECT_EQ(std::vector<std::string>({}),
              annotation.aggregate_colors({}));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 })));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 1)));
}

TEST(FastColorCompressed, get_most_frequent_colors_fast) {
    annotate::FastColorCompressed<> annotation(5);

    annotation.add_colors(0, {"Label0", "Label2", "Label8"});
    annotation.add_colors(2, {"Label1", "Label2"});
    annotation.add_colors(3, {"Label1", "Label2", "Label8"});
    annotation.add_colors(4, {"Label2", "Label8"});

    annotation.update_index();

    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }));

    EXPECT_EQ(to_set(VectorCounts({ std::make_pair("Label1", 1),
                                    std::make_pair("Label2", 1) })),
              to_set(annotation.get_most_frequent_colors({ 2 })));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_most_frequent_colors({ 0, 1, 2, 3, 4 }, 1000));
}

TEST(FastColorCompressed, aggregate_colors_fast) {
    annotate::FastColorCompressed<> annotation(5);

    annotation.add_colors(0, {"Label0", "Label2", "Label8"});
    annotation.add_colors(2, {"Label1", "Label2"});
    annotation.add_colors(3, {"Label1", "Label2", "Label8"});
    annotation.add_colors(4, {"Label2"});

    annotation.update_index();

    EXPECT_EQ(std::vector<std::string>({}),
              annotation.aggregate_colors({}));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 })));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 1)));
}

TEST(FastColorCompressed, initialize_from_ColorCompressed) {
    annotate::ColorCompressed<> annotation_init(5);

    annotation_init.add_colors(0, {"Label0", "Label2", "Label8"});
    annotation_init.add_colors(2, {"Label1", "Label2"});
    annotation_init.add_colors(3, {"Label1", "Label2", "Label8"});
    annotation_init.add_colors(4, {"Label2"});

    annotate::FastColorCompressed<> annotation(std::move(annotation_init));

    EXPECT_EQ(std::vector<std::string>({}),
              annotation.aggregate_colors({}));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 })));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 })));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.aggregate_colors({ 0, 1, 2, 3, 4 }, 1)));
}
