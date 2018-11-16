#include <random>

#include "gtest/gtest.h"

#include "annotate_row_compressed.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_row_compressed";


TEST(RowCompressed, EmptyConstructor) {
    annotate::RowCompressed<> annotation(5, false);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

TEST(RowCompressed, add_label) {
    annotate::RowCompressed<> annotation(5, false);
    annotation.add_label(0, "0");
    annotation.add_label(1, "0");
    annotation.add_label(2, "1");
    annotation.add_label(1, "2");

    EXPECT_EQ(convert_to_set({"0"}), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({"0", "2"}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({"1"}), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(4)));
}

TEST(RowCompressed, set_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));
}

TEST(RowCompressed, insert_no_empty_rows) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));

    annotation->insert_rows({});

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));
}

TEST(RowCompressed, insert_first_column_to_empty_annotation) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(0, false)
    );

    ASSERT_EQ(0u, annotation->num_labels());

    annotation->insert_rows({ 0, });
    ASSERT_EQ(0u, annotation->num_labels());

    annotation->set_labels(0, { "Label0", "Label2", "Label8" });

    ASSERT_EQ(3u, annotation->num_labels());
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
}

TEST(RowCompressed, insert_one_empty_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));

    annotation->insert_rows({ 4, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(5)));
}

TEST(RowCompressed, insert_first_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));

    annotation->insert_rows({ 0, });

    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(4)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(5)));
}

TEST(RowCompressed, insert_last_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));

    annotation->insert_rows({ 5, });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(5)));
}

TEST(RowCompressed, insert_empty_rows) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    EXPECT_EQ(5u, annotation->num_objects());
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    ASSERT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    ASSERT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    ASSERT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    ASSERT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));

    annotation->insert_rows({ 1, 2, 4, 5 });
    EXPECT_EQ(9u, annotation->num_objects());

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(4)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(5)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(6)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(7)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(8)));
}

TEST(RowCompressed, Serialization) {
    {
        annotate::RowCompressed<> annotation(5, false);
        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label2" });
        annotation.set_labels(4, { "Label8" });

        annotation.serialize(test_dump_basename + "_row_compressed");
    }
    {
        annotate::RowCompressed<> annotation(5, false);
        ASSERT_FALSE(annotation.load(test_dump_basename_vec_bad));
        ASSERT_TRUE(annotation.load(test_dump_basename_vec_good));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(RowCompressed, has_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, false)
    );
    annotation->set(0, {"Label0", "Label2", "Label8"});
    annotation->set(2, {"Label1", "Label2"});
    annotation->set(4, {"Label8"});

    EXPECT_FALSE(annotation->has_labels(0, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(0, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_TRUE(annotation->has_labels(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(0, { "Label0", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(0, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(0, {}));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(1, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(1, {}));

    EXPECT_FALSE(annotation->has_labels(2, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(2, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_labels(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(2, { "Label1", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(2, { "Label1", "Label2" }));
    EXPECT_TRUE(annotation->has_labels(2, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(2, {}));
}

TEST(RowCompressed, add_label_sequential) {
    size_t graph_half_size = 1000;
    annotate::RowCompressed<> annotation(graph_half_size * 2, false);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_label(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TEST(RowCompressed, add_label_random) {
    size_t graph_half_size = 1000;
    annotate::RowCompressed<> annotation(graph_half_size * 2, false);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, labels[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector);

TEST(RowCompressed, get_top_labels) {
    annotate::RowCompressed<> annotation(5, false);

    annotation.add_labels(0, {"Label0", "Label2", "Label8"});
    annotation.add_labels(2, {"Label1", "Label2"});
    annotation.add_labels(3, {"Label1", "Label2", "Label8"});
    annotation.add_labels(4, {"Label2", "Label8"});

    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }));

    EXPECT_EQ(to_set(VectorCounts({ std::make_pair("Label1", 1),
                                    std::make_pair("Label2", 1) })),
              to_set(annotation.get_top_labels({ 2 })));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 1000));
}

TEST(RowCompressed, get_labels) {
    annotate::RowCompressed<> annotation(5, false);

    annotation.add_labels(0, {"Label0", "Label2", "Label8"});
    annotation.add_labels(2, {"Label1", "Label2"});
    annotation.add_labels(3, {"Label1", "Label2", "Label8"});
    annotation.add_labels(4, {"Label2"});

    EXPECT_EQ(std::vector<std::string>({}),
              annotation.get_labels({}, 1));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 1)));
}


TEST(RowCompressedSparse, EmptyConstructor) {
    annotate::RowCompressed<> annotation(5, true);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

TEST(RowCompressedSparse, add_label) {
    annotate::RowCompressed<> annotation(5, true);
    annotation.add_label(0, "0");
    annotation.add_label(1, "0");
    annotation.add_label(2, "1");
    annotation.add_label(1, "2");

    EXPECT_EQ(convert_to_set({"0"}), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({"0", "2"}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({"1"}), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(4)));
}

TEST(RowCompressedSparse, set_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, true)
    );
    annotation->set_labels(0, { "Label0", "Label2", "Label8" });
    annotation->set_labels(2, { "Label1", "Label2" });
    annotation->set_labels(4, { "Label8" });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation->get_labels(4)));
}

TEST(RowCompressedSparse, Serialization) {
    {
        annotate::RowCompressed<> annotation(5, true);
        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label2" });
        annotation.set_labels(4, { "Label8" });

        annotation.serialize(test_dump_basename + "_row_compressed");
    }
    {
        annotate::RowCompressed<> annotation(5, true);
        ASSERT_FALSE(annotation.load(test_dump_basename_vec_bad));
        ASSERT_TRUE(annotation.load(test_dump_basename_vec_good));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(RowCompressedSparse, has_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::RowCompressed<>(5, true)
    );
    annotation->set(0, {"Label0", "Label2", "Label8"});
    annotation->set(2, {"Label1", "Label2"});
    annotation->set(4, {"Label8"});

    EXPECT_FALSE(annotation->has_labels(0, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(0, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_TRUE(annotation->has_labels(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(0, { "Label0", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(0, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(0, {}));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(1, { "Label0", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(1, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(1, {}));

    EXPECT_FALSE(annotation->has_labels(2, { "Label0", "Label1",
                                             "Label2", "Label4",
                                             "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_labels(2, { "Label0",
                                             "Label2", "Label4",
                                             "Label8" }));

    EXPECT_FALSE(annotation->has_labels(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_labels(2, { "Label1", "Label8" }));
    EXPECT_TRUE(annotation->has_labels(2, { "Label1", "Label2" }));
    EXPECT_TRUE(annotation->has_labels(2, { "Label2" }));
    EXPECT_TRUE(annotation->has_labels(2, {}));
}

TEST(RowCompressedSparse, add_label_sequential) {
    size_t graph_half_size = 1000;
    annotate::RowCompressed<> annotation(graph_half_size * 2, true);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_label(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TEST(RowCompressedSparse, add_label_random) {
    size_t graph_half_size = 1000;
    annotate::RowCompressed<> annotation(graph_half_size * 2, true);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, labels[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector);

TEST(RowCompressedSparse, get_top_labels) {
    annotate::RowCompressed<> annotation(5, true);

    annotation.add_labels(0, {"Label0", "Label2", "Label8"});
    annotation.add_labels(2, {"Label1", "Label2"});
    annotation.add_labels(3, {"Label1", "Label2", "Label8"});
    annotation.add_labels(4, {"Label2", "Label8"});

    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }));

    EXPECT_EQ(to_set(VectorCounts({ std::make_pair("Label1", 1),
                                    std::make_pair("Label2", 1) })),
              to_set(annotation.get_top_labels({ 2 })));

    EXPECT_EQ(VectorCounts({}),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              annotation.get_top_labels({ 0, 1, 2, 3, 4 }, 1000));
}

TEST(RowCompressedSparse, get_labels) {
    annotate::RowCompressed<> annotation(5, true);

    annotation.add_labels(0, {"Label0", "Label2", "Label8"});
    annotation.add_labels(2, {"Label1", "Label2"});
    annotation.add_labels(3, {"Label1", "Label2", "Label8"});
    annotation.add_labels(4, {"Label2"});

    EXPECT_EQ(std::vector<std::string>({}),
              annotation.get_labels({}, 1));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(annotation.get_labels({ 0, 1, 2, 3, 4 }, 1)));
}
