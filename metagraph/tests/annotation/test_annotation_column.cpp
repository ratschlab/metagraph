#include <random>

#include "gtest/gtest.h"

#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column_compressed";


TEST(ColumnCompressed, EmptyConstructor) {
    annotate::ColumnCompressed<> annotation(5);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

TEST(ColumnCompressed, add_label) {
    annotate::ColumnCompressed<> annotation(5);
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

TEST(ColumnCompressed, set_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, insert_no_empty_rows) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, insert_first_column_to_empty_annotation) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(0)
    );

    ASSERT_EQ(0u, annotation->num_labels());

    annotation->insert_rows({ 0, });
    ASSERT_EQ(0u, annotation->num_labels());

    annotation->set_labels(0, { "Label0", "Label2", "Label8" });

    ASSERT_EQ(3u, annotation->num_labels());
    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
}

TEST(ColumnCompressed, insert_one_empty_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, insert_first_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, insert_last_row) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, insert_empty_rows) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, Serialization) {
    {
        annotate::ColumnCompressed<> annotation(5);
        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label2" });
        annotation.set_labels(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good);
    }
    {
        annotate::ColumnCompressed<> annotation(5);
        ASSERT_FALSE(annotation.load(test_dump_basename_vec_bad));
        ASSERT_TRUE(annotation.load(test_dump_basename_vec_good));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(ColumnCompressed, MergeLoad) {
    {
        annotate::ColumnCompressed<> annotation(5);
        annotation.set_labels(0, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label2" });
        annotation.set_labels(4, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_1");
    }
    {
        annotate::ColumnCompressed<> annotation(5);
        annotation.set_labels(1, { "Label0", "Label2", "Label8" });
        annotation.set_labels(2, { "Label1", "Label9", "Label0" });
        annotation.set_labels(3, { "Label8" });

        annotation.serialize(test_dump_basename_vec_good + "_2");
    }
    {
        annotate::ColumnCompressed<> annotation(0);
        ASSERT_TRUE(annotation.merge_load({ test_dump_basename_vec_good + "_1",
                                            test_dump_basename_vec_good + "_2" }));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label9", "Label0" }), convert_to_set(annotation.get(2)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(3)));
        EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
    }
}

TEST(ColumnCompressed, has_labels) {
    std::unique_ptr<annotate::MultiLabelAnnotation<uint64_t, std::string>> annotation(
        new annotate::ColumnCompressed<>(5)
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

TEST(ColumnCompressed, add_label_sequential) {
    size_t graph_half_size = 1000;
    annotate::ColumnCompressed<> annotation(graph_half_size * 2);
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

TEST(ColumnCompressed, add_label_random) {
    size_t graph_half_size = 1000;
    annotate::ColumnCompressed<> annotation(graph_half_size * 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, labels[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TEST(ColumnCompressed, add_label_random_with_caching) {
    size_t graph_half_size = 1000;
    annotate::ColumnCompressed<> annotation(graph_half_size * 2, 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, labels[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get_labels(i).size());
    }
}

TEST(ColumnCompressed, get_top_labels) {
    annotate::ColumnCompressed<> annotation(5);

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

TEST(ColumnCompressed, get_labels) {
    annotate::ColumnCompressed<> annotation(5);

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

TEST(ColumnCompressed, ToRowAnnotator) {
    {
        annotate::ColumnCompressed<> annotation(0);

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels(0, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels(0, {"Label0", "Label2", "Label8"});
        annotation.add_labels(2, {"Label1", "Label2"});
        annotation.add_labels(3, {});
        annotation.add_labels(4, {"Label1", "Label2", "Label8"});
        annotation.add_labels(5, {"Label2"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorParallel) {
    {
        annotate::ColumnCompressed<> annotation(0);

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(1);
        annotation.add_labels(0, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

        for (size_t i = 0; i < 1; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    {
        annotate::ColumnCompressed<> annotation(6);
        annotation.add_labels(0, {"Label0", "Label2", "Label8"});
        annotation.add_labels(2, {"Label1", "Label2"});
        annotation.add_labels(3, {});
        annotation.add_labels(4, {"Label1", "Label2", "Label8"});
        annotation.add_labels(5, {"Label2"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

        for (size_t i = 0; i < 6; ++i) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, NoRenameColumns) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.set_labels(0, { "Label0", "Label2", "Label8" });
    annotation.set_labels(2, { "Label1", "Label2" });
    annotation.set_labels(4, { "Label8" });

    annotation.rename_columns({});

    EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, RenameColumns) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.set_labels(0, { "Label0", "Label2", "Label8" });
    annotation.set_labels(2, { "Label1", "Label2" });
    annotation.set_labels(4, { "Label8" });

    annotation.rename_columns({ { "Label2", "Label2Renamed" },
                                { "Label8", "Label8Renamed" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Label2Renamed", "Label8Renamed" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label2Renamed" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Label8Renamed" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, SwapColumns) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.set_labels(0, { "Label0", "Label2", "Label8" });
    annotation.set_labels(2, { "Label1", "Label2" });
    annotation.set_labels(4, { "Label8" });

    annotation.rename_columns({ { "Label2", "Label8" },
                                { "Label8", "Label2" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Label8", "Label2" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Label8" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Label2" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, RenameColumnsMerge) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.set_labels(0, { "Label0", "Label2", "Label8" });
    annotation.set_labels(2, { "Label1", "Label2" });
    annotation.set_labels(4, { "Label8" });

    annotation.rename_columns({ { "Label2", "Merged" },
                                { "Label8", "Merged" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, RenameColumnsMergeAll) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.set_labels(0, { "Label0", "Label2", "Label8" });
    annotation.set_labels(2, { "Label1", "Label2" });
    annotation.set_labels(4, { "Label8" });

    annotation.rename_columns({ { "Label0", "Merged" },
                                { "Label1", "Merged" },
                                { "Label2", "Merged" },
                                { "Label8", "Merged" },
                            });

    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, DumpColumns) {
    annotate::ColumnCompressed<> annotation(6);

    annotation.set_labels(1, { "Label0", "Label1" });
    annotation.set_labels(2, { "Label0", "Label1" });
    annotation.set_labels(3, { "Label1" });
    annotation.set_labels(4, { "Label1" });
    annotation.set_labels(5, { "Label0", "Label1" });

    annotation.dump_columns(test_dump_basename);

    std::vector<std::vector<uint64_t>> matrix {
        { 3, 1, 2, 5 },
        { 5, 1, 2, 3, 4, 5 }
    };

    for (uint64_t i = 0; i < 2; ++i) {
        std::ifstream in(test_dump_basename + "." + std::to_string(i)
                                            + ".raw.column.annodbg",
                         std::ios::binary);
        ASSERT_TRUE(in.good());

        std::vector<uint64_t> vect;
        vect.push_back(load_number(in));
        for (uint64_t j = 0; j < vect[0]; ++j) {
            vect.push_back(load_number(in));
        }

        EXPECT_EQ(matrix[i], vect);
    }
}
