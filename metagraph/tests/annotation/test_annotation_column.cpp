#include <random>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "annotate_column_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "threading.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_column_compressed";


TEST(ColumnCompressed, add_label_random_with_caching) {
    size_t graph_half_size = 1000;
    annotate::ColumnCompressed<> annotation(graph_half_size * 2, 2);

    std::vector<std::string> labels { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_labels({ i }, { labels[i % 2] });
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TEST(ColumnCompressed, RenameColumnsMerge) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    annotation.rename_labels({ { "Label2", "Merged" },
                               { "Label8", "Merged" } });

    EXPECT_EQ(convert_to_set({ "Label0", "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Label1", "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, RenameColumnsMergeAll) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    annotation.rename_labels({ { "Label0", "Merged" },
                               { "Label1", "Merged" },
                               { "Label2", "Merged" },
                               { "Label8", "Merged" }, });

    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(0)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(1)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(2)));
    EXPECT_EQ(convert_to_set({}), convert_to_set(annotation.get(3)));
    EXPECT_EQ(convert_to_set({ "Merged" }), convert_to_set(annotation.get(4)));
}

TEST(ColumnCompressed, GetColumn) {
    annotate::ColumnCompressed<> annotation(5);
    annotation.add_labels({ 0 }, { "Label0", "Label2", "Label8" });
    annotation.add_labels({ 2 }, { "Label1", "Label2" });
    annotation.add_labels({ 4 }, { "Label8" });

    EXPECT_EQ(1u, annotation.get_column("Label0").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label0")[0]);

    EXPECT_EQ(1u, annotation.get_column("Label1").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label1")[2]);

    EXPECT_EQ(2u, annotation.get_column("Label2").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label2")[0]);
    EXPECT_TRUE(annotation.get_column("Label2")[2]);

    EXPECT_EQ(2u, annotation.get_column("Label8").num_set_bits());
    EXPECT_TRUE(annotation.get_column("Label8")[0]);
    EXPECT_TRUE(annotation.get_column("Label8")[4]);
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
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

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

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

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

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator);

        for (size_t i = 0; i < num_rows; i += 1000) {
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
        annotation.add_labels({ 0 }, {"Label0", "Label2", "Label8"});

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

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

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

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

        annotate::RowCompressed<> row_annotator(0);
        annotation.convert_to_row_annotator(&row_annotator, 10);

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreaming) {
    {
        annotate::ColumnCompressed<> annotation(0);

        annotation.convert_to_row_annotator(test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotation.convert_to_row_annotator(test_dump_basename);

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

        annotation.convert_to_row_annotator(test_dump_basename);

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

        annotation.convert_to_row_annotator(test_dump_basename);

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
        annotation.convert_to_row_annotator(test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
}

TEST(ColumnCompressed, ToRowAnnotatorStreamingParallel) {
    set_num_threads(10);
    {
        annotate::ColumnCompressed<> annotation(0);

        annotation.convert_to_row_annotator(test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        ASSERT_EQ(0u, row_annotator.num_objects());
        ASSERT_EQ(0u, row_annotator.num_labels());
    }
    {
        annotate::ColumnCompressed<> annotation(1);

        annotation.convert_to_row_annotator(test_dump_basename);

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

        annotation.convert_to_row_annotator(test_dump_basename);

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

        annotation.convert_to_row_annotator(test_dump_basename);

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
        annotation.convert_to_row_annotator(test_dump_basename);

        annotate::RowCompressed<> row_annotator(0);
        ASSERT_TRUE(row_annotator.load(test_dump_basename));

        for (size_t i = 0; i < num_rows; i += 1000) {
            EXPECT_EQ(convert_to_set(annotation.get(i)),
                      convert_to_set(row_annotator.get(i)));
        }
    }
    set_num_threads(1);
}
