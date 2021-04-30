#include "gtest/gtest.h"

#include "common/vector_map.hpp"

#include "../test_helpers.hpp"
#include "test_annotation.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

using mtg::annot::MultiLabelEncoded;


std::vector<std::string> get_labels(const MultiLabelEncoded<std::string> &annotator,
                                    const std::vector<uint64_t> &indices,
                                    double min_label_frequency = 0.0) {
    const auto& label_encoder = annotator.get_label_encoder();
    const size_t min_count = std::max(1.0,
                                      std::ceil(min_label_frequency * indices.size()));

    std::vector<std::pair<std::string, size_t>> label_counts;
    for (size_t j = 0; j < label_encoder.size(); ++j) {
        label_counts.emplace_back(label_encoder.decode(j), 0);
    }

    for (auto i : indices) {
        for (auto j : annotator.get_matrix().get_row(i)) {
            label_counts[j].second++;
        }
    }

    std::vector<std::string> labels;
    for (auto&& pair : label_counts) {
        if (pair.second >= min_count)
            labels.emplace_back(std::move(pair.first));
    }

    return labels;
}

TYPED_TEST(AnnotatorPresetTest, call_rows_get_labels) {
    EXPECT_EQ(std::vector<std::string>({}),
              get_labels(*this->annotation, {}, 1));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels(*this->annotation, { 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 1)));
}

std::vector<std::string> get_labels_by_label(const MultiLabelEncoded<std::string> &annotator,
                                             const std::vector<uint64_t> &indices,
                                             double min_label_frequency = 0.0) {
    VectorMap<uint64_t, size_t> index_counts;
    for (auto i : indices) {
        index_counts[i] = 1;
    }

    const size_t min_count = std::max(1.0,
                                      std::ceil(min_label_frequency * indices.size()));

    auto code_counts = annotator.get_matrix().sum_rows(index_counts.values_container(), min_count, min_count);

    std::vector<std::string> labels;
    labels.reserve(code_counts.size());

    const auto &label_encoder = annotator.get_label_encoder();

    for (auto&& pair : code_counts) {
        assert(pair.second >= min_count);
        labels.push_back(label_encoder.decode(pair.first));
    }

    return labels;
}

TYPED_TEST(AnnotatorPresetTest, call_rows_get_labels_by_label) {
    EXPECT_EQ(std::vector<std::string>({}),
              get_labels(*this->annotation, {}, 1));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2, 4 }, 0.5)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2, 4 }, 0.501)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 2, 4 }, 1)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 1)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0)));

    EXPECT_EQ(convert_to_set({"Label0", "Label1", "Label2", "Label8"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.2)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.201)));

    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.4)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.401)));

    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.8)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0.801)));

    EXPECT_EQ(convert_to_set({}),
              convert_to_set(get_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 1)));
}


std::vector<std::pair<std::string, size_t>>
get_top_labels(const MultiLabelEncoded<std::string> &annotator,
               const std::vector<uint64_t> &indices,
               size_t num_top_labels = static_cast<size_t>(-1),
               double min_label_frequency = 0.0) {
    const auto& label_encoder = annotator.get_label_encoder();
    const size_t min_count = std::max(1.0,
                                      std::ceil(min_label_frequency * indices.size()));

    std::vector<std::pair<std::string, size_t>> label_counts;
    for (size_t i = 0; i < label_encoder.size(); ++i) {
        label_counts.emplace_back(label_encoder.decode(i), 0);
    }

    for (auto i : indices) {
        for (auto j : annotator.get_matrix().get_row(i)) {
            label_counts[j].second++;
        }
    }

    std::sort(label_counts.begin(), label_counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    // remove labels which don't meet num_top_labels and min_label_frequency criteria
    label_counts.erase(std::find_if(label_counts.begin(),
                                    num_top_labels < label_counts.size()
                                        ? label_counts.begin() + num_top_labels
                                        : label_counts.end(),
                                    [&](const auto &pair) {
                                        return pair.second < min_count;
                                    }),
                       label_counts.end());

    return label_counts;
}

TYPED_TEST(AnnotatorPreset3Test, call_rows_get_top_labels) {
    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              get_top_labels(*this->annotation, {}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }));

    EXPECT_EQ(convert_to_set(VectorCounts({ std::make_pair("Label1", 1),
                                            std::make_pair("Label2", 1) })),
              convert_to_set(get_top_labels(*this->annotation, { 2 })));

    EXPECT_EQ(VectorCounts({}),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels(*this->annotation, { 0, 1, 2, 3, 4 }, 1000));
}

std::vector<std::pair<std::string, size_t>>
get_top_labels_by_label(const MultiLabelEncoded<std::string> &annotator,
                        const std::vector<uint64_t> &indices,
                        size_t num_top_labels = static_cast<size_t>(-1),
                        double min_label_frequency = 0.0) {
    const size_t min_count = std::max(1.0,
                                      std::ceil(min_label_frequency * indices.size()));

    VectorMap<uint64_t, size_t> index_counts;
    for (auto i : indices) {
        index_counts[i] = 1;
    }

    auto code_counts = annotator.get_matrix().sum_rows(index_counts.values_container(), min_count);

    assert(std::all_of(
        code_counts.begin(), code_counts.end(),
        [&](const auto &code_count) { return code_count.second >= min_count; }
    ));

    std::sort(code_counts.begin(), code_counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    // leave only first |num_top_labels| top labels
    if (code_counts.size() > num_top_labels)
        code_counts.erase(code_counts.begin() + num_top_labels,
                           code_counts.end());

    const auto &label_encoder = annotator.get_label_encoder();

    std::vector<std::pair<std::string, size_t>> label_counts(code_counts.size());

    for (size_t i = 0; i < code_counts.size(); ++i) {
        label_counts[i].first = label_encoder.decode(code_counts[i].first);
        label_counts[i].second = code_counts[i].second;
    }

    return label_counts;
}

TYPED_TEST(AnnotatorPreset3Test, call_rows_get_top_labels_by_label) {
    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    EXPECT_EQ(VectorCounts({}),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({}),
              get_top_labels_by_label(*this->annotation, {}));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }));

    EXPECT_EQ(convert_to_set(VectorCounts({ std::make_pair("Label1", 1),
                                            std::make_pair("Label2", 1) })),
              convert_to_set(get_top_labels_by_label(*this->annotation, { 2 })));

    EXPECT_EQ(VectorCounts({}),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 0));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 1));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 2));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 3));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 4));

    EXPECT_EQ(VectorCounts({ std::make_pair("Label2", 4),
                             std::make_pair("Label8", 3),
                             std::make_pair("Label1", 2),
                             std::make_pair("Label0", 1) }),
              get_top_labels_by_label(*this->annotation, { 0, 1, 2, 3, 4 }, 1000));
}

} // namespace
