#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <numeric>
#include <tuple>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <json/json.h>

#include "../test_helpers.hpp"
#include "../graph/all/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/annotated_dbg.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;
using namespace mtg::graph;
using mtg::annot::matrix::CSRMatrix;

// The C5 order and the cut of a search's top N (SPEC-multigraph-envelope.md section 3)

TEST(SearchCut, LowerMedian) {
    // odd n: the middle element
    EXPECT_EQ(7u, lower_median({ 100, 7, 3 }));
    EXPECT_EQ(5u, lower_median({ 5 }));
    // even n: the lower of the two middle elements, index floor((n - 1) / 2)
    EXPECT_EQ(3u, lower_median({ 100, 7, 3, 1 }));
    EXPECT_EQ(1u, lower_median({ 1, 100 }));
    // zeros excluded: positions where the label is absent are not abundances
    EXPECT_EQ(100u, lower_median({ 0, 0, 0, 100 }));
    EXPECT_EQ(1u, lower_median({ 0, 1, 0, 100, 0 }));
    EXPECT_EQ(7u, lower_median({ 0, 100, 100, 3, 0, 0, 7, 7, 7, 7, 7 }));
    // repeated values, one per position
    EXPECT_EQ(2u, lower_median({ 2, 2, 2, 9 }));
    // no value left
    EXPECT_EQ(0u, lower_median({}));
    EXPECT_EQ(0u, lower_median({ 0, 0 }));
}

TEST(SearchCut, RowOrder) {
    // score descending first
    EXPECT_TRUE(search_row_before({ 5, 0, "z" }, { 4, 100, "a" }));
    EXPECT_FALSE(search_row_before({ 4, 100, "a" }, { 5, 0, "z" }));
    // then the median descending
    EXPECT_TRUE(search_row_before({ 5, 100, "z" }, { 5, 1, "a" }));
    EXPECT_FALSE(search_row_before({ 5, 1, "a" }, { 5, 100, "z" }));
    // then the label ascending, by its bytes
    EXPECT_TRUE(search_row_before({ 5, 7, "a" }, { 5, 7, "b" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "Z" }, { 5, 7, "a" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "a" }, { 5, 7, "ab" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "NC_1" }, { 5, 7, "NZ_0" }));
    // a strict weak order: no row before itself
    EXPECT_FALSE(search_row_before({ 5, 7, "a" }, { 5, 7, "a" }));
}

// search_top_n over |labels| with |scores| (and |medians|, if any), the labels decoded counted
struct TopN {
    std::vector<std::pair<size_t, std::string>> rows;
    SearchCut cut;
    std::vector<size_t> decoded;
};

TopN top_n(const std::vector<std::string> &labels,
           const std::vector<size_t> &scores,
           size_t n,
           const std::vector<uint64_t> &medians = {}) {
    TopN out;
    std::function<std::vector<uint64_t>(const std::vector<size_t> &)> medians_of;
    if (medians.size()) {
        medians_of = [&](const std::vector<size_t> &indexes) {
            std::vector<uint64_t> values;
            for (size_t i : indexes) {
                values.push_back(medians.at(i));
            }
            return values;
        };
    }
    out.rows = search_top_n(scores, n, [&](size_t i) -> const std::string & {
        out.decoded.push_back(i);
        return labels.at(i);
    }, medians_of, &out.cut);
    // (a label is looked up whenever compared: the candidates decoded, each once)
    std::sort(out.decoded.begin(), out.decoded.end());
    out.decoded.erase(std::unique(out.decoded.begin(), out.decoded.end()), out.decoded.end());
    return out;
}

std::vector<std::string> labels_of(const TopN &top) {
    std::vector<std::string> labels;
    for (const auto &[i, label] : top.rows) {
        labels.push_back(label);
    }
    return labels;
}

TEST(SearchCut, Complete) {
    // fewer matches than N: every one returned, in the C5 order, no N-th row
    auto top = top_n({ "c", "a", "b" }, { 3, 5, 3 }, 10);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_EQ(0u, top.cut.tied_returned);
    EXPECT_FALSE(top.cut.abundances);

    // exactly N matches: nothing is cut
    top = top_n({ "c", "a", "b" }, { 3, 5, 3 }, 3);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c" }), labels_of(top));
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
}

TEST(SearchCut, Empty) {
    auto top = top_n({}, {}, 10);
    EXPECT_TRUE(top.rows.empty());
    EXPECT_EQ(0u, top.cut.returned);
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.nth_score);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_EQ(0u, top.cut.tied_returned);
    EXPECT_TRUE(top.decoded.empty());
}

TEST(SearchCut, TiesStraddleTheCut) {
    // scores 9, then four labels at 5, then 1: N = 3 cuts the tie group at 5, which
    // competes for two places by label
    std::vector<std::string> labels { "e", "low", "d", "top", "b", "c" };
    std::vector<size_t> scores {       5,   1,    5,    9,    5,   5  };
    auto top = top_n(labels, scores, 3);
    EXPECT_EQ((std::vector<std::string>{ "top", "b", "c" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(5u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(2u, top.cut.tied_returned);
    EXPECT_FALSE(top.cut.abundances);
    // the rows below the N-th score are never decoded
    EXPECT_EQ((std::vector<size_t>{ 0, 2, 3, 4, 5 }), top.decoded);

    // the returned rows point at their candidates
    for (const auto &[i, label] : top.rows) {
        EXPECT_EQ(labels[i], label);
    }

    // a cut between two scores: the N-th score has one label, returned
    top = top_n(labels, scores, 5);
    EXPECT_EQ((std::vector<std::string>{ "top", "b", "c", "d", "e" }), labels_of(top));
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(5u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(4u, top.cut.tied_returned);
    EXPECT_EQ((std::vector<size_t>{ 0, 2, 3, 4, 5 }), top.decoded);

    // N = 1: the single top row, no tie
    top = top_n(labels, scores, 1);
    EXPECT_EQ((std::vector<std::string>{ "top" }), labels_of(top));
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(9u, top.cut.nth_score);
    EXPECT_EQ(1u, top.cut.tied);
    EXPECT_EQ(1u, top.cut.tied_returned);
    EXPECT_EQ((std::vector<size_t>{ 3 }), top.decoded);
}

TEST(SearchCut, MediansOrderTheTies) {
    // a counts annotation: the tie group at the cut is ordered by median, then by label
    std::vector<std::string> labels   { "a", "b", "c", "d", "top" };
    std::vector<size_t> scores        {  12,  12,  12,  12,  20   };
    std::vector<uint64_t> medians     {   1,   7,   7, 100,  1    };
    auto top = top_n(labels, scores, 3, medians);
    EXPECT_EQ((std::vector<std::string>{ "top", "d", "b" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(12u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(2u, top.cut.tied_returned);
    EXPECT_TRUE(top.cut.abundances);
    EXPECT_TRUE(top.cut.nth_median_read);
    EXPECT_EQ(7u, top.cut.nth_median);
    EXPECT_EQ(2u, top.cut.tied_same_median);
    // the tie member of a lower median than the N-th row's is not decoded
    EXPECT_EQ((std::vector<size_t>{ 1, 2, 3, 4 }), top.decoded);

    // N = 1 of two labels of equal k-mer counts: the higher median wins over the label
    top = top_n({ "SRR0999", "SRR1001" }, { 12, 12 }, 1, { 1, 7 });
    EXPECT_EQ((std::vector<std::string>{ "SRR1001" }), labels_of(top));
    EXPECT_EQ(12u, top.cut.nth_score);
    EXPECT_EQ(2u, top.cut.tied);
    EXPECT_EQ(1u, top.cut.tied_returned);
    EXPECT_EQ(7u, top.cut.nth_median);
    EXPECT_EQ(1u, top.cut.tied_same_median);

    // returned rows of one score are ordered by median also without a cut
    top = top_n(labels, scores, 10, medians);
    EXPECT_EQ((std::vector<std::string>{ "top", "d", "b", "c", "a" }), labels_of(top));
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.nth_median);
    EXPECT_EQ(0u, top.cut.tied_same_median);
}

// Sample 14 of SPEC-multigraph-envelope.md (round 2 of the contract notice): its row's
// kmer_abundances in the server's default format (the runs, zeros omitted), its kmer_count the
// number of query positions the runs cover (a non-zero abundance each), its cut's nth_median
// the lower median of those abundances, and the cut itself what search_top_n states for the
// sample's N (its request's top_labels) over SRR1001 and the narrated SRR0999 (12 k-mers,
// median 1)
TEST(SearchCut, Sample14IsTheEnginesCut) {
    const std::string path = "../tests/data/multigraph_contract/14_search_counts_cut.json";
    std::ifstream in(path, std::ios::binary);
    ASSERT_TRUE(in.good()) << path << " not found (run from the build directory)";
    Json::Value sample;
    Json::CharReaderBuilder builder;
    std::string errors;
    ASSERT_TRUE(Json::parseFromStream(builder, in, &sample, &errors)) << errors;
    const Json::Value &request = sample["request"];
    EXPECT_FALSE(request.isMember("num_labels"));
    const size_t n = request["top_labels"].asUInt64();
    const Json::Value &query = sample["body"][0];
    ASSERT_EQ(1u, query["results"].size());
    const Json::Value &row = query["results"][0];
    // kmer_abundances as the server writes it without --verbose (SeqSearchResult::to_json):
    // the runs of one abundance, "pos=count" or "start-end=count", the zero runs omitted
    std::vector<uint64_t> abundances;
    std::stringstream runs(row["kmer_abundances"].asString());
    size_t next = 0;
    for (std::string run; std::getline(runs, run, ':'); ) {
        const size_t eq = run.find('=');
        ASSERT_NE(std::string::npos, eq) << run;
        const std::string span = run.substr(0, eq);
        const size_t dash = span.find('-');
        const size_t first = std::stoull(span.substr(0, dash));
        const size_t last = dash == std::string::npos ? first : std::stoull(span.substr(dash + 1));
        const uint64_t value = std::stoull(run.substr(eq + 1));
        EXPECT_LE(next, first) << run;
        EXPECT_LE(first, last) << run;
        EXPECT_GT(value, 0u) << run;
        abundances.insert(abundances.end(), last - first + 1, value);
        next = last + 1;
    }
    // the positions the runs cover: those where the label is present
    EXPECT_EQ(abundances.size(), row["kmer_count"].asUInt64());
    EXPECT_EQ(12u, abundances.size());
    const uint64_t median = lower_median(abundances);
    const auto top = top_n({ "SRR0999", row["sample"].asString() },
                           { 12, row["kmer_count"].asUInt64() }, n, { 1, median });
    EXPECT_EQ((std::vector<std::string>{ row["sample"].asString() }), labels_of(top));
    ASSERT_EQ(1u, query["cuts"].size());
    const Json::Value &cut = query["cuts"][0];
    EXPECT_EQ(top.cut.returned, cut["returned"].asUInt64());
    EXPECT_EQ(top.cut.complete, cut["complete"].asBool());
    EXPECT_EQ(top.cut.nth_score, cut["nth_score"].asUInt64());
    EXPECT_EQ(top.cut.tied, cut["tied"].asUInt64());
    EXPECT_EQ(top.cut.tied_returned, cut["tied_returned"].asUInt64());
    EXPECT_EQ(top.cut.nth_median, cut["nth_median"].asUInt64());
    EXPECT_EQ(top.cut.tied_same_median, cut["tied_same_median"].asUInt64());
    EXPECT_EQ(7u, median);
}

TEST(SearchCut, NoneRequested) {
    // N = 0: nothing returned, the matches cut, no N-th row
    auto top = top_n({ "a", "b" }, { 3, 2 }, 0);
    EXPECT_TRUE(top.rows.empty());
    EXPECT_EQ(0u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_TRUE(top.decoded.empty());
}


// search_top_n against the definition: every candidate sorted in the C5 order (score
// descending, median descending, label ascending, then its index), the first N taken, the
// cut counted from that order; random scores with many ties, medians or none, and labels
// held by two candidates
TEST(SearchCut, SelectionIsTheC5Order) {
    std::mt19937 rng(11);
    for (size_t trial = 0; trial < 2000; ++trial) {
        const size_t num = rng() % 40;
        const size_t n = rng() % 45;
        const bool counts = trial % 2;
        std::vector<size_t> scores(num);
        std::vector<uint64_t> medians(num);
        std::vector<std::string> labels(num);
        for (size_t i = 0; i < num; ++i) {
            scores[i] = rng() % 5;
            medians[i] = counts ? 1 + rng() % 3 : 0;
            labels[i] = "L" + std::to_string(rng() % (num + 1));
        }
        auto top = top_n(labels, scores, n, counts ? medians : std::vector<uint64_t>());

        std::vector<size_t> order(num);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return std::make_tuple(-(int64_t)scores[i], -(int64_t)medians[i], labels[i], i)
                    < std::make_tuple(-(int64_t)scores[j], -(int64_t)medians[j], labels[j], j);
        });
        const size_t returned = std::min(n, num);
        ASSERT_EQ(returned, top.rows.size()) << trial;
        for (size_t r = 0; r < returned; ++r) {
            ASSERT_EQ(order[r], top.rows[r].first) << trial << " " << r;
            ASSERT_EQ(labels[order[r]], top.rows[r].second) << trial;
        }
        ASSERT_EQ(returned, top.cut.returned);
        ASSERT_EQ(num <= n, top.cut.complete) << trial;
        if (num) {
            ASSERT_EQ(counts, top.cut.abundances);
        }
        if (top.cut.complete || !n) {
            ASSERT_EQ(0u, top.cut.tied) << trial;
            ASSERT_FALSE(top.cut.nth_median_read) << trial;
            continue;
        }
        const size_t nth = order[n - 1];
        size_t tied = 0, tied_returned = 0, same_median = 0;
        for (size_t r = 0; r < num; ++r) {
            if (scores[order[r]] != scores[nth])
                continue;
            ++tied;
            tied_returned += r < n;
            same_median += medians[order[r]] == medians[nth];
        }
        ASSERT_EQ(scores[nth], top.cut.nth_score) << trial;
        ASSERT_EQ(tied, top.cut.tied) << trial;
        ASSERT_EQ(tied_returned, top.cut.tied_returned) << trial;
        if (!counts) {
            ASSERT_FALSE(top.cut.nth_median_read);
            ASSERT_EQ(0u, top.cut.tied_same_median);
        } else if (tied > 1) {
            ASSERT_TRUE(top.cut.nth_median_read) << trial;
            ASSERT_EQ(medians[nth], top.cut.nth_median) << trial;
            ASSERT_EQ(same_median, top.cut.tied_same_median) << trial;
        } else {
            ASSERT_FALSE(top.cut.nth_median_read) << trial;
            ASSERT_EQ(1u, top.cut.tied_same_median) << trial;
        }
    }
}

// R10-SEARCH-COST-3: a tie group at the cut takes its places by a selection, O(T + N log N),
// not by sorting the whole group (twice): a million labels of one score cost a small multiple
// of a million labels of distinct scores, of which only N are decoded. Measured as a ratio in
// one process, the best of three runs each, so that a loaded machine slows both alike.
TEST(SearchCut, LargeTieGroupIsSelectedNotSorted) {
    const size_t num = 1'000'000;
    const size_t n = 10;
    std::vector<size_t> tied(num, 5);
    std::vector<size_t> distinct(num);
    std::vector<std::string> labels(num);
    for (size_t i = 0; i < num; ++i) {
        const size_t permuted = (i * 2654435761ull) % num;
        distinct[i] = permuted + 1;
        labels[i] = "SRR" + std::to_string(1'000'000'000ull + permuted);
    }
    auto best_of_three = [&](const std::vector<size_t> &scores) {
        double best = std::numeric_limits<double>::max();
        for (int run = 0; run < 3; ++run) {
            SearchCut cut;
            auto start = std::chrono::steady_clock::now();
            auto rows = search_top_n(scores, n,
                                     [&](size_t i) -> const std::string & { return labels[i]; },
                                     {}, &cut);
            best = std::min(best, std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count());
            EXPECT_EQ(n, rows.size());
        }
        return best;
    };
    const double distinct_s = best_of_three(distinct);
    const double tied_s = best_of_three(tied);
    // a full sort of the group took 35 to 80 times the distinct case, the selection 3 to 7
    EXPECT_LT(tied_s, 20 * distinct_s) << "tied " << tied_s << " s, distinct " << distinct_s << " s";

    // the tie group's first labels in byte order
    SearchCut cut;
    auto rows = search_top_n(tied, n,
                             [&](size_t i) -> const std::string & { return labels[i]; }, {},
                             &cut);
    std::vector<std::string> sorted_labels = labels;
    std::partial_sort(sorted_labels.begin(), sorted_labels.begin() + n, sorted_labels.end());
    for (size_t r = 0; r < n; ++r) {
        EXPECT_EQ(sorted_labels[r], rows[r].second);
    }
    EXPECT_EQ(num, cut.tied);
    EXPECT_EQ(n, cut.tied_returned);
}

// GPT B2: a tie group of many rows at the cut costs about a partial selection in the C5 order
// over the candidates' indexes (the labels compared by reference, never copied but for the rows
// returned), which takes the same rows; the selection of release 10 before this fix built,
// decoded and copied every label of the group (8 times this reference)
TEST(SearchCut, ATieGroupCostsAboutAPartialSelection) {
    const size_t num = 300'000;
    const size_t n = 1000;
    std::vector<std::string> labels(num);
    for (size_t i = 0; i < num; ++i) {
        labels[i] = "sample-common-prefix-" + std::to_string((i * 999'983ull) % num);
    }
    const std::vector<size_t> scores(num, 100);
    using Clock = std::chrono::steady_clock;
    double head = std::numeric_limits<double>::max();
    double reference = std::numeric_limits<double>::max();
    for (int run = 0; run < 5; ++run) {
        auto start = Clock::now();
        SearchCut cut;
        auto rows = search_top_n(scores, n,
                                 [&](size_t i) -> const std::string & { return labels[i]; },
                                 {}, &cut);
        head = std::min(head, std::chrono::duration<double>(Clock::now() - start).count());
        start = Clock::now();
        std::vector<size_t> order(num);
        std::iota(order.begin(), order.end(), 0);
        auto c5 = [&](size_t a, size_t b) {
            if (scores[a] != scores[b])
                return scores[a] > scores[b];
            if (labels[a] != labels[b])
                return labels[a] < labels[b];
            return a < b;
        };
        std::nth_element(order.begin(), order.begin() + n - 1, order.end(), c5);
        order.resize(n);
        std::sort(order.begin(), order.end(), c5);
        reference = std::min(reference,
                             std::chrono::duration<double>(Clock::now() - start).count());
        ASSERT_EQ(n, rows.size());
        for (size_t r = 0; r < n; ++r) {
            ASSERT_EQ(order[r], rows[r].first) << r;
            ASSERT_EQ(labels[order[r]], rows[r].second) << r;
        }
        EXPECT_EQ(num, cut.tied);
        EXPECT_EQ(n, cut.tied_returned);
    }
    EXPECT_LT(head, 3 * reference) << "head " << head << " s, reference " << reference << " s";
}

// R10-SEARCH-COST-2: the medians are read only for the rows a median is compared for (the tie
// group at a cut, whose N-th row's median the cut states, and rows above the cut that share a
// score), never for a row with a score of its own
TEST(SearchCut, MediansReadOnlyWhereCompared) {
    auto asked_for = [](const std::vector<size_t> &scores, size_t n) {
        std::vector<size_t> asked;
        std::vector<std::string> labels;
        for (size_t i = 0; i < scores.size(); ++i) {
            labels.push_back("L" + std::to_string(i));
        }
        search_top_n(scores, n,
                     [&](size_t i) -> const std::string & { return labels.at(i); },
                     [&](const std::vector<size_t> &indexes) {
                         asked.insert(asked.end(), indexes.begin(), indexes.end());
                         return std::vector<uint64_t>(indexes.size(), 1);
                     }, nullptr);
        std::sort(asked.begin(), asked.end());
        return asked;
    };
    // complete, every score its own: no median read at all
    EXPECT_EQ(std::vector<size_t>(), asked_for({ 9, 5, 3, 1 }, 10));
    // complete: the two rows of score 5 only
    EXPECT_EQ((std::vector<size_t>{ 1, 2 }), asked_for({ 9, 5, 5, 3 }, 10));
    // a cut at score 5 (three rows for one place): the tie group and the two rows of score 9,
    // not the row of score 7 nor the row below the cut
    EXPECT_EQ((std::vector<size_t>{ 0, 1, 3, 4, 5 }),
              asked_for({ 9, 9, 7, 5, 5, 5, 1 }, 4));
    // a cut whose N-th score has one row: no median orders anything, none is read (the cut
    // states nth_median null)
    EXPECT_EQ(std::vector<size_t>(), asked_for({ 9, 7, 5 }, 2));
    // ... but the rows above it that share a score are ordered by median
    EXPECT_EQ((std::vector<size_t>{ 0, 1 }), asked_for({ 9, 9, 7, 5 }, 3));
    // N = 0: nothing returned, nothing read
    EXPECT_EQ(std::vector<size_t>(), asked_for({ 9, 7, 5 }, 0));
}

TEST(SearchCut, MedianCounterIsTheLowerMedian) {
    std::mt19937 rng(7);
    for (size_t trial = 0; trial < 500; ++trial) {
        const size_t n = rng() % 200;
        const uint64_t max_value = 1 + rng() % (trial % 2 ? 4 : 1000);
        std::vector<uint64_t> values;
        LowerMedianCounter counter;
        for (size_t i = 0; i < n; ++i) {
            // runs of one value, as along a query, and zeros (absent)
            const uint64_t value = rng() % (max_value + 1);
            const size_t run = 1 + rng() % 3;
            const size_t positions = 1 + rng() % 2;
            for (size_t r = 0; r < run; ++r) {
                counter.add(value, positions);
                values.insert(values.end(), positions, value);
            }
        }
        ASSERT_EQ(lower_median(values), counter.lower_median()) << trial;
    }
    LowerMedianCounter none;
    none.add(0, 5);
    EXPECT_EQ(0u, none.lower_median());
}

// An integer matrix that records the rows each read asks for
class RecordingIntMatrix : public annot::matrix::IntMatrix {
  public:
    explicit RecordingIntMatrix(std::vector<RowValues> rows) : rows_(std::move(rows)) {}

    std::vector<RowValues> get_row_values(const std::vector<annot::matrix::BinaryMatrix::Row> &rows,
                                          size_t = 1) const override {
        reads.push_back(rows.size());
        rows_read.insert(rows_read.end(), rows.begin(), rows.end());
        std::vector<RowValues> values;
        for (auto row : rows) {
            values.push_back(rows_.at(row));
        }
        return values;
    }

    const annot::matrix::BinaryMatrix& get_binary_matrix() const override {
        throw std::logic_error("not read");
    }

    mutable std::vector<size_t> reads;
    mutable std::vector<annot::matrix::BinaryMatrix::Row> rows_read;

  private:
    std::vector<RowValues> rows_;
};

TEST(SearchCut, AbundanceMediansReadInBatches) {
    // row r: abundance r + 1 in column 0, 10 in column 1
    std::vector<RecordingIntMatrix::RowValues> rows;
    for (uint64_t r = 0; r < 5; ++r) {
        rows.push_back({ { 0, r + 1 }, { 1, 10 } });
    }
    RecordingIntMatrix matrix(rows);
    // row 2 stands for three query positions: column 0 holds 1, 2, 3, 3, 3, 4, 5. The rows
    // in the query's order, not the matrix's
    std::vector<std::pair<uint64_t, size_t>> row_counts { { 3, 1 }, { 0, 1 }, { 4, 1 },
                                                          { 2, 3 }, { 1, 1 } };
    EXPECT_EQ((std::vector<uint64_t>{ 10, 3, 0 }),
              abundance_medians(matrix, row_counts, { 1, 0, 2 }, 3, 2));
    // read two rows at a time, never all of the query's rows at once
    EXPECT_EQ((std::vector<size_t>{ 2, 2, 1 }), matrix.reads);
    // in the matrix's order (its rank structures read in order, as a batched query's slice
    // reads them): a matrix other than a row-diff one
    EXPECT_EQ((std::vector<annot::matrix::BinaryMatrix::Row>{ 0, 1, 2, 3, 4 }),
              matrix.rows_read);
}

// GPT B1: a k-mer repeated at 10000 query positions costs one (abundance, positions) entry per
// label, never 10000 values. A counter holds its distinct values (a run along the query is one
// entry, merged whenever the entries doubled), and a query's rows are read collapsed, each with
// the number of positions it stands for
TEST(SearchCut, ARepeatedKmerIsOneEntryPerLabel) {
    LowerMedianCounter once;
    once.add(5, 10'000);
    EXPECT_EQ(1u, once.size());
    once.add(9, 1);
    EXPECT_EQ(2u, once.size());
    EXPECT_EQ(5u, once.lower_median());
    // one value per position, as the counts mode passes them: the entries stay a few times the
    // distinct values (two here), never one per position
    LowerMedianCounter positions;
    size_t largest = 0;
    for (size_t i = 0; i < 10'000; ++i) {
        positions.add(1 + i % 2);
        largest = std::max(largest, positions.size());
    }
    EXPECT_LE(largest, 2 * 2 + 16u);
    EXPECT_EQ(1u, positions.lower_median());
    // one row standing for 10000 positions is read once, for every label asked
    std::vector<RecordingIntMatrix::RowValues> rows(2);
    std::vector<uint32_t> columns(1000);
    for (uint32_t j = 0; j < columns.size(); ++j) {
        columns[j] = j;
        rows[0].emplace_back(j, 3);
        rows[1].emplace_back(j, 9);
    }
    RecordingIntMatrix matrix(rows);
    EXPECT_EQ(std::vector<uint64_t>(columns.size(), 3),
              abundance_medians(matrix, { { 0, 10'000 }, { 1, 1 } }, columns, columns.size(),
                                2'500));
    EXPECT_EQ((std::vector<annot::matrix::BinaryMatrix::Row>{ 0, 1 }), matrix.rows_read);
}

// An index whose annotation holds k-mer counts: one sequence, every k-mer of it annotated
// with |abundances[j]| in column j
std::unique_ptr<AnnotatedDBG> counts_index(const std::string &sequence, size_t k,
                                           const std::vector<std::string> &labels,
                                           const std::vector<uint64_t> &abundances) {
    auto graph = build_graph_batch<DBGHashOrdered>(k, { sequence });
    Vector<CSRMatrix::RowValues> rows(graph->max_index());
    graph->map_to_nodes(sequence, [&](DeBruijnGraph::node_index node) {
        auto &row = rows[AnnotatedDBG::graph_to_anno_index(node)];
        row.clear();
        for (size_t j = 0; j < labels.size(); ++j) {
            if (abundances[j])
                row.emplace_back(j, abundances[j]);
        }
    });
    annot::LabelEncoder<std::string> encoder;
    for (const auto &label : labels) {
        encoder.insert_and_encode(label);
    }
    return std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annot::IntRowAnnotator>(
            std::make_unique<CSRMatrix>(std::move(rows), labels.size()), encoder));
}

TEST(SearchCut, CountsTieBrokenByMedianAbundance) {
    const size_t k = 5;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAAC";
    const size_t num_kmers = sequence.size() - k + 1;
    // two labels on the same k-mers, abundance 1 and 100: equal k-mer counts. The label
    // order alone would return "a_low"; the median abundance returns "b_high".
    auto index = counts_index(sequence, k, { "a_low", "b_high" }, { 1, 100 });
    ASSERT_TRUE(AnnotatedDBG::holds_kmer_abundances(index->get_annotator()));
    ASSERT_TRUE(index->has_kmer_abundances());

    auto check_cut = [&](const SearchCut &cut) {
        EXPECT_EQ(1u, cut.returned);
        EXPECT_FALSE(cut.complete);
        EXPECT_EQ(num_kmers, cut.nth_score);
        EXPECT_EQ(2u, cut.tied);
        EXPECT_EQ(1u, cut.tied_returned);
        EXPECT_TRUE(cut.abundances);
        EXPECT_TRUE(cut.nth_median_read);
        EXPECT_EQ(100u, cut.nth_median);
        EXPECT_EQ(1u, cut.tied_same_median);
    };

    SearchCut cut;
    auto matches = index->get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ(num_kmers, matches[0].second);
    check_cut(cut);

    auto counts = index->get_kmer_counts(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ("b_high", std::get<0>(counts[0]));
    EXPECT_EQ(num_kmers, std::get<1>(counts[0]));
    EXPECT_EQ(std::vector<size_t>(num_kmers, 100), std::get<2>(counts[0]));
    check_cut(cut);

    auto signatures = index->get_top_label_signatures(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, signatures.size());
    EXPECT_EQ("b_high", std::get<0>(signatures[0]));
    check_cut(cut);

    // without a cut both are returned, the higher median first
    matches = index->get_top_labels(sequence, 2, 0.0, 0.0, false, &cut);
    ASSERT_EQ(2u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ("a_low", matches[1].first);
    EXPECT_TRUE(cut.complete);
    EXPECT_EQ(2u, cut.returned);
    EXPECT_EQ(0u, cut.tied);

    // the abundance sums differ: no tie, no median read
    matches = index->get_top_labels(sequence, 1, 0.0, 0.0, true, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ(100 * num_kmers, matches[0].second);
    EXPECT_EQ(100 * num_kmers, cut.nth_score);
    EXPECT_EQ(1u, cut.tied);
    EXPECT_EQ(1u, cut.tied_same_median);
    EXPECT_TRUE(cut.abundances);
    EXPECT_FALSE(cut.nth_median_read);

    // no k-mer of the query in the graph: nothing returned, complete
    matches = index->get_top_labels("TTTTTTTTTT", 1, 0.0, 0.0, false, &cut);
    EXPECT_TRUE(matches.empty());
    EXPECT_EQ(0u, cut.returned);
    EXPECT_TRUE(cut.complete);
    EXPECT_TRUE(cut.abundances);
}

TEST(SearchCut, BinaryTieOrderedByLabel) {
    // a binary annotation: ties by label (byte order), not by column
    const size_t k = 5;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAAC";
    auto index = build_anno_graph<DBGHashOrdered, annot::ColumnCompressed<>>(
            k, { sequence, sequence, sequence }, { "b", "c", "a" });
    ASSERT_FALSE(index->has_kmer_abundances());

    SearchCut cut;
    auto matches = index->get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a", matches[0].first);
    EXPECT_EQ(1u, cut.returned);
    EXPECT_FALSE(cut.complete);
    EXPECT_EQ(3u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);
    EXPECT_FALSE(cut.abundances);

    matches = index->get_top_labels(sequence, 3, 0.0, 0.0, false, &cut);
    ASSERT_EQ(3u, matches.size());
    EXPECT_EQ("a", matches[0].first);
    EXPECT_EQ("b", matches[1].first);
    EXPECT_EQ("c", matches[2].first);
    EXPECT_TRUE(cut.complete);
}

TEST(SearchCut, CoordinateRecordsTiedByHeader) {
    // a coordinate annotation with its record mapping (.seqs): two records of one column,
    // the same sequence, so every k-mer matches both. The records are ranked by header; the
    // coordinates are positions, not abundances, and order no tie.
    const size_t k = 11;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAACGGATCCATGCAAGT";
    const uint64_t n = sequence.size() - k + 1;
    auto coords = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, { sequence, sequence }, { "F", "F" }, DeBruijnGraph::BASIC, true, { 0, n });
    ASSERT_FALSE(AnnotatedDBG::holds_kmer_abundances(coords->get_annotator()));

    std::filesystem::path dir = utils::create_temp_dir("", "test_search_cut");
    const std::string prefix = dir / "coords";
    coords->get_annotator().serialize(prefix);
    auto annotator = std::make_unique<annot::ColumnCoordAnnotator>();
    ASSERT_TRUE(annotator->load(prefix));
    utils::remove_temp_dir(dir);

    // the first record's header sorts after the second's
    auto mapping = std::make_shared<annot::CoordToHeader>(
            std::vector<std::vector<std::string>>{ { "z_record", "a_record" } },
            std::vector<std::vector<uint64_t>>{ { n, n } });
    auto graph = std::dynamic_pointer_cast<DeBruijnGraph>(
            std::const_pointer_cast<SequenceGraph>(coords->get_graph_ptr()));
    ASSERT_TRUE(graph);
    AnnotatedDBG index(graph, std::move(annotator), false, mapping);
    ASSERT_FALSE(index.has_kmer_abundances());

    auto check_cut = [&](const SearchCut &cut) {
        EXPECT_EQ(1u, cut.returned);
        EXPECT_FALSE(cut.complete);
        EXPECT_EQ(2u, cut.tied);
        EXPECT_EQ(1u, cut.tied_returned);
        EXPECT_FALSE(cut.abundances);
    };

    SearchCut cut;
    auto records = index.get_kmer_coordinates(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, records.size());
    EXPECT_EQ("a_record", std::get<0>(records[0]));
    EXPECT_EQ(n, std::get<1>(records[0]));
    EXPECT_EQ(n, cut.nth_score);
    check_cut(cut);

    auto matches = index.get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a_record", matches[0].first);
    check_cut(cut);

    // the coordinate counts: ranked by the score reported
    matches = index.get_top_labels(sequence, 1, 0.0, 0.0, true, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a_record", matches[0].first);
    EXPECT_EQ(n, matches[0].second);
    check_cut(cut);

    auto counts = index.get_kmer_counts(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ("a_record", std::get<0>(counts[0]));
    check_cut(cut);

    // both records: in header order
    records = index.get_kmer_coordinates(sequence, 2, 0.0, 0.0, &cut);
    ASSERT_EQ(2u, records.size());
    EXPECT_EQ("a_record", std::get<0>(records[0]));
    EXPECT_EQ("z_record", std::get<0>(records[1]));
    EXPECT_TRUE(cut.complete);
    EXPECT_EQ(2u, cut.returned);
}

} // namespace
