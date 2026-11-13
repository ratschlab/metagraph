#include <algorithm>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "cli/query.hpp"
#include "cli/search_answer.hpp"
#include "cli/server_checks.hpp"
#include "graph/search_budget.hpp"

#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__GLIBC__)
#include <malloc.h>
#endif

#if USE_JEMALLOC
extern "C" int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen);
#endif


// The /search answer budget's answer account (SPEC-multigraph-envelope.md section 3, lane
// T5b-Y4): each FASTA record's JSON object is priced before it is built (its tree as to_json and
// the server build it, the per-position strings measured without writing them), the text at
// its exact length before it is written (the bytes Json::writeString writes), and a request
// past the bound is refused 400 search_answer_too_large, stated by GET /capabilities' `search`.

namespace mtg::test {
namespace {

using cli::SeqSearchResult;
using cli::QuerySequence;
using graph::PresenceRuns;
using graph::CountRuns;
using graph::PositionCoords;
using graph::SearchAccount;
using graph::SearchAnswerTooLarge;
using graph::SearchBudget;

constexpr uint64_t kLarge = uint64_t(1) << 40;

std::string compact(const Json::Value &value) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, value);
}

// the allocator's live bytes, where they can be read
std::optional<int64_t> live_heap_bytes() {
#if USE_JEMALLOC
    uint64_t allocated = 0, deallocated = 0;
    size_t sz = sizeof(uint64_t);
    if (!mallctl("thread.allocated", &allocated, &sz, nullptr, 0)
            && !mallctl("thread.deallocated", &deallocated, &sz, nullptr, 0)
            && (allocated || deallocated)) {
        return int64_t(allocated) - int64_t(deallocated);
    }
#endif
#if defined(__APPLE__)
    malloc_statistics_t stats;
    malloc_zone_statistics(nullptr, &stats);
    return int64_t(stats.size_in_use);
#elif defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    return int64_t(mallinfo2().uordblks);
#else
    return std::nullopt;
#endif
}

// random per-position values of |size| positions, about |fill| of them set, in runs
PresenceRuns random_mask(std::mt19937 &rng, size_t size, double fill) {
    PresenceRuns mask(size);
    for (size_t i = 0; i < size; ++i) {
        if (rng() % 1000 < fill * 1000)
            mask.add(i);
    }
    return mask;
}
CountRuns random_counts(std::mt19937 &rng, size_t size, double fill) {
    CountRuns counts(size);
    for (size_t i = 0; i < size; ) {
        const size_t run = 1 + rng() % 5;
        if (rng() % 1000 < fill * 1000)
            counts.add_run(i, std::min(size, i + run), 1 + rng() % 30000);
        i += run;
    }
    return counts;
}
PositionCoords random_coords(std::mt19937 &rng, size_t size, double fill, bool consecutive) {
    PositionCoords coords(size);
    uint64_t next = rng() % 100000000;
    for (size_t i = 0; i < size; ++i) {
        if (rng() % 1000 >= fill * 1000)
            continue;
        const size_t n = 1 + rng() % 3;
        std::vector<uint64_t> list;
        for (size_t j = 0; j < n; ++j) {
            list.push_back(consecutive ? next + j * 1000000 + i : rng() % 1000000000000);
        }
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
        for (uint64_t c : list) {
            coords.add(i, c);
        }
    }
    return coords;
}

std::string random_label(std::mt19937 &rng, size_t length, bool properties) {
    std::string label = "NC_" + std::string(length, 'x');
    for (char &c : label) {
        if (c == 'x')
            c = "ACGTacgt_.|"[rng() % 11];
    }
    if (properties)
        label += ";taxid=" + std::to_string(rng() % 100000) + ";name=abc def;score=nan";
    return label;
}

// one record's result of each mode, |rows| rows over |size| positions
SeqSearchResult record(std::mt19937 &rng, int mode, size_t rows, size_t size, size_t label,
                       bool properties, double fill, bool alignment = false) {
    QuerySequence sequence { 0, "read_" + std::string(label, 'q'), std::string(size + 30, 'A') };
    SeqSearchResult::result_type result;
    switch (mode) {
        case 0: {
            SeqSearchResult::LabelCountVec v;
            for (size_t r = 0; r < rows; ++r) {
                v.emplace_back(random_label(rng, label, properties), rng() % 1000);
            }
            result = std::move(v);
            break;
        }
        case 1: {
            SeqSearchResult::LabelSigVec v;
            for (size_t r = 0; r < rows; ++r) {
                PresenceRuns mask = random_mask(rng, size, fill);
                v.emplace_back(random_label(rng, label, properties), mask.count(), std::move(mask));
            }
            result = std::move(v);
            break;
        }
        case 2: {
            SeqSearchResult::LabelCountAbundancesVec v;
            for (size_t r = 0; r < rows; ++r) {
                CountRuns counts = random_counts(rng, size, fill);
                if (!counts.count())
                    counts.add(0, 1);
                v.emplace_back(random_label(rng, label, properties), counts.count(),
                               std::move(counts));
            }
            result = std::move(v);
            break;
        }
        default: {
            SeqSearchResult::LabelCountCoordsVec v;
            for (size_t r = 0; r < rows; ++r) {
                PositionCoords coords = random_coords(rng, size, fill, r % 2);
                v.emplace_back(random_label(rng, label, properties), coords.count(),
                               std::move(coords), r % 3 ? 5000 + r : 0);
            }
            result = std::move(v);
            break;
        }
    }
    graph::SearchCut cut;
    cut.returned = rows;
    cut.complete = false;
    cut.nth_score = 3;
    cut.tied = 2;
    cut.abundances = mode == 2;
    std::optional<cli::Alignment> aligned;
    if (alignment)
        aligned = cli::Alignment { 10, 20, std::string(40, 'M'), true };
    return SeqSearchResult(std::move(sequence), std::move(result), std::move(aligned), cut);
}

// The record's object as the server builds it (process_search_request): to_json, the pair's
// tags on every row and on its cut, the cut appended
Json::Value server_object(const SeqSearchResult &result, bool verbose, bool tagged,
                          const std::string &graph, const std::string &graph_path,
                          const std::string &annotation_path, const std::string &index_fp) {
    Json::Value json = result.to_json(verbose, 31);
    if (tagged) {
        for (Json::Value &row : json["results"]) {
            row["graph"] = graph;
            row["graph_path"] = graph_path;
            row["annotation_path"] = annotation_path;
        }
    }
    Json::Value cut = result.cut_to_json();
    cut["graph"] = tagged ? Json::Value(graph) : Json::Value();
    cut["graph_path"] = tagged ? Json::Value(graph_path) : Json::Value();
    cut["annotation_path"] = tagged ? Json::Value(annotation_path) : Json::Value();
    cut["index_fp"] = index_fp.empty() ? Json::Value() : Json::Value(index_fp);
    json["cuts"] = Json::Value(Json::arrayValue);
    json["cuts"].append(std::move(cut));
    return json;
}

// The per-position strings measured without writing them are the strings to_json writes
TEST(SearchAnswer, TheStringsAreMeasuredAsWritten) {
    std::mt19937 rng(31);
    for (size_t round = 0; round < 60; ++round) {
        const size_t size = 1 + rng() % (round < 30 ? 50 : 3000);
        const double fill = (rng() % 100) / 100.0;
        for (bool verbose : { false, true }) {
            SCOPED_TRACE(::testing::Message() << round << " size " << size << " verbose "
                         << verbose);
            const PresenceRuns mask = random_mask(rng, size, fill);
            const CountRuns counts = random_counts(rng, size, fill);
            const PositionCoords coords = random_coords(rng, size, fill, round % 2);
            SeqSearchResult::LabelSigVec sig;
            sig.emplace_back("a", mask.count(), mask);
            SeqSearchResult::LabelCountAbundancesVec abund;
            if (counts.count())
                abund.emplace_back("a", counts.count(), counts);
            SeqSearchResult::LabelCountCoordsVec coord;
            coord.emplace_back("a", coords.count(), coords, 7);
            const Json::Value s = SeqSearchResult(QuerySequence { 0, "q", "" }, std::move(sig))
                    .to_json(verbose, 31);
            EXPECT_EQ(s["results"][0]["signature"]["presence_mask"].asString().size(),
                      cli::presence_mask_text_length(mask, verbose));
            if (counts.count()) {
                const Json::Value c = SeqSearchResult(QuerySequence { 0, "q", "" },
                                                      std::move(abund)).to_json(verbose, 31);
                EXPECT_EQ(c["results"][0]["kmer_abundances"].asString().size(),
                          cli::count_runs_text_length(counts, verbose));
            }
            const Json::Value x = SeqSearchResult(QuerySequence { 0, "q", "" }, std::move(coord))
                    .to_json(verbose, 31);
            EXPECT_EQ(x["results"][0]["kmer_coords"].asString().size(),
                      cli::coords_text_length(coords, verbose));
            if (!verbose) {
                // and the text joined by ':' and its strings
                const std::vector<std::string> ranges = cli::collapse_coord_ranges(coords);
                const cli::CoordRangesText text = cli::collapsed_coords_text(coords);
                EXPECT_EQ(ranges.size(), text.ranges);
                uint64_t long_bytes = 0;
                for (const std::string &range : ranges) {
                    if (range.size() > 15)
                        long_bytes += range.size() + 1;
                }
                EXPECT_EQ(long_bytes, text.long_ranges_bytes);
            }
        }
    }
}

// The price of a record's object covers the tree the server builds of it (the allocator's live
// bytes where they can be read; the model is the deployed build's, whose nodes are the same
// size or larger than libc++'s), on shapes of many short rows, of long run strings, of labels
// with properties, tagged and untagged, verbose or not
TEST(SearchAnswer, ThePriceCoversTheTree) {
    std::mt19937 rng(37);
    size_t measured = 0;
    const std::string graph = "bact_pseudomonadota";
    const std::string graph_path = "/scratch/nvme8/refseq_r237/bact_pseudomonadota/graph_k31.dbg";
    const std::string annotation_path
            = "/scratch/nvme8/refseq_r237/bact_pseudomonadota/annotation.row_diff_brwt_coord.annodbg";
    const std::string index_fp = "0123456789abcdef";
    for (int mode : { 0, 1, 2, 3 }) {
        for (auto [rows, size, label, fill] : { std::make_tuple(2000, 40, 12, 0.9),
                                                std::make_tuple(20, 20000, 40, 0.5),
                                                std::make_tuple(300, 2000, 100, 0.1),
                                                std::make_tuple(5, 50000, 30, 0.97) }) {
            for (bool tagged : { false, true }) {
                for (bool verbose : { false, true }) {
                    const bool properties = rows == 300;
                    const SeqSearchResult result = record(rng, mode, rows, size, label,
                                                          properties, fill, rows == 5);
                    const std::string what = "mode " + std::to_string(mode) + " rows "
                            + std::to_string(rows) + " size " + std::to_string(size)
                            + (tagged ? " tagged" : "") + (verbose ? " verbose" : "");
                    const cli::SearchRecordPrice price = cli::search_record_price(
                            result, verbose, tagged, graph, graph_path, annotation_path,
                            tagged ? index_fp : "");
                    Json::Value answer(Json::arrayValue);
                    const std::optional<int64_t> before = live_heap_bytes();
                    answer.append(server_object(result, verbose, tagged, graph, graph_path,
                                                annotation_path, tagged ? index_fp : ""));
                    const std::optional<int64_t> after = live_heap_bytes();
                    if (before && after && *after > *before) {
                        ++measured;
                        EXPECT_LE(uint64_t(*after - *before), price.tree)
                                << what << ": priced " << price.tree;
                        // a model, not a gross over-count
                        EXPECT_LE(price.tree, 4 * uint64_t(*after - *before) + (64 << 10))
                                << what << ": priced " << price.tree << ", held "
                                << *after - *before;
                    }
                    // the per-position strings at least once beside the tree
                    if (mode)
                        EXPECT_LT(0u, price.transient) << what;
                }
            }
        }
    }
    if (!measured)
        GTEST_SKIP() << "no allocator counters to compare the prices with";
}

// The text is written, byte for byte, as Json::writeString writes it, into a string of its
// exact length
TEST(SearchAnswer, TheTextIsJsonWriteStrings) {
    std::mt19937 rng(41);
    for (size_t round = 0; round < 20; ++round) {
        Json::Value answer(Json::arrayValue);
        for (size_t r = 0; r < 1 + round % 4; ++r) {
            answer.append(server_object(record(rng, round % 4, 1 + rng() % 30, 1 + rng() % 500,
                                               10, round % 3 == 0, 0.5, round % 5 == 0),
                                        round % 2, round % 2, "g\"\\é", "p\n", "a\t", ""));
        }
        answer[0]["double"] = 0.1 + round;
        answer[0]["control"] = std::string("\x01\x1f", 2);
        const Json::StreamWriterBuilder builder;
        const std::string expected = Json::writeString(builder, answer);
        EXPECT_EQ(expected.size(), cli::json_text_length(answer, builder));
        const std::string got = cli::json_text_of_length(answer, builder, expected.size());
        EXPECT_EQ(expected, got);
        EXPECT_GE(got.capacity(), got.size());

        SearchBudget budget(kLarge);
        EXPECT_EQ(expected, cli::search_answer_text(answer, budget.answer(), false));
        // the text, beside it its copy for the transport (and the head and gzip's state)
        EXPECT_LE(2 * expected.size(), budget.answer().held());
        SearchBudget cached(kLarge);
        cli::search_answer_text(answer, cached.answer(), true);
        EXPECT_LE(3 * expected.size(), cached.answer().held());
        // a byte short: refused before the text is written
        SearchBudget short_one(budget.answer().held() - 1);
        EXPECT_THROW(cli::search_answer_text(answer, short_one.answer(), false),
                     SearchAnswerTooLarge);
    }
}

// The rows a record's result holds are charged to its pair's values account
TEST(SearchAnswer, TheRowsOfARecordArePriced) {
    std::mt19937 rng(43);
    for (int mode : { 0, 1, 2, 3 }) {
        const SeqSearchResult small = record(rng, mode, 1, 40, 10, false, 0.5);
        const SeqSearchResult large = record(rng, mode, 1000, 40, 30, false, 0.5);
        EXPECT_LT(cli::search_record_rows_bytes(small), cli::search_record_rows_bytes(large));
        // at least each row's label beyond the inline bytes, and its entry
        EXPECT_LE(1000u * (32 + 33), cli::search_record_rows_bytes(large)) << mode;
    }
}

// The refusal, written error first, and GET /capabilities' block
TEST(SearchAnswer, TheRefusalAndTheCapabilityBlock) {
    EXPECT_EQ("{\"error\":\"search: the answer would hold more than --search-max-answer-memory-mb "
              "4096 MiB: nothing partial is sent; ask for fewer top_labels, a larger "
              "discovery_fraction, fewer or shorter FASTA records, or no per-position values "
              "(with_signature, query_counts, query_coords)\",\"code\":\"search_answer_too_large\","
              "\"max_memory_mb\":4096}",
              cli::error_first_json_text(cli::search_answer_too_large_json(4096)));
    EXPECT_EQ(std::string(SearchAnswerTooLarge(uint64_t(1) << 20, "answer").what()),
              cli::search_answer_too_large_json(1)["error"].asString());
    EXPECT_EQ("{\"max_memory_mb\":4096,\"rule\":\"SPEC-multigraph-envelope.md section 3, the "
              "/search answer budget\",\"too_large\":{\"code\":\"search_answer_too_large\","
              "\"members\":[\"code\",\"error\",\"max_memory_mb\"],\"status\":400}}",
              compact(cli::search_capabilities_json(4096)));
    const Json::Value none = cli::search_capabilities_json(0);
    EXPECT_TRUE(none["max_memory_mb"].isNull());
    EXPECT_EQ(compact(cli::search_capabilities_json(4096)["too_large"]),
              compact(none["too_large"]));
    // the members a refusal states are the body's
    Json::Value members(Json::arrayValue);
    for (const std::string &name : cli::search_answer_too_large_json(4096).getMemberNames()) {
        members.append(name);
    }
    EXPECT_EQ(compact(members), compact(none["too_large"]["members"]));
}

} // namespace
} // namespace mtg::test
