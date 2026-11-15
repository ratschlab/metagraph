#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "annotation/binary_matrix/base/decode_budget.hpp"
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
using annot::matrix::alloc_bytes;

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
                    if (mode) {
                        EXPECT_LT(0u, price.transient) << what;
                    }
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

// The refusal, written error first, and GET /capabilities' block (lane T7-P1: the 400's
// `account`, the units refused alone, the concurrency the model total depends on)
TEST(SearchAnswer, TheRefusalAndTheCapabilityBlock) {
    for (const char *account : { "answer", "values" }) {
        EXPECT_EQ("{\"error\":\"search: the answer would hold more than "
                  "--search-max-answer-memory-mb 4096 MiB: nothing partial is sent; ask for "
                  "fewer top_labels, a larger discovery_fraction, fewer or shorter FASTA "
                  "records, or no per-position values (with_signature, query_counts, "
                  "query_coords)\",\"account\":\"" + std::string(account) + "\","
                  "\"code\":\"search_answer_too_large\",\"max_memory_mb\":4096}",
                  cli::error_first_json_text(cli::search_answer_too_large_json(4096, account)));
    }
    EXPECT_EQ(std::string(SearchAnswerTooLarge(uint64_t(1) << 20, "answer").what()),
              cli::search_answer_too_large_json(1, "answer")["error"].asString());
    EXPECT_EQ("{\"handler_threads\":2,\"max_memory_mb\":4096,\"pair_searches_at_once\":2,"
              "\"records_at_once\":2,\"rule\":\"SPEC-multigraph-envelope.md section 3, the "
              "/search answer budget\",\"too_large\":{\"accounts\":[\"answer\",\"values\"],"
              "\"code\":\"search_answer_too_large\",\"members\":[\"account\",\"code\",\"error\","
              "\"max_memory_mb\"],\"status\":400},\"unit_refused\":{\"accounts\":[\"answer\","
              "\"values\"],\"code\":\"search_answer_too_large\",\"field\":\"partial_refusal\","
              "\"members\":[\"account\",\"code\",\"max_memory_mb\",\"priced_mb\"]}}",
              compact(cli::search_capabilities_json(4096, 2, 2, 2)));
    const Json::Value none = cli::search_capabilities_json(0, 1, 3, 5);
    EXPECT_TRUE(none["max_memory_mb"].isNull());
    // the concurrency is stated without a bound too
    EXPECT_EQ(1u, none["pair_searches_at_once"].asUInt64());
    EXPECT_EQ(3u, none["records_at_once"].asUInt64());
    EXPECT_EQ(5u, none["handler_threads"].asUInt64());
    EXPECT_EQ(compact(cli::search_capabilities_json(4096, 2, 2, 2)["too_large"]),
              compact(none["too_large"]));
    EXPECT_EQ(compact(cli::search_capabilities_json(4096, 2, 2, 2)["unit_refused"]),
              compact(none["unit_refused"]));
    // the members a refusal states are the body's
    Json::Value members(Json::arrayValue);
    for (const std::string &name
            : cli::search_answer_too_large_json(4096, "values").getMemberNames()) {
        members.append(name);
    }
    EXPECT_EQ(compact(members), compact(none["too_large"]["members"]));
    // and a refused cut entry's `refused` states the members of unit_refused (priced_mb with
    // the answer account only)
    for (auto account : { graph::SearchCut::Refused::VALUES, graph::SearchCut::Refused::ANSWER }) {
        SeqSearchResult stub(QuerySequence { 0, "q", "ACGT" }, SeqSearchResult::LabelCountVec());
        stub.refuse(account, 4096, 5000);
        const Json::Value refused = stub.cut_to_json()["refused"];
        for (const Json::Value &name : none["unit_refused"]["members"]) {
            EXPECT_EQ(name.asString() != "priced_mb"
                              || account == graph::SearchCut::Refused::ANSWER,
                      refused.isMember(name.asString())) << name;
        }
        EXPECT_EQ(account == graph::SearchCut::Refused::ANSWER ? 4u : 3u, refused.size());
    }
}

// Test 6 (lane T7-P1): a refused unit's cut entry, exactly, with and without the medians'
// members (the computed cut's), for each account; its rows are dropped, its alignment and
// sequence kept
TEST(SearchAnswer, ARefusedUnitsCutEntry) {
    std::mt19937 rng(53);
    for (bool abundances : { false, true }) {
        for (bool alignment : { false, true }) {
            SCOPED_TRACE(::testing::Message() << "abundances " << abundances << " alignment "
                         << alignment);
            SeqSearchResult result = record(rng, abundances ? 2 : 3, 40, 300, 12, false, 0.5,
                                            alignment);
            graph::SearchCut cut = result.get_cut();
            ASSERT_EQ(abundances, cut.abundances);
            const Json::Value object = result.to_json(false, 31);
            SeqSearchResult values = std::move(result);
            values.refuse(graph::SearchCut::Refused::VALUES, 4096);
            std::visit([](const auto &rows) { EXPECT_TRUE(rows.empty()); }, values.get_result());
            EXPECT_EQ(abundances
                          ? "{\"complete\":false,\"nth_median\":null,\"nth_score\":null,"
                            "\"refused\":{\"account\":\"values\",\"code\":"
                            "\"search_answer_too_large\",\"max_memory_mb\":4096},\"returned\":0,"
                            "\"tied\":0,\"tied_returned\":0,\"tied_same_median\":0}"
                          : "{\"complete\":false,\"nth_score\":null,\"refused\":{\"account\":"
                            "\"values\",\"code\":\"search_answer_too_large\",\"max_memory_mb\":"
                            "4096},\"returned\":0,\"tied\":0,\"tied_returned\":0}",
                      compact(values.cut_to_json()));
            // the object: the record's members as before, no row
            Json::Value expected = object;
            expected["results"] = Json::Value(Json::arrayValue);
            EXPECT_EQ(compact(expected), compact(values.to_json(false, 31)));
            EXPECT_EQ(alignment, bool(values.get_alignment()));
            // the answer account states the unit's own price in MiB
            SeqSearchResult answer = record(rng, abundances ? 2 : 3, 4, 30, 12, false, 0.5);
            answer.refuse(graph::SearchCut::Refused::ANSWER, 1, 7);
            EXPECT_EQ(abundances
                          ? "{\"complete\":false,\"nth_median\":null,\"nth_score\":null,"
                            "\"refused\":{\"account\":\"answer\",\"code\":"
                            "\"search_answer_too_large\",\"max_memory_mb\":1,\"priced_mb\":7},"
                            "\"returned\":0,\"tied\":0,\"tied_returned\":0,"
                            "\"tied_same_median\":0}"
                          : "{\"complete\":false,\"nth_score\":null,\"refused\":{\"account\":"
                            "\"answer\",\"code\":\"search_answer_too_large\",\"max_memory_mb\":1,"
                            "\"priced_mb\":7},\"returned\":0,\"tied\":0,\"tied_returned\":0}",
                      compact(answer.cut_to_json()));
        }
    }
    // an answered unit's entry has no `refused` member
    EXPECT_FALSE(record(rng, 0, 3, 30, 12, false, 0.5).cut_to_json().isMember("refused"));
}

// Test 7 (lane T7-P1): the price of a refused stub covers the tree the server builds of it
// (allocator-measured where it can be), tagged and untagged, each account
TEST(SearchAnswer, TheStubsPriceCoversItsTree) {
    std::mt19937 rng(59);
    size_t measured = 0;
    const std::string graph = "vo_sauropsida";
    const std::string graph_path = "/scratch/nvme9/refseq_r237/vo_sauropsida/graph_k31.dbg";
    const std::string annotation_path
            = "/scratch/nvme9/refseq_r237/vo_sauropsida/annotation.row_diff_brwt_coord.annodbg";
    for (int mode : { 0, 1, 2, 3 }) {
        for (bool tagged : { false, true }) {
            for (bool alignment : { false, true }) {
                for (auto account : { graph::SearchCut::Refused::VALUES,
                                      graph::SearchCut::Refused::ANSWER }) {
                    SeqSearchResult stub = record(rng, mode, 30, 300, 20, true, 0.5, alignment);
                    const cli::SearchRecordPrice full = cli::search_record_price(
                            stub, false, tagged, graph, graph_path, annotation_path,
                            "0123456789abcdef");
                    stub.refuse(account, 4096, 123456);
                    const cli::SearchRecordPrice price = cli::search_record_price(
                            stub, false, tagged, graph, graph_path, annotation_path,
                            "0123456789abcdef");
                    // no row: a fraction of the unit's, nothing transient, no scratch
                    EXPECT_LT(price.tree, full.tree);
                    EXPECT_EQ(0u, price.transient);
                    EXPECT_EQ(0u, cli::search_record_scratch(stub, false));
                    EXPECT_EQ(0u, cli::search_record_rows_bytes(stub)
                                  - 3 * sizeof(SeqSearchResult)
                                  - (stub.get_sequence().name.size() > 15
                                         ? alloc_bytes(stub.get_sequence().name.size() + 1) : 0)
                                  - (stub.get_sequence().sequence.size() > 15
                                         ? alloc_bytes(stub.get_sequence().sequence.size() + 1)
                                         : 0)
                                  - (alignment ? alloc_bytes(stub.get_alignment()->cigar.size()
                                                             + 1) : 0));
                    Json::Value answer(Json::arrayValue);
                    const std::optional<int64_t> before = live_heap_bytes();
                    answer.append(server_object(stub, false, tagged, graph, graph_path,
                                                annotation_path, "0123456789abcdef"));
                    const std::optional<int64_t> after = live_heap_bytes();
                    if (before && after && *after > *before) {
                        ++measured;
                        EXPECT_LE(uint64_t(*after - *before), price.tree)
                                << "mode " << mode << " tagged " << tagged;
                        EXPECT_LE(price.tree, 4 * uint64_t(*after - *before) + (64 << 10));
                    }
                }
            }
        }
    }
    if (!measured)
        GTEST_SKIP() << "no allocator counters to compare the prices with";
}

// Tests 3 and 13 (lane T7-P1): a unit's object check, charge_search_unit: an object priced
// exactly the bound passes it and a byte more does not; with the field the unit is then
// refused alone, its stub states the unit's own price in MiB rounded up and the stub's price is
// charged; without it the request is refused on the answer side and nothing is charged
TEST(SearchAnswer, TheObjectCheckOfAUnit) {
    std::mt19937 rng(61);
    const SeqSearchResult unit = record(rng, 3, 200, 400, 30, false, 0.6);
    const cli::SearchRecordPrice price = cli::search_record_price(unit, false, true, "g", "g.dbg",
                                                                  "g.annodbg", "");
    const uint64_t own = price.tree + price.transient;
    const uint64_t scratch = cli::search_record_scratch(unit, false);
    for (bool partial : { false, true }) {
        SCOPED_TRACE(partial);
        // exactly the bound: answered, charged its tree and transient
        {
            SearchBudget budget(own, partial);
            SeqSearchResult copy = unit;
            uint64_t transient_charged = 0;
            EXPECT_TRUE(cli::charge_search_unit(budget, &transient_charged, &copy, false, true,
                                                "g", "g.dbg", "g.annodbg", "", "pair (g)"));
            EXPECT_EQ(graph::SearchCut::Refused::NO, copy.get_cut().refused);
            EXPECT_EQ(own, budget.answer().held());
            EXPECT_EQ(1u, budget.units());
            EXPECT_FALSE(budget.passed());
        }
        // a byte below its price: the unit's object alone passes the bound
        SearchBudget budget(own - 1, partial);
        SeqSearchResult copy = unit;
        uint64_t transient_charged = 0;
        const bool built = cli::charge_search_unit(budget, &transient_charged, &copy, false, true,
                                                   "g", "g.dbg", "g.annodbg", "", "pair (g)");
        if (!partial) {
            EXPECT_FALSE(built);
            EXPECT_TRUE(budget.passed());
            EXPECT_EQ("answer", std::string(budget.account()));
            EXPECT_EQ("the answer account of record 0 on pair (g)", budget.first_passed());
            EXPECT_EQ(0u, budget.first_held());
            EXPECT_EQ(own, budget.first_charge());
            // only the scratch of its measurement was charged
            EXPECT_EQ(scratch, budget.answer().held());
            EXPECT_EQ(0u, budget.units());
            continue;
        }
        EXPECT_TRUE(built);
        EXPECT_FALSE(budget.passed());
        EXPECT_EQ(graph::SearchCut::Refused::ANSWER, copy.get_cut().refused);
        EXPECT_EQ(graph::mib_rounded_up(own), copy.get_cut().refused_priced_mb);
        EXPECT_EQ((own - 1) >> 20, copy.get_cut().refused_max_mb);
        ASSERT_EQ(1u, budget.unit_refusals().size());
        EXPECT_EQ("record 0 on pair (g)", budget.unit_refusals()[0].where);
        EXPECT_EQ("answer", std::string(budget.unit_refusals()[0].account));
        EXPECT_EQ(0u, budget.unit_refusals()[0].held);
        EXPECT_EQ(own, budget.unit_refusals()[0].charge);
        // charged: the scratch of the measurement (the pair's transient), then the stub's tree
        const cli::SearchRecordPrice stub = cli::search_record_price(copy, false, true, "g",
                                                                     "g.dbg", "g.annodbg", "");
        EXPECT_EQ(scratch + stub.tree, budget.answer().held());
        EXPECT_EQ(scratch, transient_charged);
        EXPECT_EQ(1u, budget.units());
        // a stub refused on its values is charged its price and never checked again
        SeqSearchResult values = unit;
        values.refuse(graph::SearchCut::Refused::VALUES, 1);
        SearchBudget tiny(stub.tree + 1000, true);
        transient_charged = 0;
        EXPECT_TRUE(cli::charge_search_unit(tiny, &transient_charged, &values, false, true, "g",
                                            "g.dbg", "g.annodbg", "", "pair (g)"));
        EXPECT_EQ(graph::SearchCut::Refused::VALUES, values.get_cut().refused);
        EXPECT_TRUE(tiny.unit_refusals().empty());
    }
    // a price past the range saturates: priced_mb 2^44
    EXPECT_EQ(uint64_t(1) << 44, graph::mib_rounded_up(graph::saturating_add(
            std::numeric_limits<uint64_t>::max(), 1)));
}

// An allocator that counts what its buffers hold (each at its size class, alloc_bytes) and the
// most they held at once
struct Held {
    uint64_t now = 0;
    uint64_t most = 0;
};
template <class T>
struct CountingAllocator {
    using value_type = T;
    Held *held;
    explicit CountingAllocator(Held *held) : held(held) {}
    template <class U>
    CountingAllocator(const CountingAllocator<U> &other) : held(other.held) {}
    T* allocate(size_t n) {
        held->now += alloc_bytes(n * sizeof(T));
        held->most = std::max(held->most, held->now);
        return std::allocator<T>().allocate(n);
    }
    void deallocate(T *p, size_t n) {
        held->now -= alloc_bytes(n * sizeof(T));
        std::allocator<T>().deallocate(p, n);
    }
    template <class U>
    bool operator==(const CountingAllocator<U> &other) const { return held == other.held; }
    template <class U>
    bool operator!=(const CountingAllocator<U> &other) const { return held != other.held; }
};

// The most the two vectors of open ranges of collapse_coord_ranges (cli/query.cpp,
// call_coord_ranges) hold at once for |coords|: each position empties the next ones
// (resize(0)), adds one range per coordinate (extended or new) and swaps them with the open
// ones; a gap before a position ends every open range (a position without coordinates)
uint64_t open_ranges_peak(const PositionCoords &coords) {
    using Range = std::tuple<uint64_t, uint64_t, uint64_t>;
    Held held;
    std::vector<Range, CountingAllocator<Range>> ranges { CountingAllocator<Range>(&held) };
    std::vector<Range, CountingAllocator<Range>> next_ranges { CountingAllocator<Range>(&held) };
    auto position = [&](size_t coordinates) {
        next_ranges.resize(0);
        for (size_t c = 0; c < coordinates; ++c) {
            next_ranges.emplace_back(c, c, c);
        }
        ranges.swap(next_ranges);
    };
    size_t next = 0;
    coords.call_positions([&](size_t i, const PositionCoords::Coords &list) {
        if (i > next)
            position(0);
        position(list.size());
        next = i + 1;
    });
    return held.most;
}

// REVIEW 9726cb420 [P2], lane T6-W2: the ranges a coordinate row's collapse keeps open are
// priced at their peak, coord_ranges_bytes of the row's longest list: two vectors grown by
// doubling, the old buffer of one beside its new one at a growth (2 x buffer_bytes(open, 24)
// before, which missed both)
TEST(SearchAnswer, TheOpenRangesArePricedAtTheirPeak) {
    EXPECT_EQ(0u, cli::coord_ranges_bytes(0));
    std::mt19937 rng(47);
    for (size_t round = 0; round < 200; ++round) {
        const size_t size = 1 + rng() % 40;
        PositionCoords coords(size);
        uint64_t longest = 0;
        for (size_t i = 0; i < size; ++i) {
            if (rng() % 3 == 0)
                continue;
            const size_t n = 1 + rng() % (round < 100 ? 70 : 5000);
            longest = std::max<uint64_t>(longest, n);
            for (size_t c = 0; c < n; ++c) {
                coords.add(i, 1000 + c * (1 + rng() % 3) + i);
            }
        }
        const uint64_t peak = open_ranges_peak(coords);
        EXPECT_LE(peak, cli::coord_ranges_bytes(longest)) << round;
    }
    // the worst case is the bound: two positions with as many coordinates, just past a power
    // of two (the REVIEW's 30,000 too)
    for (size_t n : { size_t(2), size_t(3), size_t(17), size_t(1025), size_t(30000),
                      size_t(65537) }) {
        PositionCoords coords(2);
        for (size_t i = 0; i < 2; ++i) {
            for (size_t c = 0; c < n; ++c) {
                coords.add(i, 1000 + c + i);
            }
        }
        EXPECT_EQ(open_ranges_peak(coords), cli::coord_ranges_bytes(n)) << n;
    }
    EXPECT_EQ(2 * 786'432u + 393'216u, cli::coord_ranges_bytes(30000));
}

// REVIEW 9726cb420 [P2], lane T6-W2: what measuring a record's price holds (the open ranges
// of collapsed_coords_text for its coordinate rows) is charged before it is measured, as the
// record's transient: the REVIEW's record, two successive positions of 30,000 coordinates,
// fits a values account of 1 MiB, while measuring it would hold 1,966,080 B of ranges; it is
// refused with nothing charged
TEST(SearchAnswer, TheScratchOfThePriceIsChargedBeforeItIsMeasured) {
    SearchBudget budget(1 << 20);
    SearchAccount values(&budget, "values");
    PositionCoords coords(2, &values);
    for (size_t i = 0; i < 2; ++i) {
        for (uint64_t c = 0; c < 30000; ++c) {
            coords.add(i, 5000 + c + i);
        }
    }
    coords.shrink();
    EXPECT_FALSE(budget.passed());
    EXPECT_GE(uint64_t(1) << 20, values.held());
    auto record = [&](bool second) {
        SeqSearchResult::LabelCountCoordsVec rows;
        rows.emplace_back(second ? "b" : "a", 2, coords, 0);
        return SeqSearchResult(QuerySequence { size_t(second), "q", std::string(32, 'A') },
                               std::move(rows));
    };
    const SeqSearchResult first = record(false);
    const uint64_t scratch = cli::search_record_scratch(first, false);
    EXPECT_EQ(cli::coord_ranges_bytes(30000), scratch);
    EXPECT_LT(uint64_t(1) << 20, scratch);
    // verbose: the dense text is measured without allocating
    EXPECT_EQ(0u, cli::search_record_scratch(first, true));
    // other modes: nothing
    EXPECT_EQ(0u, cli::search_record_scratch(
            SeqSearchResult(QuerySequence { 0, "q", "" }, SeqSearchResult::LabelVec { "a" }),
            false));

    uint64_t transient_charged = 0;
    EXPECT_THROW(cli::charge_search_record(budget.answer(), &transient_charged, first, false,
                                           false, "", "", "", ""),
                 SearchAnswerTooLarge);
    EXPECT_EQ(0u, budget.answer().held());
    EXPECT_EQ(0u, transient_charged);
    EXPECT_TRUE(budget.passed());

    // a bound that fits: the record's tree and its transient (which covers the scratch), each
    // record's transient as the most one took
    SearchBudget large(kLarge);
    transient_charged = 0;
    const cli::SearchRecordPrice price = cli::charge_search_record(
            large.answer(), &transient_charged, first, false, false, "", "", "", "");
    EXPECT_LE(scratch, price.transient);
    EXPECT_EQ(price.transient, transient_charged);
    EXPECT_EQ(price.tree + price.transient, large.answer().held());
    const SeqSearchResult second = record(true);
    const cli::SearchRecordPrice price2 = cli::charge_search_record(
            large.answer(), &transient_charged, second, false, false, "", "", "", "");
    EXPECT_EQ(price.tree + price2.tree + std::max(price.transient, price2.transient),
              large.answer().held());
    // a bound of exactly that admits it, one byte less refuses it before its scratch
    SearchBudget exact(price.tree + price.transient);
    transient_charged = 0;
    cli::charge_search_record(exact.answer(), &transient_charged, first, false, false, "", "",
                              "", "");
    SearchBudget less(scratch - 1);
    transient_charged = 0;
    EXPECT_THROW(cli::charge_search_record(less.answer(), &transient_charged, first, false,
                                           false, "", "", "", ""),
                 SearchAnswerTooLarge);
    EXPECT_EQ(0u, less.answer().held());
}

} // namespace
} // namespace mtg::test
