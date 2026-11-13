#include <algorithm>
#include <atomic>
#include <limits>
#include <random>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "coordinate_test_graph.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "cli/query.hpp"
#include "graph/search_budget.hpp"


// The /search answer budget's values account (SPEC-multigraph-envelope.md section 3, lane
// T5b-Y4): the values a search's rows hold per query position (PresenceRuns, CountRuns,
// PositionCoords) charge every growth of their buffers to the pair's account BEFORE it is made,
// so that past the bound the query throws with nothing more held; with an account the values
// and the rows are those without one. The account only grows, so whether a request passes the
// bound does not depend on the order its threads charge in.

namespace mtg::test {
namespace {

using graph::PresenceRuns;
using graph::CountRuns;
using graph::PositionCoords;
using graph::SearchAccount;
using graph::SearchAnswerTooLarge;
using graph::SearchBudget;
using annot::matrix::alloc_bytes;

constexpr uint64_t kLarge = uint64_t(1) << 40;

std::string random_sequence(std::mt19937 &rng, size_t length) {
    std::string s(length, 'A');
    for (char &c : s) {
        c = "ACGT"[rng() % 4];
    }
    return s;
}

TEST(SearchBudget, ChargesUpToTheBound) {
    SearchBudget budget(2 << 20);
    SearchAccount values(&budget, "values", "pair (a, a.dbg, a.annodbg)");
    values.charge(1 << 20);
    budget.answer().charge(1 << 20);
    values.charge(1 << 20);
    EXPECT_EQ(uint64_t(2) << 20, values.held());
    EXPECT_FALSE(budget.passed());
    try {
        values.charge(1);
        FAIL() << "not refused";
    } catch (const SearchAnswerTooLarge &e) {
        EXPECT_EQ("search: the answer would hold more than --search-max-answer-memory-mb 2 MiB: "
                  "nothing partial is sent; ask for fewer top_labels, a larger "
                  "discovery_fraction, fewer or shorter FASTA records, or no per-position "
                  "values (with_signature, query_counts, query_coords)", std::string(e.what()));
        EXPECT_EQ(2u, e.max_mb());
        EXPECT_EQ("values", std::string(e.account()));
    }
    // nothing more held; the request is marked, and every further charge of it refused, on
    // every account (the work still running stops at its next charge)
    EXPECT_EQ(uint64_t(2) << 20, values.held());
    EXPECT_TRUE(budget.passed());
    EXPECT_EQ("the values account of pair (a, a.dbg, a.annodbg)", budget.first_passed());
    EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
    EXPECT_EQ(uint64_t(1) << 20, budget.answer().held());
    // a charge larger than the bound itself, and the largest charge
    SearchBudget small(10);
    EXPECT_THROW(small.answer().charge(11), SearchAnswerTooLarge);
    SearchBudget other(10);
    EXPECT_THROW(other.answer().charge(std::numeric_limits<uint64_t>::max()),
                 SearchAnswerTooLarge);
    EXPECT_EQ(0u, other.answer().held());
    // none charged without an account
    graph::charge_search(nullptr, 1 << 30);
}

// Charged from many threads in any order, the accounts refuse exactly when the charges sum to
// more than the bound (each account only grows), and hold their sum otherwise
TEST(SearchBudget, TheOutcomeIsTheSameInEveryOrder) {
    std::mt19937 rng(7);
    for (size_t round = 0; round < 40; ++round) {
        std::vector<uint64_t> charges(200 + rng() % 200);
        uint64_t total = 0;
        for (uint64_t &c : charges) {
            c = rng() % 1000;
            total += c;
        }
        // the bound around the total: below, at, above
        const uint64_t bound = total + (round % 3) - 1;
        SearchBudget budget(bound);
        std::atomic<size_t> refused = 0;
        std::vector<std::thread> threads;
        for (size_t t = 0; t < 4; ++t) {
            threads.emplace_back([&, t]() {
                for (size_t i = t; i < charges.size(); i += 4) {
                    try {
                        budget.answer().charge(charges[i]);
                    } catch (const SearchAnswerTooLarge &) {
                        ++refused;
                    }
                }
            });
        }
        for (auto &thread : threads) {
            thread.join();
        }
        EXPECT_EQ(total > bound, budget.passed()) << round;
        EXPECT_EQ(total > bound, refused > 0) << round;
        if (total <= bound) {
            EXPECT_EQ(total, budget.answer().held());
        }
        EXPECT_LE(budget.answer().held(), bound);
    }
}

// Built with an account, each value is the one built without (in either form), and the
// account holds at least the bytes of its buffers (exactly the model's bytes of them while the
// account is attached: the growths' net bytes add up to the last buffer's)
TEST(SearchBudget, ValuesChargeTheirBuffers) {
    std::mt19937 rng(11);
    for (size_t size : { 1, 5, 63, 64, 65, 300, 1000, 5000 }) {
        for (size_t step : { 1, 2, 3, 7, 64, 200 }) {
            SCOPED_TRACE(::testing::Message() << size << " step " << step);
            SearchBudget budget(kLarge);
            SearchAccount account(&budget, "values");

            PresenceRuns plain_mask(size);
            PresenceRuns mask(size, &account);
            for (size_t i = 0; i < size; i += step) {
                plain_mask.add(i);
                mask.add(i);
                EXPECT_LE(alloc_bytes(mask.heap_bytes()), account.held());
            }
            EXPECT_EQ(plain_mask, mask);
            EXPECT_EQ(cli::encode_presence_mask(plain_mask), cli::encode_presence_mask(mask));
            EXPECT_EQ(mask.is_dense()
                          ? alloc_bytes(PresenceRuns::dense_bytes(size) + sizeof(uint64_t))
                          : alloc_bytes(mask.heap_bytes()),
                      account.held());

            SearchAccount count_account(&budget, "values");
            CountRuns plain_counts(size);
            CountRuns counts(size, &count_account);
            for (size_t i = 0; i < size; i += step) {
                const size_t value = 1 + rng() % 3;
                plain_counts.add(i, value);
                counts.add(i, value);
                EXPECT_LE(alloc_bytes(counts.heap_bytes()), count_account.held());
            }
            EXPECT_EQ(plain_counts, counts);
            EXPECT_EQ(cli::count_run_strings(plain_counts), cli::count_run_strings(counts));

            SearchAccount coords_account(&budget, "values");
            PositionCoords plain_coords(size);
            PositionCoords coords(size, &coords_account);
            uint64_t lists = 0;
            for (size_t i = 0; i < size; i += step) {
                const size_t n = 1 + rng() % 5;
                for (size_t j = 0; j < n; ++j) {
                    plain_coords.add(i, 1000 * i + j);
                    coords.add(i, 1000 * i + j);
                }
            }
            coords.call_positions([&](size_t, const PositionCoords::Coords &list) {
                lists += PositionCoords::list_bytes(list);
            });
            EXPECT_EQ(plain_coords, coords);
            EXPECT_EQ(cli::collapse_coord_ranges(plain_coords), cli::collapse_coord_ranges(coords));
            // (heap_bytes() sums two buffers, the runs' and the lists', each charged as its own
            // allocation: the account holds at least their bytes)
            EXPECT_LE(coords.heap_bytes() + lists, coords_account.held());
        }
    }
    // long coordinate lists, grown many times: the account holds each list's buffer as the
    // model states it from the list's capacity
    for (size_t n : { 3, 5, 9, 17, 100, 1000, 20000 }) {
        SCOPED_TRACE(::testing::Message() << n << " coordinates");
        SearchBudget budget(kLarge);
        SearchAccount account(&budget, "values");
        PositionCoords coords(3, &account);
        for (size_t j = 0; j < n; ++j) {
            coords.add(1, j);
        }
        uint64_t lists = 0;
        coords.call_positions([&](size_t, const PositionCoords::Coords &list) {
            lists += PositionCoords::list_bytes(list);
        });
        EXPECT_LE(coords.heap_bytes() + lists, account.held());
        EXPECT_GE(alloc_bytes(coords.heap_bytes()) + lists + 64, account.held());
    }
    // detached by shrink(): a later growth (of it or of a copy) charges nothing
    SearchBudget budget(kLarge);
    SearchAccount account(&budget, "values");
    PresenceRuns mask(100, &account);
    mask.add(0);
    mask.shrink();
    const uint64_t held = account.held();
    PresenceRuns copy = mask;
    mask.add(5);
    mask.add(7);
    copy.add(9);
    copy.add(11);
    EXPECT_EQ(held, account.held());
}

// An account a growth does not fit refuses it before it is made: the value then holds no more
// than the account was charged with, and stays valid (its outputs as far as it was built)
TEST(SearchBudget, AGrowthThatDoesNotFitIsRefusedBeforeItIsMade) {
    for (size_t size : { 100, 1000, 100000 }) {
        for (uint64_t bound : { 0, 16, 100, 1000, 5000 }) {
            SCOPED_TRACE(::testing::Message() << size << " bound " << bound);
            {
                SearchBudget budget(bound);
                SearchAccount account(&budget, "values");
                PresenceRuns mask(size, &account);
                size_t added = 0;
                try {
                    for (size_t i = 0; i < size; i += 2) {
                        mask.add(i);
                        ++added;
                    }
                } catch (const SearchAnswerTooLarge &) {
                    EXPECT_TRUE(budget.passed());
                }
                EXPECT_LE(account.held(), bound);
                EXPECT_LE(alloc_bytes(mask.heap_bytes()), account.held());
                EXPECT_EQ(added, mask.count());
            }
            {
                SearchBudget budget(bound);
                SearchAccount account(&budget, "values");
                CountRuns counts(size, &account);
                try {
                    for (size_t i = 0; i < size; ++i) {
                        counts.add(i, 1 + i % 2);
                    }
                } catch (const SearchAnswerTooLarge &) {}
                EXPECT_LE(account.held(), bound);
                EXPECT_LE(alloc_bytes(counts.heap_bytes()), account.held());
            }
            {
                SearchBudget budget(bound);
                SearchAccount account(&budget, "values");
                PositionCoords coords(size, &account);
                try {
                    for (size_t i = 0; i < size; ++i) {
                        for (size_t j = 0; j < 4; ++j) {
                            coords.add(i, 10 * i + j);
                        }
                    }
                } catch (const SearchAnswerTooLarge &) {}
                uint64_t lists = 0;
                coords.call_positions([&](size_t, const PositionCoords::Coords &list) {
                    lists += PositionCoords::list_bytes(list);
                });
                EXPECT_LE(account.held(), bound);
                EXPECT_LE(coords.heap_bytes() + lists, account.held());
            }
        }
    }
}

using Signatures = std::vector<std::tuple<std::string, size_t, PresenceRuns>>;
using Counts = std::vector<std::tuple<std::string, size_t, CountRuns>>;
using Coordinates = std::vector<std::tuple<std::string, size_t, PositionCoords, size_t>>;

// The queries of every per-position mode with an account give the rows they give without one
// (header and column paths, discovery 0 and 1); the account's charge is the same every time;
// an account a byte smaller refuses the query
TEST(SearchBudget, QueriesWithAnAccountAnswerAsWithout) {
    std::mt19937 rng(23);
    const size_t k = 5;
    std::vector<CoordinateRecord> records;
    for (size_t i = 0; i < 30; ++i) {
        records.push_back({ random_sequence(rng, 30 + rng() % 40), "F" + std::to_string(i % 3),
                            "r" + std::to_string(i) });
    }
    const std::string query = records[0].sequence + records[3].sequence + records[3].sequence
                              + random_sequence(rng, 20) + records[7].sequence;
    const size_t top = 1000;
    for (bool by_record : { true, false }) {
        for (bool row_diff : { false, true }) {
            auto index = coordinate_test_graph(k, records, row_diff, false, false, by_record);
            for (double discovery : { 0.0, 1.0 }) {
                SCOPED_TRACE(::testing::Message() << "by_record " << by_record << " row_diff "
                             << row_diff << " discovery " << discovery);
                auto check = [&](const auto &query_with) {
                    const auto plain = query_with(nullptr);
                    SearchBudget budget(kLarge);
                    SearchAccount account(&budget, "values");
                    const auto got = query_with(&account);
                    EXPECT_EQ(plain, got);
                    SearchBudget again(kLarge);
                    SearchAccount account_again(&again, "values");
                    query_with(&account_again);
                    EXPECT_EQ(account.held(), account_again.held());
                    if (account.held()) {
                        SearchBudget short_one(account.held() - 1);
                        SearchAccount refused(&short_one, "values");
                        EXPECT_THROW(query_with(&refused), SearchAnswerTooLarge);
                        EXPECT_TRUE(short_one.passed());
                    }
                    return account.held();
                };
                const uint64_t signatures = check([&](SearchAccount *a) -> Signatures {
                    return index->get_top_label_signatures(query, top, discovery, 0, nullptr, a);
                });
                EXPECT_LT(0u, signatures);
                if (discovery == 1.0)
                    continue;
                EXPECT_LT(0u, check([&](SearchAccount *a) -> Counts {
                    return index->get_kmer_counts(query, top, discovery, 0, nullptr, a);
                }));
                EXPECT_LT(0u, check([&](SearchAccount *a) -> Coordinates {
                    return index->get_kmer_coordinates(query, top, discovery, 0, nullptr, a);
                }));
            }
        }
    }
}

} // namespace
} // namespace mtg::test
