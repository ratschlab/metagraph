#include <algorithm>
#include <atomic>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "coordinate_test_graph.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/binary_matrix/base/hopscotch_layout.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "cli/query.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
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
using graph::SearchUnitRefused;
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
    // for the log: what that account held before the charge refused, and the charge
    EXPECT_EQ(uint64_t(2) << 20, budget.first_held());
    EXPECT_EQ(1u, budget.first_charge());
    EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
    EXPECT_EQ(uint64_t(1) << 20, budget.answer().held());
    // the first account to pass the bound is the one kept
    EXPECT_EQ("the values account of pair (a, a.dbg, a.annodbg)", budget.first_passed());
    EXPECT_EQ(uint64_t(2) << 20, budget.first_held());
    EXPECT_EQ(1u, budget.first_charge());
    // a charge larger than the bound itself, and the largest charge
    SearchBudget small(10);
    small.answer().charge(3);
    EXPECT_THROW(small.answer().charge(11), SearchAnswerTooLarge);
    EXPECT_EQ("the answer account", small.first_passed());
    EXPECT_EQ(3u, small.first_held());
    EXPECT_EQ(11u, small.first_charge());
    SearchBudget other(10);
    EXPECT_THROW(other.answer().charge(std::numeric_limits<uint64_t>::max()),
                 SearchAnswerTooLarge);
    EXPECT_EQ(0u, other.answer().held());
    EXPECT_EQ(0u, other.first_held());
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(), other.first_charge());
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

// REVIEW 9726cb420 [P2], lane T6-W2: the scratch one task of an account's work holds while it
// runs (a record's store of its rows) is charged as the most one task took, times the tasks
// that run at once, raised before the scratch is allocated: a freed scratch is not charged again
// for each task, and the total is the same in every order
TEST(SearchBudget, AScratchIsTheMostOneTaskTookTimesItsThreads) {
    SearchBudget budget(kLarge);
    SearchAccount values(&budget, "values", "the graph", 3);
    EXPECT_EQ(3u, values.threads());
    values.charge(10);
    values.charge_scratch(100);
    EXPECT_EQ(310u, values.held());
    values.charge_scratch(50);
    values.charge_scratch(100);
    EXPECT_EQ(310u, values.held());
    values.charge_scratch(150);
    EXPECT_EQ(460u, values.held());
    EXPECT_EQ(150u, values.scratch());
    graph::charge_search_scratch(&values, 160);
    EXPECT_EQ(490u, values.held());
    // none charged without an account, and none for nothing
    graph::charge_search_scratch(nullptr, 1 << 30);
    graph::charge_search_scratch(&values, 0);
    EXPECT_EQ(490u, values.held());
    // one thread by default; at least one
    SearchAccount one(&budget, "values");
    one.charge_scratch(7);
    EXPECT_EQ(7u, one.held());
    SearchAccount none(&budget, "values", "", 0);
    EXPECT_EQ(1u, none.threads());

    // past the bound: refused before anything is held, the request marked; a product past the
    // range of the bytes is past every bound
    SearchBudget small(1000);
    SearchAccount two(&small, "values", "", 2);
    two.charge_scratch(400);
    EXPECT_EQ(800u, two.held());
    EXPECT_THROW(two.charge_scratch(501), SearchAnswerTooLarge);
    EXPECT_EQ(800u, two.held());
    EXPECT_TRUE(small.passed());
    SearchBudget huge(kLarge);
    SearchAccount many(&huge, "values", "", 1000);
    EXPECT_THROW(many.charge_scratch(std::numeric_limits<uint64_t>::max() / 100),
                 SearchAnswerTooLarge);
    EXPECT_EQ(0u, many.held());

    // tasks on threads, in any order: the account holds the threads times the largest scratch,
    // and the bound refuses the same requests
    std::mt19937 rng(29);
    for (size_t round = 0; round < 30; ++round) {
        std::vector<uint64_t> scratch(100 + rng() % 100);
        uint64_t most = 0;
        for (uint64_t &bytes : scratch) {
            bytes = rng() % 10000;
            most = std::max(most, bytes);
        }
        uint64_t plain = 0;
        for (size_t i = 0; i < scratch.size(); i += 7) {
            plain += 25;
        }
        // the bound around the total: below, at, above
        const uint64_t bound = plain + 4 * most + (round % 3) - 1;
        SearchBudget b(bound);
        SearchAccount account(&b, "values", "", 4);
        std::vector<std::thread> threads;
        for (size_t t = 0; t < 4; ++t) {
            threads.emplace_back([&, t]() {
                for (size_t i = t; i < scratch.size(); i += 4) {
                    try {
                        account.charge(i % 7 == 0 ? 25 : 0);
                        account.charge_scratch(scratch[i]);
                    } catch (const SearchAnswerTooLarge &) {}
                }
            });
        }
        for (auto &thread : threads) {
            thread.join();
        }
        EXPECT_EQ(plain + 4 * most > bound, b.passed()) << round;
        if (!b.passed()) {
            EXPECT_EQ(plain + 4 * most, account.held()) << round;
        }
    }
}

// A column-label search's store of its rows' values (utils::ValueStore) is priced as what it
// allocates: below one million labels its dense index over the labels and the values' buffer,
// from one million on the hash table tsl reserves for the rows kept (both sides of the
// threshold, each value type)
template <class Values>
void check_store_prices() {
    using Store = utils::ValueStore<uint32_t, Values>;
    for (size_t num_labels : { size_t(1), size_t(999'999), size_t(1'000'000), size_t(5'000'000) }) {
        for (size_t rows : { 1, 2, 3, 7, 14, 15, 100, 1000, 7373, 10000 }) {
            SCOPED_TRACE(::testing::Message() << num_labels << " labels, " << rows << " rows");
            if (rows > num_labels)
                continue;
            Store store(num_labels, rows);
            EXPECT_EQ(num_labels < 1'000'000, store.is_dense());
            EXPECT_EQ(store.is_dense(), Store::dense(num_labels));
            // spread ids, as a search's labels are
            for (size_t r = 0; r < rows; ++r) {
                store.emplace(uint32_t((r * 7919) % num_labels), Values(10));
            }
            const uint64_t priced = graph::search_store_bytes<Values>(num_labels, rows);
            if (store.is_dense()) {
                EXPECT_EQ(alloc_bytes(num_labels * sizeof(size_t))
                              + alloc_bytes(rows * sizeof(Values)), priced);
            } else {
                // the table as it is, its buckets and its overflow
                EXPECT_EQ(annot::matrix::hopscotch_held_bytes(store.sparse_values()), priced);
                EXPECT_EQ(0u, store.sparse_values().overflow_size());
            }
        }
    }
    // the dense index is the labels', whatever the rows kept: 900,000 labels and one row
    EXPECT_LE(7'200'000u, graph::search_store_bytes<Values>(900'000, 1));
}

TEST(SearchBudget, TheStoreOfAColumnSearchIsPricedAsItIsAllocated) {
    check_store_prices<PresenceRuns>();
    check_store_prices<CountRuns>();
    check_store_prices<PositionCoords>();
}

// An annotation of |num_labels| labels (row-major, so that the labels cost no column each) on
// which the label "L0" covers every k-mer of |sequence| and no other label covers any
std::unique_ptr<graph::AnnotatedDBG> many_labels(const std::string &sequence, size_t k,
                                                 size_t num_labels) {
    auto graph = std::make_shared<graph::DBGHashOrdered>(k);
    graph->add_sequence(sequence);
    Vector<annot::matrix::BinaryMatrix::SetBitPositions> rows(graph->max_index());
    graph->call_nodes([&](auto node) { rows[node - 1].push_back(0); });
    std::vector<std::string> labels;
    labels.reserve(num_labels);
    for (size_t j = 0; j < num_labels; ++j) {
        labels.push_back("L" + std::to_string(j));
    }
    return std::make_unique<graph::AnnotatedDBG>(
            graph, std::make_unique<annot::RowCompressed<>>(std::move(rows), labels));
}

// REVIEW 9726cb420 [P2], lane T6-W2: a column-label search on fewer than one million labels
// builds a dense index of them for its store whatever the rows it keeps, and its values account
// holds it (the base charged the hash table's estimate of the one row kept, a few hundred bytes,
// and answered under a bound of 1 MiB); from one million labels on, the hash table it builds
TEST(SearchBudget, AColumnSearchChargesTheStoreItBuilds) {
    std::mt19937 rng(31);
    const std::string sequence = random_sequence(rng, 60);
    for (size_t num_labels : { size_t(900'000), size_t(1'000'000) }) {
        SCOPED_TRACE(::testing::Message() << num_labels << " labels");
        auto index = many_labels(sequence, 11, num_labels);
        // what the store allocates: the dense index, or a table reserved for one row
        uint64_t store;
        if (num_labels < 1'000'000) {
            store = 900'000 * sizeof(size_t);
        } else {
            tsl::hopscotch_map<uint32_t, PresenceRuns> table;
            table.reserve(1);
            store = annot::matrix::hopscotch_held_bytes(table);
        }
        SearchBudget budget(kLarge);
        SearchAccount account(&budget, "values");
        const auto rows = index->get_top_label_signatures(sequence, 1, 0, 0, nullptr, &account);
        ASSERT_EQ(1u, rows.size());
        EXPECT_EQ("L0", std::get<0>(rows[0]));
        EXPECT_LE(store, account.held());
        // the same rows as without an account
        EXPECT_EQ(index->get_top_label_signatures(sequence, 1, 0, 0), rows);
        if (num_labels < 1'000'000) {
            // a bound of 1 MiB refuses it before the store is built
            SearchBudget small(1 << 20);
            SearchAccount refused(&small, "values");
            EXPECT_THROW(index->get_top_label_signatures(sequence, 1, 0, 0, nullptr, &refused),
                         SearchAnswerTooLarge);
            EXPECT_TRUE(small.passed());
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Lane T7-P1: the /search partial refusal per unit (a FASTA record on a pair;
// DESIGN-search-partial-refusal.md sections 3 and 4). Each unit has its own values account; with
// the request's `partial_refusal` a unit past the bound is refused alone (SearchUnitRefused, the
// request open), without it the request is refused, and the 400 names "values" whenever a
// unit's values would pass, else "answer".

constexpr char kPair[] = "pair (a, a.dbg, a.annodbg)";

// Test 1 (scope): a unit's account refuses its unit alone with the field, and the request with
// it; once the request is refused every charge throws as its state says
TEST(SearchBudget, AUnitIsRefusedAloneWithTheField) {
    SearchBudget budget(100, true);
    EXPECT_TRUE(budget.partial_refusal());
    SearchAccount unit(&budget, 3, kPair);
    EXPECT_TRUE(unit.unit());
    EXPECT_EQ(3u, unit.record());
    EXPECT_EQ("record 3 on pair (a, a.dbg, a.annodbg)", unit.where());
    EXPECT_EQ("values", std::string(unit.name()));
    EXPECT_EQ(1u, unit.threads());
    unit.charge(60);
    try {
        unit.charge(41);
        FAIL() << "not refused";
    } catch (const SearchUnitRefused &e) {
        EXPECT_EQ("values", std::string(e.account()));
    }
    // the unit alone: nothing more held, the request open, its other accounts go on
    EXPECT_EQ(60u, unit.held());
    EXPECT_TRUE(unit.refused());
    EXPECT_FALSE(budget.passed());
    SearchAccount other(&budget, 4, kPair);
    other.charge(100);
    budget.answer().charge(100);
    // sticky: the refused unit refuses again, whatever the charge, and is noted once
    EXPECT_THROW(unit.charge(1), SearchUnitRefused);
    EXPECT_THROW(unit.charge_scratch(1), SearchUnitRefused);
    ASSERT_EQ(1u, budget.unit_refusals().size());
    const auto refusal = budget.unit_refusals()[0];
    EXPECT_EQ(3u, refusal.record);
    EXPECT_EQ("record 3 on pair (a, a.dbg, a.annodbg)", refusal.where);
    EXPECT_EQ("values", std::string(refusal.account));
    EXPECT_EQ(60u, refusal.held);
    EXPECT_EQ(41u, refusal.charge);
    // the request's answer past the bound refuses the request ("answer"), and every account
    // stops at its next charge, the units' too
    EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
    EXPECT_TRUE(budget.passed());
    EXPECT_EQ("answer", std::string(budget.account()));
    EXPECT_EQ("the answer account", budget.first_passed());
    EXPECT_EQ(100u, budget.first_held());
    EXPECT_EQ(1u, budget.first_charge());
    SearchAccount late(&budget, 5, kPair);
    EXPECT_THROW(late.charge(1), SearchAnswerTooLarge);
    EXPECT_THROW(other.charge(0 + 1), SearchAnswerTooLarge);
    EXPECT_EQ(0u, late.held());
    // a unit's values past the bound do not change a request refused on its answer
    SearchAccount unit2(&budget, 6, kPair);
    EXPECT_THROW(unit2.charge(1000), SearchAnswerTooLarge);
    EXPECT_EQ("answer", std::string(budget.account()));
}

TEST(SearchBudget, AUnitRefusesTheRequestWithoutTheField) {
    SearchBudget budget(100);
    EXPECT_FALSE(budget.partial_refusal());
    SearchAccount unit(&budget, 7, "the graph");
    unit.charge(100);
    try {
        unit.charge(1);
        FAIL() << "not refused";
    } catch (const SearchAnswerTooLarge &e) {
        EXPECT_EQ("values", std::string(e.account()));
    }
    EXPECT_FALSE(unit.refused());
    EXPECT_TRUE(budget.passed());
    EXPECT_EQ("values", std::string(budget.account()));
    EXPECT_EQ("the values account of record 7 on the graph", budget.first_passed());
    EXPECT_EQ(100u, budget.first_held());
    EXPECT_EQ(1u, budget.first_charge());
    EXPECT_TRUE(budget.unit_refusals().empty());
    // every account stops at its next charge (the early stop of train 5b)
    SearchAccount other(&budget, 8, "the graph");
    EXPECT_THROW(other.charge(1), SearchAnswerTooLarge);
    EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
}

// Without the field, a refusal on the answer side leaves the account pending: no object is
// charged any more, the units' values accounts go on, and a unit's values past the bound then
// name the 400 "values"; otherwise it names "answer"
TEST(SearchBudget, AnAnswerSideRefusalLeavesTheAccountPending) {
    for (bool values_pass : { false, true }) {
        SCOPED_TRACE(values_pass);
        SearchBudget budget(100);
        SearchAccount unit(&budget, size_t(0), kPair);
        unit.charge(50);
        budget.answer().charge(90);
        EXPECT_THROW(budget.answer().charge(11), SearchAnswerTooLarge);
        EXPECT_TRUE(budget.passed());
        EXPECT_EQ("answer", std::string(budget.account()));
        // pending: the answer refuses, the units' values go on
        EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
        unit.charge(50);
        SearchAccount later(&budget, 1, kPair);
        later.charge(100);
        if (values_pass) {
            EXPECT_THROW(later.charge(1), SearchAnswerTooLarge);
            EXPECT_EQ("values", std::string(budget.account()));
            EXPECT_EQ("the values account of record 1 on pair (a, a.dbg, a.annodbg)",
                      budget.first_passed());
            EXPECT_EQ(100u, budget.first_held());
            EXPECT_EQ(1u, budget.first_charge());
            // and now every charge stops
            EXPECT_THROW(unit.charge(1), SearchAnswerTooLarge);
        } else {
            EXPECT_EQ("answer", std::string(budget.account()));
            EXPECT_EQ("the answer account", budget.first_passed());
            EXPECT_EQ(90u, budget.first_held());
            EXPECT_EQ(11u, budget.first_charge());
        }
    }
    // the object check of a unit (its own object past the bound): the answer side too
    SearchBudget budget(100);
    budget.note_object_passed("record 2 on the graph", 101);
    EXPECT_TRUE(budget.passed());
    EXPECT_EQ("answer", std::string(budget.account()));
    EXPECT_EQ("the answer account of record 2 on the graph", budget.first_passed());
    EXPECT_EQ(0u, budget.first_held());
    EXPECT_EQ(101u, budget.first_charge());
    SearchAccount unit(&budget, 3, "the graph");
    unit.charge(100);
    EXPECT_THROW(budget.answer().charge(1), SearchAnswerTooLarge);
}

// Test 3 (edges): a unit charged exactly the bound is answered and a byte more refused; the same
// for the request's account; a charge of UINT64_MAX saturates (nothing held); the saturating
// helpers
TEST(SearchBudget, TheEdgesOfTheUnitAndTheRequest) {
    for (bool partial : { false, true }) {
        SCOPED_TRACE(partial);
        const uint64_t bound = 1 << 20;
        SearchBudget budget(bound, partial);
        SearchAccount exact(&budget, size_t(0), kPair);
        exact.charge(bound - 7);
        exact.charge(7);
        EXPECT_FALSE(budget.passed());
        SearchAccount over(&budget, 1, kPair);
        over.charge(bound);
        if (partial) {
            EXPECT_THROW(over.charge(1), SearchUnitRefused);
            EXPECT_FALSE(budget.passed());
        } else {
            EXPECT_THROW(over.charge(1), SearchAnswerTooLarge);
            EXPECT_TRUE(budget.passed());
        }
        SearchBudget request(bound, partial);
        request.answer().charge(bound);
        EXPECT_FALSE(request.passed());
        EXPECT_THROW(request.answer().charge(1), SearchAnswerTooLarge);
        EXPECT_TRUE(request.passed());
        EXPECT_EQ("answer", std::string(request.account()));
        // the largest charge, on a unit and as a scratch
        SearchBudget huge(kLarge, partial);
        SearchAccount unit(&huge, 2, kPair);
        unit.charge(5);
        if (partial) {
            EXPECT_THROW(unit.charge(std::numeric_limits<uint64_t>::max()), SearchUnitRefused);
            ASSERT_EQ(1u, huge.unit_refusals().size());
            EXPECT_EQ(5u, huge.unit_refusals()[0].held);
            EXPECT_EQ(std::numeric_limits<uint64_t>::max(), huge.unit_refusals()[0].charge);
        } else {
            EXPECT_THROW(unit.charge(std::numeric_limits<uint64_t>::max()), SearchAnswerTooLarge);
            EXPECT_EQ(std::numeric_limits<uint64_t>::max(), huge.first_charge());
        }
        EXPECT_EQ(5u, unit.held());
    }
    EXPECT_EQ(3u, graph::saturating_add(1, 2));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              graph::saturating_add(std::numeric_limits<uint64_t>::max(), 1));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              graph::saturating_add(std::numeric_limits<uint64_t>::max() - 3, 4));
    EXPECT_EQ(0u, graph::mib_rounded_up(0));
    EXPECT_EQ(1u, graph::mib_rounded_up(1));
    EXPECT_EQ(1u, graph::mib_rounded_up(1 << 20));
    EXPECT_EQ(2u, graph::mib_rounded_up((1 << 20) + 1));
    EXPECT_EQ(uint64_t(1) << 44, graph::mib_rounded_up(std::numeric_limits<uint64_t>::max()));
}

// Test 2 (determinism): units with fixed charge sequences and object prices, a request's text,
// run as the server runs them (a unit's values charged as its search goes; at its end, unless
// the request is refused, its object checked against the bound and charged, or its stub's;
// after every unit, the text) in 500 orders drawn from a seed, single-threaded, in both modes:
// the status, the 400's account and the units refused (with their accounts and priced_mb) are
// the same in every order; and the same on 8 real threads
struct Unit {
    std::vector<uint64_t> charges;
    uint64_t object;  // its object's own price (tree and transient)
};
struct Outcome {
    bool refused = false;
    std::string account;
    std::map<size_t, std::string> units;  // the units refused alone: their accounts
    bool operator==(const Outcome &o) const {
        return std::tie(refused, account, units) == std::tie(o.refused, o.account, o.units);
    }
};
constexpr uint64_t kStub = 13;

// one unit's end, as process_search_request's callback: the object check, then its charge
void finish_unit(SearchBudget &budget, std::mutex &mutex, size_t u, const Unit &unit,
                 bool values_refused, std::map<size_t, std::string> *refused) {
    std::lock_guard<std::mutex> lock(mutex);
    if (budget.passed())
        return;
    uint64_t price = unit.object;
    if (values_refused) {
        (*refused)[u] = "values";
        price = kStub;
    } else if (unit.object > budget.max_bytes()) {
        const std::string where = "record " + std::to_string(u) + " on the graph";
        if (!budget.partial_refusal()) {
            budget.note_object_passed(where, unit.object);
            return;
        }
        budget.note_unit_refused({ u, where, "answer", 0, unit.object });
        (*refused)[u] = "answer " + std::to_string(graph::mib_rounded_up(unit.object));
        price = kStub;
    }
    try {
        budget.answer().charge(price);
    } catch (const SearchAnswerTooLarge &) {}
}

Outcome outcome_of(SearchBudget &budget, const std::map<size_t, std::string> &refused,
                   uint64_t text) {
    if (!budget.passed()) {
        try {
            budget.answer().charge(text);
        } catch (const SearchAnswerTooLarge &) {}
    }
    Outcome outcome;
    outcome.refused = budget.passed();
    if (outcome.refused) {
        outcome.account = budget.account();
    } else {
        outcome.units = refused;
    }
    return outcome;
}

Outcome run_in_order(const std::vector<Unit> &units, uint64_t bound, bool partial,
                     uint64_t text, std::mt19937 *order) {
    SearchBudget budget(bound, partial);
    std::mutex mutex;
    std::map<size_t, std::string> refused;
    std::vector<std::unique_ptr<SearchAccount>> accounts;
    std::vector<size_t> next(units.size(), 0);
    std::vector<bool> done(units.size(), false), values_refused(units.size(), false);
    for (size_t u = 0; u < units.size(); ++u) {
        accounts.push_back(std::make_unique<SearchAccount>(&budget, u, "the graph"));
    }
    std::vector<size_t> active(units.size());
    std::iota(active.begin(), active.end(), 0);
    while (!active.empty()) {
        const size_t i = order ? (*order)() % active.size() : 0;
        const size_t u = active[i];
        bool end = next[u] == units[u].charges.size();
        if (!end) {
            try {
                accounts[u]->charge(units[u].charges[next[u]++]);
            } catch (const SearchUnitRefused &) {
                values_refused[u] = true;
                end = true;
            } catch (const SearchAnswerTooLarge &) {
                // the request refused: this unit stops
                active.erase(active.begin() + i);
                continue;
            }
        }
        if (end) {
            finish_unit(budget, mutex, u, units[u], values_refused[u], &refused);
            active.erase(active.begin() + i);
        }
    }
    return outcome_of(budget, refused, text);
}

Outcome run_on_threads(const std::vector<Unit> &units, uint64_t bound, bool partial,
                       uint64_t text, size_t num_threads) {
    SearchBudget budget(bound, partial);
    std::mutex mutex;
    std::map<size_t, std::string> refused;
    std::atomic<size_t> next_unit = 0;
    std::vector<std::thread> threads;
    for (size_t t = 0; t < num_threads; ++t) {
        threads.emplace_back([&]() {
            for (size_t u; (u = next_unit++) < units.size(); ) {
                SearchAccount account(&budget, u, "the graph");
                bool values_refused = false;
                try {
                    for (uint64_t c : units[u].charges) {
                        account.charge(c);
                    }
                } catch (const SearchUnitRefused &) {
                    values_refused = true;
                } catch (const SearchAnswerTooLarge &) {
                    continue;
                }
                finish_unit(budget, mutex, u, units[u], values_refused, &refused);
            }
        });
    }
    for (auto &thread : threads) {
        thread.join();
    }
    return outcome_of(budget, refused, text);
}

TEST(SearchBudget, TheUnitsRefusedAndTheAccountAreTheSameInEveryOrder) {
    std::mt19937 rng(71);
    size_t refused_requests[2] = { 0, 0 }, partly[2] = { 0, 0 }, values_named = 0,
           answer_named = 0;
    for (size_t instance = 0; instance < 60; ++instance) {
        const uint64_t bound = 1000;
        std::vector<Unit> units(2 + rng() % 10);
        uint64_t objects = 0;
        for (Unit &unit : units) {
            unit.charges.resize(1 + rng() % 12);
            for (uint64_t &c : unit.charges) {
                c = rng() % (instance % 3 == 0 ? 300 : 120);
            }
            // now and then an object past the bound alone
            unit.object = rng() % 7 == 0 ? bound + 1 + rng() % 500 : rng() % 250;
            objects += unit.object;
        }
        const uint64_t text = instance % 4 == 0 ? rng() % 900 : rng() % 100;
        for (bool partial : { false, true }) {
            SCOPED_TRACE(::testing::Message() << "instance " << instance << " partial "
                         << partial);
            const Outcome expected = run_in_order(units, bound, partial, text, nullptr);
            std::mt19937 order(1000 + instance);
            for (size_t schedule = 0; schedule < 500; ++schedule) {
                ASSERT_EQ(expected, run_in_order(units, bound, partial, text, &order))
                        << "schedule " << schedule;
            }
            for (size_t round = 0; round < 5; ++round) {
                ASSERT_EQ(expected, run_on_threads(units, bound, partial, text, 8))
                        << "threads, round " << round;
            }
            refused_requests[partial] += expected.refused;
            partly[partial] += !expected.units.empty();
            if (expected.refused && !partial)
                (expected.account == "values" ? values_named : answer_named)++;
        }
    }
    // the instances reach every outcome
    EXPECT_LT(0u, refused_requests[0]);
    EXPECT_LT(0u, refused_requests[1]);
    EXPECT_LT(0u, partly[1]);
    EXPECT_EQ(0u, partly[0]);
    EXPECT_LT(0u, values_named);
    EXPECT_LT(0u, answer_named);
}

} // namespace
} // namespace mtg::test
