#include <algorithm>
#include <atomic>
#include <limits>
#include <map>
#include <random>
#include <span>
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

// The fixture's matrix without row-diff, recording the rows its batched read is asked for
class CountingRows : public annot::matrix::TupleCSCMatrix<annot::matrix::ColumnMajor> {
  public:
    using TupleCSCMatrix::TupleCSCMatrix;

    std::vector<RowTuples> get_row_tuples(std::span<const Row> rows,
                                          size_t num_threads = 1) const override {
        requested.insert(requested.end(), rows.begin(), rows.end());
        return TupleCSCMatrix::get_row_tuples(rows, num_threads);
    }

    mutable std::vector<Row> requested;
};

// Train 7 P1's R8 (a), lane T7-M5. A query by record folds its rows into the records it
// meets; past the fold's budget (kFoldValuesBudget: 16 MiB of runs in a release build, none
// in a debug one), with more records still in the race than the cut returns, the fold keeps
// only the counts and reads the rows again after the cut for the returned records (master's
// #765). The runs that read builds are values of the returned records: they are charged to
// the pair's account (as the task's scratch) before they grow, the values built from them as
// they grow. Here 1,200 records share a stretch of 1,000 k-mers, which the query holds each at
// a position of its own: 1.2 M runs before the cut. With an account the answer is the one
// without; the account's scratch holds at least the runs read again (16 B each: position,
// length, value); a bound below what the query needs, and one below those runs alone, refuse
// it with SearchAnswerTooLarge.
TEST(SearchBudget, AQueryByRecordReadAgainAfterTheCutIsCharged) {
    std::mt19937 rng(765);
    const size_t k = 11;
    const size_t num_shared = 1000;
    const std::string shared = random_sequence(rng, num_shared + k - 1);
    std::vector<CoordinateRecord> records;
    for (size_t i = 0; i < 1200; ++i) {
        records.push_back({ random_sequence(rng, 10) + shared + random_sequence(rng, 10),
                            "F" + std::to_string(i % 2), "r" + std::to_string(i) });
    }
    std::string query;
    for (size_t i = 0; i < num_shared; ++i) {
        query += shared.substr(i, k) + "N";
    }
    auto index = coordinate_test_graph<CountingRows>(k, records, false, false, false, true);
    const auto &matrix = dynamic_cast<const CountingRows &>(index->get_annotator().get_matrix());
    const size_t top = 10;
    auto check = [&](const std::string &mode, const auto &query_with) {
        SCOPED_TRACE(mode);
        matrix.requested.clear();
        const auto plain = query_with(nullptr);
        ASSERT_EQ(top, plain.size());
        // every row of the query was read twice: the fold kept only the counts
        std::map<annot::matrix::BinaryMatrix::Row, size_t> requests;
        for (auto row : matrix.requested) {
            ++requests[row];
        }
        std::map<annot::matrix::BinaryMatrix::Row, size_t> twice;
        for (auto node : index->get_graph().node_range(query)) {
            if (node)
                twice[graph::AnnotatedDBG::graph_to_anno_index(node)] = 2;
        }
        ASSERT_LT(num_shared * 9 / 10, twice.size());
        EXPECT_EQ(twice, requests);
        SearchBudget budget(kLarge);
        SearchAccount account(&budget, "values");
        EXPECT_EQ(plain, query_with(&account));
        const uint64_t runs_read_again = top * num_shared * 16;
        EXPECT_LE(runs_read_again, account.scratch());
        EXPECT_LT(account.scratch(), account.held());
        for (uint64_t bound : { account.held() - 1, runs_read_again / 2 }) {
            SearchBudget short_one(bound);
            SearchAccount refused(&short_one, "values");
            EXPECT_THROW(query_with(&refused), SearchAnswerTooLarge) << bound;
            EXPECT_TRUE(short_one.passed());
            EXPECT_GE(bound, refused.held());
        }
    };
    check("signatures", [&](SearchAccount *a) -> Signatures {
        return index->get_top_label_signatures(query, top, 0, 0, nullptr, a);
    });
    check("counts", [&](SearchAccount *a) -> Counts {
        return index->get_kmer_counts(query, top, 0, 0, nullptr, a);
    });
    check("coordinates", [&](SearchAccount *a) -> Coordinates {
        return index->get_kmer_coordinates(query, top, 0, 0, nullptr, a);
    });
}

// The fold of a query by record past its budget (a release build's 16 MiB of runs) frees the
// values of the records that can no longer be returned, and holds none for the records first
// met too late to be returned. The answers are those of a fold that frees nothing, with three
// of the five records tied at the cut, in the C5 order. Without a budget, the fold keeps only
// the counts and reads the rows again. (Review M2, lane T7-M5: a release build reaches every
// branch of the fold's pruning.)
TEST(SearchBudget, TheFoldFreesOnlyTheValuesOfRecordsThatCannotBeReturned) {
    std::mt19937 rng(764);
    const size_t k = 15;
    const size_t num_tied = 5, num_low = 1100, num_late = 10;
    const size_t num_repeats = 300, num_shared = 1000, num_late_kmers = 100;
    const std::string repeated = random_sequence(rng, k);
    const std::string shared = random_sequence(rng, num_shared + k - 1);
    const std::string late = random_sequence(rng, num_late_kmers + k - 1);
    std::vector<CoordinateRecord> records;
    auto add = [&](const std::string &sequence, const std::string &header) {
        records.push_back({ random_sequence(rng, 10) + sequence + random_sequence(rng, 10),
                            "F" + std::to_string(records.size() % 2), header });
    };
    for (size_t i = 0; i < num_tied; ++i) {
        add(repeated + random_sequence(rng, 10) + shared, "h" + std::to_string(i));
    }
    for (size_t i = 0; i < num_low; ++i) {
        add(shared, "l" + std::to_string(i));
    }
    for (size_t i = 0; i < num_late; ++i) {
        add(late, "z" + std::to_string(i));
    }
    // |repeated| at many positions (one row, counted at its first position, values held for
    // the returned records only), then the shared and the late k-mers one by one: a run each
    std::string query;
    for (size_t i = 0; i < num_repeats; ++i) {
        query += repeated + "N";
    }
    for (size_t i = 0; i < num_shared; ++i) {
        query += shared.substr(i, k) + "N";
    }
    for (size_t i = 0; i < num_late_kmers; ++i) {
        query += late.substr(i, k) + "N";
    }
    // The runs pass the budget at the 949th shared k-mer. The records of the shared k-mers
    // only then match 949 positions of the 151 left, the tied ones 300 + 949: they are freed,
    // and the late ones, which can match 100 positions, are held without values.
    ASSERT_LT(FoldBudget::kRelease, (num_tied + num_low) * num_shared * 16);
    auto index = coordinate_test_graph(k, records, false, false, false, true);
    const size_t top = 3;
    auto check = [&](const std::string &mode, const auto &query_by) {
        SCOPED_TRACE(mode);
        graph::reset_record_fold_counts();
        const auto unbounded = [&]() {
            const FoldBudget fold_budget(std::numeric_limits<size_t>::max());
            return query_by();
        }();
        EXPECT_EQ(0u, graph::record_fold_counts().freed);
        EXPECT_EQ(0u, graph::record_fold_counts().met_late);
        EXPECT_EQ(0u, graph::record_fold_counts().counts_only);
        ASSERT_EQ(top, unbounded.size());
        for (size_t r = 0; r < top; ++r) {
            EXPECT_EQ("h" + std::to_string(r), std::get<0>(unbounded[r]));
            EXPECT_EQ(num_repeats + num_shared, std::get<1>(unbounded[r]));
        }
        {
            const FoldBudget fold_budget(FoldBudget::kRelease);
            EXPECT_EQ(unbounded, query_by());
            EXPECT_EQ(num_low, graph::record_fold_counts().freed);
            EXPECT_EQ(num_late, graph::record_fold_counts().met_late);
            EXPECT_EQ(0u, graph::record_fold_counts().counts_only);
        }
        graph::reset_record_fold_counts();
        {
            const FoldBudget fold_budget(0);
            EXPECT_EQ(unbounded, query_by());
            EXPECT_EQ(1u, graph::record_fold_counts().counts_only);
        }
    };
    check("signatures", [&]() -> Signatures {
        return index->get_top_label_signatures(query, top, 0, 0);
    });
    check("counts", [&]() -> Counts {
        return index->get_kmer_counts(query, top, 0, 0);
    });
    check("coordinates", [&]() -> Coordinates {
        return index->get_kmer_coordinates(query, top, 0, 0);
    });
}

// A query by record for coordinates at a homopolymer charges the coordinates of the returned
// records at its one repeated row to the values account's scratch: each record's buffer, grown
// by doubling from one, 16 B per coordinate (not per coordinate and position). (Review m1, lane
// T7-M5.)
TEST(SearchBudget, AQueryByRecordChargesTheValuesOfItsRepeatedRows) {
    std::mt19937 rng(16);
    const size_t k = 9;
    std::vector<CoordinateRecord> records;
    for (size_t i = 0; i < 6; ++i) {
        records.push_back({ random_sequence(rng, 12) + std::string(k + 20 * i + 3, 'A')
                                + random_sequence(rng, 12),
                            "F" + std::to_string(i % 2), "r" + std::to_string(i) });
    }
    const std::string query(k + 500, 'A');
    auto index = coordinate_test_graph(k, records, false, false, false, true);
    const size_t top = 4;
    SearchBudget budget(kLarge);
    SearchAccount account(&budget, "values");
    const Coordinates result = index->get_kmer_coordinates(query, top, 0, 0, nullptr, &account);
    ASSERT_EQ(top, result.size());
    uint64_t scratch = 0;
    uint64_t spread = 0;
    for (const auto &[label, count, positions, length] : result) {
        SCOPED_TRACE(label);
        EXPECT_EQ(501u, count);
        // the record's coordinates at the row of A^k, the same at every position
        const auto dense = positions.to_dense();
        const size_t num_coords = dense[0].size();
        ASSERT_LT(2u, num_coords);
        EXPECT_EQ(num_coords, dense[500].size());
        size_t capacity = 1;
        while (capacity < num_coords) {
            capacity *= 2;
        }
        scratch += annot::matrix::buffer_bytes(capacity, 16);
        spread += num_coords * 501 * 16;
    }
    EXPECT_EQ(scratch, account.scratch());
    EXPECT_LT(account.scratch() * 100, spread);
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

} // namespace
} // namespace mtg::test
