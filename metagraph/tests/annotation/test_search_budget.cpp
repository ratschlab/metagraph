#include <algorithm>
#include <atomic>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
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
// values of the records that can no longer be returned and keeps those of a record at the
// edge: its count and every row left reach the cutoff exactly, and it then ties at the cut and
// wins by its header's bytes. It holds no values for the records first met too late to be
// returned. The answers are those of a fold that frees nothing; the counts at the cut are
// distinct. Without a budget, the fold keeps only the counts and reads the rows again.
// (Reviews M2 and A-m2, lane T7-M5: a release build reaches every branch of the fold's pruning
// and its edge; freeing at the cutoff, or cutting at the top count instead of the N-th, fails.)
TEST(SearchBudget, TheFoldFreesOnlyTheValuesOfRecordsThatCannotBeReturned) {
    std::mt19937 rng(764);
    const size_t k = 15;
    const size_t num_low = 1100, num_late = 10, num_shared = 1000, num_late_kmers = 100;
    const size_t top = 3;
    // The shared k-mers make a run each for the low records and for a, b and t: the runs
    // pass the budget at the shared k-mer |at_check|, the fold's one check.
    const size_t at_check = FoldBudget::kRelease / ((num_low + 3) * 16) + 1;
    ASSERT_LT(at_check, num_shared);
    const size_t remaining = num_shared - at_check + num_late_kmers;
    // Repeated k-mers, each counted at its first position for a, b, t and e. t matches the
    // shared k-mers up to the check and nothing after: its count there, the third, is the
    // cutoff. e matches only its repeated k-mer before the check and every row after it.
    const size_t w_a = 500, w_b = 300, w_t = 200;
    const size_t cutoff = w_t + at_check;
    const size_t w_e = cutoff - remaining;
    ASSERT_LT(num_shared + num_late_kmers, cutoff) << "the low records are freed";
    ASSERT_LT(cutoff, w_b + num_shared) << "b stays above t";
    const std::string x_a = random_sequence(rng, k), x_b = random_sequence(rng, k),
                      x_t = random_sequence(rng, k), x_e = random_sequence(rng, k);
    const std::string shared = random_sequence(rng, num_shared + k - 1);
    const std::string late = random_sequence(rng, num_late_kmers + k - 1);
    std::vector<CoordinateRecord> records;
    auto add = [&](const std::string &sequence, const std::string &header) {
        records.push_back({ random_sequence(rng, 10) + sequence + random_sequence(rng, 10),
                            "F" + std::to_string(records.size() % 2), header });
    };
    // a base other than |c|: t's and e's pieces of |shared| end and start with it, so that
    // the k-mer across the end or the start is no shared k-mer of the query
    auto other = [](char c) { return std::string(1, c == 'A' ? 'C' : 'A'); };
    add(x_a + random_sequence(rng, 10) + shared, "a");
    add(x_b + random_sequence(rng, 10) + shared, "b");
    add(x_t + random_sequence(rng, 10) + shared.substr(0, at_check + k - 1)
            + other(shared[at_check + k - 1]), "t");
    add(x_e + random_sequence(rng, 10) + other(shared[at_check - 1]) + shared.substr(at_check)
            + random_sequence(rng, 10) + late, "e");
    for (size_t i = 0; i < num_low; ++i) {
        add(shared, "l" + std::to_string(i));
    }
    for (size_t i = 0; i < num_late; ++i) {
        add(late, "z" + std::to_string(i));
    }
    // the repeated k-mers (one row each, counted at its first position, values held for the
    // returned records only), then the shared and the late k-mers one by one: a run each
    std::string query;
    for (const auto &[kmer, weight] : std::vector<std::pair<std::string, size_t>> {
                { x_a, w_a }, { x_b, w_b }, { x_t, w_t }, { x_e, w_e } }) {
        for (size_t i = 0; i < weight; ++i) {
            query += kmer + "N";
        }
    }
    for (size_t i = 0; i < num_shared; ++i) {
        query += shared.substr(i, k) + "N";
    }
    for (size_t i = 0; i < num_late_kmers; ++i) {
        query += late.substr(i, k) + "N";
    }
    auto index = coordinate_test_graph(k, records, false, false, false, true);
    // a, b, then e and t tied (e first by its header's bytes): e at the edge is returned
    const std::vector<std::pair<std::string, size_t>> expected {
        { "a", w_a + num_shared }, { "b", w_b + num_shared }, { "e", cutoff } };
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
            EXPECT_EQ(expected[r].first, std::get<0>(unbounded[r]));
            EXPECT_EQ(expected[r].second, std::get<1>(unbounded[r]));
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

// A record first met when the weight of the rows left equals the cutoff of the fold's last
// check keeps its values: it matches every row from there on, ties the N-th record at the cut
// and wins by its header's bytes. (Review M5-n1, lane T7-M5: the met-late edge; holding no
// values for a record met at the cutoff, `>` for `>=`, fails here.) A small budget makes the
// check come early, with more rows left than the cutoff.
TEST(SearchBudget, ARecordFirstMetAtTheCutoffKeepsItsValues) {
    std::mt19937 rng(1235);
    const size_t k = 15;
    const size_t budget = 16 << 10;
    const size_t top = 2;
    // The shared k-mers make a run each for a and t: the runs pass the budget at the shared
    // k-mer |at_check|. t's count there, the second, is the cutoff; t matches nothing after.
    const size_t at_check = budget / (2 * 16) + 1;
    const size_t w_a = 300, w_t = 100;
    const size_t cutoff = w_t + at_check;
    // m is first met at the shared k-mer |first_m|, where the rows left (from it on) are
    // the cutoff, and matches all of them; a matches the shared k-mers before it
    const size_t first_m = at_check + 174;
    const size_t num_shared = first_m + cutoff;
    // no second check: the runs stay within twice those at the check
    ASSERT_LE(num_shared + at_check, 4 * at_check);
    const std::string x_a = random_sequence(rng, k), x_t = random_sequence(rng, k);
    const std::string shared = random_sequence(rng, num_shared + k - 1);
    auto other = [](char c) { return std::string(1, c == 'A' ? 'C' : 'A'); };
    std::vector<CoordinateRecord> records {
        { random_sequence(rng, 10) + x_a + random_sequence(rng, 10)
              + shared.substr(0, first_m + k - 1) + other(shared[first_m + k - 1]), "F0", "a" },
        { random_sequence(rng, 10) + x_t + random_sequence(rng, 10)
              + shared.substr(0, at_check + k - 1) + other(shared[at_check + k - 1]), "F1", "t" },
        { random_sequence(rng, 10) + other(shared[first_m - 1]) + shared.substr(first_m)
              + random_sequence(rng, 10), "F0", "m" },
    };
    std::string query;
    for (const auto &[kmer, weight] : std::vector<std::pair<std::string, size_t>> {
                { x_a, w_a }, { x_t, w_t } }) {
        for (size_t i = 0; i < weight; ++i) {
            query += kmer + "N";
        }
    }
    for (size_t i = 0; i < num_shared; ++i) {
        query += shared.substr(i, k) + "N";
    }
    auto index = coordinate_test_graph(k, records, false, false, false, true);
    // a, then m and t tied at the cutoff (m first by its header's bytes)
    const std::vector<std::pair<std::string, size_t>> expected {
        { "a", w_a + first_m }, { "m", cutoff } };
    auto check = [&](const std::string &mode, const auto &query_by) {
        SCOPED_TRACE(mode);
        const auto unbounded = [&]() {
            const FoldBudget fold_budget(std::numeric_limits<size_t>::max());
            return query_by();
        }();
        ASSERT_EQ(top, unbounded.size());
        for (size_t r = 0; r < top; ++r) {
            EXPECT_EQ(expected[r].first, std::get<0>(unbounded[r]));
            EXPECT_EQ(expected[r].second, std::get<1>(unbounded[r]));
        }
        graph::reset_record_fold_counts();
        const FoldBudget fold_budget(budget);
        EXPECT_EQ(unbounded, query_by());
        EXPECT_EQ(0u, graph::record_fold_counts().freed);
        EXPECT_EQ(0u, graph::record_fold_counts().met_late);
        EXPECT_EQ(0u, graph::record_fold_counts().counts_only);
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

// A query by record for coordinates charges each returned record's coordinates at the positions
// the fold holds to the values account's scratch, as the (position, coordinate) pairs it sorts:
// one record at a time, 16 B a pair (lane T7-M5's fix round: the coordinate path's |held|)
TEST(SearchBudget, AQueryByRecordChargesTheCoordinatesItSorts) {
    std::mt19937 rng(17);
    const size_t k = 11;
    std::vector<CoordinateRecord> records;
    for (size_t i = 0; i < 8; ++i) {
        records.push_back({ random_sequence(rng, 40 + 30 * i), "F" + std::to_string(i % 2),
                            "r" + std::to_string(i) });
    }
    const std::string query = records[7].sequence + records[2].sequence
                              + records[5].sequence.substr(0, 60);
    auto index = coordinate_test_graph(k, records, false, false, false, true);
    // no k-mer repeated in the query: the fold holds the coordinates of every position
    std::set<graph::AnnotatedDBG::row_index> rows;
    size_t present = 0;
    for (auto node : index->get_graph().node_range(query)) {
        if (node) {
            rows.insert(graph::AnnotatedDBG::graph_to_anno_index(node));
            ++present;
        }
    }
    ASSERT_EQ(present, rows.size());
    const size_t top = 3;
    SearchBudget budget(kLarge);
    SearchAccount account(&budget, "values");
    const Coordinates result = index->get_kmer_coordinates(query, top, 0, 0, nullptr, &account);
    ASSERT_EQ(top, result.size());
    uint64_t scratch = 0;
    for (const auto &[label, count, positions, length] : result) {
        size_t num_coords = 0;
        for (const auto &coords : positions.to_dense()) {
            num_coords += coords.size();
        }
        scratch = std::max(scratch, annot::matrix::buffer_bytes(num_coords, 16));
    }
    EXPECT_LT(0u, scratch);
    EXPECT_EQ(scratch, account.scratch());
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

// An annotation of |num_labels| labels named with seven digits (so that their names sort as their
// ids, the order a column-label search puts rows of one count into its store) on which each
// label of |ids| covers every k-mer of |sequence|, and no other label any
std::unique_ptr<graph::AnnotatedDBG> labels_at(const std::string &sequence, size_t k,
                                               size_t num_labels,
                                               const std::vector<uint64_t> &ids) {
    auto graph = std::make_shared<graph::DBGHashOrdered>(k);
    graph->add_sequence(sequence);
    Vector<annot::matrix::BinaryMatrix::SetBitPositions> rows(graph->max_index());
    graph->call_nodes([&](auto node) {
        for (uint64_t id : ids) {
            rows[node - 1].push_back(id);
        }
    });
    std::vector<std::string> labels;
    labels.reserve(num_labels);
    for (size_t j = 0; j < num_labels; ++j) {
        std::string digits = std::to_string(j);
        labels.push_back("L" + std::string(7 - std::min<size_t>(7, digits.size()), '0') + digits);
    }
    return std::make_unique<graph::AnnotatedDBG>(
            graph, std::make_unique<annot::RowCompressed<>>(std::move(rows), labels));
}

// What a hash table's allocations hold as it is used (each at its size class, as the values
// account counts it): at their most, and, under |bound| (0: none), the first allocation that
// takes them past it and what they held before it
struct TableHeap {
    uint64_t bound = 0;
    uint64_t held = 0;
    uint64_t most = 0;
    uint64_t before_passing = 0;
    uint64_t passing = 0;
};

template <class T>
class TableHeapAllocator {
  public:
    using value_type = T;
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::false_type;

    explicit TableHeapAllocator(TableHeap *heap) noexcept : heap_(heap) {}
    template <class U>
    TableHeapAllocator(const TableHeapAllocator<U> &other) noexcept : heap_(other.heap()) {}

    T* allocate(size_t n) {
        const uint64_t bytes = alloc_bytes(n * sizeof(T));
        if (heap_->bound && !heap_->passing && heap_->held + bytes > heap_->bound) {
            heap_->before_passing = heap_->held;
            heap_->passing = bytes;
        }
        heap_->held += bytes;
        heap_->most = std::max(heap_->most, heap_->held);
        return std::allocator<T>().allocate(n);
    }
    void deallocate(T *p, size_t n) noexcept {
        heap_->held -= alloc_bytes(n * sizeof(T));
        std::allocator<T>().deallocate(p, n);
    }

    TableHeap* heap() const noexcept { return heap_; }
    template <class U>
    bool operator==(const TableHeapAllocator<U> &other) const noexcept {
        return heap_ == other.heap();
    }
    template <class U>
    bool operator!=(const TableHeapAllocator<U> &other) const noexcept {
        return heap_ != other.heap();
    }

  private:
    TableHeap *heap_;
};

// The heap of a column-label search's table (a store of one million labels or more) as the ids
// |ids| are put into it in their order, reserved for them, each with a value of |Values| over
// |num_kmers| positions: tsl's own table of the store's key and value
template <class Values>
TableHeap store_table_heap(const std::vector<uint64_t> &ids, size_t num_kmers,
                           uint64_t bound = 0) {
    using Column = annot::matrix::BinaryMatrix::Column;
    using Allocator = TableHeapAllocator<std::pair<Column, Values>>;
    TableHeap heap;
    heap.bound = bound;
    {
        Allocator allocator(&heap);
        tsl::hopscotch_map<Column, Values, std::hash<Column>, std::equal_to<Column>, Allocator>
                table(allocator);
        table.reserve(ids.size());
        for (uint64_t id : ids) {
            table.emplace(id, Values(num_kmers));
        }
    }
    return heap;
}

// GPT final review of trains 5a-7, M1: a column-label search's table (from one million labels on)
// whose label ids collide in its hash grows past its reserve as its keys are put in. Ids 256
// apart double its 4,096 buckets twice, 1.75 MiB, and 2.6 MiB at the second doubling, the old
// buckets beside the new; ids 2,048 apart leave it at its 512 buckets and add an overflow node an
// id. Each allocation is charged before it is made: the values account holds the most the table
// held at once, and a bound the growth passes refuses the search at the allocation that passes
// it, charged what the table held then and that allocation; the rows are those without an
// account. Before, the table was charged once its keys were in, after it had grown (a unit
// refused there had held 1.75 MiB under a bound of 1 MiB), and its doubling's peak never
TEST(SearchBudget, AColumnSearchChargesItsTableBeforeItGrows) {
    std::mt19937 rng(37);
    const std::string sequence = random_sequence(rng, 40);
    const size_t k = 11;
    const size_t num_kmers = sequence.size() - k + 1;
    for (const auto &[stride, num_rows] : { std::pair<uint64_t, size_t>{ 256, 3500 },
                                            std::pair<uint64_t, size_t>{ 2048, 400 } }) {
        SCOPED_TRACE(::testing::Message() << "ids " << stride << " apart, " << num_rows << " rows");
        std::vector<uint64_t> ids;
        for (size_t r = 0; r < num_rows; ++r) {
            ids.push_back(stride * r);
        }
        auto index = labels_at(sequence, k, 1'000'000, ids);
        const uint64_t reserve = graph::search_store_bytes<PresenceRuns>(1'000'000, num_rows);
        const TableHeap table = store_table_heap<PresenceRuns>(ids, num_kmers);
        // the table grows past its reserve
        ASSERT_LT(reserve, table.most);

        SearchBudget budget(kLarge);
        SearchAccount account(&budget, "values");
        const auto rows = index->get_top_label_signatures(sequence, num_rows, 0, 0, nullptr,
                                                          &account);
        ASSERT_EQ(num_rows, rows.size());
        EXPECT_EQ(index->get_top_label_signatures(sequence, num_rows, 0, 0), rows);
        // the store's scratch: the most the table held at once
        EXPECT_EQ(table.most, account.scratch());

        // a bound halfway from the reserve to that: refused at the allocation that passes it
        const uint64_t bound = (reserve + table.most) / 2;
        const TableHeap passing = store_table_heap<PresenceRuns>(ids, num_kmers, bound);
        ASSERT_LT(0u, passing.passing);
        SearchBudget small(bound);
        SearchAccount refused(&small, "values");
        EXPECT_THROW(index->get_top_label_signatures(sequence, num_rows, 0, 0, nullptr, &refused),
                     SearchAnswerTooLarge);
        ASSERT_TRUE(small.passed());
        EXPECT_LE(small.first_held(), bound);
        EXPECT_EQ(passing.before_passing + passing.passing,
                  small.first_held() + small.first_charge());

        if (stride == 256) {
            // a bound the grown table and the rows' values fit, but not the second doubling (the
            // old buckets beside the new): refused (before, answered)
            using Layout = annot::matrix::HopscotchLayout<utils::ValueStore<
                    annot::matrix::BinaryMatrix::Column, PresenceRuns>::SparseMap>;
            const uint64_t grown = Layout::array_bytes(16'384);
            ASSERT_LT(grown + num_rows * alloc_bytes(2 * sizeof(PresenceRuns::Run)), table.most);
            SearchBudget peak(table.most - 1);
            SearchAccount at_peak(&peak, "values");
            EXPECT_THROW(index->get_top_label_signatures(sequence, num_rows, 0, 0, nullptr,
                                                         &at_peak),
                         SearchAnswerTooLarge);
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

// Review m3 of ae3865df7: the units refused alone are listed (and logged) in the answer's
// order, by record and then by the pair's place in the answer (the order of the record's `cuts`
// entries), not by the pairs' names: a pair "zz" answered before a pair "aa" comes first
TEST(SearchBudget, TheUnitsRefusedAreInTheAnswersOrder) {
    SearchBudget budget(100, true);
    const std::string zz = "pair (zz, z.dbg, z.annodbg)";
    const std::string aa = "pair (aa, a.dbg, a.annodbg)";
    // refused in another order than the answer's: record 1 first, then aa before zz
    SearchAccount r1_zz(&budget, 1, zz, 0);
    SearchAccount r0_aa(&budget, size_t(0), aa, 1);
    SearchAccount r0_zz(&budget, size_t(0), zz, 0);
    EXPECT_EQ(1u, r0_aa.pair_index());
    for (SearchAccount *account : { &r1_zz, &r0_aa, &r0_zz }) {
        EXPECT_THROW(account->charge(101), SearchUnitRefused);
    }
    // an object refused on the answer side, noted with its pair's place (the same unit noted
    // twice, which the engine never does, orders by its account: the order is total)
    budget.note_unit_refused({ 0, "record 0 on " + aa, "answer", 0, 500, 1 });
    std::vector<std::tuple<size_t, size_t, std::string, std::string>> listed;
    for (const auto &refusal : budget.unit_refusals()) {
        listed.emplace_back(refusal.record, refusal.pair, refusal.where, refusal.account);
    }
    using Row = std::tuple<size_t, size_t, std::string, std::string>;
    EXPECT_EQ((std::vector<Row> {
                  Row { 0, 0, "record 0 on " + zz, "values" },
                  Row { 0, 1, "record 0 on " + aa, "answer" },
                  Row { 0, 1, "record 0 on " + aa, "values" },
                  Row { 1, 0, "record 1 on " + zz, "values" } }),
              listed);
    // a unit made without a pair's place (a single-graph server) is at place 0
    SearchAccount single(&budget, 2, "the graph");
    EXPECT_EQ(0u, single.pair_index());
}

// Review untested claim 8 of ae3865df7: real threads charging units and the request's answer
// at once (the server's pairs charge the answer account concurrently; a pair's own objects one
// at a time): the units refused alone are exactly those whose own charges pass the bound, the
// answer account holds exactly what was charged, whatever the interleaving; without the field
// the request is refused "values" exactly when some unit passes
TEST(SearchBudget, ManyThreadsChargeUnitsAndTheAnswerAtOnce) {
    std::mt19937 rng(83);
    const uint64_t bound = 10'000;
    for (size_t round = 0; round < 20; ++round) {
        std::vector<std::vector<uint64_t>> charges(500 + rng() % 1500);
        std::set<size_t> expected;
        for (size_t u = 0; u < charges.size(); ++u) {
            charges[u].resize(1 + rng() % 40);
            uint64_t sum = 0;
            bool passes = false;
            for (uint64_t &c : charges[u]) {
                c = rng() % (round % 2 ? 600 : 300);
                if (!passes && sum + c > bound)
                    passes = true;
                sum += c;
            }
            if (passes)
                expected.insert(u);
        }
        for (bool partial : { true, false }) {
            SCOPED_TRACE(::testing::Message() << "round " << round << " partial " << partial);
            // an answer bound the objects never reach: it holds exactly their sum
            SearchBudget budget(partial ? bound : std::numeric_limits<uint64_t>::max() / 2,
                                partial);
            const uint64_t unit_bound = budget.max_bytes();
            std::atomic<size_t> next = 0;
            std::atomic<uint64_t> answered = 0;
            std::vector<std::thread> threads;
            std::mutex mutex;
            std::set<size_t> refused;
            for (size_t t = 0; t < 16; ++t) {
                threads.emplace_back([&, t]() {
                    for (size_t u; (u = next++) < charges.size(); ) {
                        SearchAccount account(&budget, u, "the graph", t % 3);
                        try {
                            for (uint64_t c : charges[u]) {
                                account.charge(c);
                                if (c % 7 == 0)
                                    std::this_thread::yield();
                            }
                        } catch (const SearchUnitRefused &) {
                            std::lock_guard<std::mutex> lock(mutex);
                            refused.insert(u);
                            continue;
                        } catch (const SearchAnswerTooLarge &) {
                            continue;
                        }
                        // the unit's object, on the shared answer account (no mutex; at most
                        // 3 x 2,000 B in all, under the bound)
                        if (partial) {
                            const uint64_t object = 1 + u % 3;
                            try {
                                budget.answer().charge(object);
                                answered += object;
                            } catch (const SearchAnswerTooLarge &) {
                                answered += uint64_t(1) << 40;  // fails the test below
                            }
                        }
                    }
                });
            }
            for (auto &thread : threads) {
                thread.join();
            }
            if (partial) {
                EXPECT_FALSE(budget.passed());
                EXPECT_EQ(expected, refused);
                EXPECT_EQ(expected.size(), budget.unit_refusals().size());
                EXPECT_EQ(answered.load(), budget.answer().held());
                size_t previous = 0;
                for (const auto &refusal : budget.unit_refusals()) {
                    EXPECT_LE(previous, refusal.record);
                    previous = refusal.record;
                    EXPECT_LE(refusal.held, unit_bound);
                    EXPECT_LT(unit_bound - refusal.held, refusal.charge);
                }
            } else {
                // the bound is out of reach: nothing refused
                EXPECT_FALSE(budget.passed());
                EXPECT_TRUE(refused.empty());
            }
        }
        // without the field at the units' bound: refused "values" exactly when one passes
        SearchBudget strict(bound, false);
        std::atomic<size_t> next = 0;
        std::vector<std::thread> threads;
        for (size_t t = 0; t < 16; ++t) {
            threads.emplace_back([&]() {
                for (size_t u; (u = next++) < charges.size(); ) {
                    SearchAccount account(&strict, u, "the graph");
                    try {
                        for (uint64_t c : charges[u]) {
                            account.charge(c);
                        }
                    } catch (const SearchAnswerTooLarge &) {}
                }
            });
        }
        for (auto &thread : threads) {
            thread.join();
        }
        EXPECT_EQ(!expected.empty(), strict.passed()) << round;
        if (strict.passed()) {
            EXPECT_EQ("values", std::string(strict.account()));
        }
    }
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
