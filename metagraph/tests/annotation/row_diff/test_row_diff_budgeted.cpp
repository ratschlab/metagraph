// The budget-aware row-diff decode (DESIGN-traverse-graphlet.md §14, stage 3 of §14.1):
// IRowDiff::decode_rows / decode_row_tuples against the default decode, the costs it
// states per row, its refusals, and its byte model against the allocator.
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

#include <gtest/gtest.h>

#include "../test_annotated_dbg_helpers.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"
#include "annotation/int_matrix/row_diff/tuple_row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#if USE_JEMALLOC
extern "C" int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen);
#endif


namespace {

using namespace mtg;
using namespace mtg::annot;
using namespace mtg::annot::matrix;
using Row = BinaryMatrix::Row;
using Column = BinaryMatrix::Column;
using SetBits = BinaryMatrix::SetBitPositions;
using RowTuples = MultiIntMatrix::RowTuples;
using ColRank = std::pair<Column, uint64_t>;

std::string random_seq(size_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char &c : s) {
        c = "ACGT"[rng() % 4];
    }
    return s;
}

// a column of |rows| bits, each set with probability |density|
std::unique_ptr<bit_vector> random_column(uint64_t rows, double density, std::mt19937 &rng) {
    sdsl::bit_vector bv(rows, 0);
    std::bernoulli_distribution set(density);
    for (uint64_t i = 0; i < rows; ++i) {
        bv[i] = set(rng);
    }
    return std::make_unique<bit_vector_sd>(bv);
}

std::vector<std::unique_ptr<bit_vector>> copy_columns(const std::vector<std::unique_ptr<bit_vector>> &cols) {
    std::vector<std::unique_ptr<bit_vector>> out;
    for (const auto &c : cols) {
        out.push_back(c->copy());
    }
    return out;
}

template <class T>
T sorted(T v) {
    std::sort(v.begin(), v.end());
    return v;
}

// L1: the single-row descents return each row's set bits (and column ranks) as the batched
// reads do, for random matrices of several arities, including empty and full rows
TEST(RowDiffBudgetedDecode, DescentEqualsBatchedRead) {
    std::mt19937 rng(42);
    for (size_t arity : { 2, 3, 5 }) {
        for (double density : { 0.0, 0.05, 0.5, 1.0 }) {
            const uint64_t num_rows = 300;
            std::vector<std::unique_ptr<bit_vector>> cols;
            for (size_t j = 0; j < 37; ++j) {
                cols.push_back(random_column(num_rows, density, rng));
            }
            ColumnMajor column_major(copy_columns(cols));
            BRWT brwt = BRWTBottomUpBuilder::build(std::move(cols),
                                                   BRWTBottomUpBuilder::get_basic_partitioner(arity));
            std::vector<Row> all(num_rows);
            std::iota(all.begin(), all.end(), 0);
            const auto rows = brwt.get_rows(all);
            const auto ranks = brwt.get_column_ranks(all, 1);
            const auto cm_ranks = column_major.get_column_ranks(all, 1);
            DecodeBudget budget;
            ChargedBuffer<Column> scratch;
            ChargedBuffer<ColRank> rank_scratch;
            for (Row r = 0; r < num_rows; ++r) {
                const std::string what = "arity " + std::to_string(arity) + " density "
                        + std::to_string(density) + " row " + std::to_string(r);
                scratch.data.clear();
                ASSERT_TRUE(brwt.row_columns(r, &scratch, budget)) << what;
                EXPECT_EQ(sorted(std::vector<Column>(rows[r].begin(), rows[r].end())),
                          sorted(std::vector<Column>(scratch.data.begin(), scratch.data.end()))) << what;
                rank_scratch.data.clear();
                ASSERT_TRUE(brwt.row_column_ranks(r, &rank_scratch, budget)) << what;
                EXPECT_EQ(sorted(std::vector<ColRank>(ranks[r].begin(), ranks[r].end())),
                          sorted(std::vector<ColRank>(rank_scratch.data.begin(), rank_scratch.data.end()))) << what;
                scratch.data.clear();
                ASSERT_TRUE(column_major.row_columns(r, &scratch, budget)) << what;
                EXPECT_EQ(sorted(std::vector<Column>(rows[r].begin(), rows[r].end())),
                          std::vector<Column>(scratch.data.begin(), scratch.data.end())) << what;
                rank_scratch.data.clear();
                ASSERT_TRUE(column_major.row_column_ranks(r, &rank_scratch, budget)) << what;
                EXPECT_EQ(std::vector<ColRank>(cm_ranks[r].begin(), cm_ranks[r].end()),
                          std::vector<ColRank>(rank_scratch.data.begin(), rank_scratch.data.end())) << what;
            }
            // the scratch is all that was charged, and only as it grew
            EXPECT_EQ(scratch.charged() + rank_scratch.charged(), budget.held());
        }
    }
}

// The row-diff annotations the decode is checked on: one built by the transform pipeline
// (its anchors and fork successors), in all four budget-aware formats — over ColumnMajor
// and over a BRWT made of the same diff columns, with and without coordinates.
struct Annotation {
    std::string name;
    std::shared_ptr<graph::AnnotatedDBG> anno;     // owns the graph (and the ColumnMajor one)
    std::unique_ptr<BinaryMatrix> owned;            // the BRWT variants
    const IRowDiff *rd = nullptr;
    const BinaryMatrix *binary = nullptr;
    const MultiIntMatrix *tuples = nullptr;
};

void copy_support(const IRowDiff &from, IRowDiff *to) {
    utils::TempFile anchors, fork_succ;
    {
        std::ofstream a(anchors.name(), std::ios::binary);
        from.anchor().serialize(a);
        std::ofstream f(fork_succ.name(), std::ios::binary);
        from.fork_succ().serialize(f);
    }
    to->load_anchor(anchors.name());
    to->load_fork_succ(fork_succ.name());
}

std::vector<Annotation> annotations() {
    std::vector<std::string> seqs, labels;
    // labels sharing long stretches, so that rows carry several labels and diffs cancel
    std::vector<std::string> base { random_seq(160, 1), random_seq(160, 2), random_seq(160, 3) };
    for (uint32_t i = 0; i < 18; ++i) {
        seqs.push_back(random_seq(20, 100 + i) + base[i % 3].substr(i % 7, 120) + random_seq(25, 200 + i));
        labels.push_back("L" + std::to_string(i));
    }
    std::vector<Annotation> out;
    for (bool coords : { false, true }) {
        std::shared_ptr<graph::AnnotatedDBG> anno = test::build_anno_graph<graph::DBGSuccinct,
                RowDiffColumnAnnotator>(9, seqs, labels, graph::DeBruijnGraph::BASIC, coords);
        const BinaryMatrix &m = anno->get_annotator().get_matrix();
        Annotation cm;
        cm.name = coords ? "row_diff_coord" : "row_diff";
        cm.anno = anno;
        cm.rd = dynamic_cast<const IRowDiff*>(&m);
        cm.binary = &m;
        cm.tuples = dynamic_cast<const MultiIntMatrix*>(&m);
        EXPECT_TRUE(cm.rd && cm.rd->supports_budgeted_decode()) << cm.name;
        Annotation br;
        br.name = coords ? "row_diff_brwt_coord" : "row_diff_brwt";
        br.anno = anno;
        if (!coords) {
            const auto &rd = dynamic_cast<const RowDiff<ColumnMajor>&>(m);
            BRWT brwt = BRWTBottomUpBuilder::build(copy_columns(rd.diffs().data()),
                                                   BRWTBottomUpBuilder::get_basic_partitioner(3));
            auto rdb = std::make_unique<RowDiff<BRWT>>(rd.graph(), std::move(brwt));
            copy_support(rd, rdb.get());
            br.rd = rdb.get();
            br.binary = rdb.get();
            br.owned = std::move(rdb);
        } else {
            const auto &rd = dynamic_cast<const TupleRowDiff<TupleCSCMatrix<ColumnMajor>>&>(m);
            BRWT brwt = BRWTBottomUpBuilder::build(
                    copy_columns(rd.diffs().get_binary_matrix().data()),
                    BRWTBottomUpBuilder::get_basic_partitioner(2));
            TupleCSCMatrix<BRWT> tuples(std::move(brwt));
            std::stringstream ss;
            rd.diffs().serialize_tuples(ss);
            EXPECT_TRUE(tuples.load_tuples(ss));
            auto rdb = std::make_unique<TupleRowDiff<TupleCSCMatrix<BRWT>>>(rd.graph(), std::move(tuples));
            copy_support(rd, rdb.get());
            br.rd = rdb.get();
            br.binary = rdb.get();
            br.tuples = rdb.get();
            br.owned = std::move(rdb);
        }
        out.push_back(std::move(cm));
        out.push_back(std::move(br));
    }
    return out;
}

// the rows of the graph's nodes (a row of a dummy node has no row-diff path)
std::vector<Row> rows_of(const Annotation &a) {
    std::vector<Row> rows;
    a.rd->graph()->call_nodes([&](auto node) {
        rows.push_back(graph::AnnotatedDBG::graph_to_anno_index(node));
    });
    std::sort(rows.begin(), rows.end());
    return rows;
}

// The batches the decode is checked on: consecutive and scattered rows, duplicates, an
// empty batch and single rows
std::vector<std::vector<Row>> batches(const Annotation &a, uint32_t seed) {
    const std::vector<Row> valid = rows_of(a);
    std::mt19937 rng(seed);
    std::vector<std::vector<Row>> out { {} };
    for (size_t size : { 1, 2, 7, 33, 200 }) {
        for (size_t rep = 0; rep < 4; ++rep) {
            std::vector<Row> b;
            if (rep % 2) {
                const size_t from = rng() % valid.size();
                for (size_t i = 0; i < size; ++i) {
                    b.push_back(valid[(from + i) % valid.size()]);
                }
            } else {
                for (size_t i = 0; i < size; ++i) {
                    b.push_back(valid[rng() % valid.size()]);
                }
            }
            if (size > 2)
                b.push_back(b[size / 2]);     // a duplicate
            out.push_back(std::move(b));
        }
    }
    return out;
}

bool same_costs(const RowCost &a, const RowCost &b) {
    return a.dependency_rows == b.dependency_rows && a.dependency_entries == b.dependency_entries
        && a.demand == b.demand;
}

template <class RowT>
DecodeStatus decode(const Annotation &a, const std::vector<Row> &rows, DecodeBudget &budget,
                    std::vector<RowT> *out, std::vector<RowCost> *costs, std::vector<uint64_t> *held) {
    if constexpr(std::is_same_v<RowT, SetBits>) {
        return a.rd->decode_rows(rows, budget, out, costs, held);
    } else {
        return a.rd->decode_row_tuples(rows, budget, out, costs, held);
    }
}

template <class RowT>
std::vector<RowT> reference(const Annotation &a, const std::vector<Row> &rows) {
    if constexpr(std::is_same_v<RowT, SetBits>) {
        return a.binary->get_rows(rows);
    } else {
        return a.tuples->get_row_tuples(rows);
    }
}

// L2: with an unlimited budget the decode returns what the default path returns, and holds
// exactly its outputs afterwards
template <class RowT>
void check_equivalence(const Annotation &a) {
    for (const auto &rows : batches(a, 7)) {
        DecodeBudget budget;
        std::vector<RowT> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, budget, &out, &costs, &held)) << a.name;
        ASSERT_EQ(rows.size(), out.size());
        ASSERT_EQ(rows.size(), costs.size());
        EXPECT_EQ(reference<RowT>(a, rows), out) << a.name << " batch of " << rows.size();
        uint64_t kept = IRowDiff::output_bytes<RowT>(rows.size());
        for (uint64_t b : held) {
            kept += b;
        }
        EXPECT_EQ(kept, budget.held()) << a.name;
        for (size_t i = 0; i < rows.size(); ++i) {
            const bool anchor = a.rd->anchor()[rows[i]];
            EXPECT_EQ(anchor, costs[i].dependency_rows == 0) << a.name << " row " << rows[i];
            EXPECT_LE(budget.held() - kept + held[i], costs[i].demand) << a.name;
        }
    }
}

TEST(RowDiffBudgetedDecode, EqualsTheDefaultDecode) {
    for (const Annotation &a : annotations()) {
        check_equivalence<SetBits>(a);
        if (a.tuples)
            check_equivalence<RowTuples>(a);
    }
}

// L3: a row's costs are the same in every batch, decoding it alone stays within the demand
// it states, and a budget of exactly that peak admits it while one byte less refuses it
template <class RowT>
void check_batch_independence(const Annotation &a) {
    std::map<Row, RowCost> alone;
    size_t refused = 0;
    const std::vector<Row> valid = rows_of(a);
    for (size_t vi = 0; vi < valid.size(); vi += 3) {
        const Row r = valid[vi];
        DecodeBudget budget;
        std::vector<RowT> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, { r }, budget, &out, &costs, &held));
        EXPECT_LE(budget.peak(), costs[0].demand) << a.name << " row " << r;
        alone[r] = costs[0];
        // the exact admission
        DecodeBudget exact(budget.peak());
        out.clear(); costs.clear(); held.clear();
        EXPECT_EQ(DecodeStatus::OK, decode<RowT>(a, { r }, exact, &out, &costs, &held)) << r;
        DecodeBudget less(budget.peak() - 1);
        out.clear(); costs.clear(); held.clear();
        EXPECT_EQ(DecodeStatus::REFUSED, decode<RowT>(a, { r }, less, &out, &costs, &held)) << r;
        EXPECT_EQ(0u, less.held());
        refused++;
    }
    EXPECT_GT(refused, 10u);
    for (const auto &rows : batches(a, 11)) {
        DecodeBudget budget;
        std::vector<RowT> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, budget, &out, &costs, &held));
        for (size_t i = 0; i < rows.size(); ++i) {
            auto it = alone.find(rows[i]);
            if (it == alone.end()) {
                DecodeBudget b;
                std::vector<RowT> o;
                std::vector<RowCost> c;
                std::vector<uint64_t> h;
                ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, { rows[i] }, b, &o, &c, &h));
                it = alone.emplace(rows[i], c[0]).first;
            }
            EXPECT_TRUE(same_costs(it->second, costs[i])) << a.name << " row " << rows[i];
        }
    }
    // Batches [p, r] where p is on r's row-diff path: r's path then stops at p, decoded and
    // shared earlier in the call, so r continues from a copy of p's row. With an empty diff
    // after p the demand once took the copy's exact size instead of p's aggregate (review
    // F4: a key admitted on a warm cache and refused on a cold one).
    size_t pairs = 0;
    const graph::DeBruijnGraph &graph = *a.rd->graph();
    for (size_t vi = 0; vi < valid.size(); vi += 2) {
        const Row r = valid[vi];
        auto node = graph::AnnotatedDBG::anno_to_graph_index(r);
        for (size_t step = 0; step < 6 && !a.rd->anchor()[graph::AnnotatedDBG::graph_to_anno_index(node)]; ++step) {
            node = row_diff_successor(graph, node, a.rd->fork_succ());
            const Row p = graph::AnnotatedDBG::graph_to_anno_index(node);
            DecodeBudget budget;
            std::vector<RowT> out;
            std::vector<RowCost> costs;
            std::vector<uint64_t> held;
            ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, { p, r }, budget, &out, &costs, &held));
            auto it = alone.find(r);
            if (it == alone.end()) {
                DecodeBudget b;
                std::vector<RowT> o;
                std::vector<RowCost> c;
                std::vector<uint64_t> h;
                ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, { r }, b, &o, &c, &h));
                it = alone.emplace(r, c[0]).first;
            }
            EXPECT_TRUE(same_costs(it->second, costs[1]))
                << a.name << " row " << r << " after its path row " << p << ": demand "
                << costs[1].demand << " alone " << it->second.demand;
            pairs++;
        }
    }
    EXPECT_GT(pairs, 20u) << a.name;
}

TEST(RowDiffBudgetedDecode, CostsAreBatchIndependent) {
    for (const Annotation &a : annotations()) {
        check_batch_independence<SetBits>(a);
        if (a.tuples)
            check_batch_independence<RowTuples>(a);
    }
}

// L4: a refusal at any charge returns nothing, restores the budget and never lets the held
// bytes pass the maximum; the same budget then succeeds
template <class RowT>
void check_denials(const Annotation &a) {
    std::vector<Row> rows;
    const std::vector<Row> valid = rows_of(a);
    for (size_t i = 1; i < valid.size(); i += 17) {
        rows.push_back(valid[i]);
        rows.push_back(valid[i - 1]);
    }
    DecodeBudget probe;
    probe.charge(100);      // something held on entry
    std::vector<RowT> out;
    std::vector<RowCost> costs;
    std::vector<uint64_t> held;
    ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, probe, &out, &costs, &held));
    const uint64_t charges = probe.charges() - 1;
    ASSERT_GT(charges, 20u);
    const std::vector<RowT> expected = out;
    for (uint64_t n = 0; n < charges; n += std::max<uint64_t>(1, charges / 60)) {
        DecodeBudget budget(probe.peak() + 1000);
        budget.charge(100);
        budget.deny = [&](uint64_t ordinal) { return ordinal == n + 1; };
        std::vector<RowT> sentinel(1);
        std::vector<RowCost> c;
        std::vector<uint64_t> h;
        EXPECT_EQ(DecodeStatus::REFUSED, decode<RowT>(a, rows, budget, &sentinel, &c, &h))
                << a.name << " charge " << n;
        EXPECT_EQ(1u, sentinel.size());
        EXPECT_TRUE(c.empty() && h.empty());
        EXPECT_EQ(100u, budget.held()) << a.name << " charge " << n;
        EXPECT_LE(budget.peak(), budget.max_bytes());
        budget.deny = nullptr;
        std::vector<RowT> again;
        ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, budget, &again, &c, &h));
        EXPECT_EQ(expected, again);
    }
    // and every budget below what the batch needs refuses it, never holding more than given
    for (uint64_t max = 100; max < probe.peak(); max += std::max<uint64_t>(1, probe.peak() / 40)) {
        DecodeBudget budget(max);
        budget.charge(100);
        std::vector<RowT> o;
        std::vector<RowCost> c;
        std::vector<uint64_t> h;
        EXPECT_EQ(DecodeStatus::REFUSED, decode<RowT>(a, rows, budget, &o, &c, &h)) << max;
        EXPECT_EQ(100u, budget.held());
        EXPECT_LE(budget.peak(), max);
    }
}

TEST(RowDiffBudgetedDecode, RefusalsHaveNoSideEffects) {
    for (const Annotation &a : annotations()) {
        check_denials<SetBits>(a);
        if (a.tuples)
            check_denials<RowTuples>(a);
    }
}

// Formats without the budget-aware path say so, charge nothing and touch nothing
TEST(RowDiffBudgetedDecode, UnsupportedFormatsChargeNothing) {
    std::vector<std::string> seqs { random_seq(60, 5), random_seq(60, 6) };
    auto anno = test::build_anno_graph<graph::DBGSuccinct, RowDiffDiskAnnotator>(
            7, seqs, { "A", "B" }, graph::DeBruijnGraph::BASIC);
    const auto *rd = dynamic_cast<const IRowDiff*>(&anno->get_annotator().get_matrix());
    ASSERT_TRUE(rd);
    EXPECT_FALSE(rd->supports_budgeted_decode());
    DecodeBudget budget;
    std::vector<SetBits> out;
    std::vector<RowCost> costs;
    std::vector<uint64_t> held;
    EXPECT_EQ(DecodeStatus::UNSUPPORTED, rd->decode_rows({ 0, 1 }, budget, &out, &costs, &held));
    EXPECT_TRUE(out.empty());
    EXPECT_EQ(0u, budget.held());
    EXPECT_EQ(0u, budget.charges());
}

// what a compact copy of a row holds
uint64_t copy_bytes_of(const RowTuples &row) {
    uint64_t bytes = buffer_bytes(row.size(), sizeof(row[0]));
    for (const auto &entry : row) {
        bytes += small_vector_bytes(entry.second.size(), sizeof(uint64_t));
    }
    return bytes;
}

// L5 (§14 freeze gate, "a row-diff row whose dependencies are dense but whose result is
// tiny"): an anchor of 300 columns (or 6,000 coordinates) and a diff on the path that
// cancels all but one of them. The default decode returns the one-column row; the
// budget-aware decode states the dependencies, and a budget between the result's bytes
// and what decoding it holds refuses it whole.
TEST(RowDiffBudgetedDecode, DenseDependenciesTinyResult) {
    graph::DBGSuccinct graph(4);
    graph.add_sequence("ACTAGCTAGCTAGCTAGCTAGC");
    graph.add_sequence("ACTCTAG");
    const uint64_t num_rows = graph.max_index();
    // a row and the path to its anchor (the last row of the path)
    const Row requested = 3;
    std::vector<Row> path;
    sdsl::bit_vector anchors_bv(num_rows, 0);
    {
        graph::DeBruijnGraph::node_index node = graph::AnnotatedDBG::anno_to_graph_index(requested);
        bit_vector_small no_fork_succ;
        for (size_t step = 0; step < 3; ++step) {
            path.push_back(graph::AnnotatedDBG::graph_to_anno_index(node));
            node = row_diff_successor(graph, node, no_fork_succ);
        }
        anchors_bv[path.back()] = 1;
        // every other row is an anchor too, so that every path ends
        for (Row r = 0; r < num_rows; ++r) {
            if (std::find(path.begin(), path.end(), r) == path.end())
                anchors_bv[r] = 1;
        }
    }
    utils::TempFile anchors_file;
    {
        std::ofstream f(anchors_file.name(), std::ios::binary);
        IRowDiff::anchor_bv_type(anchors_bv).serialize(f);
    }
    const size_t num_columns = 300;
    // the anchor carries every column, the requested row's diff all but column 7: the row
    // is {7}, its dependencies 300 + 299 entries
    std::vector<std::unique_ptr<bit_vector>> cols;
    for (size_t j = 0; j < num_columns; ++j) {
        sdsl::bit_vector bv(num_rows, 0);
        bv[path.back()] = 1;
        if (j != 7)
            bv[path.front()] = 1;
        cols.push_back(std::make_unique<bit_vector_sd>(bv));
    }
    std::vector<std::unique_ptr<IRowDiff>> matrices;
    {
        auto rd = std::make_unique<RowDiff<ColumnMajor>>(&graph, ColumnMajor(copy_columns(cols)));
        rd->load_anchor(anchors_file.name());
        matrices.push_back(std::move(rd));
        auto rdb = std::make_unique<RowDiff<BRWT>>(&graph, BRWTBottomUpBuilder::build(copy_columns(cols)));
        rdb->load_anchor(anchors_file.name());
        matrices.push_back(std::move(rdb));
    }
    for (const auto &rd : matrices) {
        const auto &m = dynamic_cast<const BinaryMatrix&>(*rd);
        ASSERT_EQ(SetBits({ 7 }), m.get_rows({ requested })[0]);
        DecodeBudget budget;
        std::vector<SetBits> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, rd->decode_rows({ requested }, budget, &out, &costs, &held));
        EXPECT_EQ(SetBits({ 7 }), out[0]);
        EXPECT_EQ(path.size() - 1, costs[0].dependency_rows);
        EXPECT_EQ(num_columns, costs[0].dependency_entries);   // the anchor's; the diff is the row's own
        // a budget that holds the one-column row (inline in its SetBitPositions) is far
        // below what decoding it holds: the anchor's 300 columns and the diff's 299
        const uint64_t result = IRowDiff::output_bytes<SetBits>(1);
        EXPECT_GT(budget.peak(), 20 * result);
        EXPECT_GT(costs[0].demand, budget.peak() - 1);
        for (uint64_t max : { result, (result + budget.peak()) / 2, budget.peak() - 1 }) {
            DecodeBudget small(max);
            std::vector<SetBits> o;
            std::vector<RowCost> c;
            std::vector<uint64_t> h;
            EXPECT_EQ(DecodeStatus::REFUSED, rd->decode_rows({ requested }, small, &o, &c, &h)) << max;
            EXPECT_TRUE(o.empty());
            EXPECT_EQ(0u, small.held());
        }
    }
    // the coordinate variant: 20 coordinates per column on the anchor, the diff cancelling
    // all of them but column 7's first
    {
        std::vector<bit_vector_smart> delimiters;
        std::vector<sdsl::int_vector<>> values;
        for (size_t j = 0; j < num_columns; ++j) {
            // rows (front < back or not): the tuples are stored in row order of the column
            std::vector<std::pair<Row, std::vector<uint64_t>>> entries;
            std::vector<uint64_t> anchor_coords, diff_coords;
            for (uint64_t c = 0; c < 20; ++c) {
                // the anchor's coordinates, shifted by one at each of the two steps, meet
                // the requested row's diff after the first
                anchor_coords.push_back(1000 + 10 * c + 2);
                diff_coords.push_back(1000 + 10 * c + 1);
            }
            if (j == 7)
                diff_coords.erase(diff_coords.begin());
            entries.emplace_back(path.back(), anchor_coords);
            entries.emplace_back(path.front(), diff_coords);
            std::sort(entries.begin(), entries.end());
            std::vector<bool> delims { 1 };
            std::vector<uint64_t> vals;
            for (const auto &[row, coords] : entries) {
                for (uint64_t c : coords) {
                    vals.push_back(c);
                    delims.push_back(0);
                }
                delims.push_back(1);
            }
            sdsl::bit_vector d(delims.size());
            for (size_t i = 0; i < delims.size(); ++i) d[i] = delims[i];
            delimiters.emplace_back(d);
            sdsl::int_vector<> v(vals.size(), 0, 64);
            for (size_t i = 0; i < vals.size(); ++i) v[i] = vals[i];
            values.push_back(std::move(v));
        }
        // here every column is on both rows: the diff of column 7 keeps one coordinate
        std::vector<std::unique_ptr<bit_vector>> coord_cols;
        for (size_t j = 0; j < num_columns; ++j) {
            sdsl::bit_vector bv(num_rows, 0);
            bv[path.back()] = 1;
            bv[path.front()] = 1;
            coord_cols.push_back(std::make_unique<bit_vector_sd>(bv));
        }
        auto rd = std::make_unique<TupleRowDiff<TupleCSCMatrix<ColumnMajor>>>(
                &graph, TupleCSCMatrix<ColumnMajor>(ColumnMajor(std::move(coord_cols)),
                                                    std::move(delimiters), std::move(values)));
        rd->load_anchor(anchors_file.name());
        const auto expected = rd->get_row_tuples({ requested });
        ASSERT_EQ(1u, expected[0].size());
        EXPECT_EQ(7u, expected[0][0].first);
        EXPECT_EQ(1u, expected[0][0].second.size());
        DecodeBudget budget;
        std::vector<RowTuples> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, rd->decode_row_tuples({ requested }, budget, &out, &costs, &held));
        EXPECT_EQ(expected, out);
        EXPECT_EQ(num_columns * 21, costs[0].dependency_entries);
        const uint64_t result = copy_bytes_of(out[0]) + IRowDiff::output_bytes<RowTuples>(1);
        EXPECT_GT(budget.peak(), 20 * result);
        DecodeBudget small((result + budget.peak()) / 2);
        std::vector<RowTuples> o;
        std::vector<RowCost> c;
        std::vector<uint64_t> h;
        EXPECT_EQ(DecodeStatus::REFUSED, rd->decode_row_tuples({ requested }, small, &o, &c, &h));
        EXPECT_EQ(0u, small.held());
    }
}

// L6: the byte model bounds the heap: jemalloc's peak of a decode call never exceeds the
// peak the budget charged (with jemalloc as the process allocator)
TEST(RowDiffBudgetedDecode, ModelBoundsTheAllocator) {
#if USE_JEMALLOC
    size_t measured = 0;
    double worst = 0;
    for (const Annotation &a : annotations()) {
        for (const auto &rows : batches(a, 13)) {
            for (bool tuples : { false, true }) {
                if (tuples && !a.tuples)
                    continue;
                std::vector<SetBits> bits;
                std::vector<RowTuples> tups;
                std::vector<RowCost> costs;
                std::vector<uint64_t> held;
                DecodeBudget budget;
                uint64_t allocated = 0, deallocated = 0;
                size_t sz = sizeof(uint64_t);
                if (mallctl("thread.peak.reset", nullptr, nullptr, nullptr, 0))
                    GTEST_SKIP() << "jemalloc without thread.peak";
                mallctl("thread.allocated", &allocated, &sz, nullptr, 0);
                mallctl("thread.deallocated", &deallocated, &sz, nullptr, 0);
                const DecodeStatus status = tuples
                    ? a.rd->decode_row_tuples(rows, budget, &tups, &costs, &held)
                    : a.rd->decode_rows(rows, budget, &bits, &costs, &held);
                uint64_t peak = 0;
                mallctl("thread.peak.read", &peak, &sz, nullptr, 0);
                ASSERT_EQ(DecodeStatus::OK, status);
                if (budget.peak() > 4096 && !peak)
                    GTEST_SKIP() << "jemalloc is not the process allocator";
                EXPECT_LE(peak, budget.peak()) << a.name << " batch of " << rows.size();
                if (budget.peak())
                    worst = std::max(worst, double(peak) / budget.peak());
                measured++;
            }
        }
    }
    EXPECT_GT(measured, 50u);
    std::cerr << "jemalloc peak / model peak, worst: " << worst << std::endl;
#else
    GTEST_SKIP() << "needs jemalloc";
#endif
}

// ---- The row-diff path cache (row_diff_cache.hpp; the efficiency pass)

template <class RowT>
std::vector<RowT> cached_rows(const Annotation &a, const std::vector<Row> &rows,
                              RowDiffPathCache &cache) {
    if constexpr(std::is_same_v<RowT, SetBits>) {
        return a.rd->get_rows_cached(rows, cache);
    } else {
        return a.rd->get_row_tuples_cached(rows, cache);
    }
}

// The batches of a walk: consecutive rows along the graph in both directions, then the
// batches() mix, so that later paths meet the rows earlier calls reconstructed
std::vector<std::vector<Row>> walk_batches(const Annotation &a, uint32_t seed) {
    std::vector<std::vector<Row>> out;
    const std::vector<Row> valid = rows_of(a);
    for (size_t from = 0; from + 4 < valid.size(); from += 37) {
        for (size_t i = from; i < from + 4; ++i) {
            out.push_back({ valid[i] });
            out.push_back({ valid[i + 1], valid[i] });
        }
    }
    for (auto &b : batches(a, seed)) {
        out.push_back(std::move(b));
    }
    return out;
}

// P1: the default decode with the cache returns the rows of the default decode, call after
// call, under bounds that keep everything, evict often, and keep nothing; the cache stays
// within its bound (and a shared room) and paths stop at cached rows
template <class RowT>
void check_cached_rows(const Annotation &a) {
    for (uint64_t bound : { uint64_t(64) << 20, uint64_t(24) << 10, uint64_t(1) }) {
        RowDiffPathCache cache;
        uint64_t room = std::numeric_limits<uint64_t>::max();
        cache.set_bound(bound, [&]() { return room; });
        size_t calls = 0;
        for (const auto &rows : walk_batches(a, 5)) {
            // a shared bound that shrinks and grows (LabelOracle::make_room)
            room = calls % 7 == 3 ? bound / 3 : std::numeric_limits<uint64_t>::max();
            EXPECT_EQ(reference<RowT>(a, rows), cached_rows<RowT>(a, rows, cache))
                << a.name << " bound " << bound << " call " << calls;
            EXPECT_LE(cache.bytes(), std::min(bound, room)) << a.name << " call " << calls;
            calls++;
        }
        if (bound > (uint64_t(1) << 20)) {
            EXPECT_GT(cache.hits(), calls / 2) << a.name;
        } else if (bound == 1) {
            EXPECT_EQ(0u, cache.hits()) << a.name;
            EXPECT_EQ(0u, cache.bytes());
        }
    }
}

TEST(RowDiffPathCache, CachedDecodeEqualsTheDefaultDecode) {
    for (const Annotation &a : annotations()) {
        ASSERT_TRUE(a.rd->supports_path_cache()) << a.name;
        check_cached_rows<SetBits>(a);
        if (a.tuples)
            check_cached_rows<RowTuples>(a);
    }
}

// P2: the budget-aware decode with the cache returns the same rows with the same costs and
// held bytes as without it, call after call (the costs of a path cut at a cached row are its
// whole path's: what is admitted and charged does not depend on the cache) and at most the
// held bytes; a row decoded alone stays within its demand and its peak without the cache; rows the default decode cached (no path aggregates) are not
// used by it until it reconstructed them itself; and it refuses as the decode without the
// cache does when not even its own (smaller) decode fits
template <class RowT>
void check_budgeted_with_cache(const Annotation &a) {
    for (bool warm_by_default : { false, true }) {
        RowDiffPathCache cache;
        cache.set_bound(uint64_t(64) << 20);
        size_t calls = 0;
        for (const auto &rows : walk_batches(a, 9)) {
            if (warm_by_default && calls % 3 == 0)
                cached_rows<RowT>(a, rows, cache);     // entries without aggregates
            DecodeBudget plain_budget, cached_budget;
            std::vector<RowT> plain_out, cached_out;
            std::vector<RowCost> plain_costs, cached_costs;
            std::vector<uint64_t> plain_held, cached_held;
            ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, plain_budget, &plain_out,
                                                     &plain_costs, &plain_held));
            DecodeStatus status;
            if constexpr(std::is_same_v<RowT, SetBits>) {
                status = a.rd->decode_rows(rows, cached_budget, &cached_out, &cached_costs,
                                           &cached_held, &cache);
            } else {
                status = a.rd->decode_row_tuples(rows, cached_budget, &cached_out,
                                                 &cached_costs, &cached_held, &cache);
            }
            ASSERT_EQ(DecodeStatus::OK, status) << a.name;
            EXPECT_EQ(plain_out, cached_out) << a.name << " call " << calls;
            ASSERT_EQ(plain_costs.size(), cached_costs.size());
            for (size_t i = 0; i < rows.size(); ++i) {
                EXPECT_TRUE(same_costs(plain_costs[i], cached_costs[i]))
                    << a.name << " call " << calls << " row " << rows[i] << ": demand "
                    << cached_costs[i].demand << " without the cache " << plain_costs[i].demand;
            }
            // a returned row continued from a cached row's exact copy holds at most what it
            // holds reconstructed along its whole path (what is admitted is its demand)
            ASSERT_EQ(plain_held.size(), cached_held.size());
            for (size_t i = 0; i < rows.size(); ++i) {
                EXPECT_LE(cached_held[i], plain_held[i]) << a.name << " row " << rows[i];
            }
            EXPECT_LE(cached_budget.held(), plain_budget.held()) << a.name;
            // A row read alone holds at most what it holds decoded to its anchor (several rows
            // can hold more at once: the copies of the cached rows their paths stop at, where
            // the decode to the anchors held their small diffs; charged all the same, and a
            // run that does not fit is retried in halves, down to single rows)
            if (rows.size() == 1) {
                EXPECT_LE(cached_budget.peak(), plain_budget.peak()) << a.name << " call " << calls;
                EXPECT_LE(cached_budget.peak(), cached_costs[0].demand) << a.name;
            }
            calls++;
        }
        EXPECT_GT(cache.hits(), calls / 2) << a.name;
    }
    // a budget the cached decode does not fit refuses it whole: nothing returned, the budget
    // as on entry (both decodes on copies of one cache state: the first caches what it reads)
    RowDiffPathCache cache;
    cache.set_bound(uint64_t(64) << 20);
    const std::vector<Row> valid = rows_of(a);
    size_t refused = 0;
    for (size_t vi = 0; vi + 1 < valid.size(); vi += 11) {
        std::vector<RowT> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        DecodeBudget probe;
        RowDiffPathCache probe_cache = cache, less_cache = cache;
        DecodeStatus status;
        if constexpr(std::is_same_v<RowT, SetBits>) {
            status = a.rd->decode_rows({ valid[vi] }, probe, &out, &costs, &held, &probe_cache);
        } else {
            status = a.rd->decode_row_tuples({ valid[vi] }, probe, &out, &costs, &held,
                                             &probe_cache);
        }
        ASSERT_EQ(DecodeStatus::OK, status);
        if (probe.peak()) {
            DecodeBudget less(probe.peak() - 1);
            std::vector<RowT> sentinel(1);
            std::vector<RowCost> c;
            std::vector<uint64_t> h;
            if constexpr(std::is_same_v<RowT, SetBits>) {
                status = a.rd->decode_rows({ valid[vi] }, less, &sentinel, &c, &h, &less_cache);
            } else {
                status = a.rd->decode_row_tuples({ valid[vi] }, less, &sentinel, &c, &h,
                                                 &less_cache);
            }
            EXPECT_EQ(DecodeStatus::REFUSED, status) << a.name << " row " << valid[vi];
            EXPECT_EQ(1u, sentinel.size());
            EXPECT_EQ(0u, less.held());
            refused++;
        }
        // the cache warms along the rows (the next rows' paths meet these)
        cache = probe_cache;
    }
    EXPECT_GT(refused, 5u) << a.name;
}

TEST(RowDiffPathCache, BudgetedCostsDoNotDependOnTheCache) {
    for (const Annotation &a : annotations()) {
        check_budgeted_with_cache<SetBits>(a);
        if (a.tuples)
            check_budgeted_with_cache<RowTuples>(a);
    }
}

// P3: the cache's own bookkeeping: two generations within the bound, a hit in the older one
// moves the row to the current one, a row larger than half the bound is not kept, the
// aggregates are added to a row cached without them, trim() and clear()
TEST(RowDiffPathCache, Generations) {
    using Cache = RowDiffCache<SetBits>;
    Cache cache;
    EXPECT_FALSE(cache.enabled());
    cache.insert(1, SetBits{ 1, 2, 3 });
    EXPECT_EQ(0u, cache.size());     // off: nothing kept
    const uint64_t entry = Cache::kEntryBytes + row_copy_bytes(SetBits{ 1, 2, 3 });
    cache.set_bound(8 * entry);
    for (Row r = 0; r < 4; ++r) {
        cache.insert(r, SetBits{ 1, 2, 3 });
    }
    EXPECT_EQ(4u, cache.size());
    EXPECT_EQ(4 * entry, cache.bytes());
    cache.insert(4, SetBits{ 1, 2, 3 });     // the current generation is full: rotated
    EXPECT_EQ(5u, cache.size());
    // the counts move with their tables (review of the fixes: a rotation that dropped the
    // current generation's count with the older table counted 1 entry of the 5 held)
    EXPECT_EQ(5 * entry, cache.bytes());
    ASSERT_TRUE(cache.find(0, false));       // moved to the current generation
    EXPECT_EQ(5 * entry, cache.bytes());
    for (Row r = 5; r < 8; ++r) {
        cache.insert(r, SetBits{ 1, 2, 3 });
    }
    EXPECT_LE(cache.bytes(), 8 * entry);
    EXPECT_EQ(cache.size() * entry, cache.bytes());
    EXPECT_TRUE(cache.find(0, false));       // kept by the move
    EXPECT_FALSE(cache.find(1, false));      // dropped with the older generation
    EXPECT_FALSE(cache.find(0, true));       // no aggregates
    PathAggregates agg;
    agg.length = 3;
    cache.insert(0, SetBits{ 1, 2, 3 }, 0, &agg);
    ASSERT_TRUE(cache.find(0, true));
    EXPECT_EQ(3u, cache.find(0, true)->path.length);
    SetBits wide(4 * entry / sizeof(Column));
    std::iota(wide.begin(), wide.end(), 0);
    cache.insert(100, wide);
    EXPECT_FALSE(cache.find(100, false));    // more than half the bound
    cache.trim(entry);
    EXPECT_LE(cache.bytes(), entry);
    EXPECT_EQ(cache.size() * entry, cache.bytes());
    cache.clear();
    EXPECT_EQ(0u, cache.size());
    EXPECT_EQ(0u, cache.bytes());

    // A shared bound that shrinks (LabelOracle::make_room trims to it) keeps what both
    // generations hold within it: with the older generation's count lost at its rotation, the
    // trim saw 3 of the 7 entries held and kept all 7 beside a label cache that had grown
    uint64_t room = 8 * entry;
    cache.set_bound(8 * entry, [&room]() { return room; });
    for (Row r = 0; r < 7; ++r) {
        cache.insert(r, SetBits{ 1, 2, 3 });
    }
    EXPECT_EQ(7u, cache.size());
    EXPECT_EQ(7 * entry, cache.bytes());
    room = 4 * entry;
    cache.trim(cache.limit());
    EXPECT_LE(cache.size() * entry, 4 * entry);
    EXPECT_EQ(cache.size() * entry, cache.bytes());
    cache.insert(7, SetBits{ 1, 2, 3 });
    EXPECT_LE(cache.size() * entry, 4 * entry);
    EXPECT_EQ(cache.size() * entry, cache.bytes());
}

// Review of the efficiency pass, finding 4: a generation dropped releases its table. A cleared
// hopscotch map keeps its bucket array, sized for the most entries it held, so after many
// narrow rows a few wide ones held the narrow phase's arrays beside their bound (74.7 MiB of
// heap for 63.9 MiB accounted at a 64 MiB bound); the tables now hold no more than the share
// of the entries the bound counts (kEntryBytes), and trim() and clear() free them
TEST(RowDiffPathCache, DroppedGenerationsReleaseTheirTables) {
    using Cache = RowDiffCache<SetBits>;
    // a bucket: its key and entry, and the neighbourhood bitmap
    const uint64_t bucket = sizeof(std::pair<Row, Cache::Entry>) + sizeof(uint64_t);
    Cache cache;
    cache.set_bound(uint64_t(4) << 20);
    // narrow rows: thousands of entries a generation, rotated many times
    for (Row r = 0; r < 100000; ++r) {
        cache.insert(r, SetBits{ 1 });
    }
    const size_t narrow = cache.bucket_count();
    EXPECT_LE(narrow * bucket, cache.size() * Cache::kEntryBytes);
    // wide rows: a few hundred a generation
    SetBits wide(1000);
    std::iota(wide.begin(), wide.end(), 0);
    for (Row r = 200000; r < 202000; ++r) {
        cache.insert(r, wide);
    }
    EXPECT_GT(cache.size(), 100u);
    EXPECT_LE(cache.bytes(), uint64_t(4) << 20);
    EXPECT_LE(cache.bucket_count() * bucket, cache.size() * Cache::kEntryBytes)
        << cache.bucket_count() << " buckets for " << cache.size() << " entries";
    EXPECT_LT(cache.bucket_count() * 4, narrow);
    cache.trim(0);
    EXPECT_EQ(0u, cache.bucket_count());
    cache.insert(1, SetBits{ 1 });
    EXPECT_GT(cache.bucket_count(), 0u);
    cache.clear();
    EXPECT_EQ(0u, cache.bucket_count());
}

// Review of the efficiency pass, finding 5: with the cache a budget-aware read holds less (its
// paths are cut at cached rows), so the lookahead under a memory budget, which reads ahead
// until a run does not fit, read more with the cache than without it (UHGG 16S, annotate mode,
// 8 MiB: 75k rows warmed against 52-55k, 13-27% more instructions). Admitted as without the
// cache (admit_as_uncached), a read also charges, until its stored rows are read, the rows
// beyond the cached rows of its paths — per anchor its longest cached path's — and their
// stored rows: a row read alone holds then what decoding it to its anchor held at that point
// (that decode's peak, less at most its scratch and its reconstruction's peak beyond the
// stored rows). Its rows, costs and held bytes are those of the read without the mode
template <class RowT>
void check_admitted_as_uncached(const Annotation &a) {
    RowDiffPathCache cache;
    cache.set_bound(uint64_t(64) << 20);
    // every row cached with its aggregates
    for (const auto &rows : walk_batches(a, 3)) {
        DecodeBudget budget;
        std::vector<RowT> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        if constexpr(std::is_same_v<RowT, SetBits>) {
            ASSERT_EQ(DecodeStatus::OK, a.rd->decode_rows(rows, budget, &out, &costs, &held, &cache));
        } else {
            ASSERT_EQ(DecodeStatus::OK, a.rd->decode_row_tuples(rows, budget, &out, &costs, &held,
                                                                &cache));
        }
    }
    auto &typed = [&]() -> RowDiffCache<RowT>& {
        if constexpr(std::is_same_v<RowT, SetBits>) {
            return cache.rows;
        } else {
            return cache.tuples;
        }
    }();
    size_t spared = 0;
    for (Row row : rows_of(a)) {
        const auto *entry = typed.find(row, true);
        if (!entry || entry->path.length < 2)
            continue;
        const PathAggregates path = entry->path;
        DecodeBudget plain_budget, cut_budget, admitted_budget;
        std::vector<RowT> plain_out, cut_out, admitted_out;
        std::vector<RowCost> plain_costs, cut_costs, admitted_costs;
        std::vector<uint64_t> plain_held, cut_held, admitted_held;
        ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, { row }, plain_budget, &plain_out,
                                                 &plain_costs, &plain_held));
        for (bool admit : { false, true }) {
            cache.set_admit_as_uncached(admit);
            DecodeBudget &budget = admit ? admitted_budget : cut_budget;
            auto *out = admit ? &admitted_out : &cut_out;
            auto *costs = admit ? &admitted_costs : &cut_costs;
            auto *held = admit ? &admitted_held : &cut_held;
            if constexpr(std::is_same_v<RowT, SetBits>) {
                ASSERT_EQ(DecodeStatus::OK, a.rd->decode_rows({ row }, budget, out, costs, held,
                                                              &cache));
            } else {
                ASSERT_EQ(DecodeStatus::OK, a.rd->decode_row_tuples({ row }, budget, out, costs,
                                                                    held, &cache));
            }
        }
        cache.set_admit_as_uncached(false);
        EXPECT_EQ(plain_out, admitted_out) << a.name;
        EXPECT_TRUE(same_costs(plain_costs[0], admitted_costs[0])) << a.name << " row " << row;
        EXPECT_EQ(cut_held, admitted_held) << a.name;
        // the spared rows charged: at least their stored rows beside the cached row's copy, and
        // within what the decode to the anchor held at most, its demand
        EXPECT_GE(admitted_budget.peak(), cut_budget.peak() + path.path_stored) << a.name;
        EXPECT_LE(admitted_budget.peak(), plain_costs[0].demand + plain_costs[0].demand / 2)
            << a.name << " row " << row;
        EXPECT_GE(admitted_budget.peak() + path.recon_peak + path.full
                      + 64 * (path.max_scratch + 1) + 4096,
                  plain_budget.peak()) << a.name << " row " << row;
        spared++;
    }
    EXPECT_GT(spared, 10u) << a.name;
}

TEST(RowDiffPathCache, LookaheadAdmittedAsWithoutTheCache) {
    for (const Annotation &a : annotations()) {
        // a coordinate annotation caches its tuple rows, from which its rows are made too
        if (a.tuples) {
            check_admitted_as_uncached<RowTuples>(a);
        } else {
            check_admitted_as_uncached<SetBits>(a);
        }
    }
}

// R10 (review of feature level 4): the retention rule (RowDiffCache::keeps). Whatever it
// keeps — every row (checkpoint 1), the requested rows and checkpoints, almost nothing — the
// default and the budget-aware decodes return the rows of the decode without the cache, and
// the budget-aware one their costs; a call copies at most n x (successors + 3) + s /
// checkpoint rows (n rows asked for, s stored rows read) into the cache when no row counts as
// narrow, where keeping every row copied all s; and the narrow rows are all kept
template <class RowT>
void check_retention(const Annotation &a) {
    struct Rule { uint32_t checkpoint, successors; uint64_t narrow; };
    uint64_t total_all = 0, total_rule = 0;
    for (const Rule rule : { Rule{ 1, 0, 0 }, Rule{ 4, 1, 0 }, Rule{ 16, 8, 0 },
                             Rule{ uint32_t(1) << 30, 0, 0 },
                             Rule{ uint32_t(1) << 30, 0, uint64_t(1) << 40 } }) {
        for (bool budgeted : { false, true }) {
            RowDiffPathCache cache;
            cache.set_bound(uint64_t(64) << 20);
            cache.set_retention(rule.checkpoint, rule.successors, rule.narrow);
            auto &typed = [&]() -> RowDiffCache<RowT>& {
                if constexpr(std::is_same_v<RowT, SetBits>) {
                    return cache.rows;
                } else {
                    return cache.tuples;
                }
            }();
            size_t calls = 0;
            for (const auto &rows : walk_batches(a, 13)) {
                const uint64_t kept_before = typed.rows_inserted;
                const uint64_t stored_before = typed.stored_rows_read;
                const uint64_t bytes_before = typed.bytes_inserted;
                if (!budgeted) {
                    EXPECT_EQ(reference<RowT>(a, rows), cached_rows<RowT>(a, rows, cache))
                        << a.name << " rule " << rule.checkpoint << " call " << calls;
                } else {
                    DecodeBudget plain_budget, cached_budget;
                    std::vector<RowT> plain_out, cached_out;
                    std::vector<RowCost> plain_costs, cached_costs;
                    std::vector<uint64_t> plain_held, cached_held;
                    ASSERT_EQ(DecodeStatus::OK, decode<RowT>(a, rows, plain_budget, &plain_out,
                                                             &plain_costs, &plain_held));
                    DecodeStatus status;
                    if constexpr(std::is_same_v<RowT, SetBits>) {
                        status = a.rd->decode_rows(rows, cached_budget, &cached_out,
                                                   &cached_costs, &cached_held, &cache);
                    } else {
                        status = a.rd->decode_row_tuples(rows, cached_budget, &cached_out,
                                                         &cached_costs, &cached_held, &cache);
                    }
                    ASSERT_EQ(DecodeStatus::OK, status);
                    EXPECT_EQ(plain_out, cached_out) << a.name << " call " << calls;
                    ASSERT_EQ(plain_costs.size(), cached_costs.size());
                    for (size_t i = 0; i < rows.size(); ++i) {
                        EXPECT_TRUE(same_costs(plain_costs[i], cached_costs[i]))
                            << a.name << " rule " << rule.checkpoint << " call " << calls;
                    }
                }
                const uint64_t kept = typed.rows_inserted - kept_before;
                const uint64_t stored = typed.stored_rows_read - stored_before;
                if (rule.checkpoint == 1 || rule.narrow) {
                    // every row the call reconstructed (each stored row read is one of them)
                    EXPECT_EQ(stored, kept) << a.name << " call " << calls;
                } else {
                    EXPECT_LE(kept, rows.size() * (rule.successors + 3) + stored / rule.checkpoint)
                        << a.name << " rule " << rule.checkpoint << " call " << calls;
                }
                // (a kept row's copy can hold no heap: up to two columns are inline)
                if (!kept) {
                    EXPECT_EQ(bytes_before, typed.bytes_inserted) << a.name;
                }
                calls++;
            }
            EXPECT_LE(typed.peak_bytes, uint64_t(64) << 20);
            if (!budgeted && rule.checkpoint == 1)
                total_all = typed.rows_inserted;
            if (!budgeted && rule.checkpoint == 16)
                total_rule = typed.rows_inserted;
        }
    }
    // on these walks (paths of up to 100 rows) the rule copies fewer rows than keeping all
    EXPECT_LT(total_rule, total_all) << a.name;
}

TEST(RowDiffPathCache, RetentionKeepsRowsAndCostsAndBoundsTheCopies) {
    for (const Annotation &a : annotations()) {
        if (a.tuples) {
            check_retention<RowTuples>(a);
        } else {
            check_retention<SetBits>(a);
        }
    }
}

// The flat storage of a cached tuple row (StoredRow<RowTuples>) gives back the row exactly,
// holds less than a RowTuples copy of a row whose columns carry more than two coordinates,
// and states what that copy holds (Entry::copy_bytes, what a hit is charged)
TEST(RowDiffPathCache, FlatTupleRowsRoundTrip) {
    RowTuples row;
    for (uint32_t c = 0; c < 300; ++c) {
        RowTuples::value_type::second_type coords;
        for (uint64_t i = 0; i < 1 + c % 5; ++i) {
            coords.push_back(1000 * c + 7 * i);
        }
        row.emplace_back(3 * c, coords);
    }
    StoredRow<RowTuples> stored;
    stored.store(row);
    RowTuples back;
    stored.load(&back);
    EXPECT_EQ(row, back);
    EXPECT_LT(stored.bytes(), row_copy_bytes(row));
    RowDiffCache<RowTuples> cache;
    cache.set_bound(uint64_t(1) << 20);
    cache.insert(5, row, 3);
    const auto *entry = cache.find(5, false);
    ASSERT_TRUE(entry);
    EXPECT_EQ(row_copy_bytes(row), entry->copy_bytes);
    EXPECT_EQ(3u, entry->depth);
    EXPECT_EQ(stored.bytes() + RowDiffCache<RowTuples>::kEntryBytes, cache.bytes());
    EXPECT_EQ(1u, cache.rows_inserted);
    EXPECT_EQ(stored.bytes(), cache.bytes_inserted);
    RowTuples empty;
    stored.store(empty);
    stored.load(&back);
    EXPECT_TRUE(back.empty());
}

} // namespace
