// The budget-aware row-diff decode (DESIGN-traverse-graphlet.md §14, §14.1):
// IRowDiff::decode_rows / decode_row_tuples against the default decode, the costs it
// states per row, its refusals, and its byte model against the allocator.
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <type_traits>

#include <gtest/gtest.h>

#include "../test_annotated_dbg_helpers.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_disk/row_disk.hpp"
#include "annotation/binary_matrix/row_flat/flat_matrix.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/int_matrix/row_diff/tuple_row_diff.hpp"
#include "annotation/int_matrix/row_disk/int_row_disk.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/binary_matrix/base/hopscotch_layout.hpp"
#include "common/hashers/seeded_hash.hpp"

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

// The single-row descents return each row's set bits (and column ranks) as the batched reads
// do, for random matrices of several arities, including empty and full rows
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
    std::shared_ptr<utils::TempFile> file;          // the file a RowDisk reads its rows from
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

// With an unlimited budget the decode returns what the default path returns, and holds
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

// A row's costs are the same in every batch, decoding it alone stays within the demand it
// states, and a budget of exactly that peak admits it while one byte less refuses it
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
    // after p the demand must take p's aggregate, not the copy's exact size (otherwise a key
    // is admitted on a warm cache and refused on a cold one).
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

// A refusal at any charge returns nothing, restores the budget and never lets the held bytes
// pass the maximum; the same budget then succeeds
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

// ---- The deployed forms beyond BRWT and ColumnMajor (E14): the row-major index matrices
// (row_diff_flat, row_diff_sparse, row_diff_disk) and the counts (row_diff_int_brwt)

// every row of |m| in row order, sorted, as the row-major constructors take them
std::function<void(const BinaryMatrix::RowCallback&)> rows_in_order(const BinaryMatrix &m) {
    return [&m](const BinaryMatrix::RowCallback &callback) {
        std::vector<Row> all(m.num_rows());
        std::iota(all.begin(), all.end(), 0);
        for (const auto &row : m.get_rows(all)) {
            callback(sorted(row));
        }
    };
}

// |m| written as a RowDisk into |file| and loaded from it into |disk| (in place: a RowDisk
// reads its rows from the file, which must outlive it, and maps its boundary)
void load_row_disk(const BinaryMatrix &m, const utils::TempFile &file, RowDisk *disk) {
    RowDisk::serialize(file.name(), rows_in_order(m), m.num_columns(), m.num_rows(),
                       m.num_relations());
    utils::named_ifstream in(file.name());
    ASSERT_TRUE(disk->load(in));
}

// The single-row reads of the row-major forms return each row's set bits, ascending, as the
// batched reads do, for random matrices including empty and full rows; the scratch is all
// they keep charged (RowDisk's read alone opens a reader, charged around the read)
TEST(RowDiffBudgetedDecode, RowMajorReadsEqualTheBatchedRead) {
    std::mt19937 rng(43);
    for (double density : { 0.0, 0.05, 0.5, 1.0 }) {
        const uint64_t num_rows = 300;
        std::vector<std::unique_ptr<bit_vector>> cols;
        for (size_t j = 0; j < 37; ++j) {
            cols.push_back(random_column(num_rows, density, rng));
        }
        ColumnMajor column_major(std::move(cols));
        const uint64_t set_bits = column_major.num_relations();
        RowFlat<> flat(rows_in_order(column_major), 37, num_rows, set_bits);
        RowSparse sparse(rows_in_order(column_major), 37, num_rows, set_bits);
        utils::TempFile file;
        RowDisk disk;
        load_row_disk(column_major, file, &disk);
        std::vector<Row> all(num_rows);
        std::iota(all.begin(), all.end(), 0);
        const auto rows = column_major.get_rows(all);
        ASSERT_EQ(rows, disk.get_rows(all));
        for (const BinaryMatrix *m : { (const BinaryMatrix*)&flat, (const BinaryMatrix*)&sparse,
                                       (const BinaryMatrix*)&disk }) {
            DecodeBudget budget;
            ChargedBuffer<Column> scratch;
            for (Row r = 0; r < num_rows; ++r) {
                const std::string what = "density " + std::to_string(density) + " row "
                        + std::to_string(r);
                scratch.data.clear();
                bool ok = false;
                if (m == &flat) {
                    ok = flat.row_columns(r, &scratch, budget);
                } else if (m == &sparse) {
                    ok = sparse.row_columns(r, &scratch, budget);
                } else {
                    ok = disk.row_columns(r, &scratch, budget);
                }
                ASSERT_TRUE(ok) << what;
                EXPECT_EQ(std::vector<Column>(rows[r].begin(), rows[r].end()),
                          std::vector<Column>(scratch.data.begin(), scratch.data.end())) << what;
            }
            // the scratch is all that was charged, and only as it grew
            EXPECT_EQ(scratch.charged(), budget.held());
            if (m == &disk && density > 0) {
                // the reader was charged while it was open
                EXPECT_GE(budget.peak(), disk.reader_bytes() + scratch.charged());
            }
        }
        // a budget below what a reader holds refuses RowDisk's read alone, keeping nothing
        DecodeBudget small(disk.reader_bytes() - 1);
        ChargedBuffer<Column> scratch;
        EXPECT_FALSE(disk.row_columns(0, &scratch, small));
        EXPECT_EQ(0u, small.held());
        EXPECT_TRUE(scratch.data.empty());
        // and a growth that does not fit refuses a row-major read
        if (density == 1.0) {
            DecodeBudget none(0);
            EXPECT_FALSE(flat.row_columns(0, &scratch, none));
            EXPECT_FALSE(sparse.row_columns(0, &scratch, none));
            EXPECT_EQ(0u, none.held());
        }
    }
}

// The count annotation over the row-diff annotation |base|'s graph, anchors and fork
// successors: column j of row r counts 1 + (r * 31 + j * 7) % 5 where |base| has the bit, so
// that counts change along the row-diff paths where the labels do not. Its stored rows as
// IntRowDiff stores them: an anchor its counts, any other row the non-zero changes from its
// row-diff successor (row(r) = row(succ) + diff(r)), each encoded (encode_diff)
Annotation count_annotation(const Annotation &base) {
    const graph::DeBruijnGraph &graph = *base.rd->graph();
    const uint64_t num_rows = base.binary->num_rows();
    const uint64_t num_columns = base.binary->num_columns();
    const auto presence = base.binary->get_rows(rows_of(base));
    std::vector<std::vector<std::pair<Column, int64_t>>> counts(num_rows);
    {
        const std::vector<Row> valid = rows_of(base);
        for (size_t i = 0; i < valid.size(); ++i) {
            for (Column j : presence[i]) {
                counts[valid[i]].emplace_back(j, 1 + (valid[i] * 31 + j * 7) % 5);
            }
        }
    }
    std::vector<std::vector<std::pair<Row, uint64_t>>> stored(num_columns);
    for (Row r : rows_of(base)) {
        std::map<Column, int64_t> diff;
        for (const auto &[j, c] : counts[r]) {
            diff[j] += c;
        }
        if (!base.rd->anchor()[r]) {
            const Row succ = graph::AnnotatedDBG::graph_to_anno_index(row_diff_successor(
                    graph, graph::AnnotatedDBG::anno_to_graph_index(r), base.rd->fork_succ()));
            for (const auto &[j, c] : counts[succ]) {
                diff[j] -= c;
            }
        }
        for (const auto &[j, d] : diff) {
            if (d)
                stored[j].emplace_back(r, encode_diff(d));
        }
    }
    std::vector<std::unique_ptr<bit_vector>> cols;
    std::vector<CountsVector> values;
    for (Column j = 0; j < num_columns; ++j) {
        std::sort(stored[j].begin(), stored[j].end());
        sdsl::bit_vector bv(num_rows, 0);
        sdsl::int_vector<> v(stored[j].size(), 0, 64);
        for (size_t i = 0; i < stored[j].size(); ++i) {
            bv[stored[j][i].first] = 1;
            v[i] = stored[j][i].second;
        }
        cols.push_back(std::make_unique<bit_vector_sd>(bv));
        values.emplace_back(std::move(v));
    }
    using Counts = CSCMatrix<BRWT, CountsVector>;
    auto rd = std::make_unique<IntRowDiff<Counts>>(
            &graph, Counts(BRWTBottomUpBuilder::build(std::move(cols),
                                                      BRWTBottomUpBuilder::get_basic_partitioner(2)),
                           std::move(values)));
    copy_support(*base.rd, rd.get());
    Annotation a;
    a.name = "row_diff_int_brwt";
    a.anno = base.anno;
    a.rd = rd.get();
    a.binary = rd.get();
    a.owned = std::move(rd);
    return a;
}

// The forms of this section, over the diffs of the binary row_diff annotation of annotations()
// (the same stored rows, so the same rows, dependencies and, for flat and sparse, the same
// demands as row_diff; row_diff_disk's demands add its reader), and the counts
std::vector<Annotation> deployed_forms() {
    std::vector<Annotation> base = annotations();
    const Annotation &rd = base[0];
    EXPECT_EQ("row_diff", rd.name);
    const auto &diffs = dynamic_cast<const RowDiff<ColumnMajor>&>(*rd.binary).diffs();
    const graph::DeBruijnGraph *graph = rd.rd->graph();
    std::vector<Annotation> out;
    {
        auto flat = std::make_unique<RowDiff<RowFlat<>>>(graph, rows_in_order(diffs),
                diffs.num_columns(), diffs.num_rows(), diffs.num_relations());
        copy_support(*rd.rd, flat.get());
        Annotation a;
        a.name = "row_diff_flat";
        a.anno = rd.anno;
        a.rd = flat.get();
        a.binary = flat.get();
        a.owned = std::move(flat);
        out.push_back(std::move(a));
    }
    {
        auto sparse = std::make_unique<RowDiff<RowSparse>>(graph, rows_in_order(diffs),
                diffs.num_columns(), diffs.num_rows(), diffs.num_relations());
        copy_support(*rd.rd, sparse.get());
        Annotation a;
        a.name = "row_diff_sparse";
        a.anno = rd.anno;
        a.rd = sparse.get();
        a.binary = sparse.get();
        a.owned = std::move(sparse);
        out.push_back(std::move(a));
    }
    {
        Annotation a;
        a.file = std::make_shared<utils::TempFile>();
        auto disk = std::make_unique<RowDiff<RowDisk>>(graph);
        load_row_disk(diffs, *a.file, &disk->diffs());
        copy_support(*rd.rd, disk.get());
        a.name = "row_diff_disk";
        a.anno = rd.anno;
        a.rd = disk.get();
        a.binary = disk.get();
        a.owned = std::move(disk);
        out.push_back(std::move(a));
    }
    out.push_back(count_annotation(rd));
    for (const Annotation &a : out) {
        EXPECT_TRUE(a.rd->supports_budgeted_decode()) << a.name;
    }
    return out;
}

TEST(RowDiffBudgetedDecode, EveryDeployedFormEqualsTheDefaultDecode) {
    for (const Annotation &a : deployed_forms()) {
        check_equivalence<SetBits>(a);
    }
}

TEST(RowDiffBudgetedDecode, EveryDeployedFormCostsAreBatchIndependent) {
    for (const Annotation &a : deployed_forms()) {
        check_batch_independence<SetBits>(a);
    }
}

TEST(RowDiffBudgetedDecode, EveryDeployedFormRefusalsHaveNoSideEffects) {
    for (const Annotation &a : deployed_forms()) {
        check_denials<SetBits>(a);
    }
}

// The row-major forms decode the stored rows row_diff (ColumnMajor) decodes: the same rows and
// costs, row_diff_disk's demand that of row_diff plus its reader (opened once per call, in
// every row's demand), which a budget below it refuses whole. Its peak is row_diff's at most
// plus the reader, which is held from the first stored row on: where row_diff's peak is the
// trace's growth (before any stored row is read: the old and new tables of a growing
// container, which a call of narrow rows may never reach again), less than that
TEST(RowDiffBudgetedDecode, RowMajorFormsChargeAsTheColumnForm) {
    std::vector<Annotation> base = annotations();
    const Annotation &column = base[0];
    for (const Annotation &a : deployed_forms()) {
        if (a.name == "row_diff_int_brwt")
            continue;
        const uint64_t reader = a.name == "row_diff_disk"
            ? dynamic_cast<const RowDiff<RowDisk>&>(*a.binary).diffs().reader_bytes() : 0;
        for (const auto &rows : batches(a, 17)) {
            DecodeBudget b1, b2;
            std::vector<SetBits> o1, o2;
            std::vector<RowCost> c1, c2;
            std::vector<uint64_t> h1, h2;
            ASSERT_EQ(DecodeStatus::OK, decode<SetBits>(column, rows, b1, &o1, &c1, &h1));
            ASSERT_EQ(DecodeStatus::OK, decode<SetBits>(a, rows, b2, &o2, &c2, &h2));
            EXPECT_EQ(o1, o2) << a.name;
            EXPECT_EQ(h1, h2) << a.name;
            ASSERT_EQ(c1.size(), c2.size());
            for (size_t i = 0; i < rows.size(); ++i) {
                EXPECT_EQ(c1[i].dependency_rows, c2[i].dependency_rows) << a.name;
                EXPECT_EQ(c1[i].dependency_entries, c2[i].dependency_entries) << a.name;
                EXPECT_EQ(c1[i].demand + reader, c2[i].demand) << a.name;
            }
            EXPECT_EQ(b1.held(), b2.held()) << a.name;
            if (!rows.empty()) {
                EXPECT_LE(b1.peak(), b2.peak()) << a.name << " batch of " << rows.size();
                EXPECT_LE(b2.peak(), b1.peak() + reader)
                        << a.name << " batch of " << rows.size();
                if (!reader) {
                    EXPECT_EQ(b1.peak(), b2.peak()) << a.name << " batch of " << rows.size();
                }
            }
            if (reader && !rows.empty()) {
                DecodeBudget small(reader - 1);
                std::vector<SetBits> o;
                std::vector<RowCost> c;
                std::vector<uint64_t> h;
                EXPECT_EQ(DecodeStatus::REFUSED, decode<SetBits>(a, rows, small, &o, &c, &h));
                EXPECT_EQ(0u, small.held());
            }
        }
    }
}

// Counts are read as presence, and presence needs the counts: a column stays in a row when a
// diff changes its count without cancelling it. Decoded as bits (a RowDiff over the same
// binary diff matrix, symmetric differences), those rows would lose the column; the
// budget-aware decode of the counts returns what get_rows() returns
TEST(RowDiffBudgetedDecode, CountsAreReadAsPresence) {
    std::vector<Annotation> base = annotations();
    const Annotation counts = count_annotation(base[0]);
    const auto &rd = dynamic_cast<const IntRowDiff<CSCMatrix<BRWT, CountsVector>>&>(*counts.binary);
    const std::vector<Row> valid = rows_of(counts);
    // the counts annotation's labels are those of the base (the fixture is right)
    EXPECT_EQ(base[0].binary->get_rows(valid), counts.binary->get_rows(valid));
    // the same stored entries decoded as bits
    std::vector<std::unique_ptr<bit_vector>> cols;
    for (Column j = 0; j < rd.num_columns(); ++j) {
        sdsl::bit_vector bv(rd.num_rows(), 0);
        for (Row r : rd.diffs().get_column(j)) {
            bv[r] = 1;
        }
        cols.push_back(std::make_unique<bit_vector_sd>(bv));
    }
    RowDiff<ColumnMajor> as_bits(rd.graph(), ColumnMajor(std::move(cols)));
    copy_support(rd, &as_bits);
    size_t differ = 0;
    const auto bits = as_bits.get_rows(valid);
    const auto presence = counts.binary->get_rows(valid);
    for (size_t i = 0; i < valid.size(); ++i) {
        differ += bits[i] != presence[i];
    }
    EXPECT_GT(differ, valid.size() / 10);
    DecodeBudget budget;
    std::vector<SetBits> out;
    std::vector<RowCost> costs;
    std::vector<uint64_t> held;
    ASSERT_EQ(DecodeStatus::OK, rd.decode_rows(valid, budget, &out, &costs, &held));
    EXPECT_EQ(presence, out);
    // decode_row_tuples is not a count annotation's
    std::vector<RowTuples> tuples;
    EXPECT_EQ(DecodeStatus::UNSUPPORTED,
              rd.decode_row_tuples(valid, budget, &tuples, &costs, &held));
}

// ---- Index matrices read without a row-diff (BudgetedRows): a plain BRWT, a plain ColumnMajor

// The rows of a plain matrix, each read alone: what get_rows() returns, no dependency rows, a
// demand that a row read alone stays within and that a budget of exactly its peak admits while
// one byte less refuses it; costs independent of the batch; a refused call keeps nothing
TEST(BudgetedRowsAlone, PlainBrwtAndColumnMajor) {
    std::mt19937 rng(44);
    for (double density : { 0.0, 0.05, 0.5, 1.0 }) {
        const uint64_t num_rows = 200;
        std::vector<std::unique_ptr<bit_vector>> cols;
        for (size_t j = 0; j < 23; ++j) {
            cols.push_back(random_column(num_rows, density, rng));
        }
        ColumnMajor column_major(copy_columns(cols));
        BRWT brwt = BRWTBottomUpBuilder::build(std::move(cols),
                                               BRWTBottomUpBuilder::get_basic_partitioner(3));
        RowFlat<> flat(rows_in_order(column_major), 23, num_rows, column_major.num_relations());
        // a row-major or row-diff matrix has none (a row-diff one has IRowDiff::decode_rows)
        EXPECT_FALSE(budgeted_rows(flat));
        for (const BinaryMatrix *m : { (const BinaryMatrix*)&column_major,
                                       (const BinaryMatrix*)&brwt }) {
            const auto reader = budgeted_rows(*m);
            ASSERT_TRUE(reader);
            std::vector<Row> all(num_rows);
            std::iota(all.begin(), all.end(), 0);
            std::vector<Row> rows(all.rbegin(), all.rend());
            rows.push_back(5);      // a duplicate
            DecodeBudget budget;
            budget.charge(100);     // something held on entry
            std::vector<SetBits> out;
            std::vector<RowCost> costs;
            std::vector<uint64_t> held;
            ASSERT_EQ(DecodeStatus::OK, reader->decode_rows(rows, budget, &out, &costs, &held));
            EXPECT_EQ(m->get_rows(rows), out);
            uint64_t kept = 100 + IRowDiff::output_bytes<SetBits>(rows.size());
            for (size_t i = 0; i < rows.size(); ++i) {
                kept += held[i];
                EXPECT_EQ(0u, costs[i].dependency_rows);
                EXPECT_EQ(0u, costs[i].dependency_entries);
                EXPECT_EQ(small_vector_bytes(out[i].size(), sizeof(Column)), held[i]);
            }
            EXPECT_EQ(kept, budget.held());
            size_t refused = 0;
            for (Row r = 0; r < num_rows; r += 7) {
                DecodeBudget alone;
                std::vector<SetBits> o;
                std::vector<RowCost> c;
                std::vector<uint64_t> h;
                ASSERT_EQ(DecodeStatus::OK, reader->decode_rows({ r }, alone, &o, &c, &h));
                EXPECT_LE(alone.peak(), c[0].demand) << r;
                EXPECT_EQ(costs[num_rows - 1 - r].demand, c[0].demand) << r;
                DecodeBudget exact(alone.peak());
                EXPECT_EQ(DecodeStatus::OK, reader->decode_rows({ r }, exact, &o, &c, &h));
                DecodeBudget less(alone.peak() - 1);
                std::vector<SetBits> sentinel(1);
                EXPECT_EQ(DecodeStatus::REFUSED,
                          reader->decode_rows({ r }, less, &sentinel, &c, &h));
                EXPECT_EQ(1u, sentinel.size());
                EXPECT_EQ(0u, less.held());
                refused++;
            }
            EXPECT_GT(refused, 20u);
            // a refusal at any charge keeps nothing, and the same budget then succeeds
            const uint64_t charges = budget.charges() - 1;
            for (uint64_t n = 0; n < charges; n += std::max<uint64_t>(1, charges / 40)) {
                DecodeBudget b(budget.peak() + 1000);
                b.charge(100);
                b.deny = [&](uint64_t ordinal) { return ordinal == n + 1; };
                std::vector<SetBits> sentinel(1);
                std::vector<RowCost> c;
                std::vector<uint64_t> h;
                EXPECT_EQ(DecodeStatus::REFUSED, reader->decode_rows(rows, b, &sentinel, &c, &h));
                EXPECT_EQ(1u, sentinel.size());
                EXPECT_EQ(100u, b.held());
                b.deny = nullptr;
                std::vector<SetBits> again;
                ASSERT_EQ(DecodeStatus::OK, reader->decode_rows(rows, b, &again, &c, &h));
                EXPECT_EQ(out, again);
            }
        }
    }
}

// A plain BRWT whose column groups are not ranges of columns (as transform_anno --greedy
// builds them): its batched read returns a row's columns in the tree's order, not sorted, and
// the budget-aware read of each row alone returns exactly that -- the same columns in the same
// order, so a caller that keeps a row's first columns (the annotate recorder's
// max_labels_per_node) or numbers them in row order keeps and numbers the same ones under a
// budget as without one. The charges are those of any row read alone: the scratch and the copy
// at the row's size.
TEST(BudgetedRowsAlone, PlainBrwtKeepsTheOrderOfItsBatchedRead) {
    std::mt19937 rng(45);
    const uint64_t num_rows = 300;
    const size_t num_columns = 12;
    // the leaves in groups {1, 4, 7, 10}, {2, 5, 8, 11}, {0, 3, 6, 9} (in this order), then
    // the three under the root
    const BRWTBottomUpBuilder::Partitioner interleaved = [](const auto &columns) {
        BRWTBottomUpBuilder::Partition groups(columns.size() <= 3 ? 1 : 3);
        for (uint64_t j = 0; j < columns.size(); ++j) {
            groups[j % groups.size()].push_back(j);
        }
        std::rotate(groups.begin(), groups.begin() + 1, groups.end());
        return groups;
    };
    for (double density : { 0.3, 0.7, 1.0 }) {
        std::vector<std::unique_ptr<bit_vector>> cols;
        for (size_t j = 0; j < num_columns; ++j) {
            cols.push_back(random_column(num_rows, density, rng));
        }
        BRWT brwt = BRWTBottomUpBuilder::build(std::move(cols), interleaved);
        const auto reader = budgeted_rows(brwt);
        ASSERT_TRUE(reader);
        std::vector<Row> rows(num_rows);
        std::iota(rows.begin(), rows.end(), 0);
        std::shuffle(rows.begin(), rows.end(), rng);
        const std::vector<SetBits> expected = brwt.get_rows(rows);
        // the fixture is what it is for: rows the tree does not return in column order
        size_t unsorted = 0;
        for (const SetBits &row : expected) {
            unsorted += !std::is_sorted(row.begin(), row.end());
        }
        ASSERT_GT(unsorted, num_rows / 4) << density;
        DecodeBudget budget;
        std::vector<SetBits> out;
        std::vector<RowCost> costs;
        std::vector<uint64_t> held;
        ASSERT_EQ(DecodeStatus::OK, reader->decode_rows(rows, budget, &out, &costs, &held));
        EXPECT_EQ(expected, out) << density;
        for (size_t i = 0; i < rows.size(); ++i) {
            EXPECT_EQ(small_vector_bytes(out[i].size(), sizeof(Column)), held[i]);
            // read alone: the demand is what the read holds at its peak, exactly
            DecodeBudget alone;
            std::vector<SetBits> o;
            std::vector<RowCost> c;
            std::vector<uint64_t> h;
            ASSERT_EQ(DecodeStatus::OK, reader->decode_rows({ rows[i] }, alone, &o, &c, &h));
            EXPECT_EQ(expected[i], o[0]);
            EXPECT_EQ(costs[i].demand, c[0].demand);
            EXPECT_LE(alone.peak(), c[0].demand);
        }
    }
}

// Formats without the budget-aware path say so, charge nothing and touch nothing (a count
// row-diff annotation on disk, row_diff_int_disk: not a deployed form)
TEST(RowDiffBudgetedDecode, UnsupportedFormatsChargeNothing) {
    IntRowDiff<IntRowDisk> rd;
    EXPECT_FALSE(rd.supports_budgeted_decode());
    DecodeBudget budget;
    std::vector<SetBits> out;
    std::vector<RowCost> costs;
    std::vector<uint64_t> held;
    EXPECT_EQ(DecodeStatus::UNSUPPORTED, rd.decode_rows({ 0, 1 }, budget, &out, &costs, &held));
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

// "A row-diff row whose dependencies are dense but whose result is tiny": an anchor of 300
// columns (or 6,000 coordinates) and a diff on the path that cancels all but one of them. The
// default decode returns the one-column row; the budget-aware decode states the dependencies,
// and a budget between the result's bytes and what decoding it holds refuses it whole.
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
        ASSERT_EQ(SetBits({ 7 }), m.get_rows(std::vector<Row>{ requested })[0]);
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
        const auto expected = rd->get_row_tuples(std::vector<Row>{ requested });
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

// The byte model bounds the heap: jemalloc's peak of a decode call never exceeds the peak the
// budget charged (with jemalloc as the process allocator)
TEST(RowDiffBudgetedDecode, ModelBoundsTheAllocator) {
#if USE_JEMALLOC
    size_t measured = 0;
    double worst = 0;
    std::vector<Annotation> all = annotations();
    for (Annotation &a : deployed_forms()) {
        all.push_back(std::move(a));
    }
    for (const Annotation &a : all) {
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

// ---- The row-diff path cache (row_diff_cache.hpp)

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

// The default decode with the cache returns the rows of the default decode, call after call,
// under bounds that keep everything, evict often, and keep nothing; the cache stays within its
// bound (and a shared room) and paths stop at cached rows
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
    // the row-major forms (the counts have no path cache)
    for (const Annotation &a : deployed_forms()) {
        if (a.name == "row_diff_int_brwt") {
            EXPECT_FALSE(a.rd->supports_path_cache());
            continue;
        }
        ASSERT_TRUE(a.rd->supports_path_cache()) << a.name;
        check_cached_rows<SetBits>(a);
    }
}

// The budget-aware decode with the cache returns the same rows with the same costs and held
// bytes as without it, call after call (the costs of a path cut at a cached row are its whole
// path's: what is admitted and charged does not depend on the cache) and at most the held
// bytes; a row decoded alone stays within its demand and its peak without the cache; rows the
// default decode cached (no path aggregates) are not used by it until it reconstructed them
// itself; and it refuses as the decode without the cache does when not even its own (smaller)
// decode fits
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
    for (const Annotation &a : deployed_forms()) {
        if (a.name != "row_diff_int_brwt")
            check_budgeted_with_cache<SetBits>(a);
    }
}

// The cache's own bookkeeping: two generations within the bound, its rows and its tables
// counted as they are (bytes() = the rows' stored bytes + table_bytes()), an insert made only
// with room for what it can add to its table (its doubling where it can double it); a hit in
// the older generation moves the row to the current one; a row larger than half the bound is
// not kept; the aggregates are added to a row cached without them; trim() and clear()
TEST(RowDiffPathCache, Generations) {
    using Cache = RowDiffCache<SetBits>;
    using Layout = HopscotchLayout<Cache::Table>;
    Cache cache;
    EXPECT_FALSE(cache.enabled());
    cache.insert(1, SetBits{ 1, 2, 3 });
    EXPECT_EQ(0u, cache.size());     // off: nothing kept
    SetBits row(1000);
    std::iota(row.begin(), row.end(), 0);
    const uint64_t r = row_copy_bytes(row);
    // what the rows and the tables hold, as the cache counts them
    auto counted = [&]() { return cache.size() * r + cache.table_bytes(); };
    // tsl doubles a table at 90% load: 2 buckets for one entry, 4 for two or three, 8 for
    // four to seven. An insert keeps room for its table's doubling only where it can double
    // it: at that threshold (the 2nd, 4th and 8th row) or from 80% load (the 7th row, into 8
    // buckets: kHopscotchGrowthLoad); else for one overflow node. Half the bound is the
    // fourth row's insert (three rows, 4 buckets and the 8 beside them, the fourth row): a
    // generation takes six rows in 8 buckets (the fifth and the sixth fit where the 4
    // buckets were), not the seventh (its doubling to 16 buckets). Room for a doubling at
    // every insert rotated it at the fifth row (R4 review, S1)
    const uint64_t half = 4 * r + Layout::array_bytes(4) + Layout::array_bytes(8)
                        + Layout::node_bytes();
    ASSERT_LE(2 * r, Layout::array_bytes(4));
    const uint64_t bound = 2 * half;
    cache.set_bound(bound);
    for (Row i = 0; i < 6; ++i) {
        cache.insert(i, row);
        EXPECT_EQ(counted(), cache.bytes());
    }
    EXPECT_EQ(6u, cache.size());
    EXPECT_EQ(Layout::array_bytes(8), cache.table_bytes());
    EXPECT_EQ(6 * r + Layout::array_bytes(8), cache.bytes());
    cache.insert(6, row);     // the current generation is full: rotated
    EXPECT_EQ(7u, cache.size());
    // the counts move with their tables (a rotation that dropped the current generation's
    // count with the older table would count 1 entry of the 7 held)
    EXPECT_EQ(7 * r + Layout::array_bytes(8) + Layout::array_bytes(2), cache.bytes());
    ASSERT_TRUE(cache.find(0, false));       // moved to the current generation
    // the row's bytes moved; the current table doubled for it, the older one kept its buckets
    EXPECT_EQ(7 * r + Layout::array_bytes(8) + Layout::array_bytes(4), cache.bytes());
    EXPECT_EQ(counted(), cache.bytes());
    for (Row i = 7; i < 12; ++i) {
        cache.insert(i, row);
        EXPECT_EQ(counted(), cache.bytes());
        EXPECT_LE(cache.bytes(), bound);
    }
    EXPECT_TRUE(cache.find(0, false));       // kept by the move
    EXPECT_FALSE(cache.find(1, false));      // dropped with the older generation
    EXPECT_FALSE(cache.find(0, true));       // no aggregates
    PathAggregates agg;
    agg.length = 3;
    cache.insert(0, row, 0, &agg);
    ASSERT_TRUE(cache.find(0, true));
    EXPECT_EQ(3u, cache.find(0, true)->path.length);
    EXPECT_LE(cache.peak_bytes, bound);
    EXPECT_GE(cache.peak_bytes, cache.bytes());
    SetBits wide(bound / 2 / sizeof(Column));
    std::iota(wide.begin(), wide.end(), 0);
    cache.insert(100, wide);
    EXPECT_FALSE(cache.find(100, false));    // more than half the bound
    EXPECT_EQ(1u, cache.rows_refused);
    cache.trim(half);
    EXPECT_LE(cache.bytes(), half);
    EXPECT_EQ(counted(), cache.bytes());
    cache.clear();
    EXPECT_EQ(0u, cache.size());
    EXPECT_EQ(0u, cache.bytes());

    // A shared bound that shrinks (LabelOracle::make_room trims to it) keeps what both
    // generations hold within it: with the older generation's count lost at its rotation, the
    // trim saw 3 of the 7 entries held and kept all 7 beside a label cache that had grown
    uint64_t room = bound;
    cache.set_bound(bound, [&room]() { return room; });
    for (Row i = 0; i < 7; ++i) {
        cache.insert(i, row);
        EXPECT_LE(cache.bytes(), bound);
    }
    EXPECT_EQ(7u, cache.size());
    EXPECT_EQ(counted(), cache.bytes());
    room = half;
    cache.trim(cache.limit());
    EXPECT_LE(cache.bytes(), half);
    EXPECT_EQ(counted(), cache.bytes());
    cache.insert(7, row);
    EXPECT_LE(cache.bytes(), half);
    EXPECT_EQ(counted(), cache.bytes());
    // a move that has no room for the current table's doubling is not made: the row is
    // found where it is
    room = cache.bytes();
    const size_t held = cache.size();
    for (Row i = 0; i < 8; ++i) {
        cache.find(i, false);
        EXPECT_LE(cache.bytes(), room);
    }
    EXPECT_EQ(held, cache.size());
}

// The bound holds whatever the rows' ids (GPT review 7, finding 8). Rows hashed by their id
// (std::hash: the identity) take a table's bucket from the id's low bits. Ids bit-reversed (the
// row i is bitreverse_20(i)) share their low bits within any run of consecutive i, so tsl finds
// a neighbourhood full and doubles a table down to 10% load: the tables held far more than the
// old per-entry charge (640 B, the share at 40% load) allowed for. Ids that share their low 20
// bits all fall on one bucket, which no doubling separates: they go to the overflow list.
// Counted as they are, with room for the doubling of each insert that can double its table
// (tsl's threshold, or from 80% load, where a full neighbourhood doubles a seeded table), the
// cache holds at most its bound after every insert, with the identity hash as with the seeded
// one. Its peak stays within the bound with the seeded hash; with the identity hash the
// bit-reversed ids, which share a neighbourhood, double a table below that load or twice in one
// insert, which passes the bound while it runs (the residual SPEC-labeled-traversal-core.md
// section 8.4 states). The seeded tables stay at the occupancy of random ids (their buckets
// within 3 per entry / 0.8, the padding aside), so the bound is spent on rows
uint64_t bitreverse_20(uint64_t i) {
    uint64_t r = 0;
    for (int b = 0; b < 20; ++b) {
        r |= ((i >> b) & 1) << (19 - b);
    }
    return r;
}

template <class Hash>
void check_bound_whatever_the_ids(const char *name) {
    using Cache = RowDiffCache<SetBits, Hash>;
    using Layout = HopscotchLayout<typename Cache::Table>;
    const bool seeded = std::is_same<Hash, common::SeededHash>::value;
    const uint64_t padding = 2 * Layout::array_bytes(1);
    struct Ids {
        const char *pattern;
        std::function<Row(Row)> id;
        Row rows;
        std::vector<uint64_t> bounds;
    };
    // (the overflow lists are searched linearly: fewer rows for them)
    const std::vector<Ids> patterns = {
        { "bit-reversed", [](Row i) { return bitreverse_20(i); }, 200000,
          { uint64_t(1) << 20, uint64_t(8) << 20 } },
        { "i << 20", [](Row i) { return i << 20; }, 30000, { uint64_t(1) << 20 } },
    };
    for (const Ids &ids : patterns) {
        const std::string pattern = ids.pattern;
        // the most the tables held per entry beyond their padding, under any of the bounds
        double most_per_entry = 0;
        for (uint64_t bound : ids.bounds) {
            Cache cache;
            cache.set_bound(bound);
            uint64_t most = 0;
            double most_here = 0;
            for (Row i = 0; i < ids.rows; ++i) {
                cache.insert(ids.id(i), SetBits{ 1 });
                most = std::max(most, cache.bytes());
                ASSERT_LE(cache.bytes(), bound) << name << ", " << pattern << ", row " << i;
                ASSERT_EQ(cache.size() * row_copy_bytes(SetBits{ 1 }) + cache.table_bytes(),
                          cache.bytes()) << name << ", " << pattern;
                if (cache.size() >= 200 && cache.table_bytes() > padding) {
                    most_here = std::max(most_here,
                            double(cache.table_bytes() - padding) / cache.size());
                }
            }
            most_per_entry = std::max(most_per_entry, most_here);
            EXPECT_GE(cache.peak_bytes, most) << name << ", " << pattern;
            EXPECT_GT(cache.size(), 0u) << name << ", " << pattern;
            if (seeded) {
                EXPECT_LE(cache.peak_bytes, bound) << name << ", " << pattern;
                EXPECT_LE(cache.table_bytes(),
                          cache.size() * 3 * Layout::kBucketBytes * 5 / 4 + padding)
                    << name << ", " << pattern << ": " << cache.bucket_count()
                    << " buckets for " << cache.size();
                // the rows the old charge kept, about: the bound spent on rows, not buckets
                EXPECT_GT(cache.size() * 4, bound / (Layout::kBucketBytes * 3 * 5 / 4 + 8))
                    << name << ", " << pattern;
                EXPECT_LT(most_here, 640) << name << ", " << pattern;
            } else if (pattern == "i << 20") {
                // one bucket's neighbourhood holds 62 of them: the tables stop at 128 buckets
                // (90% load of the 62), the other rows in the overflow lists
                EXPECT_LE(cache.bucket_count(), 2u * 128) << name << ", " << pattern;
                EXPECT_GT(cache.size(), 1000u) << name << ", " << pattern;
            }
            std::cerr << name << ", ids " << pattern << " under " << (bound >> 20) << " MiB: "
                      << cache.size() << " rows held, tables " << cache.table_bytes() << " B ("
                      << cache.bucket_count() << " buckets), at most " << int(most_here)
                      << " B of table per entry, peak " << cache.peak_bytes << " B" << std::endl;
        }
        if (!seeded && pattern == "bit-reversed") {
            // what the identity hash made of them (under the smaller bound, where a table
            // holds a short run of rows): more than the old charge
            EXPECT_GT(most_per_entry, 640) << name << ", " << pattern;
        }
    }
}

TEST(RowDiffPathCache, TheBoundHoldsWhateverTheRowIds) {
    check_bound_whatever_the_ids<std::hash<Row>>("identity hash");
    check_bound_whatever_the_ids<common::SeededHash>("seeded hash");
}

// An insert can double a table at tsl's threshold (90% load) and, when the neighbourhood of its
// key's bucket is full, before it: the seeded tables did so from 82.5% load at the lowest
// (hopscotch_layout.hpp, kHopscotchGrowthLoad 0.8). Every growth of a seeded table, over random and
// consecutive keys and many seeds, is one hopscotch_insert_can_grow foresaw, and what an insert
// held while it ran is within the room hopscotch_insert_bytes states (R4 review, S1)
TEST(RowDiffPathCache, EveryGrowthOfASeededTableIsForeseen) {
    using Table = tsl::hopscotch_map<uint64_t, uint64_t, common::SeededHash>;
    using Layout = HopscotchLayout<Table>;
    // 8 buckets: no growth before the 6th entry (80% load), tsl's threshold at the 7th
    Table eight(8);
    ASSERT_EQ(8u, eight.bucket_count());
    for (uint64_t key = 0; key < 7; ++key) {
        EXPECT_EQ(key >= 6, hopscotch_insert_can_grow(eight)) << key;
        EXPECT_EQ(key >= 6 ? hopscotch_growth_bytes(eight)
                           : Layout::node_bytes(),
                  hopscotch_insert_bytes(eight)) << key;
        eight.emplace(key, key);
    }
    EXPECT_TRUE(hopscotch_insert_can_grow(Table()));
    uint64_t growths = 0;
    for (uint64_t seed = 1; seed <= 40; ++seed) {
        Table table(0, common::SeededHash{ seed });
        std::mt19937_64 rng(seed);
        for (uint64_t i = 0; i < (1u << 14); ++i) {
            const bool can = hopscotch_insert_can_grow(table);
            const uint64_t room = hopscotch_insert_bytes(table);
            const uint64_t before = hopscotch_held_bytes(table);
            const size_t buckets = table.bucket_count();
            table.emplace(seed % 2 ? rng() : i, i);
            if (table.bucket_count() == buckets) {
                EXPECT_LE(hopscotch_held_bytes(table), before + room);
                continue;
            }
            growths++;
            ASSERT_TRUE(can) << "seed " << seed << ": " << table.size() - 1 << " entries in "
                             << buckets << " buckets";
            ASSERT_EQ(Layout::next_bucket_count(buckets), table.bucket_count());
            // the old and the new bucket array at once
            EXPECT_LE(Layout::array_bytes(buckets) + Layout::array_bytes(table.bucket_count()),
                      before + room);
        }
    }
    EXPECT_GE(growths, 40u * 13);
}

// A generation fills its table to the load at which an insert can double it (80%) before it
// rotates. Room for a doubling at every insert left a generation of narrow rows at 45% load of
// a table whose doubling did not fit half the bound: 13-40% fewer narrow tuple rows under
// 128 MiB than the old 640 B charge kept (R4 review, S1; most rows held under 128 MiB over
// 1.5 M random ids, old charge / room at every insert / now: tuple rows of 48 B 195,082 /
// 117,966 / 209,714, of 352 B 135,300 / 117,966 / 209,714)
template <class RowT>
void check_a_generation_fills_its_table(const RowT &row, const char *name) {
    using Cache = RowDiffCache<RowT>;
    using Layout = HopscotchLayout<typename Cache::Table>;
    const uint64_t bound = uint64_t(8) << 20;
    const uint64_t b = StoredRow<RowT>::bytes_of(row);
    // the largest table a generation affords: its doubling from half as many buckets, made at
    // 90% load of those at the latest, fits half the bound
    uint64_t top = 2;
    while (Layout::array_bytes(top) + Layout::array_bytes(2 * top) + Layout::node_bytes()
               + (top * 9 / 10 + 1) * b <= bound / 2) {
        top *= 2;
    }
    Cache cache;
    cache.set_bound(bound);
    std::mt19937_64 rng(7);
    size_t most = 0;
    for (uint64_t i = 0; i < 6 * top; ++i) {
        cache.insert(rng(), row);
        most = std::max(most, cache.size());
        ASSERT_LE(cache.bytes(), bound) << name;
    }
    // both generations at 80% load of the largest table
    EXPECT_GE(most, 2 * size_t(float(top) * kHopscotchGrowthLoad))
        << name << ": tables of " << top << " buckets";
    EXPECT_LE(cache.peak_bytes, bound) << name;
}

TEST(RowDiffPathCache, AGenerationFillsItsTable) {
    check_a_generation_fills_its_table(SetBits{ 1 }, "binary rows");
    check_a_generation_fills_its_table(RowTuples{ { 0, { 1 } } }, "tuple rows");
}

// A generation dropped releases its table. A cleared hopscotch map keeps its bucket array,
// sized for the most entries it held, so after many narrow rows a few wide ones would hold the
// narrow phase's arrays beside their bound (74.7 MiB of heap for 63.9 MiB accounted at a
// 64 MiB bound); the tables are counted as they are (table_bytes(), within the bound), and
// trim() and clear() free them
TEST(RowDiffPathCache, DroppedGenerationsReleaseTheirTables) {
    using Cache = RowDiffCache<SetBits>;
    using Layout = HopscotchLayout<Cache::Table>;
    Cache cache;
    const uint64_t bound = uint64_t(4) << 20;
    cache.set_bound(bound);
    // narrow rows: thousands of entries a generation, rotated many times
    for (Row r = 0; r < 100000; ++r) {
        cache.insert(r, SetBits{ 1 });
    }
    const size_t narrow = cache.bucket_count();
    EXPECT_EQ(Layout::array_bytes(cache.bucket_count()) > 0, cache.table_bytes() > 0);
    EXPECT_LE(cache.bytes(), bound);
    // the seeded tables hold about what random ids take: 3 buckets / 0.8 per entry at most
    EXPECT_LE(cache.table_bytes(), cache.size() * 3 * Layout::kBucketBytes * 5 / 4
                                       + 2 * Layout::array_bytes(1));
    // wide rows: a few hundred a generation
    SetBits wide(1000);
    std::iota(wide.begin(), wide.end(), 0);
    for (Row r = 200000; r < 202000; ++r) {
        cache.insert(r, wide);
    }
    EXPECT_GT(cache.size(), 100u);
    EXPECT_LE(cache.bytes(), bound);
    EXPECT_EQ(cache.size() * row_copy_bytes(wide) + cache.table_bytes(), cache.bytes());
    EXPECT_LT(cache.bucket_count() * 4, narrow);
    cache.trim(0);
    EXPECT_EQ(0u, cache.bucket_count());
    EXPECT_EQ(0u, cache.table_bytes());
    cache.insert(1, SetBits{ 1 });
    EXPECT_GT(cache.bucket_count(), 0u);
    cache.clear();
    EXPECT_EQ(0u, cache.bucket_count());
}

// With the cache a budget-aware read holds less (its paths are cut at cached rows), so the
// lookahead under a memory budget, which reads ahead until a run does not fit, would read
// more with the cache than without it (UHGG 16S, annotate mode, 8 MiB: 75k rows warmed
// against 52-55k, 13-27% more instructions). Admitted as without the cache
// (admit_as_uncached), a read also charges, until its stored rows are read, the rows beyond
// the cached rows of its paths — per anchor its longest cached path's — and their stored rows:
// a row read alone holds then what decoding it to its anchor held at that point (that decode's
// peak, less at most its scratch and its reconstruction's peak beyond the stored rows). Its
// rows, costs and held bytes are those of the read without the mode
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

// The retention rule (RowDiffCache::keeps). Whatever it keeps — every row (checkpoint 1), the
// requested rows and checkpoints, almost nothing — the default and the budget-aware decodes
// return the rows of the decode without the cache, and the budget-aware one their costs; a
// call copies at most n x (successors + 3) + s / checkpoint rows (n rows asked for, s stored
// rows read) into the cache when no row counts as narrow, where keeping every row would copy
// all s; and the narrow rows are all kept
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
    EXPECT_EQ(stored.bytes() + cache.table_bytes(), cache.bytes());
    // a tuple row's bucket: its key and entry, and the neighbourhood bitmap (176 B), and the
    // one-entry table's array: 2 buckets and the padding of a neighbourhood
    EXPECT_EQ(HopscotchLayout<RowDiffCache<RowTuples>::Table>::array_bytes(2),
              cache.table_bytes());
    EXPECT_EQ(1u, cache.rows_inserted);
    EXPECT_EQ(stored.bytes(), cache.bytes_inserted);
    RowTuples empty;
    stored.store(empty);
    stored.load(&back);
    EXPECT_TRUE(back.empty());
}

#if USE_JEMALLOC
// the probe block of thread_allocated(), escaped so that the compiler cannot elide its
// allocation (a new/delete pair whose memory is never read may be removed)
char *volatile allocation_probe = nullptr;
#endif

// What the calling thread allocated so far (jemalloc's thread.allocated), or nullopt where
// jemalloc is not the process allocator (the counter does not move for a 1 MiB block)
std::optional<uint64_t> thread_allocated() {
#if USE_JEMALLOC
    auto read = []() {
        uint64_t allocated = 0;
        size_t sz = sizeof(allocated);
        return mallctl("thread.allocated", &allocated, &sz, nullptr, 0) ? 0 : allocated;
    };
    const uint64_t before = read();
    {
        std::unique_ptr<char[]> probe(new char[1 << 20]);
        allocation_probe = probe.get();
        if (read() < before + (1 << 20))
            return std::nullopt;
    }
    allocation_probe = nullptr;
    return read();
#else
    return std::nullopt;
#endif
}

// The cache admits a row before it copies it. A cache bounded at 1 KiB, its shared room 1 KiB
// too, given a row of 1,048,576 columns would otherwise allocate 4,210,688 bytes for a copy it
// then refuses, its bytes, peak and inserted rows all stating 0. A refused row is never
// copied: nothing is allocated (by jemalloc's counters, where it is the allocator, also at the
// moment the cache reads its room), and the refusal is counted (rows_refused); the same for a
// tuple row of 65,536 columns
TEST(RowDiffPathCache, ARefusedRowIsNeverCopied) {
    {
        const SetBits full(1u << 20, 42);
        RowDiffCache<SetBits> cache;
        std::optional<uint64_t> at_room;
        cache.set_bound(1024, [&]() {
            at_room = thread_allocated();
            return uint64_t(1024);
        });
        const std::optional<uint64_t> before = thread_allocated();
        cache.insert(1, full);
        const std::optional<uint64_t> after = thread_allocated();
        EXPECT_EQ(0u, cache.bytes());
        EXPECT_EQ(0u, cache.size());
        EXPECT_EQ(0u, cache.peak_bytes);
        EXPECT_EQ(0u, cache.rows_inserted);
        EXPECT_EQ(0u, cache.bytes_inserted);
        EXPECT_EQ(1u, cache.rows_refused);
        EXPECT_FALSE(cache.find(1, false));
        if (before && after && at_room) {
            // the measurements' own 1 MiB probe blocks excluded: thread_allocated() adds one
            // to the counter before it reads it
            const uint64_t probes = uint64_t(1) << 20;
            EXPECT_EQ(*before + 2 * probes, *after) << "the insert allocated "
                                                    << *after - *before - 2 * probes << " B";
            EXPECT_EQ(*before + probes, *at_room);
            std::cerr << "refused row of 1,048,576 columns under a 1 KiB bound: "
                      << *after - *before - 2 * probes << " B allocated" << std::endl;
        } else {
            std::cerr << "jemalloc is not the process allocator: the cache's counters only"
                      << std::endl;
        }
    }
    {
        RowTuples full;
        for (Column c = 0; c < (1u << 16); ++c) {
            full.emplace_back(c, RowTuples::value_type::second_type{ 1, 2, 3 });
        }
        RowDiffCache<RowTuples> cache;
        cache.set_bound(1024);
        const std::optional<uint64_t> before = thread_allocated();
        cache.insert(1, full);
        const std::optional<uint64_t> after = thread_allocated();
        EXPECT_EQ(0u, cache.bytes());
        EXPECT_EQ(0u, cache.peak_bytes);
        EXPECT_EQ(0u, cache.rows_inserted);
        EXPECT_EQ(1u, cache.rows_refused);
        // braced: gtest's EXPECT_* expands to an if/else, which GCC's -Wdangling-else rejects
        if (before && after) {
            EXPECT_EQ(*before + (uint64_t(1) << 20), *after);
        }
        // what a stored copy would hold is stated without making it
        StoredRow<RowTuples> stored;
        stored.store(full);
        EXPECT_EQ(stored.bytes(), StoredRow<RowTuples>::bytes_of(full));
    }
}

// ... and a kept row is copied after the evictions that make room for it: copied before
// them, a third row of 1.5 MiB would sit beside two held under a 4 MiB bound (4.5 MiB at
// once) until the oldest is dropped; the cache's heap never exceeds its bound during an
// insert (with its tables), measured by jemalloc's peak where it is the allocator, and the
// peak it states is what it held
TEST(RowDiffPathCache, AKeptRowIsCopiedWithinTheBound) {
    const uint64_t bound = uint64_t(4) << 20;
    SetBits wide(((uint64_t(3) << 19) - 64) / sizeof(Column));
    std::iota(wide.begin(), wide.end(), 0);
    const uint64_t b = row_copy_bytes(wide);
    // each generation's table: one entry, 2 buckets and the padding
    const uint64_t table = HopscotchLayout<
            RowDiffCache<SetBits>::Table>::array_bytes(2);
    ASSERT_LE(2 * (b + table), bound);
    ASSERT_GT(3 * b, bound);
    RowDiffCache<SetBits> cache;
    cache.set_bound(bound);
    cache.insert(1, wide);
    cache.insert(2, wide);
    EXPECT_EQ(2 * (b + table), cache.bytes());
#if USE_JEMALLOC
    size_t sz = sizeof(uint64_t);
    const bool peak_known = !mallctl("thread.peak.reset", nullptr, nullptr, nullptr, 0);
#else
    const bool peak_known = false;
#endif
    const uint64_t held_before = cache.bytes();
    // the third: the current generation (row 2) becomes the older one, row 1 is dropped
    cache.insert(3, wide);
    uint64_t peak = 0;
#if USE_JEMALLOC
    if (peak_known)
        mallctl("thread.peak.read", &peak, &sz, nullptr, 0);
#endif
    EXPECT_EQ(2 * (b + table), cache.bytes());
    EXPECT_FALSE(cache.find(1, false));
    EXPECT_EQ(3u, cache.rows_inserted);
    EXPECT_EQ(0u, cache.rows_refused);
    EXPECT_EQ(2 * (b + table), cache.peak_bytes);
    EXPECT_LE(cache.peak_bytes, bound);
    if (peak_known && thread_allocated()) {
        // what the cache held before the insert, and the most the insert added to it
        EXPECT_LE(held_before + peak, bound) << "the insert added up to " << peak << " B";
        EXPECT_LT(peak, b / 2);
        std::cerr << "kept row of " << b << " B under a " << bound << " B bound holding "
                  << held_before << " B: the insert added at most " << peak << " B" << std::endl;
    }
}

} // namespace
