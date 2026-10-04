#pragma once

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/binary_matrix/row_diff/row_diff_cache.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

graph::DeBruijnGraph::node_index row_diff_successor(const graph::DeBruijnGraph &graph,
                                                    graph::DeBruijnGraph::node_index node,
                                                    const bit_vector &rd_succ);

namespace matrix {

const std::string kRowDiffAnchorExt = ".anchors";
const std::string kRowDiffForkSuccExt = ".rd_succ";

// Index matrices with the single-row reads the budget-aware decode path needs
// (BRWT, ColumnMajor: row_columns(); TupleCSCMatrix over them: row_tuples())
template <class M, class = void>
struct HasRowColumns : std::false_type {};
template <class M>
struct HasRowColumns<M, std::void_t<decltype(std::declval<const M&>().row_columns(
        BinaryMatrix::Row(), static_cast<ChargedBuffer<BinaryMatrix::Column>*>(nullptr),
        std::declval<DecodeBudget&>()))>> : std::true_type {};

template <class M, class = void>
struct HasRowTuples : std::false_type {};
template <class M>
struct HasRowTuples<M, std::void_t<decltype(std::declval<const M&>().get_binary_matrix().row_column_ranks(
        BinaryMatrix::Row(),
        static_cast<ChargedBuffer<std::pair<BinaryMatrix::Column, uint64_t>>*>(nullptr),
        std::declval<DecodeBudget&>())), decltype(&M::row_tuples)>> : std::true_type {};

// How the budget-aware decode path (IRowDiff::decode_rows) reads one STORED row (a diff,
// or an anchor's full row) of the index matrix: into the empty |*row|, every byte charged
// to |budget| before it is allocated, |*bytes| set to what is charged for |*row|. The
// scratch it reuses across the rows of one call is charged as it grows (scratch_charged(),
// released by the decoder at the end of the call); scratch_peak(n) is what it grows to for
// a row of n scratch entries (scratch_entries()) read alone, for the row's demand.
template <class RowT>
class RowFetcher {
  public:
    virtual ~RowFetcher() {}
    virtual bool fetch(BinaryMatrix::Row r, RowT *row, uint64_t *bytes,
                       DecodeBudget &budget) = 0;
    virtual uint64_t scratch_charged() const = 0;
    virtual uint64_t scratch_peak(uint64_t entries) const = 0;
    virtual uint64_t scratch_entries(const RowT &row) const = 0;
};

class IRowDiff {
  public:
    typedef bit_vector_small anchor_bv_type;
    typedef bit_vector_small fork_succ_bv_type;

    virtual ~IRowDiff() {}

    const graph::DeBruijnGraph* graph() const { return graph_; }
    void set_graph(const graph::DeBruijnGraph *graph) { graph_ = graph; }

    void load_fork_succ(const std::string &filename);
    void load_anchor(const std::string &filename);

    const anchor_bv_type& anchor() const { return anchor_; }

    const fork_succ_bv_type& fork_succ() const { return fork_succ_; }

    /**
     * The budget-aware decode path (decode_budget.hpp; DESIGN-traverse-graphlet.md §14,
     * stage 3 of §14.1), opt-in: the default functions (get_rows, get_row_tuples,
     * get_rows_dict, ...) are untouched by it and keep their speed and threading.
     *
     * decode_rows() / decode_row_tuples() return the rows of |rows| like get_rows() /
     * get_row_tuples(), single-threaded (call them outside any OpenMP region), charging
     * every buffer to |budget| before it is allocated: the row-diff trace, each dependency
     * row read alone by the index matrix's single-row descent, each coordinate tuple, the
     * reconstruction's buffers. A charge that does not fit refuses the call: REFUSED, with
     * |*out|, |*costs|, |*held| untouched and |budget|.held() as at entry. OK: |*costs|
     * gives each row's RowCost (its whole row-diff path, batch-independent), |*held| the
     * bytes charged for each returned row, and |budget|.held() grew by those plus
     * output_bytes<RowT>(rows.size()), the three vectors' buffers: the caller releases
     * them when it frees them. UNSUPPORTED: this matrix has no such path (its index matrix
     * has no single-row descent); nothing was charged. On OK the output vectors' previous
     * contents are replaced (and were never the budget's).
     */
    //
    // |cache| (optional): the request's path cache (row_diff_cache.hpp), as for
    // get_rows_cached(): a path also stops at a cached row whose path aggregates are known,
    // which then gives the costs of its whole path, so that every RowCost — and what the
    // caller admits and charges by it — is the same with and without the cache; only the
    // decoding (and the charges of its transient buffers) shrinks. A cached row is copied
    // into the call, charged as its copy; every row the call reconstructs is cached with its
    // aggregates within the cache's own bound, which is not this budget's (the caller's
    // allotment): also by a call that is refused afterwards (the rows are the same either way).
    // With the cache's admit_as_uncached set, the call also charges, from its trace until its
    // stored rows are read, what decoding the rows the cache spared it would have held: per
    // anchor the rows of its longest cached path beyond the cached row (in the trace's
    // containers and per-visit arrays) and that path's stored rows — so that a read that reads
    // ahead until a run does not fit (the lookahead) reads about what it read without the cache.
    virtual bool supports_budgeted_decode() const { return false; }
    virtual DecodeStatus decode_rows(const std::vector<BinaryMatrix::Row> &rows,
                                     DecodeBudget &budget,
                                     std::vector<BinaryMatrix::SetBitPositions> *out,
                                     std::vector<RowCost> *costs,
                                     std::vector<uint64_t> *held,
                                     RowDiffPathCache *cache = nullptr) const {
        return DecodeStatus::UNSUPPORTED;
    }
    virtual DecodeStatus decode_row_tuples(const std::vector<BinaryMatrix::Row> &rows,
                                           DecodeBudget &budget,
                                           std::vector<MultiIntMatrix::RowTuples> *out,
                                           std::vector<RowCost> *costs,
                                           std::vector<uint64_t> *held,
                                           RowDiffPathCache *cache = nullptr) const {
        return DecodeStatus::UNSUPPORTED;
    }

    /**
     * The default decode with the request's path cache (row_diff_cache.hpp; the efficiency
     * pass): the rows get_rows() / get_row_tuples() return, serial. Each row-diff path
     * stops at an anchor, at a row visited before in the call (as get_rd_ids), or at a row
     * |cache| holds, whose full row then starts the reconstruction as an anchor's stored row
     * does; every row the call reconstructs — the requested rows, the rows on their paths and
     * the anchors read — is cached within the cache's bound. Unsupported (supports_path_cache
     * false): std::logic_error.
     */
    virtual bool supports_path_cache() const { return false; }
    virtual std::vector<BinaryMatrix::SetBitPositions>
    get_rows_cached(const std::vector<BinaryMatrix::Row> &rows, RowDiffPathCache &cache) const {
        throw std::logic_error("this annotation has no row-diff path cache");
    }
    virtual std::vector<MultiIntMatrix::RowTuples>
    get_row_tuples_cached(const std::vector<BinaryMatrix::Row> &rows,
                          RowDiffPathCache &cache) const {
        throw std::logic_error("this annotation has no row-diff path cache");
    }
    // the buffers of a successful call's three output vectors for |n| rows
    template <class RowT>
    static uint64_t output_bytes(size_t n) {
        return buffer_bytes(n, sizeof(RowT)) + buffer_bytes(n, sizeof(RowCost))
                + buffer_bytes(n, sizeof(uint64_t));
    }

  protected:
    // The body of decode_rows() / decode_row_tuples() (row_diff_budgeted.cpp, instantiated
    // for SetBitPositions and RowTuples): serial, mirroring get_rd_ids() and call_rows()
    // with one group. Kept out of this header so that a translation unit using RowDiff does
    // not compile it.
    template <class RowT>
    DecodeStatus decode_budgeted(const std::vector<BinaryMatrix::Row> &rows,
                                 DecodeBudget &budget, RowFetcher<RowT> &fetcher,
                                 std::vector<RowT> *out, std::vector<RowCost> *costs,
                                 std::vector<uint64_t> *held,
                                 RowDiffCache<RowT> *cache = nullptr) const;

    // get row-diff paths starting at |row_ids|
    // Returns: (rd_ids, rd_paths_trunc, times_traversed, groups)
    // groups[g] records indices of row_ids paths traced in group g.
    // Rows from different groups access disjoint rd_rows entries, enabling
    // parallel reconstruction grouped by thread.
    std::tuple<std::vector<BinaryMatrix::Row>,
               std::vector<std::vector<size_t>>,
               std::vector<size_t>,
               std::vector<std::vector<size_t>>>
    get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids, size_t num_threads = 1) const;

    template <class F, class G, class H, class Callback>
    void call_rows(const std::vector<BinaryMatrix::Row> &row_ids, F call_rd_rows, G add_diff,
                   H decode_diffs, Callback call_row, size_t num_threads) const;

    // call_rows() with one thread and the path cache (get_rows_cached): |call_row| receives
    // each requested row as an rvalue
    template <class RowT, class F, class G, class H, class Callback>
    void call_rows_cached(const std::vector<BinaryMatrix::Row> &row_ids,
                          RowDiffCache<RowT> &cache, F call_rd_rows, G add_diff,
                          H decode_diffs, Callback call_row) const;

    const graph::DeBruijnGraph *graph_ = nullptr;
    anchor_bv_type anchor_;
    fork_succ_bv_type fork_succ_;
};

/**
 * Sparsified representation of the underlying #BinaryMatrix that stores diffs between
 * successive nodes, rather than the full annotation.
 * The successor of a node (that is the node to diff against) is determined by a path in
 * an external graph structure, a #graph::DBGSuccinct, which has the property that rows
 * on a path are likely identical or very similar.
 *
 * RowDiff sparsification can be applied to any BinaryMatrix instance.
 * The row-diff binary matrix is defined by three data structures:
 *   1. #diffs_ the underlying sparsified (diffed) #BinaryMatrix
 *   2. #anchor_ rows marked as anchors are stored in full
 *   3. #graph_ the graph that was used to determine adjacent rows for sparsification
 * Retrieving data from RowDiff requires the associated #graph_. In order to get the
 * annotation for  i-th row, we start traversing the node corresponding to i in #graph_
 * and accumulate the values in #diffs until we hit an anchor node, which is stored in
 * full.
 */
//NOTE: Clang aggressively abuses the clause in the C++ standard (14.7.1/11) that allows
// virtual methods in template classes to not be instantiated if unused and mistakenly
// does not instantiate the virtual methods in this class, so I had to move definitions
// to the header (gcc works fine)
template <class BaseMatrix>
class RowDiff : public IRowDiff, public BinaryMatrix {
  public:
    using base_matrix_type = BaseMatrix;

    template <typename... Args>
    RowDiff(const graph::DeBruijnGraph *graph = nullptr, Args&&... args)
        : diffs_(std::forward<Args>(args)...) { graph_ = graph; }

    /**
     * Returns the number of set bits in the row-diff transformed matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override;

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const override;

    /**
     * Return rows (in arbitrary order) and update the row indexes in |rows|
     * to point to their respective rows in the vector returned.
     *
     * In contrast to get_rows_dict in most other classes, the rows here
     * are not deduplicated. Benchmarks on real queries showed that merging
     * identical reconstructed rows rarely shrinks the batch size by much
     * (only a few percent on dense annotation rows; occasionally on the order
     * of ~40% when duplication was high), while hashing requires a critical
     * section that dominates the query time. Hence, we skip deduplication.
     */
    std::vector<SetBitPositions>
    get_rows_dict(std::vector<Row> *rows, size_t num_threads) const override;

    bool load(std::istream &f) override;
    void serialize(std::ostream &f) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    // the budget-aware decode path (IRowDiff), over BRWT and ColumnMajor
    bool supports_budgeted_decode() const override { return HasRowColumns<BaseMatrix>::value; }
    DecodeStatus decode_rows(const std::vector<Row> &rows, DecodeBudget &budget,
                             std::vector<SetBitPositions> *out, std::vector<RowCost> *costs,
                             std::vector<uint64_t> *held,
                             RowDiffPathCache *cache = nullptr) const override;

    // the default decode with the path cache (IRowDiff), over any index matrix
    bool supports_path_cache() const override { return true; }
    std::vector<SetBitPositions> get_rows_cached(const std::vector<Row> &rows,
                                                 RowDiffPathCache &cache) const override;

  private:
    static void add_diff(const SetBitPositions &diff, SetBitPositions *row);

    BaseMatrix diffs_;
};

// The stored rows of a binary row-diff matrix for the budget-aware decode path: the index
// matrix's single-row descent into the scratch, sorted, then copied into the row (charged
// at its exact size before the copy)
template <class BaseMatrix>
class BudgetedRowFetcher : public RowFetcher<BinaryMatrix::SetBitPositions> {
  public:
    explicit BudgetedRowFetcher(const BaseMatrix &matrix) : matrix_(matrix) {}

    bool fetch(BinaryMatrix::Row r, BinaryMatrix::SetBitPositions *row, uint64_t *bytes,
               DecodeBudget &budget) override {
        scratch_.data.clear();
        *bytes = 0;
        if (!matrix_.row_columns(r, &scratch_, budget))
            return false;
        std::sort(scratch_.data.begin(), scratch_.data.end());
        const uint64_t b = small_vector_bytes(scratch_.data.size(), sizeof(BinaryMatrix::Column));
        if (!budget.charge(b))
            return false;
        *bytes = b;
        row->assign(scratch_.data.begin(), scratch_.data.end());
        return true;
    }
    uint64_t scratch_charged() const override { return scratch_.charged(); }
    uint64_t scratch_peak(uint64_t entries) const override {
        return ChargedBuffer<BinaryMatrix::Column>::peak_for(entries);
    }
    uint64_t scratch_entries(const BinaryMatrix::SetBitPositions &row) const override {
        return row.size();
    }

  private:
    const BaseMatrix &matrix_;
    ChargedBuffer<BinaryMatrix::Column> scratch_;
};

template <class BaseMatrix>
DecodeStatus RowDiff<BaseMatrix>::decode_rows(const std::vector<Row> &rows, DecodeBudget &budget,
                                              std::vector<SetBitPositions> *out,
                                              std::vector<RowCost> *costs,
                                              std::vector<uint64_t> *held,
                                              RowDiffPathCache *cache) const {
    if constexpr(HasRowColumns<BaseMatrix>::value) {
        BudgetedRowFetcher<BaseMatrix> fetcher(diffs_);
        return decode_budgeted<SetBitPositions>(rows, budget, fetcher, out, costs, held,
                                                cache ? &cache->rows : nullptr);
    } else {
        return DecodeStatus::UNSUPPORTED;
    }
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows_cached(const std::vector<Row> &row_ids,
                                     RowDiffPathCache &cache) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows_cached<SetBitPositions>(row_ids, cache.rows,
        [this](const std::vector<Row> &rd_ids) { return diffs_.get_rows(rd_ids, 1); },
        add_diff, [](SetBitPositions *) {},
        [&](size_t i, SetBitPositions &&row) { rows[i] = std::move(row); }
    );
    return rows;
}


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    assert(graph_ && "graph must be loaded");
    assert(diffs_.num_rows() == graph_->max_index());
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    std::vector<Row> result;
    // TODO: implement a more efficient algorithm
    graph_->call_nodes([&](auto node) {
        auto row = graph::AnnotatedDBG::graph_to_anno_index(node);
        SetBitPositions set_bits = get_rows({ row })[0];
        if (std::binary_search(set_bits.begin(), set_bits.end(), column))
            result.push_back(row);
    });
    return result;
}

template <class F, class G, class H, class Callback>
void IRowDiff::call_rows(const std::vector<BinaryMatrix::Row> &row_ids, F call_rd_rows,
                         G add_diff, H decode_diffs, Callback call_row, size_t num_threads) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == graph_->max_index() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    if (row_ids.empty())
        return;

    // No sorting in order not to break the topological order for row-diff annotation

    // get row-diff paths
    Timer timer;
    // Unique row-diff row IDs to fetch and decode
    std::vector<BinaryMatrix::Row> rd_ids;
    // Truncated reconstruction paths (indices into rd_ids) per queried row
    std::vector<std::vector<size_t>> rd_paths_trunc;
    // Multiplicity for each queried row path returned by get_rd_ids()
    std::vector<size_t> times_traversed;
    // Independent groups of row paths that can be reconstructed in parallel
    std::vector<std::vector<size_t>> groups;
    std::tie(rd_ids, rd_paths_trunc, times_traversed, groups) = get_rd_ids(row_ids, num_threads);
    double rd_traversal_time = timer.elapsed();
    timer.reset();

    auto rd_rows = call_rd_rows(rd_ids, num_threads);
    double call_rd_rows_time = timer.elapsed();
    timer.reset();
    std::vector<BinaryMatrix::Row>().swap(rd_ids);

    size_t num_rd_bits = 0;
    size_t total_capacity = 0;
    // 200 rows per task, to make the task dispatch time negligible
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 200) \
        reduction(+:num_rd_bits,total_capacity)
    for (size_t i = 0; i < rd_rows.size(); ++i) {
        decode_diffs(&rd_rows[i]);
        std::sort(rd_rows[i].begin(), rd_rows[i].end(), utils::LessFirst());
        num_rd_bits += rd_rows[i].size();
        total_capacity += rd_rows[i].capacity();
    }

    double decode_diffs_time = timer.elapsed();
    timer.reset();

    // Reconstruct annotation rows from row-diff.
    // Since get_rd_ids skips cross-thread deduplication, rows from different
    // groups access disjoint rd_rows entries. We use that to reconstruct
    // each group in parallel, while preserving the sequential dependency
    // order within each group.
    using RowType = typename std::decay_t<decltype(rd_rows)>::value_type;

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t g = 0; g < groups.size(); ++g) {
        RowType result;
        for (size_t i : groups[g]) {
            auto it = rd_paths_trunc[i].rbegin();
            result = rd_rows[*it];
            if (!(--times_traversed[*it]))
                RowType().swap(rd_rows[*it]);
            for (++it ; it != rd_paths_trunc[i].rend(); ++it) {
                add_diff(rd_rows[*it], &result);
                if (--times_traversed[*it]) {
                    rd_rows[*it] = result;
                } else {
                    RowType().swap(rd_rows[*it]);
                }
            }
            call_row(i, result);
        }
    }

    common::logger->trace("RD query [threads: {}, rows: {} -> {} ({:.1f}x)] -- "
            "traversal: {:.2f} sec, call_rd_rows: {:.2f} sec (set bits: {}, capacity: {}), "
            "decoding: {:.2f} sec, reconstruction: {:.2f} sec",
            num_threads, row_ids.size(), rd_rows.size(), (double)rd_rows.size()/row_ids.size(),
            rd_traversal_time, call_rd_rows_time, num_rd_bits, total_capacity,
            decode_diffs_time, timer.elapsed());

    assert(times_traversed == std::vector<size_t>(rd_rows.size(), 0));
    assert(std::all_of(rd_rows.begin(), rd_rows.end(), [](const auto &v) { return v.empty(); }));
}

template <class RowT, class F, class G, class H, class Callback>
void IRowDiff::call_rows_cached(const std::vector<BinaryMatrix::Row> &row_ids,
                                RowDiffCache<RowT> &cache, F call_rd_rows, G add_diff,
                                H decode_diffs, Callback call_row) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == graph_->max_index() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);
    using Row = BinaryMatrix::Row;
    using node_index = graph::DeBruijnGraph::node_index;

    const size_t n = row_ids.size();
    if (!n)
        return;
    cache.trim(cache.limit());

    // ---- trace the paths (as get_rd_ids with one thread): each stops at a row visited
    // before in this call, at a cached row (looked up first: a cached anchor is not read
    // again) or at an anchor
    VectorSet<Row> visited;
    std::vector<uint8_t> cached;       // per visit: its full row is in |cache|
    std::vector<uint32_t> steps;       // the paths' visits, path by path
    std::vector<size_t> start(n + 1);
    for (size_t i = 0; i < n; ++i) {
        start[i] = steps.size();
        node_index node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row_ids[i]);
        while (true) {
            assert(graph_->in_graph(node));
            const Row row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);
            auto [it, is_new] = visited.emplace(row);
            steps.push_back(it - visited.begin());
            if (!is_new)
                break;
            cached.push_back(cache.find(row, false) != nullptr);
            if (cached.back() || anchor_[row])
                break;
            node = row_diff_successor(*graph_, node, fork_succ_);
        }
    }
    start[n] = steps.size();
    const auto &rows = visited.values_container();
    const size_t num_visits = rows.size();

    // ---- the stored rows of the visits not cached, read in one call in ascending row order
    // (as call_rows reads them), and the cached full rows, copied
    std::vector<uint32_t> order;
    order.reserve(num_visits);
    for (uint32_t v = 0; v < num_visits; ++v) {
        if (!cached[v])
            order.push_back(v);
    }
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return rows[a] < rows[b]; });
    std::vector<Row> rd_ids(order.size());
    for (size_t j = 0; j < order.size(); ++j) {
        rd_ids[j] = rows[order[j]];
    }
    std::vector<RowT> slot(num_visits);
    {
        auto stored = call_rd_rows(rd_ids);
        assert(stored.size() == order.size());
        for (size_t j = 0; j < order.size(); ++j) {
            RowT &row = slot[order[j]];
            row = std::move(stored[j]);
            decode_diffs(&row);
            std::sort(row.begin(), row.end(), utils::LessFirst());
        }
    }
    std::vector<uint8_t> known(num_visits, 0);   // its full row is cached or reconstructed
    for (uint32_t v = 0; v < num_visits; ++v) {
        if (cached[v]) {
            const auto *entry = cache.find(rows[v], false);
            assert(entry);
            slot[v] = entry->row;
            known[v] = 1;
            cache.hits++;
        }
    }
    std::vector<uint32_t> times(num_visits, 0);
    for (uint32_t v : steps) {
        times[v]++;
    }

    // ---- reconstruct in the requested order (as call_rows with one group), caching every
    // full row the call learns: the anchors read and each row on a path
    for (size_t i = 0; i < n; ++i) {
        const uint32_t *path = steps.data() + start[i];
        const uint32_t *path_end = steps.data() + start[i + 1];
        const uint32_t last = path_end[-1];
        RowT result;
        if (--times[last]) {
            result = slot[last];
        } else {
            result = std::move(slot[last]);
            RowT().swap(slot[last]);
        }
        if (!known[last]) {
            // an anchor read in this call: its stored row is its full row
            cache.insert(rows[last], result);
            known[last] = 1;
        }
        for (const uint32_t *p = path_end - 1; p != path; ) {
            const uint32_t y = *--p;
            add_diff(slot[y], &result);
            cache.insert(rows[y], result);
            known[y] = 1;
            if (--times[y]) {
                slot[y] = result;
            } else {
                RowT().swap(slot[y]);
            }
        }
        call_row(i, std::move(result));
    }
    assert(std::all_of(times.begin(), times.end(), [](uint32_t t) { return !t; }));
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](size_t i, const SetBitPositions &row) { rows[i] = row; },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    std::vector<SetBitPositions> rows_dict(rows->size());
    call_rows(*rows,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](size_t i, const SetBitPositions &row) {
            rows_dict[i] = row;
            (*rows)[i] = i;
        },
        num_threads
    );
    return rows_dict;
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(std::istream &f) {
    auto pos = f.tellg();
    std::string version(4, '\0');
    if (f.read(version.data(), 4) && version == "v2.0") {
        if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
            auto anchor_start = f.tellg();
            if (!anchor_.load(f) || !fork_succ_.load(f))
                return false;
            // anchor_ / fork_succ_ are accessed randomly, hint the kernel.
            utils::madvise_random_range(f, anchor_start, f.tellg() - anchor_start);
        }
    } else {
        // backward compatibility
        f.seekg(pos);
        if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
            auto anchor_start = f.tellg();
            if (!anchor_.load(f))
                return false;

            common::logger->warn(
                "Loading old version of RowDiff without a fork routing bitmap."
                " The last outgoing edges will be used as successors.");
            fork_succ_ = fork_succ_bv_type();
            // anchor_ is accessed randomly, hint the kernel.
            utils::madvise_random_range(f, anchor_start, f.tellg() - anchor_start);
        }
    }
    return diffs_.load(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(std::ostream &f) const {
    f.write("v2.0", 4);
    if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
        anchor_.serialize(f);
        fork_succ_.serialize(f);
    }
    diffs_.serialize(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::add_diff(const SetBitPositions &diff, SetBitPositions *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.empty())
        return;

    SetBitPositions result;
    result.reserve(row->size() + diff.size());
    std::set_symmetric_difference(row->begin(), row->end(),
                                  diff.begin(), diff.end(),
                                  std::back_inserter(result));
    row->swap(result);
}

} // namespace matrix
} // namespace annot
} // namespace mtg
