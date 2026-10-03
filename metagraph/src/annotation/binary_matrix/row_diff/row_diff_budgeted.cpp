// The budget-aware row-diff decode (DESIGN-traverse-graphlet.md §14, "row decoding is
// budget-aware inside the decoder"; stage 3 of §14.1): IRowDiff::decode_budgeted, the body
// of decode_rows() / decode_row_tuples(). See decode_budget.hpp for the contract.
#include "row_diff.hpp"

#include <algorithm>
#include <iterator>

#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace annot {
namespace matrix {

namespace {

using Row = BinaryMatrix::Row;
using Column = BinaryMatrix::Column;
using SetBitPositions = BinaryMatrix::SetBitPositions;
using RowTuples = MultiIntMatrix::RowTuples;
using Tuple = MultiIntMatrix::Tuple;
using node_index = graph::DeBruijnGraph::node_index;

// What a full copy of a row holds (a shared row kept for a later path, a row taken from
// one): copy construction allocates each buffer at its exact size from an empty one
uint64_t copy_bytes(const SetBitPositions &row) {
    return small_vector_bytes(row.size(), sizeof(Column));
}
uint64_t copy_bytes(const RowTuples &row) {
    uint64_t bytes = buffer_bytes(row.size(), sizeof(row[0]));
    for (const auto &entry : row) {
        bytes += small_vector_bytes(entry.second.size(), sizeof(uint64_t));
    }
    return bytes;
}

uint64_t entries_of(const SetBitPositions &row) { return row.size(); }
uint64_t entries_of(const RowTuples &row) {
    uint64_t n = row.size();
    for (const auto &entry : row) {
        n += entry.second.size();
    }
    return n;
}

// The reconstruction step of RowDiff::add_diff, with its allocation known in advance:
// the symmetric difference into a buffer reserved for |row| + |diff| (none for an empty
// diff, which changes nothing). xor_bound() is what it allocates, the new row's bytes.
uint64_t xor_bound(const SetBitPositions &row, const SetBitPositions &diff) {
    return diff.empty() ? 0 : small_vector_bytes(row.size() + diff.size(), sizeof(Column));
}
// returns the bytes of the bound freed again (none here)
uint64_t xor_into(const SetBitPositions &diff, SetBitPositions *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));
    if (diff.empty())
        return 0;
    SetBitPositions result;
    result.reserve(row->size() + diff.size());
    std::set_symmetric_difference(row->begin(), row->end(), diff.begin(), diff.end(),
                                  std::back_inserter(result));
    row->swap(result);
    return 0;
}

// The reconstruction step of TupleRowDiff::add_diff (the same result: columns merged, a
// column whose diff tuple is empty dropped, the coordinates of a column in both
// symmetric-differenced, every coordinate shifted by one), with each merged tuple reserved
// for both tuples at once instead of grown by back_inserter, so that what it allocates is
// known before: xor_bound(). A merged tuple that comes out empty is removed and its bytes
// returned (freed again).
uint64_t xor_bound(const RowTuples &row, const RowTuples &diff) {
    if (diff.empty())
        return 0;
    uint64_t bytes = buffer_bytes(row.size() + diff.size(), sizeof(row[0]));
    auto it = row.begin();
    auto it2 = diff.begin();
    while (it != row.end() && it2 != diff.end()) {
        if (it->first < it2->first) {
            bytes += small_vector_bytes(it->second.size(), sizeof(uint64_t));
            ++it;
        } else if (it->first > it2->first) {
            bytes += small_vector_bytes(it2->second.size(), sizeof(uint64_t));
            ++it2;
        } else {
            if (it2->second.size()) {
                bytes += small_vector_bytes(it->second.size() + it2->second.size(),
                                            sizeof(uint64_t));
            }
            ++it;
            ++it2;
        }
    }
    for (; it != row.end(); ++it) {
        bytes += small_vector_bytes(it->second.size(), sizeof(uint64_t));
    }
    for (; it2 != diff.end(); ++it2) {
        bytes += small_vector_bytes(it2->second.size(), sizeof(uint64_t));
    }
    return bytes;
}
uint64_t xor_into(const RowTuples &diff, RowTuples *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));
    uint64_t freed = 0;
    if (diff.size()) {
        RowTuples result;
        result.reserve(row->size() + diff.size());
        auto it = row->begin();
        auto it2 = diff.begin();
        while (it != row->end() && it2 != diff.end()) {
            if (it->first < it2->first) {
                result.push_back(*it);
                ++it;
            } else if (it->first > it2->first) {
                result.push_back(*it2);
                ++it2;
            } else {
                if (it2->second.size()) {
                    result.emplace_back(it->first, Tuple{});
                    Tuple &merged = result.back().second;
                    merged.reserve(it->second.size() + it2->second.size());
                    std::set_symmetric_difference(it->second.begin(), it->second.end(),
                                                  it2->second.begin(), it2->second.end(),
                                                  std::back_inserter(merged));
                    // as TupleRowDiff::add_diff: a column without coordinates is dropped
                    if (merged.empty()) {
                        freed += small_vector_bytes(it->second.size() + it2->second.size(),
                                                    sizeof(uint64_t));
                        result.pop_back();
                    }
                }
                ++it;
                ++it2;
            }
        }
        std::copy(it, row->end(), std::back_inserter(result));
        std::copy(it2, diff.end(), std::back_inserter(result));
        row->swap(result);
    }
    for (auto &[j, tuple] : *row) {
        for (uint64_t &c : tuple) {
            c -= 1;    // TupleRowDiff::SHIFT: coordinates increase by 1 at each edge
        }
    }
    return freed;
}

// The visited rows of the trace: an open-addressing table (row -> visit index) whose growth
// is charged before it allocates, at load at most 1/2, so that what it holds is a function
// of how many rows it indexes (peak_for()). A library hash map grows by policies (a full
// neighbourhood, an overflow list) that no model bounds as tightly.
class VisitIndex {
  public:
    static constexpr Row kEmpty = std::numeric_limits<Row>::max();

    // the visit of |row|, or UINT32_MAX
    uint32_t find(Row row) const {
        if (keys_.empty())
            return UINT32_MAX;
        for (size_t i = slot_of(row); ; i = (i + 1) & (keys_.size() - 1)) {
            if (keys_[i] == row)
                return vals_[i];
            if (keys_[i] == kEmpty)
                return UINT32_MAX;
        }
    }
    // |row| must not be indexed yet. false: the growth did not fit.
    bool insert(Row row, uint32_t visit, DecodeBudget &budget) {
        if (2 * (size_ + 1) > keys_.size()) {
            const size_t next = std::max<size_t>(16, 2 * keys_.size());
            const uint64_t bytes = table_bytes(next);
            if (!budget.charge(bytes))
                return false;
            std::vector<Row> keys(next, kEmpty);
            std::vector<uint32_t> vals(next);
            keys.swap(keys_);
            vals.swap(vals_);
            for (size_t i = 0; i < keys.size(); ++i) {
                if (keys[i] != kEmpty)
                    place(keys[i], vals[i]);
            }
            std::vector<Row>().swap(keys);
            std::vector<uint32_t>().swap(vals);
            budget.release(charged_);
            charged_ = bytes;
        }
        place(row, visit);
        ++size_;
        return true;
    }
    uint64_t charged() const { return charged_; }

    // the most indexing |n| rows from empty holds at once: the last two tables
    static uint64_t peak_for(uint64_t n) {
        size_t prev = 0, cap = 0;
        for (uint64_t size = 0; size < n; ++size) {
            if (2 * (size + 1) > cap) {
                prev = cap;
                cap = std::max<size_t>(16, 2 * cap);
            }
        }
        return table_bytes(cap) + table_bytes(prev);
    }

  private:
    static uint64_t table_bytes(size_t n) {
        return buffer_bytes(n, sizeof(Row)) + buffer_bytes(n, sizeof(uint32_t));
    }
    size_t slot_of(Row row) const {
        // Fibonacci hashing: row ids of one path are scattered, but share low bits
        return (row * 0x9E3779B97F4A7C15ULL) >> (64 - __builtin_ctzll(keys_.size()));
    }
    void place(Row row, uint32_t visit) {
        size_t i = slot_of(row);
        while (keys_[i] != kEmpty) {
            i = (i + 1) & (keys_.size() - 1);
        }
        keys_[i] = row;
        vals_[i] = visit;
    }

    std::vector<Row> keys_;
    std::vector<uint32_t> vals_;
    size_t size_ = 0;
    uint64_t charged_ = 0;
};

// One row the trace visited, with the aggregates over its whole row-diff path from it to
// the anchor (computed when it is first reconstructed, from its successor's), which make
// its costs batch-independent: every row on a requested row's path is visited in the
// same call, and a path stops only at an anchor or at a row reconstructed before.
struct Visit {
    Row row;
    uint64_t stored = 0;        // bytes charged for its stored row (a diff or a full anchor)
    uint64_t entries = 0;       // entries of its stored row (columns, coordinates)
    uint64_t scratch = 0;       // scratch entries reading it took
    uint64_t slot_bytes = 0;    // bytes charged for what its slot holds now
    uint32_t times = 0;         // paths still to use it (get_rd_ids' times_traversed)
    uint32_t length = 0;        // rows on its path, itself included (0: not reconstructed)
    uint64_t dep_entries = 0;   // entries the other rows of its path store
    uint64_t path_stored = 0;   // bytes of the stored rows of its path
    uint64_t max_scratch = 0;   // the most scratch entries a row of its path took
    uint64_t recon_peak = 0;    // the most the reconstruction of its row holds beyond them
    uint64_t full = 0;          // bytes charged for its reconstructed row
};

} // namespace

template <class RowT>
DecodeStatus IRowDiff::decode_budgeted(const std::vector<Row> &rows, DecodeBudget &budget,
                                       RowFetcher<RowT> &fetcher, std::vector<RowT> *out,
                                       std::vector<RowCost> *costs,
                                       std::vector<uint64_t> *held) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == graph_->max_index() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    const uint64_t at_entry = budget.held();
#ifndef NDEBUG
    const uint64_t peak_before = budget.begin_peak();
#endif
    auto refuse = [&]() {
        // every local container is freed by its scope: only the account is restored
        budget.restore(at_entry);
#ifndef NDEBUG
        budget.end_peak(peak_before);
#endif
        return DecodeStatus::REFUSED;
    };
    const size_t n = rows.size();

    // the outputs and the paths' offsets, at their exact sizes
    const uint64_t fixed = output_bytes<RowT>(n) + buffer_bytes(n + 1, sizeof(uint32_t));
    if (!budget.charge(fixed))
        return refuse();
    std::vector<RowT> result_rows(n);
    std::vector<RowCost> result_costs(n);
    std::vector<uint64_t> result_held(n);
    std::vector<uint32_t> start(n + 1);

    // ---- trace the row-diff paths (as get_rd_ids with one thread): each path stops at an
    // anchor or at a row visited before, whose reconstruction precedes it
    VisitIndex index;
    ChargedBuffer<Visit> visits;
    ChargedBuffer<uint32_t> steps;
    for (size_t i = 0; i < n; ++i) {
        start[i] = steps.data.size();
        node_index node = graph::AnnotatedSequenceGraph::anno_to_graph_index(rows[i]);
        while (true) {
            assert(graph_->in_graph(node));
            const Row row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);
            uint32_t v = index.find(row);
            const bool is_new = v == UINT32_MAX;
            if (is_new) {
                v = visits.data.size();
                if (!visits.reserve_one(budget) || !index.insert(row, v, budget))
                    return refuse();
                visits.data.push_back(Visit{ row });
            }
            if (!steps.reserve_one(budget))
                return refuse();
            steps.data.push_back(v);
            visits.data[v].times++;
            if (!is_new || anchor_[row])
                break;
            node = row_diff_successor(*graph_, node, fork_succ_);
        }
    }
    start[n] = steps.data.size();
    const size_t num_visits = visits.data.size();

    // ---- read every stored row, in ascending row order: the order in which the index
    // matrix's descent is as fast as a batched read (a random order is twice as slow)
    const uint64_t per_visit = buffer_bytes(num_visits, sizeof(RowT))
                             + buffer_bytes(num_visits, sizeof(uint32_t));
    if (!budget.charge(per_visit))
        return refuse();
    std::vector<RowT> slot(num_visits);
    std::vector<uint32_t> order(num_visits);
    for (uint32_t v = 0; v < num_visits; ++v) {
        order[v] = v;
    }
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return visits.data[a].row < visits.data[b].row;
    });
    for (uint32_t v : order) {
        Visit &visit = visits.data[v];
        if (!fetcher.fetch(visit.row, &slot[v], &visit.stored, budget))
            return refuse();
        visit.slot_bytes = visit.stored;
        visit.entries = entries_of(slot[v]);
        visit.scratch = fetcher.scratch_entries(slot[v]);
    }

    // ---- reconstruct in the requested order (as call_rows with one group); a stored row
    // is freed once its last path has used it, a row a later path ends at is kept whole
    auto take = [&](Visit &visit, RowT &from, RowT *row, uint64_t *row_bytes) {
        if (--visit.times) {
            const uint64_t bytes = copy_bytes(from);
            if (!budget.charge(bytes))
                return false;
            *row = from;
            *row_bytes = bytes;
        } else {
            *row = std::move(from);
            RowT().swap(from);
            *row_bytes = visit.slot_bytes;
            visit.slot_bytes = 0;
        }
        return true;
    };
    for (size_t i = 0; i < n; ++i) {
        const uint32_t *path = steps.data.data() + start[i];
        const uint32_t *path_end = steps.data.data() + start[i + 1];
        const uint32_t end = path_end[-1];
        Visit &last = visits.data[end];
        if (!last.length) {
            // an anchor first reached here: its stored row is its full row
            assert(anchor_[last.row]);
            last.length = 1;
            last.dep_entries = 0;
            last.path_stored = last.stored;
            last.max_scratch = last.scratch;
            last.recon_peak = last.stored;
            last.full = last.stored;
        }
        RowT result;
        uint64_t result_bytes = 0;
        if (!take(last, slot[end], &result, &result_bytes))
            return refuse();
        uint32_t succ = end;
        for (const uint32_t *p = path_end - 1; p != path; ) {
            const uint32_t y = *--p;
            Visit &visit = visits.data[y];
            const Visit &next = visits.data[succ];
            const uint64_t bound = xor_bound(result, slot[y]);
            if (!budget.charge(bound))
                return refuse();
            const bool changes = !slot[y].empty();
            const uint64_t freed = xor_into(slot[y], &result);
            if (changes) {
                // the new row replaced the old one, which is freed
                budget.release(result_bytes + freed);
                result_bytes = bound - freed;
            } else {
                budget.release(bound);
            }
            // the aggregates of |y| from its successor's (its first reconstruction: a row
            // is reconstructed once per call, on the first path through it)
            assert(!visit.length);
            visit.length = next.length + 1;
            visit.dep_entries = next.dep_entries + next.entries;
            visit.path_stored = next.path_stored + visit.stored;
            visit.max_scratch = std::max(next.max_scratch, visit.scratch);
            visit.recon_peak = std::max(next.recon_peak, next.full + bound);
            // From the path, not from what this call built the row from: an empty diff keeps
            // the incoming row, whose bytes are its successor's aggregate when the row is
            // decoded alone but an exact-size copy when the path continued from a row shared
            // earlier in the call — the demand must not depend on the batch (review F4)
            visit.full = changes ? result_bytes : next.full;
            // the stored diff is used up; a row a later path ends at keeps the full row
            budget.release(visit.slot_bytes);
            visit.slot_bytes = 0;
            RowT().swap(slot[y]);
            if (--visit.times) {
                const uint64_t bytes = copy_bytes(result);
                if (!budget.charge(bytes))
                    return refuse();
                slot[y] = result;
                visit.slot_bytes = bytes;
            }
            succ = y;
        }
        const Visit &first = visits.data[path[0]];
        RowCost &cost = result_costs[i];
        cost.dependency_rows = first.length - 1;
        cost.dependency_entries = first.dep_entries;
        // what decoding this row alone holds at most: the call's fixed buffers for one row,
        // the trace's containers grown to its path, the per-visit arrays, the stored rows
        // of its path, the scratch grown to its widest row, and the reconstruction's peak
        // beyond the stored rows (a row and the buffer of its next step)
        const uint64_t length = first.length;
        cost.demand = output_bytes<RowT>(1) + buffer_bytes(2, sizeof(uint32_t))
                    + VisitIndex::peak_for(length)
                    + ChargedBuffer<Visit>::peak_for(length)
                    + ChargedBuffer<uint32_t>::peak_for(length)
                    + buffer_bytes(length, sizeof(RowT)) + buffer_bytes(length, sizeof(uint32_t))
                    + first.path_stored + fetcher.scratch_peak(first.max_scratch)
                    + first.recon_peak;
        result_rows[i] = std::move(result);
        result_held[i] = result_bytes;
    }

    // ---- what the call held besides its outputs is freed with its scope
    uint64_t transient = index.charged() + visits.charged() + steps.charged() + per_visit
                       + fetcher.scratch_charged() + buffer_bytes(n + 1, sizeof(uint32_t));
    for (const Visit &visit : visits.data) {
        assert(!visit.times && !visit.slot_bytes);
        transient += visit.slot_bytes;
    }
    budget.release(transient);
#ifndef NDEBUG
    {
        uint64_t kept = output_bytes<RowT>(n);
        for (uint64_t b : result_held) {
            kept += b;
        }
        assert(budget.held() == at_entry + kept);
        // a row decoded alone stays within the demand it states (which the consumers admit
        // every returned row against)
        assert(n != 1 || budget.peak() - at_entry <= result_costs[0].demand);
        budget.end_peak(peak_before);
    }
#endif
    out->swap(result_rows);
    costs->swap(result_costs);
    held->swap(result_held);
    return DecodeStatus::OK;
}

template DecodeStatus
IRowDiff::decode_budgeted<SetBitPositions>(const std::vector<Row> &, DecodeBudget &,
                                           RowFetcher<SetBitPositions> &,
                                           std::vector<SetBitPositions> *,
                                           std::vector<RowCost> *,
                                           std::vector<uint64_t> *) const;
template DecodeStatus
IRowDiff::decode_budgeted<RowTuples>(const std::vector<Row> &, DecodeBudget &,
                                     RowFetcher<RowTuples> &, std::vector<RowTuples> *,
                                     std::vector<RowCost> *, std::vector<uint64_t> *) const;

} // namespace matrix
} // namespace annot
} // namespace mtg
