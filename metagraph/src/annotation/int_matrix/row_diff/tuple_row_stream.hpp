#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace mtg::annot::matrix {

// Merge shifted immutable spans, including modulo-2^64 sentinels and explicit
// empty-cell resets. Storage is proportional to path cells, never occurrences.
template <class Matrix>
DecodeStatus stream_tuple_row(const Matrix &matrix, BinaryMatrix::Row row,
        DecodeBudget &budget, RowCost *cost,
        const std::function<uint64_t(BinaryMatrix::Column, uint64_t)> &emit,
        const std::function<bool(uint64_t)> &stop) {
    using Row = BinaryMatrix::Row;
    using Column = BinaryMatrix::Column;
    using Values = typename std::decay_t<decltype(matrix.diffs())>::values_type;
    struct Input {
        Column column;
        const Values *values;
        uint64_t begin, end, pivot, at, shift;
        bool wrapped = false;
        uint64_t value() const { return (*values)[at] - shift; }
        bool next() {
            if (++at == end && !wrapped) { at = begin; wrapped = true; }
            return wrapped ? at < pivot : at < end;
        }
        bool seek(uint64_t wanted, uint64_t *work) {
            const uint64_t remaining = wrapped ? pivot - at : end - at + pivot - begin;
            const uint64_t tail = wrapped ? 0 : end - at;
            auto index = [&](uint64_t i) { return wrapped || i < tail ? at + i : begin + i - tail; };
            uint64_t lo = 0, hi = remaining;
            while (lo < hi) {
                ++*work;
                const auto mid = lo + (hi - lo) / 2;
                if ((*values)[index(mid)] - shift < wanted) lo = mid + 1;
                else hi = mid;
            }
            if (lo == remaining) return false;
            const uint64_t found = index(lo);
            wrapped |= lo >= tail;
            at = found;
            return true;
        }
    };
    const uint64_t entry = budget.held(), prior_peak = budget.begin_peak();
    struct Restore {
        DecodeBudget &budget; uint64_t entry, peak;
        ~Restore() { budget.restore(entry); budget.end_peak(peak); }
    } restore { budget, entry, prior_peak };
    ChargedBuffer<Row> path;
    ChargedBuffer<Input> inputs;
    ChargedBuffer<std::pair<Column, uint64_t>> ranks;
    uint64_t work = 0;
    auto halted = [&]() { return stop && stop(work); };
    if (row >= matrix.num_rows()) throw std::out_of_range("Stream row outside annotation");
    while (true) {
        if (halted()) return DecodeStatus::REFUSED;
        if (std::find(path.data.begin(), path.data.end(), row) != path.data.end())
            throw std::runtime_error("Coordinate RowDiff path has no anchor");
        if (!path.reserve_one(budget)) return DecodeStatus::REFUSED;
        path.data.push_back(row); ++work;
        if (matrix.anchor()[row]) break;
        if (!matrix.graph()) throw std::logic_error("Coordinate stream requires graph");
        const auto next = annot::row_diff_successor(*matrix.graph(), row + 1, matrix.fork_succ());
        if (!next || next - 1 >= matrix.num_rows())
            throw std::runtime_error("Invalid coordinate RowDiff successor");
        row = next - 1;
    }
    const uint64_t reset_bytes = buffer_bytes(matrix.num_columns(), sizeof(uint8_t));
    if (!budget.charge(reset_bytes)) return DecodeStatus::REFUSED;
    std::vector<uint8_t> reset(matrix.num_columns());
    for (size_t i = 0; i < path.data.size(); ++i) {
        if (halted()) return DecodeStatus::REFUSED;
        const uint64_t shift = i + (i + 1 < path.data.size());
        if (!matrix.diffs().call_row_tuple_views(path.data[i], &ranks, budget,
                [&](Column col, const Values &values, uint64_t begin, uint64_t end) {
            ++work;
            if (reset[col]) return true;
            if (begin == end) { reset[col] = true; return true; }
            if (!inputs.reserve_one(budget)) return false;
            const uint64_t pivot = std::lower_bound(values.begin() + begin,
                    values.begin() + end, shift) - values.begin();
            inputs.data.push_back({col, &values, begin, end, pivot,
                                   pivot == end ? begin : pivot, shift, pivot == end});
            return true;
        })) return DecodeStatus::REFUSED;
    }
    std::sort(inputs.data.begin(), inputs.data.end(),
              [](const Input &a, const Input &b) { return a.column < b.column; });
    using Head = std::pair<uint64_t, size_t>;
    ChargedBuffer<Head> heap;
    auto compare = std::greater<Head>();
    uint64_t next_poll = work;
    for (size_t first = 0; first < inputs.data.size();) {
        const Column column = inputs.data[first].column;
        size_t last = first;
        heap.data.clear();
        while (last < inputs.data.size() && inputs.data[last].column == column) {
            if (!heap.reserve_one(budget)) return DecodeStatus::REFUSED;
            heap.data.emplace_back(inputs.data[last].value(), last); ++last;
        }
        std::make_heap(heap.data.begin(), heap.data.end(), compare);
        while (!heap.data.empty()) {
            if (work >= next_poll) {
                if (halted()) return DecodeStatus::REFUSED;
                next_poll = work + 1024;
            }
            const uint64_t coordinate = heap.data.front().first;
            bool odd = false;
            do {
                std::pop_heap(heap.data.begin(), heap.data.end(), compare);
                const size_t id = heap.data.back().second;
                heap.data.pop_back(); odd = !odd; ++work;
                if (inputs.data[id].next()) {
                    heap.data.emplace_back(inputs.data[id].value(), id);
                    std::push_heap(heap.data.begin(), heap.data.end(), compare);
                }
            } while (!heap.data.empty() && heap.data.front().first == coordinate);
            if (!odd) continue;
            const uint64_t wanted = emit(column, coordinate);
            if (wanted <= coordinate) throw std::logic_error("Coordinate visitor must advance");
            if (wanted == std::numeric_limits<uint64_t>::max()) break;
            if (wanted == coordinate + 1) continue;
            size_t kept = 0;
            for (const Head &head : heap.data) {
                auto &input = inputs.data[head.second];
                if (input.seek(wanted, &work)) heap.data[kept++] = {input.value(), head.second};
            }
            heap.data.resize(kept);
            std::make_heap(heap.data.begin(), heap.data.end(), compare);
        }
        first = last;
    }
    *cost = {path.data.size() - 1, work, budget.peak() - entry};
    return DecodeStatus::OK;
}
} // namespace mtg::annot::matrix
