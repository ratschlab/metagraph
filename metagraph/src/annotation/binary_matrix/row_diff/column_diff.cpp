#include "column_diff.hpp"

#include "annotation/binary_matrix/column_sparse/column_major.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
ColumnDiff<BaseMatrix>::ColumnDiff(const graph::DBGSuccinct *graph,
                                   BaseMatrix &&diffs,
                                   const std::string& terminal_file)
    : graph_(graph), diffs_(std::move(diffs)), terminal_file_(terminal_file) {
    load_terminal();
}

template <class BaseMatrix>
void ColumnDiff<BaseMatrix>::load_terminal() {
    std::ifstream f(terminal_file_, ios::binary);
    if (!f.good()) {
        common::logger->error("Could not open terminal file {}", terminal_file_);
        std::exit(1);
    }
    terminal_.load(f);
    f.close();
}

template <class BaseMatrix>
void ColumnDiff<BaseMatrix>::serialize(const std::string &name) const {
    std::ofstream f(name, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool ColumnDiff<BaseMatrix>::load(const std::string &name) {
    std::ifstream f(name, ios::binary);
    bool result = load(f);
    f.close();

    load_terminal();

    return result;
}

template <class BaseMatrix>
bool ColumnDiff<BaseMatrix>::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> ColumnDiff<BaseMatrix>::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

template <class BaseMatrix>
BinaryMatrix::SetBitPositions ColumnDiff<BaseMatrix>::get_row(Row row) const {
    Vector<uint64_t> result = get_diff(row);

    uint64_t boss_edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(row));
    const graph::boss::BOSS &boss = graph_->get_boss();

    while (!terminal_[row]) {
        graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
        assert(boss_edge > 1 && w != 0);

        // fwd always selects the last outgoing edge for a given node
        boss_edge = boss.fwd(boss_edge, w % boss.alph_size);
        row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                graph_->boss_to_kmer_index(boss_edge));
        merge(&result, get_diff(row));
    };
    return result;
}

void merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2) {
    if (diff2.empty()) {
        return;
    }
    Vector<uint64_t> diff1;
    std::swap(*result, diff1);
    uint64_t idx1 = 0;
    uint64_t idx2 = 0;
    while (idx1 < diff1.size() && idx2 < diff2.size()) {
        if (diff1[idx1] == diff2[idx2]) {
            idx1++;
            idx2++;
        } else if (diff1[idx1] < diff2[idx2]) {
            result->push_back(diff1[idx1]);
            idx1++;
        } else {
            result->push_back(diff2[idx2]);
            idx2++;
        }
    }
    while (idx1 < diff1.size()) {
        result->push_back(diff1[idx1]);
        idx1++;
    }
    while (idx2 < diff2.size()) {
        result->push_back(diff2[idx2]);
        idx2++;
    }
}

template
class ColumnDiff<ColumnMajor>;

} // namespace binmat
} // namespace annot
} // namespace mtg
