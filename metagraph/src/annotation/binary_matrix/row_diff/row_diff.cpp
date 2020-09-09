#include "row_diff.hpp"

#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace annot {
namespace binmat {


bool RowDiff::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<BinaryMatrix::Row> RowDiff::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions RowDiff::get_row(Row row) const {
    Vector<uint64_t> result = get_diff(row);

    uint64_t boss_edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(row));
    const graph::boss::BOSS &boss = graph_->get_boss();

    while (!terminal()[row]) {
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

bool RowDiff::load(std::istream &f)  {
    f >> num_columns_;
    f >> num_relations_;
    diffs_.load(f);
    boundary_.load(f);
    terminal_.load(f);
    sdsl::util::init_support(sboundary_, &boundary_);
    return true;
}

void RowDiff::serialize(std::ostream &f) const  {
    f << num_columns_;
    f << num_relations_;
    diffs_.serialize(f);
    boundary_.serialize(f);
    terminal_.serialize(f);
}

void RowDiff::serialize(const std::string &name) const {
    std::ofstream f(name, ios::binary);
    serialize(f);
    f.close();
}

bool RowDiff::load(const std::string &name) {
    std::ifstream f(name, ios::binary);
    bool result = load(f);
    f.close();
    return result;
}

Vector<uint64_t> RowDiff::get_diff(uint64_t node_id) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    Vector<uint64_t> result;
    uint64_t start_idx = node_id == 0 ? 0 : sboundary_.select(node_id) + 1;

    while (boundary_[start_idx] == 0) {
        result.push_back(diffs_[start_idx - node_id]);
        start_idx++;
    }
    return result;
}

void RowDiff::merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2) {
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

} // namespace binmat
} // namespace annot
} // namespace mtg
