#include "row_diff.hpp"

namespace mtg {
namespace annot {
namespace binmat {


bool RowDiff::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    size_t i = 0;
    for (; i < set_bits.size() && column < set_bits[i]; ++i)
        ;
    return i < set_bits.size() && column == set_bits[i];
}

std::vector<BinaryMatrix::Row> RowDiff::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

uint64_t RowDiff::num_relations() const {
    uint64_t result = 0;
    for (Row row = 0; row < num_rows(); ++row) {
        result += get_row(row).size();
    }
    return result;
}

BinaryMatrix::SetBitPositions RowDiff::get_row(Row node_id) const {
    Vector<uint64_t> result = get_diff(node_id);
    while (!terminal()[node_id]) {
        uint64_t boss_edge = graph->kmer_to_boss_index(node_id);
        const graph::boss::BOSS &boss = graph->get_boss();
        graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
        assert(boss_edge > 1 && w != 0);
        graph::boss::BOSS::edge_index next_id = boss.fwd(boss_edge, w % boss.alph_size);

        node_id = graph->boss_to_kmer_index(next_id);
        merge(&result, get_diff(node_id));
    };
    return result;
}

bool RowDiff::load(std::istream &f)  {
    f >> num_columns_;
    diffs_.load(f);
    boundary_.load(f);
    terminal_.load(f);
    return true;
}

void RowDiff::serialize(std::ostream &f) const  {
    f << num_columns_;
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
    uint64_t start_idx = node_id == 0 ? 0 : sboundary.select(node_id) + 1;

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
        } else {
            if (diff1[idx1] < diff2[idx2]) {
                result->push_back(diff1[idx1]);
                idx1++;
            } else {
                result->push_back(diff2[idx2]);
                idx2++;
            }
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
