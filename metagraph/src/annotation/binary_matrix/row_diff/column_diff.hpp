#pragma once

#include <fstream>

#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <sdsl/rrr_vector.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg {
namespace annot {
namespace binmat {

/**
 * Row-major representation of annotations where rows are stored as differences vs a
 * predecessor row. The predecessor of a row is determined by a path in an external graph
 * structure, a #graph::DBGSuccinct, which has the property that rows on a path are likely
 * very similar.
 * The row-diff binary matrix is defined by three data structures:
 *   1. #diffs_ stores the concatenated diffs row by row. Since each row-diff has variable
 *      length, the boundary_ vector marks the end of each row-diff.
 *   2. boundary_ marks the end of a row diff in #diffs_
 *   3. terminal_ Rows marked as terminal are stored in full.
 * Retrieving data from ColumnDiff requires the associated #graph_. In order to get the i-th
 * row, we start traversing the node corresponding to i in #graph_ and accumulate the
 * diffs until we hit a terminal node, which is stored in full.
 */
//NOTE: Clang aggressively abuses the clause in the C++ standard (14.7.1/11) that allows
// virtual methods in template classes to not be instantiated if unused and mistakenly
// does not instantiate the virtual methods in this class, so I had to move definitions
// to the header
 template <class BaseMatrix>
class ColumnDiff : public BinaryMatrix {
  public:
    ColumnDiff() {}

    ColumnDiff(const graph::DBGSuccinct *graph,
               BaseMatrix &&diffs,
               const std::string &terminal_file)
        : graph_(graph), diffs_(std::move(diffs)), terminal_file_(terminal_file) {
        load_terminal(terminal_file_, &terminal_);
    }

    uint64_t num_columns() const override { return diffs_.num_columns(); }

    const std::string& terminal_file() const { return terminal_file_; }

    /**
     * Returns the number of set bits in the matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }

    uint64_t num_rows() const override { return diffs_.num_rows(); }
    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

    bool get(Row row, Column column) const override {
        SetBitPositions set_bits = get_row(row);
        SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
        return v != set_bits.end() && *v == column;
    }

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override  {
        std::vector<Row> result;
        for (Row row = 0; row < num_rows(); ++row) {
            if (get(row, column))
                result.push_back(row);
        }
        return result;
    }

    SetBitPositions get_row(Row row) const override  {
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

    bool load(std::istream &f) override {
        uint64_t len;
        common::logger->trace("Loading terminal nodes from {}", terminal_file_);
        f.read(reinterpret_cast<char *>(&len), sizeof(uint64_t));
        terminal_file_ = std::string(len, '\0');
        f.read(terminal_file_.data(), len);

        return diffs_.load(f);
    }

    void serialize(std::ostream &f) const override {
        uint64_t len = terminal_file_.size();
        f.write(reinterpret_cast<char *>(&len), sizeof(uint64_t));
        f.write(terminal_file_.c_str(), len);
        diffs_.serialize(f);
    };

    void serialize(const std::string &name) const;
    bool load(const std::string &name);

    using ColumnCallback = typename BaseMatrix::ColumnCallback;
    void call_columns(const  ColumnCallback &callback) {
        diffs_.call_columns(callback);
    }

    const sdsl::rrr_vector<> &terminal() const { return terminal_; }

    const BaseMatrix &diffs() const { return diffs_; }
    BaseMatrix &diffs() { return diffs_; }

    Vector<uint64_t> get_diff(uint64_t node_id) const { return diffs_.get_row(node_id); }

  private:
    static void load_terminal(const std::string &filename, sdsl::rrr_vector<> *terminal) {
        std::ifstream f(filename, ios::binary);
        if (!f.good()) {
            common::logger->error("Could not open terminal file {}", filename);
            std::exit(1);
        }
        terminal->load(f);
        f.close();
    }

    static void merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2) {
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

    const graph::DBGSuccinct *graph_;

    BaseMatrix diffs_;
    sdsl::rrr_vector<> terminal_;

    std::string terminal_file_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg
