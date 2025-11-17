#include "brwt.hpp"

#include <cstring>
#include <queue>
#include <numeric>

#include <progress_bar.hpp>

#include "common/algorithms.hpp"
#include "common/serialization.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace annot {
namespace matrix {

const size_t kNumRowsInBlock = 250'000;


bool BRWT::get(Row row, Column column) const {
    assert(row < num_rows());
    assert(column < num_columns());

    // if leaf
    if (!child_nodes_.size())
        return (*nonzero_rows_)[row];

    uint64_t rank = nonzero_rows_->conditional_rank1(row);
    // terminate if the index bit is unset
    if (!rank)
        return false;

    // Check if this rank position is all-ones (not in nonones_rows_)
    // nonones_rows_ is indexed by rank in nonzero_rows_ (0-indexed)
    if (nonones_rows_ && !(rank = nonones_rows_->conditional_rank1(rank - 1)))
        return true; // All-ones row: return true for all columns

    auto child_node = assignments_.group(column);
    return child_nodes_[child_node]->get(rank - 1, assignments_.rank(column));
}

std::vector<BRWT::SetBitPositions>
BRWT::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    Vector<Column> slice;
    // expect at least 3 relations per row
    slice.reserve(row_ids.size() * 4);

    slice_rows(row_ids, &slice);

    assert(slice.size() >= row_ids.size());

    auto row_begin = slice.begin();

    for (size_t i = 0; i < rows.size(); ++i) {
        // every row in `slice` ends with `-1`
        auto row_end = std::find(row_begin, slice.end(),
                                 std::numeric_limits<Column>::max());
        rows[i].assign(row_begin, row_end);
        row_begin = row_end + 1;
    }

    return rows;
}

void BRWT::slice_rows(Row begin, Row end, Vector<Column> *slice) const {
    // TODO: It may be faster if index columns are queried in ranges instead of with element-wise
    // access queries.
    slice_rows(utils::arange<Row>(begin, end - begin), slice);
}

void BRWT::call_rows(const std::function<void(const SetBitPositions &)> &callback,
                     bool show_progress) const {
    Vector<Column> slice;
    ProgressBar progress_bar(num_rows(), "Queried BRWT rows", std::cerr, !show_progress);

    #pragma omp parallel for ordered num_threads(get_num_threads()) schedule(dynamic) private(slice)
    for (uint64_t i = 0; i < num_rows(); i += kNumRowsInBlock) {
        uint64_t begin = i;
        uint64_t end = std::min(i + kNumRowsInBlock, num_rows());
        assert(begin <= end);

        slice.resize(0);
        slice_rows(begin, end, &slice);

        #pragma omp ordered
        {
            SetBitPositions row;
            for (auto row_begin = slice.begin(); row_begin < slice.end(); ) {
                // every row in `slice` ends with `-1`
                auto row_end = std::find(row_begin, slice.end(),
                                         std::numeric_limits<Column>::max());
                row.assign(row_begin, row_end);
                std::sort(row.begin(), row.end());
                callback(row);
                ++progress_bar;
                row_begin = row_end + 1;
            }
        }
    }
}


std::vector<Vector<std::pair<BRWT::Column, uint64_t>>>
BRWT::get_column_ranks(const std::vector<Row> &row_ids) const {
    std::vector<Vector<std::pair<Column, uint64_t>>> rows(row_ids.size());

    Vector<std::pair<Column, uint64_t>> slice;
    // expect at least 3 relations per row
    slice.reserve(row_ids.size() * 4);

    slice_rows(row_ids, &slice);

    assert(slice.size() >= row_ids.size());

    auto row_begin = slice.begin();

    for (size_t i = 0; i < rows.size(); ++i) {
        // every row in `slice` ends with `-1`
        auto row_end = row_begin;
        while (row_end->first != std::numeric_limits<Column>::max()) {
            ++row_end;
            assert(row_end != slice.end());
        }
        rows[i].assign(row_begin, row_end);
        row_begin = row_end + 1;
    }

    return rows;
}

// If T = Column
//      return positions of set bits.
// If T = std::pair<Column, uint64_t>
//      return positions of set bits with their column ranks.
// Appends to `slice`
template <typename T>
void BRWT::slice_rows(const std::vector<Row> &row_ids, Vector<T> *slice) const {
    T delim;
    if constexpr(utils::is_pair_v<T>) {
        delim = std::make_pair(std::numeric_limits<Column>::max(), 0);
    } else {
        delim = std::numeric_limits<Column>::max();
    }

    // check if this is a leaf
    if (!child_nodes_.size()) {
        assert(assignments_.size() == 1);

        for (Row i : row_ids) {
            assert(i < num_rows());

            if constexpr(utils::is_pair_v<T>) {
                if (uint64_t rank = nonzero_rows_->conditional_rank1(i)) {
                    // only a single column is stored in leaves
                    slice->emplace_back(0, rank);
                }
            } else {
                if ((*nonzero_rows_)[i]) {
                    // only a single column is stored in leaves
                    slice->push_back(0);
                }
            }
            slice->push_back(delim);
        }

        return;
    }

    // construct indexing for children and the inverse mapping
    std::vector<Row> child_row_ids;
    child_row_ids.reserve(row_ids.size());

    std::vector<bool> skip_row(row_ids.size(), true);

    for (size_t i = 0; i < row_ids.size(); ++i) {
        assert(row_ids[i] < num_rows());

        uint64_t global_offset = row_ids[i];

        // if next word contains 5 or more positions, query the whole word
        // we assume that get_int is roughly 5 times slower than operator[]
        if (i + 4 < row_ids.size()
                && row_ids[i + 4] < global_offset + 64
                && row_ids[i + 4] >= global_offset
                && global_offset + 64 <= nonzero_rows_->size()) {
            // get the word
            uint64_t word = nonzero_rows_->get_int(global_offset, 64);
            uint64_t rank = -1ULL;

            do {
                // check index
                uint8_t offset = row_ids[i] - global_offset;
                if (word & (1ULL << offset)) {
                    if (rank == -1ULL)
                        rank = global_offset > 0
                                ? nonzero_rows_->rank1(global_offset - 1)
                                : 0;

                    // map index from parent's to children's coordinate system
                    child_row_ids.push_back(rank + sdsl::bits::cnt(word & sdsl::bits::lo_set[offset + 1]) - 1);
                    skip_row[i] = false;
                }
            } while (++i < row_ids.size()
                        && row_ids[i] < global_offset + 64
                        && row_ids[i] >= global_offset);
            --i;

        } else {
            // check index
            if (uint64_t rank = nonzero_rows_->conditional_rank1(global_offset)) {
                // map index from parent's to children's coordinate system
                child_row_ids.push_back(rank - 1);
                skip_row[i] = false;
            }
        }
    }

    std::vector<bool> is_all_ones_row(child_row_ids.size(), false);
    // Now update child_row_ids using nonones_rows_ to exclude all-ones rows
    if (nonones_rows_) {
        for (size_t i = 0, j = 0; i < row_ids.size() && j < child_row_ids.size(); ++i) {
            if (skip_row[i])
                continue;
            // Check if this is an all-ones row
            uint64_t rank = nonones_rows_->conditional_rank1(child_row_ids[j]);
            if (!rank) {
                // All-ones row: mark it and remove from child_row_ids
                is_all_ones_row[j] = true;
            } else {
                // Update to adjusted rank
                child_row_ids[j] = rank - 1;
            }
            ++j;
        }
        utils::erase(&child_row_ids, is_all_ones_row);
    }

    if (!child_row_ids.size()) {
        for (size_t i = 0; i < row_ids.size(); ++i) {
            if (!skip_row[i]) {
                // Handle all-ones rows: return all columns for these rows
                for (Column col = 0; col < num_columns(); ++col) {
                    if constexpr(utils::is_pair_v<T>) {
                        throw std::runtime_error("Stripping all-ones rows from BRWT with attributes is not supported");
                    } else {
                        slice->push_back(col);
                    }
                }
            }
            slice->push_back(delim);
        }
        return;
    }

    // TODO: query by columns and merge them in the very end to avoid remapping
    //       the same column indexes many times when propagating to the root.
    // TODO: implement a cache efficient method for merging the columns.

    // query all children subtrees and get relations from them
    size_t slice_start = slice->size();

    std::vector<size_t> pos(child_nodes_.size());

    for (size_t j = 0; j < child_nodes_.size(); ++j) {
        pos[j] = slice->size();
        child_nodes_[j]->slice_rows<T>(child_row_ids, slice);

        assert(slice->size() >= pos[j] + child_row_ids.size());

        // transform column indexes
        for (size_t i = pos[j]; i < slice->size(); ++i) {
            auto &v = (*slice)[i];
            if (v != delim) {
                auto &col = utils::get_first(v);
                col = assignments_.get(j, col);
            }
        }
    }

    size_t slice_offset = slice->size();

    for (size_t i = 0, j = 0; i < row_ids.size(); ++i) {
        if (!skip_row[i]) {
            if (!is_all_ones_row[j++]) {
                // merge rows from child submatrices
                for (size_t &p : pos) {
                    while ((*slice)[p++] != delim) {
                        slice->push_back((*slice)[p - 1]);
                    }
                }
            } else {
                // Handle all-ones rows: return all columns for these rows
                for (Column col = 0; col < num_columns(); ++col) {
                    if constexpr(utils::is_pair_v<T>) {
                        throw std::runtime_error("Stripping all-ones rows from BRWT with attributes is not supported");
                    } else {
                        slice->push_back(col);
                    }
                }
            }
        }
        slice->push_back(delim);
    }

    slice->erase(slice->begin() + slice_start, slice->begin() + slice_offset);
}

std::vector<BRWT::Row> BRWT::get_column(Column column) const {
    assert(column < num_columns());

    auto num_nonzero_rows = nonzero_rows_->num_set_bits();

    // check if the column is empty
    if (!num_nonzero_rows)
        return {};

    // check whether it is a leaf
    if (!child_nodes_.size()) {
        // return the index column
        std::vector<BRWT::Row> result;
        result.reserve(num_nonzero_rows);
        nonzero_rows_->call_ones([&](auto i) { result.push_back(i); });
        return result;
    }

    // Step 1: Make the recursive call to get child row indices
    auto child_node = assignments_.group(column);
    auto rows = child_nodes_[child_node]->get_column(assignments_.rank(column));

    // check if we need to update the row indexes
    if (num_nonzero_rows == nonzero_rows_->size()
            && (!nonones_rows_ || nonones_rows_->size() == nonones_rows_->num_set_bits()))
        return rows;

    // Step 2: Map child row indices to parent ranks (accounting for all-ones rows) and merge with all-ones rows
    // Child row indices are indices in the child's coordinate system, which correspond to
    // ranks in the parent's nonones_rows_ system (0-indexed)
    if (nonones_rows_) {
        std::vector<uint64_t> nnz_rows;
        nnz_rows.reserve(rows.size() + nonones_rows_->size() - nonones_rows_->num_set_bits());

        uint64_t idx = 0;
        uint64_t child_idx = 0;
        auto it = rows.begin();
        nonones_rows_->call_ones([&](uint64_t i) {
            while (idx < i) {
                nnz_rows.push_back(idx++);  // all-ones rows
            }
            if (it != rows.end() && child_idx == *it) {
                nnz_rows.push_back(idx);
                ++it;
            }
            idx++;
            child_idx++;
        });
        assert(it == rows.end());
        while (idx < nonones_rows_->size()) {
            nnz_rows.push_back(idx++);
        }
        rows = nnz_rows;
    }

    // Step 3: Map parent ranks to actual row indices (accounting for zero rows)
    for (size_t i = 0; i < rows.size(); ++i) {
        rows[i] = nonzero_rows_->select1(rows[i] + 1);
    }
    return rows;
}

bool BRWT::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        if (!assignments_.load(in))
            return false;

        if (!nonzero_rows_->load(in))
            return false;

        // Check for version string to determine format (only present when nonones_rows_ exists)
        auto start_pos = in.tellg();
        char version_buf[9];
        in.read(version_buf, 9);
        if (in.gcount() == 9 && std::memcmp(version_buf, "brwt-v2.0", 9) == 0) {
            // New format: version string present, so nonones_rows_ exists
            nonones_rows_ = std::make_unique<bit_vector_smallrank>();
            if (!nonones_rows_->load(in))
                return false;
        } else {
            // old format: no version string, rewind to before the read attempt
            in.clear();
            in.seekg(start_pos);
            nonones_rows_.reset();
        }

        size_t num_child_nodes = load_number(in);
        child_nodes_.clear();
        child_nodes_.reserve(num_child_nodes);
        for (size_t i = 0; i < num_child_nodes; ++i) {
            child_nodes_.emplace_back(new BRWT());
            if (!child_nodes_.back()->load(in))
                return false;
        }
        return !child_nodes_.size()
                    || child_nodes_.size() == assignments_.num_groups();
    } catch (...) {
        return false;
    }
}

void BRWT::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error when dumping BRWT");

    assignments_.serialize(out);

    assert(!child_nodes_.size()
                || child_nodes_.size() == assignments_.num_groups());

    nonzero_rows_->serialize(out);

    // Serialize nonones_rows_ if present
    if (nonones_rows_) {
        // Write version string for new format (only when nonones_rows_ exists)
        const char *version = "brwt-v2.0";
        out.write(version, 9);  // 9 bytes including null terminator
        nonones_rows_->serialize(out);
    }
    // If nonones_rows_ doesn't exist, format is exactly the same as old format
    // (no version string, no sentinel - goes directly to num_child_nodes)

    serialize_number(out, child_nodes_.size());
    for (const auto &child : child_nodes_) {
        child->serialize(out);
    }
}

uint64_t BRWT::num_relations() const {
    if (!child_nodes_.size())
        return nonzero_rows_->num_set_bits();

    uint64_t num_set_bits = 0;
    for (const auto &submatrix_ptr : child_nodes_) {
        num_set_bits += submatrix_ptr->num_relations();
    }
    if (nonones_rows_)
        num_set_bits += (nonones_rows_->size() - nonones_rows_->num_set_bits()) * num_columns();

    return num_set_bits;
}

double BRWT::avg_arity() const {
    if (!child_nodes_.size())
        return 0;

    uint64_t num_nodes = 0;
    uint64_t total_num_child_nodes = 0;

    BFT([&](const BRWT &node) {
        if (node.child_nodes_.size()) {
            num_nodes++;
            total_num_child_nodes += node.child_nodes_.size();
        }
    });

    return num_nodes
            ? static_cast<double>(total_num_child_nodes) / num_nodes
            : 0;
}

uint64_t BRWT::num_nodes() const {
    uint64_t num_nodes = 0;

    BFT([&num_nodes](const BRWT &) { num_nodes++; });

    return num_nodes;
}

double BRWT::shrinking_rate() const {
    double rate_sum = 0;
    uint64_t num_nodes = 0;

    BFT([&](const BRWT &node) {
        if (node.child_nodes_.size()) {
            num_nodes++;
            rate_sum += static_cast<double>(node.nonzero_rows_->num_set_bits())
                            / node.nonzero_rows_->size();
        }
    });

    return rate_sum / num_nodes;
}

void BRWT::print_tree_structure(std::ostream &os) const {
    BFT([&os](const BRWT &node) {
        // print node and its stats
        os << &node << "," << node.nonzero_rows_->size()
                    << "," << node.nonzero_rows_->num_set_bits();
        // print all its children
        for (const auto &child : node.child_nodes_) {
            os << "," << child.get();
        }
        os << std::endl;
    });
}

void BRWT::BFT(std::function<void(const BRWT &node)> callback) const {
    std::queue<const BRWT*> nodes_queue;
    nodes_queue.push(this);

    while (!nodes_queue.empty()) {
        const auto &node = *nodes_queue.front();

        callback(node);

        for (const auto &child_node : node.child_nodes_) {
            const auto *brwt_node_ptr = dynamic_cast<const BRWT*>(child_node.get());
            if (brwt_node_ptr)
                nodes_queue.push(brwt_node_ptr);
        }
        nodes_queue.pop();
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg
