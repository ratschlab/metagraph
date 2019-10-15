#include "BRWT.hpp"

#include <queue>
#include <numeric>

#include "serialization.hpp"


bool BRWT::get(Row row, Column column) const {
    assert(row < num_rows());
    assert(column < num_columns());

    // terminate if the index bit is unset
    if (!nonzero_rows_[row])
        return false;

    // return true if this is a leaf
    if (!child_nodes_.size())
        return true;

    auto child_node = assignments_.group(column);
    return child_nodes_[child_node]->get(nonzero_rows_.rank1(row) - 1,
                                         assignments_.rank(column));
}

std::vector<BRWT::Column> BRWT::get_row(Row row) const {
    assert(row < num_rows());

    // check if the row is empty
    if (!nonzero_rows_[row])
        return {};

    // check whether it is a leaf
    if (!child_nodes_.size()) {
        assert(assignments_.size() == 1);

        // the bit is set
        return utils::arange<Column>(0, assignments_.size());
    }

    // check all child nodes
    std::vector<Column> row_set_bits;
    uint64_t index_in_child = nonzero_rows_.rank1(row) - 1;

    for (size_t i = 0; i < child_nodes_.size(); ++i) {
        const auto &child = *child_nodes_[i];

        for (auto col_id : child.get_row(index_in_child)) {
            row_set_bits.push_back(assignments_.get(i, col_id));
        }
    }
    return row_set_bits;
}

std::vector<std::vector<BRWT::Column>>
BRWT::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<std::vector<Column>> rows(row_ids.size());

    // check whether it is a leaf
    if (!child_nodes_.size()) {
        assert(assignments_.size() == 1);

        for (size_t i = 0; i < row_ids.size(); ++i) {
            assert(row_ids[i] < num_rows());

            if (nonzero_rows_[row_ids[i]])
                rows[i] = utils::arange<Column>(0, assignments_.size());
        }

        return rows;
    }

    // construct indexing for children and the inverse mapping
    std::vector<Row> child_row_ids;
    child_row_ids.reserve(row_ids.size());

    std::vector<Row> from_child_to_parent;
    from_child_to_parent.reserve(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        assert(row_ids[i] < num_rows());

        // check if the row is not empty
        if (nonzero_rows_[row_ids[i]]) {
            // map index from parent's to children's coordinate system
            child_row_ids.push_back(nonzero_rows_.rank1(row_ids[i]) - 1);
            from_child_to_parent.push_back(i);
        }
    }

    // query all children subtrees
    for (size_t j = 0; j < child_nodes_.size(); ++j) {
        auto child_rows = child_nodes_[j]->get_rows(child_row_ids);

        // push rows from children back to |rows|
        for (size_t i = 0; i < child_rows.size(); ++i) {
            auto &row = rows[from_child_to_parent[i]];
            auto &child_row = child_rows[i];

            for (auto child_col_id : child_row) {
                row.push_back(assignments_.get(j, child_col_id));
            }
        }
    }

    return rows;
}

std::vector<BRWT::Row> BRWT::get_column(Column column) const {
    assert(column < num_columns());

    auto num_nonzero_rows = nonzero_rows_.num_set_bits();

    // check if the column is empty
    if (!num_nonzero_rows)
        return {};

    // check whether it is a leaf
    if (!child_nodes_.size()) {
        // return the index column
        std::vector<BRWT::Row> result;
        result.reserve(num_nonzero_rows);
        nonzero_rows_.call_ones([&](auto i) { result.push_back(i); });
        return result;
    }

    auto child_node = assignments_.group(column);
    auto rows = child_nodes_[child_node]->get_column(assignments_.rank(column));

    // check if we need to update the row indexes
    if (num_nonzero_rows == nonzero_rows_.size())
        return rows;

    // shift indexes
    for (size_t i = 0; i < rows.size(); ++i) {
        rows[i] = nonzero_rows_.select1(rows[i] + 1);
    }
    return rows;
}

bool BRWT::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        if (!assignments_.load(in))
            return false;

        if (!nonzero_rows_.load(in))
            return false;

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

    nonzero_rows_.serialize(out);

    serialize_number(out, child_nodes_.size());
    for (const auto &child : child_nodes_) {
        child->serialize(out);
    }
}

uint64_t BRWT::num_relations() const {
    if (!child_nodes_.size())
        return nonzero_rows_.num_set_bits();

    uint64_t num_set_bits = 0;
    for (const auto &submatrix_ptr : child_nodes_) {
        num_set_bits += submatrix_ptr->num_relations();
    }

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
            rate_sum += static_cast<double>(node.nonzero_rows_.num_set_bits())
                            / node.nonzero_rows_.size();
        }
    });

    return rate_sum / num_nodes;
}

uint64_t BRWT::total_column_size() const {
    uint64_t total_size = 0;

    BFT([&](const BRWT &node) {
        total_size += node.nonzero_rows_.size();
    });

    return total_size;
}

uint64_t BRWT::total_num_set_bits() const {
    uint64_t total_num_set_bits = 0;

    BFT([&](const BRWT &node) {
        total_num_set_bits += node.nonzero_rows_.num_set_bits();
    });

    return total_num_set_bits;
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
