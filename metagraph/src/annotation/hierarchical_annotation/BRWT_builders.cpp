#include "BRWT_builders.hpp"
#include "omp.h"

#include <cmath>


BRWT
BRWTBuilder::initialize(NodeBRWT&& node, bit_vector&& nonzero_rows) {
    BRWT brwt;

    brwt.assignments_ = utils::RangePartition(node.column_arrangement,
                                              node.group_sizes);
    brwt.nonzero_rows_ = nonzero_rows.convert_to<bit_vector_rrr<>>();

    brwt.child_nodes_ = std::move(node.child_nodes);

    return brwt;
}

BRWTBottomUpBuilder::Partitioner
BRWTBottomUpBuilder::get_basic_partitioner(size_t arity) {
    assert(arity > 1u);

    return [arity](const VectorsPtr &vectors) {
        Partition partition((vectors.size() + arity - 1) / arity);
        for (size_t i = 0; i < vectors.size(); ++i) {
            partition[i / arity].push_back(i);
        }
        return partition;
    };
}

sdsl::bit_vector
compute_or(const std::vector<std::unique_ptr<bit_vector>> &columns) {
    uint64_t size = columns.at(0)->size();

    sdsl::bit_vector vector_or(size, false);

    for (const auto &col_ptr : columns) {
        assert(col_ptr.get());
        assert(col_ptr->size() == size);
        if (dynamic_cast<bit_vector_stat*>(&(*col_ptr))) {
            vector_or |= dynamic_cast<bit_vector_stat&>(*col_ptr).get();
        } else {
            col_ptr->call_ones([&vector_or](auto i) { vector_or[i] = true; });
        }
    }
    return vector_or;
}

sdsl::bit_vector generate_subindex(const sdsl::bit_vector &column,
                                   const bit_vector &reference) {
    assert(column.size() == reference.size());

    sdsl::bit_vector subindex(reference.num_set_bits(), false);

    uint64_t j = 0;
    reference.call_ones([&](auto i) {
        if (column[i])
            subindex[j] = true;

        j++;
    });

    return subindex;
}

std::pair<BRWTBottomUpBuilder::NodeBRWT, std::unique_ptr<bit_vector>>
BRWTBottomUpBuilder::merge(std::vector<NodeBRWT> &&nodes,
                           std::vector<std::unique_ptr<bit_vector>> &&index) {
    assert(nodes.size());
    assert(nodes.size() == index.size());

    if (nodes.size() == 1) {
        return { std::move(nodes[0]), std::move(index[0]) };
    }

    NodeBRWT parent;

    // build the parent aggregated column
    auto parent_index = std::make_unique<bit_vector_stat>(compute_or(index));

    // initialize child nodes
    for (size_t i = 0; i < nodes.size(); ++i) {
        // define column assignments for parent
        parent.column_arrangement.insert(parent.column_arrangement.end(),
                                         nodes[i].column_arrangement.begin(),
                                         nodes[i].column_arrangement.end());
        parent.group_sizes.push_back(nodes[i].column_arrangement.size());
        std::iota(nodes[i].column_arrangement.begin(),
                  nodes[i].column_arrangement.end(), 0);

        // shrink index column
        bit_vector_rrr<> shrinked_index(
            generate_subindex(
                index[i]->convert_to<sdsl::bit_vector>(),
                *parent_index
            )
        );
        index[i].reset();

        parent.child_nodes.emplace_back(
            new BRWT(initialize(std::move(nodes[i]),
                                std::move(shrinked_index)))
        );
    }

    return { std::move(parent), std::move(parent_index) };
}

template <typename T>
std::vector<T> subset(std::vector<T> *vector,
                      const std::vector<BRWT::Column> indexes) {
    assert(vector);

    std::vector<T> result;
    result.reserve(indexes.size());
    for (auto j : indexes) {
        result.push_back(std::move((*vector)[j]));
    }
    return result;
}

BRWT BRWTBottomUpBuilder::build(VectorsPtr&& columns,
                                Partitioner partitioner,
                                size_t num_threads) {
    if (!columns.size())
        return BRWT();

    std::vector<NodeBRWT> nodes(columns.size());

    // construct leaves
    for (size_t i = 0; i < columns.size(); ++i) {
        nodes[i].column_arrangement = { i };
        nodes[i].group_sizes = { 1 };
    }

    while (nodes.size() > 1) {
        auto groups = partitioner(columns);

        assert(groups.size() > 0);
        assert(groups.size() < columns.size());

        std::vector<NodeBRWT> parent_nodes(groups.size());
        std::vector<std::unique_ptr<bit_vector>> parent_columns(groups.size());

        #pragma omp parallel for num_threads(num_threads)
        for (size_t g = 0; g < groups.size(); ++g) {
            const auto &group = groups[g];
            assert(group.size());

            auto parent = merge(subset(&nodes, group),
                                subset(&columns, group));

            parent_nodes[g] = std::move(parent.first);
            parent_columns[g] = std::make_unique<bit_vector_rrr<>>(
                parent.second->convert_to<bit_vector_rrr<>>()
            );
        }

        nodes = std::move(parent_nodes);
        columns = std::move(parent_columns);
    }

    return initialize(std::move(nodes.at(0)),
                      std::move(*columns.at(0)));
}


void BRWTOptimizer::relax(BRWT *brwt_matrix,
                          uint64_t max_arity,
                          size_t num_threads) {
    assert(brwt_matrix);

    std::deque<BRWT*> parents;

    // traverse the tree and collect the nodes
    // order: leaves first, root last
    brwt_matrix->BFT([&](const BRWT &node) {
        if (node.child_nodes_.size())
            parents.push_front(const_cast<BRWT*>(&node));
    });

    while (!parents.empty()) {
        auto &parent = *parents.front();
        parents.pop_front();

        assert(parent.child_nodes_.size());

        NodeBRWT updated_parent;

        for (size_t g = 0; g < parent.child_nodes_.size(); ++g) {
            // we leave the column arrangement unchanged
            auto num_columns = parent.child_nodes_[g]->num_columns();
            for (size_t r = 0; r < num_columns; ++r) {
                updated_parent.column_arrangement.push_back(
                    parent.assignments_.get(g, r)
                );
            }
            // column arrangement in the parrent is set
            // now we need to initialize group_sizes and child_nodes
            add_submatrix(
                std::move(parent.child_nodes_[g]),
                &updated_parent,
                max_arity - std::min(max_arity,
                        static_cast<uint64_t>(updated_parent.group_sizes.size())
                            + parent.child_nodes_.size() - g - 1),
                num_threads
            );
        }

        parent = BRWTBuilder::initialize(std::move(updated_parent),
                                         std::move(parent.nonzero_rows_));
    }
}

void BRWTOptimizer::add_submatrix(std::unique_ptr<BinaryMatrix>&& submatrix,
                                  NodeBRWT *parent,
                                  uint64_t max_delta_arity,
                                  size_t num_threads) {
    assert(parent);

    const auto *brwt_cptr = dynamic_cast<const BRWT*>(submatrix.get());

    // we don't relax other submatrix representations
    bool prune = brwt_cptr != NULL;

    // we don't remove leaves in BRWT
    if (prune && !brwt_cptr->child_nodes_.size())
        prune = false;

    // check if the arity will not be exceeded
    if (prune && brwt_cptr->child_nodes_.size() > max_delta_arity)
        prune = false;

    // we relax only those child nodes, that have only BRWT nodes
    // because nodes of other types can't be reassigned
    if (prune) {
        for (const auto &child_node : brwt_cptr->child_nodes_) {
            if (!dynamic_cast<const BRWT *>(child_node.get())) {
                prune = false;
                break;
            }
        }
    }

    // we don't remove nodes if it doesn't reduce the space requirements
    if (!prune || pruning_delta(*brwt_cptr) > 0) {
        // leave this submatrix unchanged
        parent->group_sizes.push_back(submatrix->num_columns());
        parent->child_nodes.push_back(std::move(submatrix));
    } else {
        // at this point we are going to remove this node and
        // reassign all its children to its parent directly
        reassign(std::unique_ptr<BRWT>(dynamic_cast<BRWT*>(submatrix.release())),
                 parent,
                 num_threads);
    }
}

void BRWTOptimizer::reassign(std::unique_ptr<BRWT>&& node,
                             NodeBRWT *parent,
                             size_t num_threads) {
    assert(node);
    assert(parent);

    // TODO: assignemnts?

    const auto node_index = node->nonzero_rows_.convert_to<bit_vector_stat>();

    auto offset = parent->group_sizes.size();

    parent->group_sizes.resize(parent->group_sizes.size() + node->child_nodes_.size());
    parent->child_nodes.resize(parent->child_nodes.size() + node->child_nodes_.size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < node->child_nodes_.size(); ++i) {

        auto submatrix = std::move(node->child_nodes_[i]);

        assert(dynamic_cast<BRWT*>(submatrix.get()));

        std::unique_ptr<BRWT> grand_child {
            dynamic_cast<BRWT*>(submatrix.release())
        };
        sdsl::bit_vector subindex(node_index.size(), false);

        uint64_t child_i = 0;
        node_index.call_ones([&](auto i) {
            if (grand_child->nonzero_rows_[child_i])
                subindex[i] = true;

            child_i++;
        });

        grand_child->nonzero_rows_ = bit_vector_rrr<>(subindex);

        parent->group_sizes[offset + i] = grand_child->num_columns();
        parent->child_nodes[offset + i] = std::move(grand_child);
    }
}

double logbinomial(uint64_t n, uint64_t m) {
    return (lgamma(n + 1)
                - lgamma(m + 1)
                - lgamma(n - m + 1)) / log(2);
}

double entropy(double q) {
    assert(q >= 0);
    assert(q <= 1);

    if (q == 0 || q == 1)
        return 0;

    return q * log2(q) + (1 - q) * log2(1 - q);
}

double bv_space_taken_rrr(uint64_t size,
                          uint64_t num_set_bits,
                          uint8_t block_size) {
    return logbinomial(size, num_set_bits)
            + std::ceil(log2(block_size + 1) / block_size) * size
            + sizeof(bit_vector_rrr<>) * 8;
}

double bv_space_taken(const bit_vector &type,
                      uint64_t size,
                      uint64_t num_set_bits) {
    assert(size >= num_set_bits);

    if (dynamic_cast<const bit_vector_stat *>(&type)) {
        return size + sizeof(bit_vector_stat) * 8;

    } else if (dynamic_cast<const bit_vector_rrr<15> *>(&type)) {
        return bv_space_taken_rrr(size, num_set_bits, 15);

    } else if (dynamic_cast<const bit_vector_rrr<31> *>(&type)) {
        return bv_space_taken_rrr(size, num_set_bits, 31);

    } else if (dynamic_cast<const bit_vector_rrr<63> *>(&type)) {
        return bv_space_taken_rrr(size, num_set_bits, 63);

    } else if (dynamic_cast<const bit_vector_rrr<127> *>(&type)) {
        return bv_space_taken_rrr(size, num_set_bits, 127);

    } else if (dynamic_cast<const bit_vector_rrr<255> *>(&type)) {
        return bv_space_taken_rrr(size, num_set_bits, 255);

    } else {
        throw std::runtime_error("Error: unknown space taken for this bit_vector");
    }
}

double BRWTOptimizer::pruning_delta(const BRWT &node) {
    assert(node.child_nodes_.size());

    double delta = 0;

    for (const auto &submatrix : node.child_nodes_) {
        auto *brwt_child = dynamic_cast<const BRWT*>(submatrix.get());
        assert(brwt_child);

        assert(brwt_child->num_rows() <= node.num_rows());
        assert(brwt_child->nonzero_rows_.size() <= node.num_rows());

        // updated vector
        delta += bv_space_taken(brwt_child->nonzero_rows_,
                                node.num_rows(),
                                brwt_child->nonzero_rows_.num_set_bits());

        // old index vector
        delta -= bv_space_taken(brwt_child->nonzero_rows_,
                                brwt_child->nonzero_rows_.size(),
                                brwt_child->nonzero_rows_.num_set_bits());
    }

    // removed index vector
    delta -= bv_space_taken(node.nonzero_rows_,
                            node.nonzero_rows_.size(),
                            node.nonzero_rows_.num_set_bits());

    return delta;
}
